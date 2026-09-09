#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
abo_panel.py — shared variant-panel loader and scoring engine for the
nf-core/abotyper pipeline (v2.0.0).

This module is the ONLY place that understands the panel file format
(abo_variant_panel.yaml / .yml / .json / .csv / .tsv). All four
position-consuming scripts (pysam_haploscan.py, stats_from_pileup.py,
predict_abo_phenotype.py, aggregate_abo_reports.py) import this module
instead of hardcoding position dictionaries.

Adding a new diagnostic variant to the assay requires editing ONLY the
panel file — no script needs to change.

Supported panel file formats (auto-detected from extension):
    .yaml / .yml   - human-authored source of truth (recommended for editing)
    .json          - same schema as YAML, machine-friendly
    .csv / .tsv    - flat table, convenient for spreadsheet review/edit

Use `python3 abo_panel.py --export-csv out.csv --panel abo_variant_panel.yaml`
to regenerate the flat table from the YAML source after editing.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Dict, List, Optional, Union

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

__version__ = "2.0.0"

# ---------------------------------------------------------------------------
# Flat field list — defines both the CSV column order and the VariantMarker
# schema. Keep this list in sync with abo_variant_panel.yaml row keys.
# ---------------------------------------------------------------------------
FIELD_NAMES = [
    "id", "cdna_pos", "cdna_change", "exon", "ref_base", "alt_base",
    "variant_type", "legacy_exon6_pos", "legacy_exon7_pos",
    "ng006669_2_pos", "amplicon_pos",
    "category", "call_rule", "ref_call_label", "alt_call_label",
    "ref_threshold_pct", "alt_threshold_pct", "het_band_pct",
    "variant_threshold_pct", "min_reads", "associated_alleles",
    "interpretation_ref", "interpretation_alt", "notes", "source",
]

_INT_FIELDS = {
    "cdna_pos", "legacy_exon6_pos", "legacy_exon7_pos",
    "ng006669_2_pos", "amplicon_pos", "min_reads",
}
_FLOAT_FIELDS = {
    "ref_threshold_pct", "alt_threshold_pct", "het_band_pct",
    "variant_threshold_pct",
}


def _coerce(key: str, value):
    """Coerce a raw (possibly stringified) value to its proper type."""
    if value is None or value == "":
        return None
    if key in _INT_FIELDS:
        try:
            return int(value)
        except (ValueError, TypeError):
            return None
    if key in _FLOAT_FIELDS:
        try:
            return float(value)
        except (ValueError, TypeError):
            return None
    return str(value)


@dataclass
class VariantMarker:
    """One diagnostic position from the panel."""

    id: str
    cdna_pos: Optional[int] = None
    cdna_change: str = ""
    exon: str = ""
    ref_base: str = ""
    alt_base: str = ""
    variant_type: str = "snp"
    legacy_exon6_pos: Optional[int] = None
    legacy_exon7_pos: Optional[int] = None
    ng006669_2_pos: Optional[int] = None
    amplicon_pos: Optional[int] = None
    category: str = ""
    call_rule: str = "named_marker"
    ref_call_label: str = ""
    alt_call_label: str = ""
    ref_threshold_pct: Optional[float] = None
    alt_threshold_pct: Optional[float] = None
    het_band_pct: Optional[float] = None
    variant_threshold_pct: Optional[float] = None
    min_reads: Optional[int] = None
    associated_alleles: str = ""
    interpretation_ref: str = ""
    interpretation_alt: str = ""
    notes: str = ""
    source: str = ""

    # -- convenience -------------------------------------------------
    def ng_ref_base_matches(self, ng_seq: str) -> Optional[bool]:
        """Sanity-check helper: does ng006669_2_pos in `ng_seq` actually carry
        this row's expected ref_base? Returns None if not checkable."""
        if self.ng006669_2_pos is None or self.ref_base in ("REF", "ALT", ""):
            return None
        idx = self.ng006669_2_pos - 1
        if idx < 0 or idx >= len(ng_seq):
            return None
        return ng_seq[idx].upper() == self.ref_base.upper()

    def resolved_position(self, prefer: str = "auto") -> Optional[int]:
        """
        Return the position to use for pileup/BAM lookups in THIS run.

        prefer: "auto" (amplicon_pos if set, else ng006669_2_pos if set,
                else legacy for the matching exon, else None), "amplicon",
                "ng006669_2", "legacy", or "legacy6"/"legacy7"

        ng006669_2_pos takes priority over the legacy exon6/7 coordinates
        because it is a position in the single combined reference
        (NG_006669.2) that the single-reference pipeline topology aligns
        against directly, whereas legacy_exon6_pos/legacy_exon7_pos are
        only valid coordinates in the old, separate mini-amplicon
        references and must not be used when reads are aligned to
        NG_006669.2 (or any other combined reference).
        """
        if prefer in ("auto", "amplicon") and self.amplicon_pos is not None:
            return self.amplicon_pos
        if prefer == "amplicon":
            return None
        if prefer in ("auto", "ng006669_2") and self.ng006669_2_pos is not None:
            return self.ng006669_2_pos
        if prefer == "ng006669_2":
            return None
        if self.exon == "Exon 6" or prefer == "legacy6":
            return self.legacy_exon6_pos
        if self.exon == "Exon 7" or prefer == "legacy7":
            return self.legacy_exon7_pos
        # Positions outside exon6/7 have no legacy coordinate by definition.
        return None

    def column_label(self) -> str:
        """
        Column header used in aggregate_abo_reports.py output. Preserves
        the exact legacy naming ("Exon7_pos422") for backward compatibility
        with any LIS/Excel integration already built around it; falls back
        to a cDNA-based label for positions that don't have one yet.
        """
        exon_num = "".join(ch for ch in self.exon if ch.isdigit()) or "?"
        pos = self.resolved_position()
        if pos is not None:
            return f"Exon{exon_num}_pos{pos}"
        return f"Exon{exon_num}_c{self.cdna_pos}"

    def is_calibrated(self) -> bool:
        return self.resolved_position() is not None

    def alt_bases(self) -> List[str]:
        """
        Some rows encode more than one alt allele at the same position
        (e.g. c.700C>G / c.700C>T resolving to different subtypes). This
        splits alt_base on '_or_' into a clean list; single-allele rows
        return a 1-item list.
        """
        raw = self.alt_base or ""
        return [b.strip() for b in raw.split("_or_") if b.strip()]


@dataclass
class StructuralVariant:
    id: str
    cdna_change: str = ""
    exon: str = ""
    size_bp_approx: Optional[int] = None
    associated_alleles: str = ""
    notes: str = ""
    source: str = ""


class Panel:
    """A loaded, indexed variant panel."""

    def __init__(self, variants: List[VariantMarker],
                 structural_variants: Optional[List[StructuralVariant]] = None,
                 meta: Optional[dict] = None):
        self.variants = variants
        self.structural_variants = structural_variants or []
        self.meta = meta or {}
        self._by_id = {v.id: v for v in variants}

    # -- lookups -------------------------------------------------------
    def get(self, marker_id: str) -> Optional[VariantMarker]:
        return self._by_id.get(marker_id)

    def by_exon(self, exon_label: str) -> List[VariantMarker]:
        return [v for v in self.variants if v.exon == exon_label]

    def exons(self) -> List[str]:
        seen = []
        for v in self.variants:
            if v.exon not in seen:
                seen.append(v.exon)
        return seen

    def by_category(self, category: str) -> List[VariantMarker]:
        return [v for v in self.variants if v.category == category]

    def categories(self) -> List[str]:
        seen = []
        for v in self.variants:
            if v.category not in seen:
                seen.append(v.category)
        return seen

    def primary_markers(self) -> List[VariantMarker]:
        return self.by_category("primary")

    def calibrated(self) -> List[VariantMarker]:
        return [v for v in self.variants if v.is_calibrated()]

    def uncalibrated(self) -> List[VariantMarker]:
        return [v for v in self.variants if not v.is_calibrated()]

    def indel_diagnostic_positions(self, exon_label: str) -> List[int]:
        """Positions (resolved) at THIS exon where indels are the diagnostic
        variant itself (not sequencing noise) — used by
        stats_from_pileup.py / pysam_haploscan.py to decide whether to fold
        indel counts into the denominator."""
        out = []
        for v in self.by_exon(exon_label):
            if v.variant_type in ("deletion", "insertion", "indel",
                                   "snp_or_indel", "dup_or_del"):
                pos = v.resolved_position()
                if pos is not None:
                    out.append(pos)
        return out


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _load_yaml_or_json(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in (".yaml", ".yml"):
        if yaml is None:
            raise RuntimeError(
                "pyyaml is required to read .yaml panels. "
                "Install with: pip install pyyaml"
            )
        return yaml.safe_load(text)
    return json.loads(text)


def _load_table(path: Path) -> dict:
    """Load a flat CSV/TSV panel into the same nested dict shape the
    YAML/JSON loader produces."""
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = [dict(row) for row in reader]
    # Structural variants, if present, are distinguished by an
    # is_structural=1 column (added by export_csv). Otherwise all rows
    # are treated as scorable variants.
    variants = [r for r in rows if str(r.get("is_structural", "")).strip()
                not in ("1", "true", "True")]
    svs = [r for r in rows if str(r.get("is_structural", "")).strip()
           in ("1", "true", "True")]
    for r in variants + svs:
        r.pop("is_structural", None)
    return {"meta": {}, "variants": variants, "structural_variants": svs}


def load_panel(path: Union[str, Path]) -> Panel:
    """Load a variant panel from any supported file format."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Panel file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in (".yaml", ".yml", ".json"):
        raw = _load_yaml_or_json(path)
    elif suffix in (".csv", ".tsv"):
        raw = _load_table(path)
    else:
        raise ValueError(
            f"Unsupported panel format '{suffix}'. Use .yaml/.yml/.json/.csv/.tsv"
        )

    variants = []
    for row in raw.get("variants", []):
        clean = {k: _coerce(k, v) for k, v in row.items() if k in FIELD_NAMES}
        # Required minimum fields
        clean.setdefault("id", row.get("id", f"unnamed_{len(variants)}"))
        for f_ in fields(VariantMarker):
            if f_.name not in clean:
                clean[f_.name] = f_.default if not isinstance(f_.default, type(field())) else None
        variants.append(VariantMarker(**{k: clean.get(k) for k in
                                          [f_.name for f_ in fields(VariantMarker)]}))

    svs = []
    for row in raw.get("structural_variants", []) or []:
        clean = {k: row.get(k) for k in
                 ["id", "cdna_change", "exon", "size_bp_approx",
                  "associated_alleles", "notes", "source"]}
        if clean.get("size_bp_approx") not in (None, ""):
            try:
                clean["size_bp_approx"] = int(clean["size_bp_approx"])
            except (ValueError, TypeError):
                clean["size_bp_approx"] = None
        svs.append(StructuralVariant(**clean))

    return Panel(variants=variants, structural_variants=svs, meta=raw.get("meta", {}))


# ---------------------------------------------------------------------------
# Export: YAML/JSON -> flat CSV (for spreadsheet-based review/editing)
# ---------------------------------------------------------------------------

def export_csv(panel: Panel, out_path: Union[str, Path]) -> None:
    out_path = Path(out_path)
    columns = FIELD_NAMES + ["is_structural"]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for v in panel.variants:
            row = {k: getattr(v, k) for k in FIELD_NAMES}
            row["is_structural"] = 0
            writer.writerow({k: ("" if val is None else val) for k, val in row.items()})
        for sv in panel.structural_variants:
            row = {k: "" for k in FIELD_NAMES}
            row.update({
                "id": sv.id, "cdna_change": sv.cdna_change, "exon": sv.exon,
                "associated_alleles": sv.associated_alleles, "notes": sv.notes,
                "source": sv.source,
            })
            row["is_structural"] = 1
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Scoring engines — shared by predict_abo_phenotype.py / aggregate_abo_reports.py
# ---------------------------------------------------------------------------

def _pct_lookup(row_values: Dict[str, float], base: str) -> float:
    """
    row_values holds the percentage columns available for a pileup row:
    A, G, C, T, Del, Ins (case-insensitive keys expected). base may be a
    nucleotide, 'del', 'ins', 'dup', 'REF', or 'ALT'/'ALT' placeholders
    used by primary_biallelic rows whose alt is a mixed snp/indel event.
    """
    if base is None:
        return 0.0
    b = base.strip().lower()
    mapping = {
        "del": "del", "delg": "del", "dup": "ins", "ins": "ins",
    }
    key = mapping.get(b, base.strip().upper())
    for k, v in row_values.items():
        if k.strip().lower() == str(key).strip().lower():
            try:
                return float(v)
            except (TypeError, ValueError):
                return 0.0
    return 0.0


def call_primary_biallelic(marker: VariantMarker, row_values: Dict[str, float]) -> str:
    """
    Generic replacement for the bespoke per-position branches in the
    legacy get_type()/get_type_exon6(). Reproduces the original asymmetric
    homozygous/heterozygous thresholds, now driven by the panel row.
    """
    alt_pct = _pct_lookup(row_values, marker.alt_base)
    ref_pct = 100.0 - alt_pct if marker.ref_base in ("REF",) else _pct_lookup(
        row_values, marker.ref_base) or (100.0 - alt_pct)

    alt_thr = marker.alt_threshold_pct if marker.alt_threshold_pct is not None else 80
    ref_thr = marker.ref_threshold_pct if marker.ref_threshold_pct is not None else 80
    band = marker.het_band_pct if marker.het_band_pct is not None else 20

    if alt_pct >= alt_thr and alt_pct > ref_pct:
        return marker.alt_call_label
    if ref_pct >= ref_thr and ref_pct > alt_pct:
        return marker.ref_call_label
    if abs(alt_pct - ref_pct) <= band or (
        (100 - band) > ref_pct > band and (100 - band) > alt_pct > band
    ):
        if marker.ref_call_label and marker.alt_call_label:
            return f"{marker.ref_call_label} and {marker.alt_call_label}"
    return ""


def call_named_marker(marker: VariantMarker, row_values: Dict[str, float],
                       nreads: float = 0) -> Optional[Dict[str, object]]:
    """
    Generic replacement for scan_bw_markers()/scan_a2_markers()'s per-marker
    branches. Returns a dict describing the observed call, or None if the
    marker doesn't fire. Handles rows with multiple possible alt alleles
    (e.g. c.700 C>G vs C>T) by checking each independently.
    """
    threshold = marker.variant_threshold_pct if marker.variant_threshold_pct is not None else 15
    min_reads = marker.min_reads if marker.min_reads is not None else 0

    if min_reads and nreads and nreads < min_reads:
        return {"fired": False, "reason": "low_coverage"}

    for alt in marker.alt_bases():
        pct = _pct_lookup(row_values, alt)
        if pct >= threshold:
            confidence = "Confirmed" if pct >= 80 else "Possible"
            return {
                "fired": True,
                "alt": alt,
                "pct": pct,
                "confidence": confidence,
                "label": marker.alt_call_label or marker.category,
                "text": f"{confidence} {marker.alt_call_label or marker.category}"
                        f"({marker.cdna_change}:{alt}={pct:.0f}%)",
            }
    return {"fired": False}


# ---------------------------------------------------------------------------
# CLI: panel validation / CSV export helper
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Validate an ABO variant panel and/or export it to CSV.",
    )
    ap.add_argument("--panel", required=True, help="Path to panel (.yaml/.json/.csv/.tsv)")
    ap.add_argument("--export-csv", help="Write a flat CSV copy of the panel here")
    ap.add_argument("--export-yaml", help="Write a YAML copy of the panel here")
    ap.add_argument("--summary", action="store_true", help="Print a coverage summary")
    args = ap.parse_args()

    panel = load_panel(args.panel)
    print(f"Loaded {len(panel.variants)} variants and "
          f"{len(panel.structural_variants)} structural variants from {args.panel}")

    if args.summary or not (args.export_csv or args.export_yaml):
        print("\nBy exon:")
        for exon in panel.exons():
            print(f"  {exon}: {len(panel.by_exon(exon))}")
        print("\nBy category:")
        for cat in panel.categories():
            print(f"  {cat}: {len(panel.by_category(cat))}")
        n_cal = len(panel.calibrated())
        n_unc = len(panel.uncalibrated())
        print(f"\nCalibrated positions (usable now): {n_cal}")
        print(f"Uncalibrated positions (awaiting calibrate_panel_positions.py "
              f"or a legacy exon6/7 coordinate): {n_unc}")
        if n_unc:
            print("  " + ", ".join(v.id for v in panel.uncalibrated()))

    if args.export_csv:
        export_csv(panel, args.export_csv)
        print(f"Wrote CSV table -> {args.export_csv}")

    if args.export_yaml:
        data = {
            "meta": panel.meta,
            "variants": [
                {k: getattr(v, k) for k in FIELD_NAMES} for v in panel.variants
            ],
            "structural_variants": [
                {k: getattr(sv, k) for k in
                 ["id", "cdna_change", "exon", "size_bp_approx",
                  "associated_alleles", "notes", "source"]}
                for sv in panel.structural_variants
            ],
        }
        if yaml is None:
            raise RuntimeError("pyyaml required for --export-yaml")
        with open(args.export_yaml, "w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)
        print(f"Wrote YAML -> {args.export_yaml}")


if __name__ == "__main__":
    main()
