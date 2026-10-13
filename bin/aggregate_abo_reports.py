#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import re
import sys
import glob
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd
from xlsxwriter.utility import xl_col_to_name

from abo_panel import load_panel, Panel, VariantMarker, call_named_marker

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024-2025, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = [
    "Fredrick Mobegi",
    "Benedict Matern",
    "Mathijs Groeneweg",
    "Claude Sonnet 5",
]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


"""
ABO Blood Group Report Aggregator.

Column layout, per-position parsing, and subtype-marker scanning are driven
by an external variant panel (--panel, default abo_variant_panel.yaml) via
abo_panel.py, so adding a diagnostic position or subtype category requires
editing the panel file only, not this script.

Note: c.657C>T (legacy pos283) is NOT specific to ABO*BW.06 -- it is part of
the shared B-lineage 7-SNP backbone present in nearly all B/B3/Bweak/Bel/
cisAB/BA alleles. Bw.06's actual private variant is c.1036A>G (not yet in
the exon6/7 amplicon's position set -- add via the panel once calibrated
against the combined amplicon).

QC / failure capture (v2.1.0)
-----------------------------
Calls are NOT changed. Each sample is annotated with four extra Notes
columns -- A2Trace, EvidenceBasis, QCStatus, QCFlags -- and the run prints
(and writes to ABO_run_qc.txt) a QC summary. The checks come from an audit of
a 6,218-sample run and from reading the upstream converters:

  * A2 traceability. c.467C>T (SNP) and c.1061delC (indel) are the two
    primary A2 markers. If only one is present, the A subtype is traced as
    "partial" and flagged. Observed: 84 samples had c.467 variant without
    c.1061delC; both known A2 donors that came back as plain "A" had exactly
    this pattern. Het c.467 without the deletion falls through
    determine_a_subtype() to a bare "A"; hom c.467 without the deletion is
    reported as "A1". Both are flagged.
  * A subtype contradicted by evidence (A1 with an A2-lineage marker, A2 with
    no primary A2 marker) and A1 calls made on very few reads.
  * Missing evidence. A position absent from the report is zero-filled, and
    clair2metrics.py >= 1.1.0 writes an all-zero row for a position nothing
    covers. Before 2.1.0 such a row was typed as a call (heterozygous at a
    SNP primary, "variant" at a named marker) and ReadReliability ignored
    zero-read positions. Now a zero row gives no call at that position (so a
    missing primary makes the genotype Unknown), ReadReliability is "Very
    Low", the position is listed in QCFlags and EvidenceBasis = "no data".
  * Reference-only calls. A genotype made only because no variant was seen at
    any primary position (EvidenceBasis = "reference-only") cannot be told
    apart from a failed call. clair2metrics.py writes 100% reference with the
    sample median depth when a position has neither a variant record nor a
    gVCF reference block, so this needs the run-level checks below.
  * Marker inconsistency: c.796 and c.803 are both B-defining and must
    agree; "Unknown" genotypes get an explicit reason.
  * Replicates. Samples sharing a Sequencing_ID must agree on genotype and on
    A2 evidence; disagreements mark every replicate for review. In the audit
    85 of 110 discordant replicate groups had one A1/A1 replicate.
  * #Reads caveat. In the audited run #Reads was identical at every primary
    position of every sample (a sample-wide read count, not depth at the
    site), so ReadReliability and the low-read flag cannot confirm
    site-level evidence. The run summary warns when this is the case.

Part of the nf-core/abotyper pipeline: https://github.com/fmobegi/nf-core-abotyper
"""


# ===========================================================================
# Phase Confidence Scoring — reads Haplotypes.tsv from pysam_haploscan.py
# ===========================================================================


@dataclass
class PhaseEvidence:
    """Summary of within-amplicon phase evidence for a sample+exon."""

    total_reads: int = 0
    cis_confirmed: Dict[str, int] = field(default_factory=dict)
    trans_confirmed: Dict[str, int] = field(default_factory=dict)
    ambiguous: Dict[str, int] = field(default_factory=dict)
    haplotype_counts: Dict[str, int] = field(default_factory=dict)


def parse_haplotypes_tsv(haplo_file: str) -> Optional[PhaseEvidence]:
    """Parse a Haplotypes.tsv file produced by pysam_haploscan.py."""
    if not haplo_file or not os.path.exists(haplo_file):
        return None

    try:
        df = pd.read_csv(haplo_file, sep="\t")
        if df.empty:
            return None

        evidence = PhaseEvidence(total_reads=len(df))

        if "Haplotype" in df.columns:
            for haplo_str, count in df["Haplotype"].value_counts().items():
                evidence.haplotype_counts[str(haplo_str)] = int(count)

        return evidence

    except Exception:
        return None


def compute_phase_confidence(
    evidence_exons: List[Optional[PhaseEvidence]],
    extended_genotype: str,
) -> Tuple[str, str]:
    """
    Compute phase confidence from haplotype tables.

    Accepts a list of PhaseEvidence (one per exon covered by the assay)
    rather than a fixed exon6/exon7 pair, and picks whichever has the most
    reads. Returns (confidence_level, detail_string).
    """
    evidence_list = [e for e in evidence_exons if e and e.total_reads > 0]
    if not evidence_list:
        return "N/A", "No haplotype data available"

    if "/" in extended_genotype:
        alleles = extended_genotype.split("/")
        if len(alleles) == 2 and alleles[0] == alleles[1]:
            return "High", "Homozygous; phasing trivial"

    evidence = max(evidence_list, key=lambda e: e.total_reads)

    total = evidence.total_reads
    haplotypes = evidence.haplotype_counts

    if not haplotypes:
        return "N/A", "No haplotype calls in data"

    sorted_haplos = sorted(haplotypes.items(), key=lambda x: x[1], reverse=True)

    top2_count = sum(count for _, count in sorted_haplos[:2])
    top2_fraction = top2_count / total if total > 0 else 0

    n_distinct = len([h for h, c in sorted_haplos if c >= max(3, total * 0.02)])

    if top2_fraction >= 0.90 and n_distinct <= 3:
        detail = f"Top 2 haplotypes cover {top2_fraction*100:.0f}% of {total} reads"
        return "High", detail
    elif top2_fraction >= 0.75 and n_distinct <= 5:
        detail = f"Top 2 haplotypes cover {top2_fraction*100:.0f}% of {total} reads ({n_distinct} distinct patterns)"
        return "Moderate", detail
    else:
        detail = f"High haplotype diversity: {n_distinct} patterns, top 2 cover only {top2_fraction*100:.0f}%"
        return "Low", detail


def phase_relationship(evidence: Optional[PhaseEvidence], pos_a: int, pos_b: int) -> str:
    """
    Whether two panel positions' ALT alleles sit on the same physical
    haplotype (cis) or different ones (trans), read directly off the
    haplotype-string keys already parsed into evidence.haplotype_counts
    (each string lists every ALT-carrying panel position on that side,
    e.g. "p23952=T;p24544=del" -- see ReadHaplotype.to_string() /
    clair2haplotypes.py). Advisory only: this does not feed genotype
    scoring, it only informs warning text.

    Returns "cis", "trans", "ambiguous", or "unknown" (not enough phased
    evidence covering both positions to say anything).
    """
    if not evidence or not evidence.haplotype_counts:
        return "unknown"

    tag_a, tag_b = f"p{pos_a}=", f"p{pos_b}="
    cis_reads = a_only_reads = b_only_reads = 0
    for haplo_str, count in evidence.haplotype_counts.items():
        has_a = tag_a in haplo_str
        has_b = tag_b in haplo_str
        if has_a and has_b:
            cis_reads += count
        elif has_a:
            a_only_reads += count
        elif has_b:
            b_only_reads += count

    total = cis_reads + a_only_reads + b_only_reads
    if total == 0:
        return "unknown"
    if cis_reads / total >= 0.7:
        return "cis"
    if a_only_reads > 0 and b_only_reads > 0 and cis_reads / total < 0.2:
        return "trans"
    return "ambiguous"


# ===========================================================================
# Main ABOReportParser class
# ===========================================================================


class ABOReportParser:
    """
    Collates all ABO phenotype results from each sample into a general table
    and generates an Excel worksheet and a CSV file for export to LIS soft
    or other general purpose lab management systems. Column layout, position
    parsing, and subtype marker scanning are driven by an external variant
    panel (see abo_panel.py).
    """

    def __init__(self, input_dir, panel: Panel, default_barcode="barcode00",
                 min_depth=100, qc_report="ABO_run_qc.txt"):
        self.input_dir = input_dir
        self.panel = panel
        self.default_barcode = default_barcode
        self.results = []
        self.failed_samples = []
        self.min_depth = min_depth            # reads below which an A1 call is flagged
        self.qc_report_path = qc_report       # run-level QC text report ("" disables)
        self.qc_records = []                  # one dict per sample, for the run summary
        self.replicate_stats = {"groups": 0, "discordant_genotype": 0, "discordant_a2": 0}

        # Positions actually usable in THIS run (resolved via amplicon_pos
        # if calibrated, else legacy exon6/exon7 coordinate). Anything else
        # in the panel is silently skipped -- it will start being scored
        # automatically the moment the panel is calibrated for it.
        self.primary_markers: List[VariantMarker] = [
            m for m in panel.primary_markers() if m.resolved_position() is not None
        ]
        self.subtype_markers_by_exon: Dict[str, List[VariantMarker]] = {}
        self.exon_position_lists: Dict[str, List[int]] = {}
        for exon_label in panel.exons():
            resolved = [m for m in panel.by_exon(exon_label) if m.resolved_position() is not None]
            if resolved:
                self.exon_position_lists[exon_label] = sorted(
                    {m.resolved_position() for m in resolved}
                )
                self.subtype_markers_by_exon[exon_label] = [
                    m for m in resolved if m.category != "primary"
                ]

        self.initialize_columns()

        self.pattern_with_barcode = re.compile(
            r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[A-Z0-9]+)?(-[A-Z0-9]+)?_barcode\d+$",
            re.IGNORECASE,
        )
        self.pattern_without_barcode = re.compile(
            r"^[A-Z0-9]+(_[A-Z0-9]+)*$", re.IGNORECASE
        )

        self.processing_stats = {
            "with_barcode": 0,
            "without_barcode": 0,
            "pattern_matched": 0,
            "pattern_failed": 0,
        }

    # -----------------------------------------------------------------
    # Column layout
    # -----------------------------------------------------------------
    def initialize_columns(self):
        """Build the MultiIndex column layout dynamically from the panel."""
        header_cols = ["", ""]  # Sample: Barcode, Sequencing_ID
        header_cols += ["Result"] * 3  # Phenotype, Genotype, ExtendedGenotype

        self.position_columns: List[str] = []  # ordered list of column labels
        for exon_label, positions in self.exon_position_lists.items():
            for pos in positions:
                col_label = self._column_label_for(exon_label, pos)
                self.position_columns.append(col_label)
                header_cols += [col_label] * 10

        self.notes_columns = ["ASubtype", "ReadReliability", "PhaseConfidence", "BwSubtype",
                              "A2Trace", "EvidenceBasis", "QCStatus", "QCFlags"]
        header_cols += ["Notes"] * len(self.notes_columns)

        column_metrics = ["#Reads", "Mat", "Mis", "Ins", "Del", "A", "G", "C", "T", "Type"]
        header_rows = (
            ["Barcode", "Sequencing_ID"]
            + ["Phenotype", "Genotype", "ExtendedGenotype"]
            + column_metrics * len(self.position_columns)
            + self.notes_columns
        )

        self.columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

    def _column_label_for(self, exon_label: str, pos: int) -> str:
        exon_num = "".join(ch for ch in exon_label if ch.isdigit()) or "?"
        return f"Exon{exon_num}_pos{pos}"

    def _marker_for_column(self, col_label: str) -> Optional[VariantMarker]:
        for m in self.panel.variants:
            if m.resolved_position() is not None and m.column_label() == col_label:
                return m
        return None

    # -----------------------------------------------------------------
    # Sample/filename handling
    # -----------------------------------------------------------------
    def extract_sample_info(self, filename):
        """Extract sample name and barcode from filename."""
        if "_barcode" in filename:
            parts = filename.rsplit("_barcode", 1)
            sample_name = parts[0]
            try:
                barcode_num = parts[1]
                barcode_int = int(barcode_num)
                if 0 <= barcode_int <= 99:
                    barcode = f"barcode{barcode_num.zfill(2)}"
                else:
                    print(f"Warning: Unusual barcode number {barcode_int} for {filename}")
                    barcode = f"barcode{barcode_num}"
            except ValueError:
                print(f"Warning: Invalid barcode format in {filename}, using as-is")
                barcode = f"barcode{parts[1]}"
            return sample_name, barcode, "explicit"
        else:
            sample_name = filename
            barcode = self.default_barcode
            print(f"No barcode found in '{filename}', using default barcode: {barcode}")
            return sample_name, barcode, "default"

    # -----------------------------------------------------------------
    # Per-exon report parser
    # -----------------------------------------------------------------
    def parse_exon_report(self, filename: str, exon_label: str) -> pd.DataFrame:
        """Parse an *.ABOPhenotype.txt report section for one exon, using
        the panel's position list for that exon."""
        all_positions = self.exon_position_lists.get(exon_label, [])
        empty_cols = ["Exon", "Position", "#Reads", "Mat", "Mis", "Ins", "Del", "A", "G", "C", "T"]

        try:
            with open(filename, "r", encoding="utf-8") as f:
                lines = f.readlines()

            exon_num = "".join(ch for ch in exon_label if ch.isdigit()) or "?"
            header_token = f"{exon_label} position(1-based):"

            positions, counts = [], []
            mat_values, mis_values, ins_values, del_values = [], [], [], []
            a_values, g_values, c_values, t_values = [], [], [], []

            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if header_token in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))
                                    stats_idx = j + 2
                                    if stats_idx < len(lines) and "Mat" in lines[stats_idx - 1]:
                                        stats = lines[stats_idx].split()
                                        if len(stats) >= 8:
                                            mat, mis, ins, dele, a, g, c, t = [
                                                float(x) for x in stats
                                            ]
                                            positions.append(pos)
                                            counts.append(count)
                                            mat_values.append(mat)
                                            mis_values.append(mis)
                                            ins_values.append(ins)
                                            del_values.append(dele)
                                            a_values.append(a)
                                            g_values.append(g)
                                            c_values.append(c)
                                            t_values.append(t)
                                            break
                i += 1

            df = pd.DataFrame({
                "Exon": [exon_num] * len(positions),
                "Position": positions,
                "#Reads": counts,
                "Mat": mat_values, "Mis": mis_values, "Ins": ins_values, "Del": del_values,
                "A": a_values, "G": g_values, "C": c_values, "T": t_values,
            })

            for pos in all_positions:
                if pos not in df["Position"].values:
                    df = pd.concat([df, pd.DataFrame({
                        "Exon": [exon_num], "Position": [pos], "#Reads": [0],
                        "Mat": [0], "Mis": [0], "Ins": [0], "Del": [0],
                        "A": [0], "G": [0], "C": [0], "T": [0],
                    })], ignore_index=True)

            df = df.sort_values("Position").reset_index(drop=True)
            df["Type"] = df.apply(
                lambda row: self.get_type_generic(exon_label, int(row["Position"]), row),
                axis=1,
            )
            return df

        except Exception as e:
            print(f"Error parsing {exon_label} file {filename}: {str(e)}")
            empty_df = pd.DataFrame(columns=empty_cols)
            exon_num = "".join(ch for ch in exon_label if ch.isdigit()) or "?"
            for pos in all_positions:
                empty_df = pd.concat([empty_df, pd.DataFrame({
                    "Exon": [exon_num], "Position": [pos], "#Reads": [0],
                    "Mat": [0], "Mis": [0], "Ins": [0], "Del": [0],
                    "A": [0], "G": [0], "C": [0], "T": [0], "Type": [""],
                })], ignore_index=True)
            return empty_df

    # -----------------------------------------------------------------
    # Type caller
    # -----------------------------------------------------------------

    @staticmethod
    def _no_evidence_row(row_values: dict) -> bool:
        """True when a position carries no evidence at all: Mat, Mis, Ins and Del are all zero.
        That is what clair2metrics.py (>=1.1.0) writes for a position nothing covers, and what
        parse_exon_report() zero-fills for a position missing from the report. Typing such a row
        would invent a call: both the ref and alt fractions are 0, so the biallelic caller reads
        it as heterozygous, and the named-marker caller reads "100% non-reference" as a variant."""
        keys = ("Mat", "Mis", "Ins", "Del")
        if not all(k in row_values for k in keys):
            return False
        try:
            return sum(float(row_values.get(k) or 0) for k in keys) <= 0
        except (TypeError, ValueError):
            return False

    def get_type_generic(self, exon_label: str, pos: int, row: pd.Series) -> str:
        """Determine the blood-type label for a single position, using the
        panel's call_rule for that position."""
        marker = None
        for m in self.panel.by_exon(exon_label):
            if m.resolved_position() == pos:
                marker = m
                break
        if marker is None:
            return ""

        row_values = {
            "A": row.get("A", 0), "G": row.get("G", 0),
            "C": row.get("C", 0), "T": row.get("T", 0),
            "Del": row.get("Del", 0), "Ins": row.get("Ins", 0),
            "Mat": row.get("Mat", 0), "Mis": row.get("Mis", 0),
        }
        if self._no_evidence_row(row_values):
            return ""

        if marker.call_rule == "primary_biallelic":
            return self._call_primary_biallelic_row(marker, row_values)

        if marker.call_rule == "named_marker":
            # The per-position Type cell just wants "variant" / "" -- the
            # fully-annotated marker text is produced separately by
            # scan_category_markers() for the Notes/BwSubtype/ASubtype columns.
            threshold = marker.variant_threshold_pct if marker.variant_threshold_pct is not None else 25
            non_ref_pct = 100.0 - float(row_values.get(marker.ref_base, 0) or 0) \
                if marker.ref_base in ("A", "G", "C", "T") else \
                sum(float(row_values.get(b, 0) or 0) for b in marker.alt_bases())
            if non_ref_pct >= threshold:
                return marker.alt_call_label or "variant"
            return ""

        return ""

    def _call_primary_biallelic_row(self, marker: VariantMarker, row_values: dict) -> str:
        """Row-level primary_biallelic caller (mirrors abo_panel.call_primary_biallelic
        but works off a pandas row's already-extracted dict for speed/clarity here)."""
        def pct_of(base):
            if base is None:
                return 0.0
            b = base.strip().lower()
            key_map = {"del": "Del", "ins": "Ins", "dup": "Ins"}
            key = key_map.get(b, base.strip().upper())
            return float(row_values.get(key, 0) or 0)

        alt_pct = pct_of(marker.alt_base)
        ref_pct = pct_of(marker.ref_base) if marker.ref_base not in ("REF", "") else (100.0 - alt_pct)

        alt_thr = marker.alt_threshold_pct if marker.alt_threshold_pct is not None else 80
        ref_thr = marker.ref_threshold_pct if marker.ref_threshold_pct is not None else 80
        band = marker.het_band_pct if marker.het_band_pct is not None else 20

        if alt_pct >= alt_thr and alt_pct > ref_pct:
            return marker.alt_call_label
        if ref_pct >= ref_thr and ref_pct > alt_pct:
            return marker.ref_call_label
        if abs(alt_pct - ref_pct) <= band or (band < ref_pct < 100 - band and band < alt_pct < 100 - band):
            if marker.ref_call_label and marker.alt_call_label:
                return f"{marker.ref_call_label} and {marker.alt_call_label}"
        return ""

    def _primary_state(self, marker: VariantMarker, row_values: dict) -> str:
        """Classify a primary_biallelic marker into 'ref' | 'alt' | 'het' | 'none'.

        Uses the SAME thresholds as _call_primary_biallelic_row, but returns a
        semantic state instead of a display string. This decouples the
        phenotype decision tree from how the panel happens to word its
        ref_call_label/alt_call_label, so re-labelling the panel can never
        again silently break genotype assignment.
        """

        if self._no_evidence_row(row_values):
            return "none"

        def pct_of(base):
            if base is None:
                return 0.0
            b = base.strip().lower()
            key_map = {"del": "Del", "ins": "Ins", "dup": "Ins"}
            key = key_map.get(b, base.strip().upper())
            return float(row_values.get(key, 0) or 0)

        alt_pct = pct_of(marker.alt_base)
        ref_pct = pct_of(marker.ref_base) if marker.ref_base not in ("REF", "") else (100.0 - alt_pct)

        alt_thr = marker.alt_threshold_pct if marker.alt_threshold_pct is not None else 80
        ref_thr = marker.ref_threshold_pct if marker.ref_threshold_pct is not None else 80
        band = marker.het_band_pct if marker.het_band_pct is not None else 20

        if alt_pct >= alt_thr and alt_pct > ref_pct:
            return "alt"
        if ref_pct >= ref_thr and ref_pct > alt_pct:
            return "ref"
        if abs(alt_pct - ref_pct) <= band or (band < ref_pct < 100 - band and band < alt_pct < 100 - band):
            return "het"
        return "none"

    # -----------------------------------------------------------------
    # Subtype marker scanner
    # -----------------------------------------------------------------

    def scan_category_markers(
        self,
        category: str,
        row_lookup: Dict[int, pd.Series],
        require_b_allele: bool = False,
        require_a_allele: bool = False,
        type_exon7_422: str = "",
        type_exon7_429: str = "",
    ) -> Tuple[List[str], List[str]]:
        """
        Works for ANY panel category (bw, ba, cisab, a2, a3, am, aweak, ael,
        bweak, bel, b3, onull, ...). Returns (fired_marker_texts, warnings).
        """
        fired, warnings = [], []

        if require_b_allele:
            has_b = ("B" in type_exon7_422) or ("B" in type_exon7_429)
            if not has_b:
                return fired, warnings
        if require_a_allele:
            has_a = ("A" in type_exon7_422) or ("A" in type_exon7_429)
            if not has_a:
                return fired, warnings

        for marker in self.panel.by_category(category):
            pos = marker.resolved_position()
            if pos is None or pos not in row_lookup:
                continue
            row = row_lookup[pos]
            row_values = {
                "A": row.get("A", 0), "G": row.get("G", 0),
                "C": row.get("C", 0), "T": row.get("T", 0),
                "Del": row.get("Del", 0), "Ins": row.get("Ins", 0),
            }
            nreads = row.get("#Reads", 0)
            result = call_named_marker(marker, row_values, nreads=nreads)
            if result and result.get("fired"):
                fired.append(result["text"])
            elif result and result.get("reason") == "low_coverage":
                warnings.append(
                    f"Low coverage at {marker.cdna_change} — {marker.category} "
                    f"subtyping unreliable"
                )

        return fired, warnings

    # -----------------------------------------------------------------
    # Phenotype/genotype assignment
    # -----------------------------------------------------------------

    def assign_phenotype_genotype(
        self,
        df,
        phase_evidence_by_exon: Optional[Dict[str, PhaseEvidence]] = None,
    ):
        """Assign phenotype/genotype information with panel-driven subtype
        marker scanning and phase confidence."""
        phase_evidence_by_exon = phase_evidence_by_exon or {}
        # Same "most reads wins" selection compute_phase_confidence uses --
        # in combined mode every exon key holds the same PhaseEvidence
        # object anyway (one shared Haplotypes.tsv per sample).
        _evidence_candidates = [e for e in phase_evidence_by_exon.values() if e and e.total_reads > 0]
        resolved_phase_evidence = max(_evidence_candidates, key=lambda e: e.total_reads) if _evidence_candidates else None
        try:

            def safe_get_type(df, pos_key, default=""):
                try:
                    return df.at[0, (pos_key, "Type")]
                except (KeyError, IndexError):
                    return default

            def safe_get_reads(df, pos_key, default=0):
                try:
                    value = df.at[0, (pos_key, "#Reads")]
                    return value if pd.notna(value) else default
                except (KeyError, IndexError):
                    return default

            def safe_get_row(df, pos_key) -> pd.Series:
                try:
                    return pd.Series({
                        "#Reads": df.at[0, (pos_key, "#Reads")],
                        "A": df.at[0, (pos_key, "A")], "G": df.at[0, (pos_key, "G")],
                        "C": df.at[0, (pos_key, "C")], "T": df.at[0, (pos_key, "T")],
                        "Del": df.at[0, (pos_key, "Del")],
                        "Mat": df.at[0, (pos_key, "Mat")], "Mis": df.at[0, (pos_key, "Mis")],
                        "Ins": df.at[0, (pos_key, "Ins")],
                    })
                except (KeyError, IndexError):
                    return pd.Series({"#Reads": 0, "A": 0, "G": 0, "C": 0, "T": 0, "Del": 0,
                                      "Mat": 0, "Mis": 0, "Ins": 0})

            # ----- Resolve the primary markers by their semantic role -----
            # (column label depends on resolved_position(), so look it up
            # dynamically rather than assuming a literal string.)

            def col_for(marker_id: str) -> Optional[str]:
                m = self.panel.get(marker_id)
                if m is None or m.resolved_position() is None:
                    return None
                return m.column_label()

            col_o1 = col_for("o1_marker")
            col_796 = col_for("b_vs_ao_796")
            col_802 = col_for("o2_marker_802")
            col_803 = col_for("ao_vs_b_803")
            col_804 = col_for("o34_homopolymer_804")
            col_467 = col_for("a1_a2_467")
            col_1061 = col_for("a1_a2_1061del")

            # Normalize the primary markers into the exact legacy vocabulary
            # the decision tree below matches on. The panel's display labels
            # (ref_call_label/alt_call_label) are worded differently from the
            # tree's hardcoded literals (e.g. panel "A or B or O" vs tree
            # "O and (A or B)"; het "A or B or O and O1" vs tree
            # "O1 and (A or B or O)"), so matching the raw Type cell makes
            # every branch fail and every sample fall through to "Unknown".
            # We classify each marker's biallelic STATE and emit the legacy
            # string, so display wording and decision logic stay decoupled.

            def primary_type(col, marker_id, ref_s, alt_s, het_s):
                m = self.panel.get(marker_id)
                if m is None or not col:
                    return ""
                state = self._primary_state(m, dict(safe_get_row(df, col)))
                return {"ref": ref_s, "alt": alt_s, "het": het_s, "none": ""}[state]

            type_exon6 = primary_type(
                col_o1, "o1_marker", "A or B or O", "O1", "O1 and (A or B or O)")
            type_exon7_422 = primary_type(
                col_796, "b_vs_ao_796", "A or O", "B", "(A or O) and B")
            type_exon7_428 = primary_type(
                col_802, "o2_marker_802", "O and (A or B)", "O2", "O2 and (O or A or B)")
            type_exon7_429 = primary_type(
                col_803, "ao_vs_b_803", "A or O", "B", "(A or O) and B")
            type_exon7_431 = primary_type(
                col_804, "o34_homopolymer_804", "O and (A or B)", "O3", "O3 and (O or A or B)")
            type_exon7_93 = safe_get_type(df, col_467) if col_467 else ""

            def primary_state(col, marker_id):
                m = self.panel.get(marker_id)
                if m is None or not col:
                    return "none"
                return self._primary_state(m, dict(safe_get_row(df, col)))

            state_exon7_685 = primary_state(col_1061, "a1_a2_1061del")

            nreads6 = safe_get_reads(df, col_o1) if col_o1 else 0
            nreads_exon7_p422 = safe_get_reads(df, col_796) if col_796 else 0
            nreads_exon7_p428 = safe_get_reads(df, col_802) if col_802 else 0
            nreads_exon7_p429 = safe_get_reads(df, col_803) if col_803 else 0
            nreads_exon7_p431 = safe_get_reads(df, col_804) if col_804 else 0

            # ----- Build a position -> row lookup per exon for the generic
            # subtype scanner (covers ALL categories, not just A2/Bw) -----
            row_lookup: Dict[int, pd.Series] = {}
            for exon_label, positions in self.exon_position_lists.items():
                for pos in positions:
                    col_label = self._column_label_for(exon_label, pos)
                    row_lookup[pos] = safe_get_row(df, col_label)

            Phenotype = "Unknown"
            Genotype = "Unknown"
            ExtendedGenotype = "Unknown"
            BwSubtype = ""

            def determine_a2_subtype(markers_found: List[str]):
                has_del = "c.1061del" in markers_found
                has_907 = "c.907A" in markers_found
                has_1032 = "c.1032A" in markers_found
                has_297 = "c.297G" in markers_found
                if has_907:
                    return "A2.06", None
                if has_1032 and has_del:
                    return "A2.01", "c.1032G>A"
                if has_297 and has_del:
                    return "A2.01", "c.297A>G"
                if has_del:
                    return "A2.01", None
                if has_1032:
                    return "A2.01", "c.1032G>A"
                if has_297:
                    return "A2.01", "c.297A>G"
                return "A2", None

            def scan_a2_markers():
                """Scans c.1061del / c.907A / c.1032A / c.297G / c.266T /
                c.268C / cXXXvar A2 markers, reading thresholds/positions
                from the panel."""
                markers, warns = [], []

                # c.1061delC is heterozygous in every AO/AB sample carrying a
                # single A2.01 allele, so it reports as "A1 and A2", never
                # the bare alt label "A2" -- check the semantic state (alt OR
                # het), not the ref-only-vs-alt-only display string.
                has_del = state_exon7_685 in ("alt", "het") and nreads_exon7_p422 >= 30
                a907_marker = self.panel.get("a2p_907_a206")
                a1032_marker = self.panel.get("a2p_1032_a201")
                has_907 = False
                has_1032 = False
                if a907_marker and a907_marker.resolved_position() in row_lookup:
                    r = call_named_marker(a907_marker, dict(row_lookup[a907_marker.resolved_position()]),
                                           nreads=row_lookup[a907_marker.resolved_position()].get("#Reads", 0))
                    has_907 = bool(r and r.get("fired")) and nreads_exon7_p422 >= 30
                if a1032_marker and a1032_marker.resolved_position() in row_lookup:
                    r = call_named_marker(a1032_marker, dict(row_lookup[a1032_marker.resolved_position()]),
                                           nreads=row_lookup[a1032_marker.resolved_position()].get("#Reads", 0))
                    has_1032 = bool(r and r.get("fired")) and nreads_exon7_p422 >= 30

                a297_marker = self.panel.get("a2_marker_297_a201")
                has_297 = False
                a297_phase_note = None
                if a297_marker and a297_marker.resolved_position() in row_lookup:
                    row = row_lookup[a297_marker.resolved_position()]
                    has_297 = float(row.get("G", 0) or 0) >= 25
                    if has_297:
                        # c.297A>G is dual-purpose (panel notes): also part of
                        # the shared B-lineage backbone. If phasing shows it
                        # riding on the SAME haplotype as the B primary
                        # marker, it's likely just B-backbone noise, not
                        # independently confirmatory for A2.01 -- advisory
                        # only, does not change has_297 or any call below.
                        b_marker = self.panel.get("b_vs_ao_796")
                        if b_marker and b_marker.resolved_position() is not None:
                            rel = phase_relationship(
                                resolved_phase_evidence,
                                a297_marker.resolved_position(),
                                b_marker.resolved_position(),
                            )
                            if rel == "cis":
                                a297_phase_note = (
                                    "Phase check: c.297A>G co-occurs with the B allele on the "
                                    "same haplotype (Clair3 phasing) -- likely B-lineage backbone "
                                    "signal, not independently confirmatory for A2.01"
                                )

                if has_del:
                    markers.append("c.1061del")
                if has_907:
                    markers.append("c.907A")
                if has_1032:
                    markers.append("c.1032A")
                if has_297:
                    markers.append("c.297G")
                if a297_phase_note:
                    warns.append(a297_phase_note)

                o1v_fired = has_297 and not has_del and not has_907 and not has_1032
                if o1v_fired:
                    for m_ in ("c.297G", "c.266T", "c.268C"):
                        if m_ in markers:
                            markers.remove(m_)
                    warns.append(
                        "Note: c.297G detected without c.1061del — possible A2-derived "
                        "O allele (O1v); manual review recommended"
                    )

                if not o1v_fired:
                    m266 = self.panel.get("a2_marker_266")
                    m268 = self.panel.get("a2_marker_268")
                    for m_marker, label in [(m266, "c.266T"), (m268, "c.268C")]:
                        if m_marker and m_marker.resolved_position() in row_lookup:
                            row = row_lookup[m_marker.resolved_position()]
                            row_values = {"A": row.get("A", 0), "G": row.get("G", 0),
                                          "C": row.get("C", 0), "T": row.get("T", 0)}
                            r = call_named_marker(m_marker, row_values, nreads=row.get("#Reads", 0))
                            if r and r.get("fired"):
                                markers.append(label)

                if nreads_exon7_p422 >= 500 and not o1v_fired:
                    for m_ in self.panel.by_category("a2_panel"):
                        if m_.id in ("a2p_907_a206", "a2p_1032_a201"):
                            continue  # scored separately above
                        pos = m_.resolved_position()
                        if pos is None or pos not in row_lookup:
                            continue
                        row = row_lookup[pos]
                        row_values = {"A": row.get("A", 0), "G": row.get("G", 0),
                                      "C": row.get("C", 0), "T": row.get("T", 0)}
                        r = call_named_marker(m_, row_values, nreads=row.get("#Reads", 0))
                        if r and r.get("fired"):
                            label = "c." + str(m_.cdna_pos) + "var"
                            markers.append(label)

                return markers, warns

            def determine_a_subtype():
                """Returns (subtype_label, warning_text, zygosity).

                zygosity is "het" (one A1 allele + one subtype allele) or
                "hom" (both alleles carry the subtype) -- callers that
                construct a homozygous-looking ExtendedGenotype (the AA
                branch) MUST check this before writing f"{subtype}/{subtype}",
                otherwise a genuinely heterozygous A1/A2.01 sample gets
                silently doubled into a false A2.01/A2.01 homozygote. This
                previously wasn't tracked at all -- scan_a2_markers() only
                returned *whether* a subtype marker fired, never whether it
                fired at ~50% (het) or ~90%+ (hom) allele fraction.
                """
                markers, warns = scan_a2_markers()
                if markers:
                    clean_subtype, nt_notation = determine_a2_subtype(markers)
                    if nt_notation:
                        warns.insert(0, nt_notation)
                    warning_text = "; ".join(warns) if warns else None

                    zygosity = "hom"
                    if "c.1061del" in markers:
                        if state_exon7_685 == "het":
                            zygosity = "het"
                    else:
                        # Subtype was called from a non-primary marker
                        # (c.907A/c.1032A/c.297G) -- check ITS allele
                        # fraction directly via the same ref/alt/het
                        # classifier used for primary markers, since
                        # scan_a2_markers()/call_named_marker() only report
                        # whether the marker fired, not its zygosity.
                        driving_marker_id = {
                            "c.907A": "a2p_907_a206",
                            "c.1032A": "a2p_1032_a201",
                            "c.297G": "a2_marker_297_a201",
                        }
                        for label, marker_id in driving_marker_id.items():
                            if label not in markers:
                                continue
                            m = self.panel.get(marker_id)
                            if m and m.resolved_position() in row_lookup:
                                state = self._primary_state(
                                    m, dict(row_lookup[m.resolved_position()]))
                                if state == "het":
                                    zygosity = "het"
                            break

                    return clean_subtype, warning_text, zygosity
                warning_text = "; ".join(warns) if warns else None
                if type_exon7_93 in ("A1.02 or A2", "A1"):
                    return "A1", warning_text, "hom"
                return "", None, "hom"

            # ----- PART 1: PRIMARY PHENOTYPING LOGIC -----
            a_subtype_warning = None

            if (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O and (A or B)"):
                a_subtype, a_subtype_warning, _ = determine_a_subtype()
                Phenotype, Genotype = "A", "AO"
                ExtendedGenotype = f"{a_subtype}/O1" if a_subtype else "A/O1"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O and (A or B)"):
                a_subtype, a_subtype_warning, _ = determine_a_subtype()
                Phenotype, Genotype = "A", "AO"
                ExtendedGenotype = f"{a_subtype}/O2" if a_subtype else "A/O2"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O3 and (O or A or B)"):
                a_subtype, a_subtype_warning, _ = determine_a_subtype()
                Phenotype, Genotype = "A", "AO"
                ExtendedGenotype = f"{a_subtype}/O3" if a_subtype else "A/O3"

            elif (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "(A or O) and B"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "(A or O) and B"
                    and type_exon7_431 == "O and (A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "B", "BO", "B/O1"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "(A or O) and B"
                    and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "(A or O) and B"
                    and type_exon7_431 == "O and (A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "B", "BO", "B/O2"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "(A or O) and B"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "(A or O) and B"
                    and type_exon7_431 == "O3 and (O or A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "B", "BO", "B/O3"

            elif (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O and (A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O1/O2"

            elif (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O3 and (O or A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O1/O3"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O3 and (O or A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O2/O3"

            elif (type_exon6 == "O1" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O and (A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O1/O1"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O2" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O and (A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O2/O2"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O3"):
                Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O3/O3"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
                    and type_exon7_431 == "O and (A or B)"):
                a_subtype, a_subtype_warning, a_subtype_zygosity = determine_a_subtype()
                if a_subtype and a_subtype_zygosity == "het":
                    # One A1 allele (reference/silent at every A2 marker) +
                    # one subtype allele -- NOT two copies of the subtype.
                    Phenotype, Genotype, ExtendedGenotype = "A", "AA", f"A1/{a_subtype}"
                elif a_subtype:
                    Phenotype, Genotype, ExtendedGenotype = a_subtype, "AA", f"{a_subtype}/{a_subtype}"
                else:
                    Phenotype, Genotype, ExtendedGenotype = "A", "AA", "A/A"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "B"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "B"
                    and type_exon7_431 == "O and (A or B)"):
                Phenotype, Genotype, ExtendedGenotype = "B", "BB", "B/B"

            elif (type_exon6 == "A or B or O" and type_exon7_422 == "(A or O) and B"
                    and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "(A or O) and B"
                    and type_exon7_431 == "O and (A or B)"):
                a_subtype, a_subtype_warning, _ = determine_a_subtype()
                Phenotype, Genotype = "AB", "AB"
                ExtendedGenotype = f"{a_subtype}/B" if a_subtype else "A/B"

            else:
                Phenotype, Genotype, ExtendedGenotype = "Unknown", "Unknown", "Unknown"

            # ----- PART 2: subtype-marker scan across ALL panel categories -----
            bw_markers, bw_warnings = self.scan_category_markers(
                "bw", row_lookup, require_b_allele=True,
                type_exon7_422=type_exon7_422, type_exon7_429=type_exon7_429,
            )
            ba_markers, ba_warnings = self.scan_category_markers(
                "ba", row_lookup,
                type_exon7_422=type_exon7_422, type_exon7_429=type_exon7_429,
            )
            cisab_markers, cisab_warnings = self.scan_category_markers(
                "cisab", row_lookup,
                type_exon7_422=type_exon7_422, type_exon7_429=type_exon7_429,
            )

            other_subtype_categories = [
                c for c in self.panel.categories()
                if c not in ("primary", "a2_panel", "a2", "bw", "ba", "cisab")
            ]
            other_markers: List[str] = []
            other_warnings: List[str] = []
            for cat in other_subtype_categories:
                fired, warns = self.scan_category_markers(cat, row_lookup)
                other_markers.extend(fired)
                other_warnings.extend(warns)

            all_subtype_texts = bw_markers + ba_markers + cisab_markers + other_markers
            all_subtype_warnings = bw_warnings + ba_warnings + cisab_warnings + other_warnings

            if all_subtype_texts:
                BwSubtype = "; ".join(all_subtype_texts)
                if "B" in Phenotype or "B" in Genotype:
                    BwSubtype = "weak-B/subtype flags: " + BwSubtype
            elif all_subtype_warnings:
                BwSubtype = "; ".join(all_subtype_warnings)

            # ----- PART 3: PHASE CONFIDENCE -----
            phase_level, phase_detail = compute_phase_confidence(
                list(phase_evidence_by_exon.values()), ExtendedGenotype
            )
            PhaseConfidence = f"{phase_level}: {phase_detail}"

            # ----- PART 4: RELIABILITY SCORING -----
            read_counts = [nreads6, nreads_exon7_p422, nreads_exon7_p428,
                           nreads_exon7_p429, nreads_exon7_p431]
            valid_read_counts = []
            for count in read_counts:
                try:
                    if pd.notna(count) and float(count) > 0:
                        valid_read_counts.append(float(count))
                except (ValueError, TypeError):
                    continue

            if valid_read_counts:
                min_reads = min(valid_read_counts)
                if min_reads <= 20:
                    Reliability = "Very Low(<=20 reads)"
                elif min_reads <= 40:
                    Reliability = "Low (<=40 reads)"
                elif min_reads >= 500:
                    Reliability = "Robust(>=500 reads)"
                else:
                    Reliability = "Normal"
                if max(valid_read_counts) / min(valid_read_counts) > 5:
                    Reliability += " (Variable coverage)"
            else:
                Reliability = "Unknown (no read data)"

            # ----- PART 5: QC / FAILURE CAPTURE (v2.1.0) -----
            # Annotates the call; never changes it. See the module docstring.
            qc_review: List[str] = []
            _copies = {"ref": 0, "het": 1, "hom": 2}

            def cell(col, key):
                try:
                    v = df.at[0, (col, key)]
                    return float(v) if pd.notna(v) else 0.0
                except (KeyError, IndexError, TypeError, ValueError):
                    return 0.0

            def alt_frac(marker_id):
                m = self.panel.get(marker_id)
                if m is None or m.resolved_position() is None:
                    return None
                col = m.column_label()
                best = 0.0
                for tok in m.alt_bases():
                    t = tok.strip().lower()
                    key = "Del" if t == "del" else ("Ins" if t in ("ins", "dup") else tok.strip().upper())
                    best = max(best, cell(col, key))
                return best

            def vstate(frac):
                if frac is None:
                    return None
                return "hom" if frac >= 80 else ("het" if frac >= 20 else "ref")

            s261 = vstate(alt_frac("o1_marker"))
            s796 = vstate(alt_frac("b_vs_ao_796"))
            s802 = vstate(alt_frac("o2_marker_802"))
            s803 = vstate(alt_frac("ao_vs_b_803"))
            s467 = vstate(alt_frac("a1_a2_467"))
            s1061 = vstate(alt_frac("a1_a2_1061del"))
            primary_states = [s for s in (s261, s796, s802, s803, s467, s1061) if s is not None]

            # evidence actually present at each primary position
            no_evidence, primary_reads, unscored = [], [], []
            for pid, short in (("o1_marker", "c.261delG"), ("b_vs_ao_796", "c.796C>A"),
                               ("o2_marker_802", "c.802G>A"), ("ao_vs_b_803", "c.803G>C"),
                               ("a1_a2_467", "c.467C>T"), ("a1_a2_1061del", "c.1061delC")):
                m = self.panel.get(pid)
                if m is None or m.resolved_position() is None:
                    continue
                c = m.column_label()
                reads = cell(c, "#Reads")
                primary_reads.append(reads)
                if reads <= 0 or (cell(c, "Mat") + cell(c, "Mis") + cell(c, "Ins") + cell(c, "Del")) <= 0:
                    no_evidence.append(short)
                else:
                    # non-reference reads that are NOT the panel's alt allele (e.g. a different base at
                    # c.803). 100 - Mat is the non-reference fraction in both the pileup and the Clair3
                    # metrics formats (Mis is not: it already includes Ins/Del in the Clair3 format).
                    nonref = 100.0 - cell(c, "Mat")
                    scored = alt_frac(pid) or 0.0
                    if nonref - scored >= 20.0:
                        unscored.append(f"{short} (~{int(round(nonref - scored))}% of reads)")
            reads_uniform = len(primary_reads) > 1 and len(set(primary_reads)) == 1
            if unscored:
                qc_review.append("allele not scored by the panel at " + ", ".join(unscored))
            if no_evidence:
                qc_review.append("no evidence at primary position(s) " + ", ".join(no_evidence)
                                 + " (no call made there)")

            tokens = [t for t in ExtendedGenotype.split("/")] if ExtendedGenotype != "Unknown" else []
            has_A = any(t.startswith("A") for t in tokens)

            # --- A2 trace ---
            A2Trace = "n/a (no A allele)"
            a2_missing = [x for x in ("c.467C>T", "c.1061delC") if x in no_evidence]
            if has_A and a2_missing:
                A2Trace = "5 A2 unknown: no data at " + " and ".join(a2_missing)
                if any(t == "A" for t in tokens):
                    qc_review.append("A subtype unresolved (generic A)")
            elif has_A and s467 is not None and s1061 is not None:
                e467, e1061 = s467 != "ref", s1061 != "ref"
                if e467 and e1061:
                    A2Trace = "1 A2 confirmed (c.467 + c.1061delC)"
                elif e467:
                    A2Trace = "2 A2 partial: c.467 variant, c.1061delC not called"
                    qc_review.append("A2 partial: c.467 variant without c.1061delC (possible missed deletion or A1.02)")
                elif e1061:
                    A2Trace = "3 A2 partial: c.1061delC variant, c.467 not called"
                    qc_review.append("A2 partial: c.1061delC variant without c.467")
                else:
                    A2Trace = "4 No A2 evidence (both primary markers reference)"
                a2_tokens = [t for t in tokens if t.startswith("A2")]
                a1_tokens = [t for t in tokens if t == "A1"]
                if any(t == "A" for t in tokens):
                    qc_review.append("A subtype unresolved (generic A)")
                if a1_tokens and not a2_tokens and not A2Trace.startswith("4"):
                    qc_review.append("reported A1 despite A2-lineage marker")
                if a2_tokens and A2Trace.startswith("4"):
                    qc_review.append("reported A2 without a primary A2 marker")
                if A2Trace.startswith("4"):
                    a2_reads = [cell(self.panel.get(p).column_label(), "#Reads")
                                for p in ("a1_a2_467", "a1_a2_1061del")
                                if self.panel.get(p) is not None and self.panel.get(p).resolved_position() is not None]
                    if a2_reads and min(a2_reads) < self.min_depth:
                        kind = "same count at every primary position" if reads_uniform else "site-level"
                        qc_review.append(f"A1 call on low read count ({int(min(a2_reads))} < {self.min_depth}, {kind})")

            # --- marker consistency ---
            b_pair_disagree = s796 is not None and s803 is not None and s796 != s803
            if b_pair_disagree:
                qc_review.append(f"B markers disagree (c.796 = {s796}, c.803 = {s803})")
            if Genotype == "Unknown":
                b_copies = max(_copies.get(s796, 0), _copies.get(s803, 0))
                n_alleles = _copies.get(s261, 0) + _copies.get(s802, 0) + b_copies
                if n_alleles > 2:
                    qc_review.append(f"Unknown genotype: more than two alleles (O1 x{_copies.get(s261, 0)}, "
                                     f"O2 x{_copies.get(s802, 0)}, B x{b_copies})")
                elif not b_pair_disagree:
                    qc_review.append("Unknown genotype: primary marker pattern not recognised")

            # --- evidence basis ---
            if no_evidence:
                EvidenceBasis = "no data"
            elif unscored:
                EvidenceBasis = "unscored allele"
            elif primary_states and all(s == "ref" for s in primary_states):
                EvidenceBasis = "reference-only"
            else:
                EvidenceBasis = "variant evidence"

            QCStatus = "REVIEW" if qc_review else "PASS"
            QCFlags = "; ".join(qc_review)
            try:
                _sid = str(df.at[0, ("", "Sequencing_ID")])
            except (KeyError, IndexError):
                _sid = ""
            self.qc_records.append({
                "sample": _sid, "genotype": Genotype, "ext": ExtendedGenotype, "a2trace": A2Trace,
                "evidence": EvidenceBasis, "s261": s261, "reads_uniform": reads_uniform,
            })

            # A primary position with zero reads is missing evidence, not a position to skip.
            _zero = [r for r in primary_reads if r <= 0]
            if _zero:
                Reliability = "Very Low(<=20 reads)"

            ASubtype = a_subtype_warning or ""

            if any("cisAB" in w for w in all_subtype_warnings):
                Reliability += "; " + "; ".join(w for w in all_subtype_warnings if "cisAB" in w)

            df[("Result", "Phenotype")] = Phenotype
            df[("Result", "Genotype")] = Genotype
            df[("Result", "ExtendedGenotype")] = ExtendedGenotype
            df[("Notes", "ASubtype")] = ASubtype
            df[("Notes", "ReadReliability")] = Reliability
            df[("Notes", "PhaseConfidence")] = PhaseConfidence
            df[("Notes", "BwSubtype")] = BwSubtype
            df[("Notes", "A2Trace")] = A2Trace
            df[("Notes", "EvidenceBasis")] = EvidenceBasis
            df[("Notes", "QCStatus")] = QCStatus
            df[("Notes", "QCFlags")] = QCFlags

            return df

        except Exception as e:
            print(f"Error in assign_phenotype_genotype: {str(e)}")
            import traceback
            traceback.print_exc()
            df[("Result", "Phenotype")] = "Error"
            df[("Result", "Genotype")] = "Error"
            df[("Result", "ExtendedGenotype")] = "Error"
            df[("Notes", "ASubtype")] = ""
            df[("Notes", "ReadReliability")] = "Error processing"
            df[("Notes", "PhaseConfidence")] = "Error"
            df[("Notes", "BwSubtype")] = "Error"
            df[("Notes", "A2Trace")] = ""
            df[("Notes", "EvidenceBasis")] = "error"
            df[("Notes", "QCStatus")] = "REVIEW"
            df[("Notes", "QCFlags")] = "processing error"
            return df

    # -----------------------------------------------------------------
    # File discovery
    # -----------------------------------------------------------------
    #
    # Two on-disk layouts are supported: "combined" mode -- one
    # sample_dir/combined/*.ABOPhenotype.txt (+ *.Haplotypes.tsv) covering
    # every exon, produced against a single combined exon2-7 reference --
    # and legacy per-exon mode -- one sample_dir/<exon_dirname>/ subfolder
    # per exon, produced by the original dual-mini-amplicon topology.
    # "combined" is tried first; if absent, the legacy layout is used.
    _COMBINED_DIRNAME = "combined"

    def _find_haplotype_file(self, sample_dir: str, exon_dirname: str) -> Optional[str]:
        exon_dir = os.path.join(sample_dir, exon_dirname)
        if not os.path.isdir(exon_dir):
            return None
        haplo_files = glob.glob(os.path.join(exon_dir, "*.Haplotypes.tsv"))
        if haplo_files:
            return haplo_files[0]
        haplo_files = glob.glob(os.path.join(exon_dir, "*", "*.Haplotypes.tsv"))
        if haplo_files:
            return haplo_files[0]
        return None

    def _exon_dirname(self, exon_label: str) -> str:
        """Directory naming convention: 'Exon 6' -> 'exon6'."""
        return exon_label.lower().replace(" ", "")

    def process_file(self, filename):
        """Process a single sample directory (see File discovery above)."""
        try:
            sample_name, barcode, pattern_type = self.extract_sample_info(filename)
            sample_dir = os.path.join(self.input_dir, filename)

            result_df = pd.DataFrame(columns=self.columns)
            result_df.loc[0, ("", "Barcode")] = barcode.replace("barcode", "")
            result_df.loc[0, ("", "Sequencing_ID")] = sample_name

            any_exon_found = False
            phase_evidence_by_exon: Dict[str, PhaseEvidence] = {}

            combined_dir = os.path.join(sample_dir, self._COMBINED_DIRNAME)
            combined_phenotype_files = (
                glob.glob(os.path.join(combined_dir, "*.ABOPhenotype.txt"))
                if os.path.isdir(combined_dir) else []
            )
            combined_mode = bool(combined_phenotype_files) and os.path.getsize(combined_phenotype_files[0]) > 0
            combined_haplo_file = (
                self._find_haplotype_file(sample_dir, self._COMBINED_DIRNAME) if combined_mode else None
            )

            for exon_label in self.exon_position_lists:
                if combined_mode:
                    phenotype_file = combined_phenotype_files[0]
                    haplo_file = combined_haplo_file
                else:
                    exon_dirname = self._exon_dirname(exon_label)
                    exon_dir = os.path.join(sample_dir, exon_dirname)
                    if not os.path.exists(exon_dir):
                        continue

                    phenotype_files = glob.glob(os.path.join(exon_dir, "*.ABOPhenotype.txt"))
                    if not phenotype_files:
                        continue
                    if os.path.getsize(phenotype_files[0]) == 0:
                        print(f"Empty {exon_label} phenotype file (0 kb) for {filename}. Skipping exon.")
                        continue
                    phenotype_file = phenotype_files[0]
                    haplo_file = self._find_haplotype_file(sample_dir, exon_dirname)

                any_exon_found = True
                exon_data = self.parse_exon_report(phenotype_file, exon_label)
                if not exon_data.empty:
                    for pos in self.exon_position_lists[exon_label]:
                        pos_df = exon_data[exon_data["Position"] == pos]
                        if not pos_df.empty:
                            col_label = self._column_label_for(exon_label, pos)
                            for col in ["#Reads", "Mat", "Mis", "Ins", "Del", "A", "G", "C", "T", "Type"]:
                                if col in pos_df.columns:
                                    result_df.loc[0, (col_label, col)] = pos_df.iloc[0][col]

                phase_evidence = parse_haplotypes_tsv(haplo_file)
                if phase_evidence:
                    phase_evidence_by_exon[exon_label] = phase_evidence
                    print(f"  Phase data ({exon_label}): {phase_evidence.total_reads} reads, "
                          f"{len(phase_evidence.haplotype_counts)} distinct haplotypes")

            if not any_exon_found:
                error_msg = "No exon phenotype files found (checked combined/ and all per-exon folders)"
                print(f"Skipping file {filename}. {error_msg}.")
                self.failed_samples.append({"sample": filename, "reason": error_msg})
                return

            result_df = self.assign_phenotype_genotype(
                result_df, phase_evidence_by_exon=phase_evidence_by_exon,
            )

            self.results.append(result_df)
            print(f"Successfully processed {filename}")

        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            import traceback
            traceback.print_exc()

    def process_files(self):
        """Process all files in the input directory that match expected patterns."""
        print(f"Scanning directory: {self.input_dir}")

        valid_directories = []
        for filename in os.listdir(self.input_dir):
            if os.path.isdir(os.path.join(self.input_dir, filename)):
                match_with_barcode = self.pattern_with_barcode.match(filename)
                match_without_barcode = self.pattern_without_barcode.match(filename)
                if match_with_barcode or match_without_barcode:
                    valid_directories.append(filename)
                else:
                    self.processing_stats["pattern_failed"] += 1
                    self.failed_samples.append(
                        {"sample": filename, "reason": "Filename pattern not recognized"}
                    )

        if not valid_directories:
            print("No valid sample directories found matching expected patterns.")
            return

        print("\033[92m\n ********* Started combining samples to single file ********* \033[0m\n")

        for filename in valid_directories:
            try:
                self.processing_stats["pattern_matched"] += 1
                print(f"\nProcessing file: {filename}")

                sample_name, barcode, pattern_type = self.extract_sample_info(filename)

                if pattern_type == "explicit":
                    self.processing_stats["with_barcode"] += 1
                    print(f"Extracted Sample: {sample_name}, Barcode: {barcode}")
                else:
                    self.processing_stats["without_barcode"] += 1
                    print(f"Extracted Sample: {sample_name}, Barcode: {barcode} (default)")

                if any(char in sample_name for char in ["<", ">", ":", '"', "|", "?", "*"]):
                    print(f"Warning: Sample name '{sample_name}' contains potentially problematic characters")

                self.process_file(filename)
                print(f"Done adding Sample {sample_name} with barcode {barcode} to merged data frame")
            except Exception as e:
                print(f"\nError processing file {filename}: {e}")
                self.failed_samples.append({"sample": filename, "reason": f"Processing error: {str(e)}"})
            finally:
                print(f"Finished processing file: {filename}")

        print(f"\n--- Processing Statistics ---")
        print(f"Files with explicit barcode: {self.processing_stats['with_barcode']}")
        print(f"Files using default barcode: {self.processing_stats['without_barcode']}")
        print(f"Total pattern matches: {self.processing_stats['pattern_matched']}")
        print(f"Pattern match failures: {self.processing_stats['pattern_failed']}")
        print(f"----------------------------")

    def merge_dataframes(self):
        if not self.results:
            print("Warning: No sample results to merge.")
            return pd.DataFrame(columns=self.columns)
        final_df = pd.concat(self.results)
        final_df[("", "Barcode")] = final_df[("", "Barcode")].astype(int)
        final_df = final_df.sort_values(by=[("", "Sequencing_ID"), ("", "Barcode")], ascending=True)
        final_df = final_df.reset_index(drop=True)   # one-row frames all carried index 0
        self._apply_replicate_checks(final_df)
        return final_df

    def save_results_to_file(self, final_df):
        """Save results to text and Excel files; column count/labels are dynamic."""
        try:
            final_df.to_csv("./ABO_result.txt", sep="\t", index=False)
            print("Results saved successfully to text file.")
        except Exception as txt_err:
            print(f"Error saving to text file: {txt_err}")
            return

        try:
            read_count_cols = []
            if isinstance(final_df.columns, pd.MultiIndex):
                for i, col in enumerate(final_df.columns):
                    if col[1] == "#Reads":
                        read_count_cols.append(i)
            else:
                for i, col in enumerate(final_df.columns):
                    if "#Reads" in str(col):
                        read_count_cols.append(i)

            writer = pd.ExcelWriter("./ABO_result.xlsx", engine="xlsxwriter")

            if isinstance(final_df.columns, pd.MultiIndex):
                final_df.columns = final_df.columns.droplevel()

            final_df.to_excel(writer, sheet_name="ABO_Result", header=True, index=False, startrow=1)

            workbook = writer.book
            worksheet = writer.sheets["ABO_Result"]

            data_format = workbook.add_format({"bg_color": "white", "font_color": "black", "border": 1})
            header_format = workbook.add_format({"bold": True, "fg_color": "#007399", "border": 1, "font_color": "white"})
            red_bg_format = workbook.add_format({"bg_color": "#e2725b", "font_color": "black"})
            orange_bg_format = workbook.add_format({"bg_color": "#ff9a00", "font_color": "black"})
            header_format.set_align("center")
            header_format.set_align("vcenter")

            num_rows, num_cols = final_df.shape
            # locate by name: the Notes block now has 8 columns, so a fixed offset from the end is unsafe
            _cols = list(final_df.columns)
            reliability_col = xl_col_to_name(_cols.index("ReadReliability")) if "ReadReliability" in _cols else xl_col_to_name(num_cols - 1)

            print(f"Data has {num_rows} rows, starting at row 3 with two header rows")

            for col_idx in read_count_cols:
                col_letter = xl_col_to_name(col_idx)
                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",
                    {"type": "cell", "criteria": "<=", "value": 20, "format": red_bg_format},
                )
                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",
                    {"type": "cell", "criteria": "between", "minimum": 21, "maximum": 40, "format": orange_bg_format},
                )

            print(f"Applying read count conditional formatting to columns: "
                  f"{[xl_col_to_name(i) for i in read_count_cols]}")

            try:
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",
                    {"type": "formula", "criteria": f'=${reliability_col}3="Very Low(<=20 reads)"', "format": red_bg_format},
                )
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",
                    {"type": "formula", "criteria": f'=${reliability_col}3="Low (<=40 reads)"', "format": orange_bg_format},
                )
            except Exception as format_err:
                print(f"Warning: Could not apply row-level conditional formatting: {format_err}")

            if "QCStatus" in _cols:
                qc_letter = xl_col_to_name(_cols.index("QCStatus"))
                worksheet.conditional_format(
                    f"{qc_letter}3:{qc_letter}{num_rows + 2}",
                    {"type": "cell", "criteria": "==", "value": '"REVIEW"', "format": orange_bg_format},
                )

            for row in range(num_rows):
                for col in range(num_cols):
                    cell_value = final_df.iat[row, col]
                    if not pd.isna(cell_value):
                        worksheet.write(row + 2, col, cell_value, data_format)

            # Per-position header merge ranges
            column_start = 5
            merge_ranges = []
            for col_label in self.position_columns:
                start_col = column_start
                end_col = start_col + 9
                start_letter = xl_col_to_name(start_col)
                end_letter = xl_col_to_name(end_col)
                merge_ranges.append((f"{start_letter}1:{end_letter}1", col_label))
                column_start = end_col + 1

            notes_start = xl_col_to_name(column_start)
            notes_end = xl_col_to_name(column_start + len(self.notes_columns) - 1)

            worksheet.merge_range("A1:B1", "Sample", header_format)
            worksheet.merge_range("C1:E1", "Result", header_format)
            worksheet.merge_range(f"{notes_start}1:{notes_end}1", "Notes", header_format)

            for merge_range in merge_ranges:
                worksheet.merge_range(merge_range[0], merge_range[1], header_format)

            for col in range(num_cols):
                cell_value = final_df.columns[col]
                if not pd.isna(cell_value):
                    worksheet.write(1, col, cell_value, header_format)

            writer.close()
            print("Results saved successfully to Excel file.")
        except Exception as excel_err:
            print(f"Error saving to Excel file: {excel_err}")
            import traceback
            traceback.print_exc()

        # LIS export
        self.df_for_lis_soft = pd.DataFrame()
        self.df_for_lis_soft["Sample ID"] = final_df["Sequencing_ID"]
        self.df_for_lis_soft["Shipment Date"] = ""

        if ("Genotype" in final_df.columns and not final_df["Genotype"].isnull().all()
                and not (final_df["Genotype"] == "Unknown").all()):
            valid_genotype_mask = (final_df["Genotype"] != "Unknown") & final_df["Genotype"].notnull()
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type1"] = final_df.loc[valid_genotype_mask, "Genotype"].str[0]
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type2"] = final_df.loc[valid_genotype_mask, "Genotype"].str[1]
        else:
            self.df_for_lis_soft["ABO Geno Type1"] = ""
            self.df_for_lis_soft["ABO Geno Type2"] = ""

        plain_phenotype = final_df["Phenotype"].copy()

        if "BwSubtype" in final_df.columns:
            bw_mask = final_df["BwSubtype"].astype(str).str.len() > 0
            annotated_phenotype = final_df["Phenotype"].copy()
            annotated_phenotype.loc[bw_mask] = (
                final_df.loc[bw_mask, "Phenotype"] + " [" + final_df.loc[bw_mask, "BwSubtype"] + "]"
            )
            final_df["Phenotype"] = annotated_phenotype

        self.df_for_lis_soft["ABO Pheno Type"] = plain_phenotype
        self.df_for_lis_soft["RH"] = ""
        self.df_for_lis_soft["Blood Type"] = plain_phenotype
        self.df_for_lis_soft["ABORH Comments"] = ""

        if "ASubtype" in final_df.columns:
            self.df_for_lis_soft["A Subtype"] = final_df["ASubtype"]
        if "BwSubtype" in final_df.columns:
            self.df_for_lis_soft["Bw Subtype"] = final_df["BwSubtype"]
        if "PhaseConfidence" in final_df.columns:
            self.df_for_lis_soft["Phase Confidence"] = final_df["PhaseConfidence"]
        if "ReadReliability" in final_df.columns:
            self.df_for_lis_soft["Read Reliability"] = final_df["ReadReliability"]

        if isinstance(final_df.columns, pd.MultiIndex):
            reads_df = final_df.loc[:, (slice(None), "#Reads")]
            self.df_for_lis_soft["#Reads"] = reads_df.mean(axis=1)
        else:
            reads_columns = [col for col in final_df.columns if "#Reads" in str(col)]
            if reads_columns:
                self.df_for_lis_soft["#Reads"] = final_df[reads_columns].mean(axis=1)
            else:
                self.df_for_lis_soft["#Reads"] = 0

        for _src, _dst in (("QCStatus", "QC Status"), ("QCFlags", "QC Flags"),
                           ("A2Trace", "A2 Trace"), ("EvidenceBasis", "Evidence Basis")):
            if _src in final_df.columns:
                self.df_for_lis_soft[_dst] = final_df[_src]

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")
        print(f"LIS export file created successfully with {len(self.df_for_lis_soft)} samples")

    # -----------------------------------------------------------------
    # QC: replicate cross-check and run-level summary (v2.1.0)
    # -----------------------------------------------------------------
    def _apply_replicate_checks(self, final_df):
        """Samples sharing a Sequencing_ID (replicates across barcodes/runs)
        must agree on genotype and on A2 evidence. Disagreement marks EVERY
        member of the group for review; nothing is changed."""
        sid = final_df[("", "Sequencing_ID")].astype(str)
        stats = {"groups": 0, "discordant_genotype": 0, "discordant_a2": 0}
        for _, idx in sid.groupby(sid).groups.items():
            idx = list(idx)
            if len(idx) < 2:
                continue
            stats["groups"] += 1
            grp = final_df.loc[idx]
            genos = set(grp[("Result", "Genotype")].astype(str))
            a2 = grp[("Notes", "A2Trace")].astype(str)
            a2_has = set(not x.startswith("4") for x in a2 if "n/a" not in x and not x.startswith("5"))
            members = ", ".join(
                f"barcode{int(b):02d}={e}" for b, e in zip(grp[("", "Barcode")], grp[("Result", "ExtendedGenotype")]))
            msg = None
            if len(genos) > 1:
                stats["discordant_genotype"] += 1
                msg = f"replicates disagree on genotype ({members})"
            elif len(a2_has) > 1:
                stats["discordant_a2"] += 1
                msg = f"replicates disagree on A2 evidence ({members})"
            if msg:
                for i in idx:
                    old = str(final_df.at[i, ("Notes", "QCFlags")] or "")
                    final_df.at[i, ("Notes", "QCFlags")] = (old + "; " if old else "") + msg
                    final_df.at[i, ("Notes", "QCStatus")] = "REVIEW"
        self.replicate_stats = stats

    def report_qc(self, final_df):
        """Print (and optionally write) a run-level QC summary. Called after
        the results are saved, so final_df has single-level column names."""
        lines = ["=== QC SUMMARY ==="]
        n = len(final_df)
        if n == 0:
            return
        review = int((final_df["QCStatus"] == "REVIEW").sum()) if "QCStatus" in final_df.columns else 0
        lines.append(f"Samples: {n}   flagged REVIEW: {review} ({100.0 * review / n:.1f}%)")

        if "QCFlags" in final_df.columns:
            tally = Counter()
            for flags in final_df["QCFlags"].astype(str):
                for f in flags.split("; "):
                    f = f.strip()
                    if not f or f == "nan":
                        continue
                    f = re.sub(r"\(barcode.*?\)$", "", f).strip()     # drop per-group member lists
                    f = re.sub(r"\(\d+ < \d+,", "(<N,", f)           # collapse read counts
                    f = re.sub(r"\(O1 x\d, O2 x\d, B x\d\)", "", f).strip()
                    tally[f] += 1
            if tally:
                lines.append("Flags (a sample can carry more than one):")
                for f, c in tally.most_common():
                    lines.append(f"  {c:6d}  {f}")

        if "EvidenceBasis" in final_df.columns:
            lines.append("Evidence basis: " + ", ".join(
                f"{k}={v}" for k, v in final_df["EvidenceBasis"].value_counts().items()))
        if "A2Trace" in final_df.columns:
            lines.append("A2 trace:")
            for k, v in final_df["A2Trace"].value_counts().sort_index().items():
                lines.append(f"  {v:6d}  {k}")
        rs = self.replicate_stats
        if rs["groups"]:
            bad = rs["discordant_genotype"] + rs["discordant_a2"]
            lines.append(f"Replicate groups: {rs['groups']}   disagreeing: {bad} "
                         f"(genotype {rs['discordant_genotype']}, A2 evidence {rs['discordant_a2']})")

        notes = []
        recs = self.qc_records
        if recs:
            uni = sum(1 for r in recs if r["reads_uniform"]) / len(recs)
            if uni >= 0.9:
                notes.append(
                    f"#Reads is identical at every primary position for {100 * uni:.0f}% of samples, i.e. it is a "
                    "sample-wide read count, not depth at the site. ReadReliability and the low-read flag therefore "
                    "cannot confirm site-level evidence. Emit per-position depth upstream (AlignmentStatistics 'Depth') "
                    "and write it into the report's 'Aligned Read Count'.")
            ref_only = sum(1 for r in recs if r["evidence"] == "reference-only") / len(recs)
            if ref_only >= 0.5:
                notes.append(
                    f"{100 * ref_only:.0f}% of samples are reference-only (no variant at any primary position). "
                    "These cannot be told apart from failed calls; confirm a subset with independent typing.")
            c = Counter(r["s261"] for r in recs if r["s261"] is not None)
            tot = sum(c.values())
            if tot >= 200:
                q = (c["het"] + 2 * c["hom"]) / (2 * tot)
                exp_het = 2 * q * (1 - q) * tot
                if exp_het > 0:
                    F = 1 - c["het"] / exp_het
                    if F > 0.30:
                        notes.append(
                            f"c.261delG shows a heterozygote deficit (observed {c['het']}, expected {exp_het:.0f} under "
                            f"Hardy-Weinberg; F={F:.2f}). Possible under selection or relatedness; also what missed "
                            "heterozygous calls look like.")
            aa = sum(1 for r in recs if r["genotype"] == "AA")
            ao = sum(1 for r in recs if r["genotype"] == "AO")
            if aa + ao >= 100 and aa / (aa + ao) > 0.5:
                notes.append(
                    f"A/A is {100 * aa / (aa + ao):.0f}% of group A (A/O {100 * ao / (aa + ao):.0f}%). In most populations "
                    "A/O carriers outnumber A/A; check against independently typed samples.")
        if notes:
            lines.append("Notes (informational):")
            for t in notes:
                lines.append("  - " + t)

        text = "\n".join(lines)
        print("\n" + text)
        if self.qc_report_path:
            try:
                with open(self.qc_report_path, "w", encoding="utf-8") as fh:
                    fh.write(text + "\n")
                print(f"QC report written to {self.qc_report_path}")
            except OSError as e:
                print(f"Could not write QC report: {e}")

    def run(self):
        """Run the ABOReportParser."""
        self.process_files()
        final_df = self.merge_dataframes()
        if final_df.empty:
            print("No results to report. Exiting.")
            return
        print("\n\nFinal Results:")
        print("-" * 336)
        print(final_df.to_string(index=False))
        print("-" * 336)
        self.save_results_to_file(final_df)

        if self.failed_samples:
            print("\n\nFailed Samples Summary:")
            print("-" * 80)
            for sample in self.failed_samples:
                print(f"Sample: {sample['sample']} - Reason: {sample['reason']}")
            print("-" * 80)
            print(f"Total failed samples: {len(self.failed_samples)}")
        else:
            print("\nAll samples processed successfully.")

        print(f"\n=== PROCESSING SUMMARY ===")
        print(f"Total directories scanned: "
              f"{self.processing_stats['pattern_matched'] + self.processing_stats['pattern_failed']}")
        print(f"Successfully processed: {len(self.results)}")
        print(f"Samples with explicit barcode: {self.processing_stats['with_barcode']}")
        print(f"Samples with default barcode ({self.default_barcode}): {self.processing_stats['without_barcode']}")
        print(f"Failed samples: {len(self.failed_samples)}")
        print(f"Pattern recognition failures: {self.processing_stats['pattern_failed']}")

        if len(self.results) > 0:
            try:
                phenotype_counts = final_df["Phenotype"].value_counts()
                print(f"\nPhenotype Distribution:")
                for phenotype, count in phenotype_counts.items():
                    print(f"  {phenotype}: {count}")

                try:
                    extended_genotype_counts = final_df["ExtendedGenotype"].value_counts()
                    print(f"\nExtended Genotype Distribution (with A subtypes):")
                    a_subtypes, ab_subtypes, other_genotypes = {}, {}, {}
                    for genotype, count in extended_genotype_counts.items():
                        if "A1" in str(genotype) or "A2" in str(genotype):
                            if "B" in str(genotype):
                                ab_subtypes[genotype] = count
                            else:
                                a_subtypes[genotype] = count
                        else:
                            other_genotypes[genotype] = count
                    if a_subtypes:
                        print(f"  A Subtypes:")
                        for genotype, count in sorted(a_subtypes.items()):
                            print(f"    {genotype}: {count}")
                    if ab_subtypes:
                        print(f"  AB Subtypes:")
                        for genotype, count in sorted(ab_subtypes.items()):
                            print(f"    {genotype}: {count}")
                    if other_genotypes:
                        print(f"  Other Genotypes:")
                        for genotype, count in sorted(other_genotypes.items()):
                            print(f"    {genotype}: {count}")
                except Exception as e:
                    print(f"Could not analyze extended genotype distribution: {e}")

                try:
                    bw_col = final_df.get("BwSubtype")
                    if bw_col is not None:
                        bw_detected = bw_col[bw_col.astype(str).str.len() > 0]
                        if not bw_detected.empty:
                            print(f"\nSubtype Flags ({len(bw_detected)} samples):")
                            for val, count in bw_detected.value_counts().items():
                                print(f"    {val}: {count}")
                except Exception as e:
                    print(f"Could not analyze subtype flag distribution: {e}")

                try:
                    pc_col = final_df.get("PhaseConfidence")
                    if pc_col is not None:
                        print(f"\nPhase Confidence Distribution:")
                        for val, count in pc_col.value_counts().head(10).items():
                            print(f"    {val}: {count}")
                except Exception as e:
                    print(f"Could not analyze phase confidence: {e}")

            except Exception as e:
                print(f"Could not analyze phenotype distribution: {e}")

        print(f"============================")
        self.report_qc(final_df)


# ===========================================================================
# CLI entry point
# ===========================================================================

def create_argument_parser():
    parser = argparse.ArgumentParser(
        prog="aggregate_abo_reports.py",
        description="ABO Blood Group Report Aggregator — Combines individual sample ABO phenotype results into unified reports.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  python aggregate_abo_reports.py /path/to/results
  python aggregate_abo_reports.py /path/to/results --panel abo_variant_panel.yaml
  python aggregate_abo_reports.py /path/to/results --default-barcode barcode99 --verbose

EXPECTED DIRECTORY STRUCTURE (one subfolder per exon PRESENT IN THE PANEL,
not just exon6/exon7 -- e.g. once calibrated, exon2/exon3/exon4/exon5
subfolders are picked up automatically):
  input_directory/
  |-- SAMPLE1_barcode01/
  |   |-- exon6/
  |   |   |-- *.ABOPhenotype.txt
  |   |   `-- *.Haplotypes.tsv
  |   `-- exon7/
  |       |-- *.ABOPhenotype.txt
  |       `-- *.Haplotypes.tsv
  |-- SAMPLE2_barcode02/
  `-- SAMPLE3/  (will use default barcode)

OUTPUT FILES:
  - ABO_result.txt / ABO_result.xlsx / final_export.csv

Column layout and position parsing are driven by --panel; add a new
diagnostic position by editing the panel file, not this script.

For more information, see: https://github.com/fmobegi/nf-core-abotyper
        """,
    )

    parser.add_argument("input_directory",
                         help="Path to directory containing sample subdirectories with ABO phenotype results")
    parser.add_argument("--panel", default="abo_variant_panel.yaml",
                         help="Path to the variant panel file (.yaml/.json/.csv/.tsv). Default: %(default)s")
    parser.add_argument("--default-barcode", "-b", default="barcode00", metavar="BARCODE",
                         help="Default barcode for samples without explicit barcode suffix (default: %(default)s)")
    parser.add_argument("--min-depth", type=int, default=100, metavar="N",
                         help="Flag A1 calls made on fewer than N reads at c.467/c.1061delC (default: %(default)s)")
    parser.add_argument("--qc-report", default="ABO_run_qc.txt", metavar="FILE",
                         help="Write the run-level QC summary here; use '' to disable (default: %(default)s)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output for debugging")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    return parser


def validate_arguments(args):
    if not os.path.exists(args.input_directory):
        print(f"Error: Input directory '{args.input_directory}' does not exist.")
        sys.exit(1)
    if not os.path.isdir(args.input_directory):
        print(f"Error: '{args.input_directory}' is not a directory.")
        sys.exit(1)
    if not re.match(r"^barcode\d{1,2}$", args.default_barcode):
        print(f"Warning: Unusual barcode format '{args.default_barcode}'. Expected format: barcodeXX")
    return True


if __name__ == "__main__":
    parser = create_argument_parser()
    args = parser.parse_args()

    validate_arguments(args)

    if args.verbose:
        print(f"Verbose mode enabled")
        print(f"Input directory: {args.input_directory}")
        print(f"Default barcode: {args.default_barcode}")
        print(f"Panel: {args.panel}")
        print(f"Script version: {__version__}")

    try:
        panel = load_panel(args.panel)
    except Exception as exc:
        print(f"! CRITICAL: could not load variant panel '{args.panel}': {exc}")
        sys.exit(1)

    print(f"Using default barcode: {args.default_barcode}")
    print(f"Loaded panel v{panel.meta.get('panel_version', '?')} with {len(panel.variants)} variants "
          f"({len(panel.calibrated())} usable now, {len(panel.uncalibrated())} awaiting calibration)")

    parser_instance = ABOReportParser(args.input_directory, panel, args.default_barcode,
                                      min_depth=args.min_depth, qc_report=args.qc_report)
    parser_instance.run()
    print("All done!\n")
