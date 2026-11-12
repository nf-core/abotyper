#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genotype calling: primary-marker states, the genotype table, A-subtype and subtype-marker scans, and per-haplotype allele naming.

Part of aggregate_abo_reports.py (nf-core/abotyper); split out of that script so each part can be read and
tested on its own. This file is part of nf-core/abotyper and is released under the MIT License (see LICENSE).
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

import pandas as pd

from abo_panel import VariantMarker, call_named_marker
from abo_naming import _unpack_fired, encode_allele, encode_extended_genotype, primary_carriers
from abo_phase import PhaseEvidence, compute_phase_confidence, phase_relationship



# State ("ref", "het" or "alt") of the five primary markers, in this order:
#   c.261delG (role o1), c.796C>A (b_796), c.802G>A (o2_802), c.803G>C (b_803), c.804 homopolymer (o34_804)
# -> (Phenotype, Genotype, ExtendedGenotype, kind).
# kind "fixed": as written. "a_allele": "{a}" is the A allele named by the A-subtype step ("A" when it names none).
# "aa": both alleles are A, and the A-subtype step's zygosity decides the names (see resolve_primary_genotype).
# Any other combination, or a marker without evidence ("none"), is Unknown.
GENOTYPE_TABLE = {
    ("het", "ref", "ref", "ref", "ref"): ("A", "AO", "{a}/O1", "a_allele"),
    ("ref", "ref", "het", "ref", "ref"): ("A", "AO", "{a}/O2", "a_allele"),
    ("ref", "ref", "ref", "ref", "het"): ("A", "AO", "{a}/O3", "a_allele"),
    ("het", "het", "ref", "het", "ref"): ("B", "BO", "B/O1", "fixed"),
    ("ref", "het", "het", "het", "ref"): ("B", "BO", "B/O2", "fixed"),
    ("ref", "het", "ref", "het", "het"): ("B", "BO", "B/O3", "fixed"),
    ("het", "ref", "het", "ref", "ref"): ("O", "OO", "O1/O2", "fixed"),
    ("het", "ref", "ref", "ref", "het"): ("O", "OO", "O1/O3", "fixed"),
    ("ref", "ref", "het", "ref", "het"): ("O", "OO", "O2/O3", "fixed"),
    ("alt", "ref", "ref", "ref", "ref"): ("O", "OO", "O1/O1", "fixed"),
    ("ref", "ref", "alt", "ref", "ref"): ("O", "OO", "O2/O2", "fixed"),
    ("ref", "ref", "ref", "ref", "alt"): ("O", "OO", "O3/O3", "fixed"),
    ("ref", "ref", "ref", "ref", "ref"): ("A", "AA", None, "aa"),
    ("ref", "alt", "ref", "alt", "ref"): ("B", "BB", "B/B", "fixed"),
    ("ref", "het", "ref", "het", "ref"): ("AB", "AB", "{a}/B", "a_allele"),
}


PRIMARY_GENOTYPE_ROLES = ("o1", "b_796", "o2_802", "b_803", "o34_804")


def resolve_primary_genotype(states, determine_a_subtype):
    """(Phenotype, Genotype, ExtendedGenotype, A-subtype warning) from the five primary-marker states.

    determine_a_subtype() -> (subtype, warning, zygosity) is called only for genotypes with an A allele, because the
    A2 scan it runs records the markers that later name the ISBT alleles.
    """
    rule = GENOTYPE_TABLE.get(tuple(states))
    if rule is None:
        return "Unknown", "Unknown", "Unknown", None
    phenotype, genotype, extended, kind = rule
    if kind == "fixed":
        return phenotype, genotype, extended, None
    a_subtype, warning, zygosity = determine_a_subtype()
    if kind == "a_allele":
        return phenotype, genotype, extended.format(a=a_subtype or "A"), warning
    if a_subtype and zygosity == "het":      # one A1 allele and one subtype allele, not two copies of the subtype
        return "A", "AA", f"A1/{a_subtype}", warning
    if a_subtype:
        return a_subtype, "AA", f"{a_subtype}/{a_subtype}", warning
    return "A", "AA", "A/A", warning


def missing_required_markers(panel) -> List[str]:
    """Typing roles (panel meta.roles, or the defaults) whose marker id the panel does not define."""
    return panel.missing_roles()


class CallingMixin:
    """Calling methods of ABOReportParser (they use the parser's panel, columns and per-sample state)."""

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

    def scan_category_markers(
        self,
        category: str,
        row_lookup: Dict[int, pd.Series],
        require_b_allele: bool = False,
        require_a_allele: bool = False,
        b_marker_states: Tuple[str, ...] = (),
    ) -> Tuple[List[str], List[str]]:
        """
        Works for ANY panel category (bw, ba, cisab, a2, a3, am, aweak, ael,
        bweak, bel, b3, onull, ...). Returns (fired_marker_texts, warnings).
        """
        fired, warnings = [], []

        if require_b_allele:
            has_b = any(st in ("alt", "het") for st in b_marker_states)      # c.796 / c.803 carry a B allele
            if not has_b:
                return fired, warnings
        if require_a_allele:
            has_a = any(st in ("ref", "het") for st in b_marker_states)      # an allele without B (A or O)
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
                self._fired_markers.append((marker, result.get("alt")))
                self._fired_pct[marker.id] = result.get("pct")
            elif result and not result.get("reason") and nreads and nreads >= (marker.min_reads or 0):
                self._absent_markers.append(marker)      # read with enough depth and not carried
            elif result and result.get("reason") == "low_coverage":
                warnings.append(
                    f"Low coverage at {marker.cdna_change} — {marker.category} "
                    f"subtyping unreliable"
                )

        return fired, warnings

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

            col_o1 = col_for(self.panel.role_id("o1"))
            col_796 = col_for(self.panel.role_id("b_796"))
            col_802 = col_for(self.panel.role_id("o2_802"))
            col_803 = col_for(self.panel.role_id("b_803"))
            col_804 = col_for(self.panel.role_id("o34_804"))
            col_467 = col_for(self.panel.role_id("a1a2_467"))
            col_1061 = col_for(self.panel.role_id("a2_1061del"))

            def primary_state(col, marker_id):
                m = self.panel.get(marker_id)
                if m is None or not col:
                    return "none"
                return self._primary_state(m, dict(safe_get_row(df, col)))

            # the five primary markers' states key the genotype table (GENOTYPE_TABLE)
            state_o1 = primary_state(col_o1, self.panel.role_id("o1"))
            state_796 = primary_state(col_796, self.panel.role_id("b_796"))
            state_802 = primary_state(col_802, self.panel.role_id("o2_802"))
            state_803 = primary_state(col_803, self.panel.role_id("b_803"))
            state_804 = primary_state(col_804, self.panel.role_id("o34_804"))
            type_exon7_93 = safe_get_type(df, col_467) if col_467 else ""

            state_exon7_685 = primary_state(col_1061, self.panel.role_id("a2_1061del"))
            state_exon7_93 = primary_state(col_467, self.panel.role_id("a1a2_467"))
            nreads_exon7_p685 = safe_get_reads(df, col_1061) if col_1061 else 0
            background_notes: List[str] = []   # panel changes that are not part of any ISBT allele definition

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
            self._fired_markers = []   # non-primary markers that fired; refine ISBT allele names
            self._absent_markers = []  # markers scanned with enough coverage that did not fire; rule alleles out
            self._fired_pct = {}

            def determine_a2_subtype(markers_found: List[str]):
                has_del = "c.1061del" in markers_found
                has_297 = "c.297G" in markers_found
                # ISBT: ABO*A2.01 = c.467C>T + c.1061delC; ABO*A2.06 = c.1061delC alone. c.907G>A and c.1032G>A
                # are pipeline background markers, not part of any ISBT definition, so they never name an allele.
                if has_del and state_exon7_93 == "ref" and nreads_exon7_p422 >= self.panel.threshold("a2_min_reads"):
                    return "A2.06", None
                if has_297 and has_del:
                    return "A2.01", "c.297A>G"
                if has_del:
                    return "A2.01", None
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
                has_del = state_exon7_685 in ("alt", "het") and nreads_exon7_p422 >= self.panel.threshold("a2_min_reads")
                a907_marker = self.panel.role("background_907")
                a1032_marker = self.panel.role("background_1032")
                has_907 = False
                has_1032 = False
                if a907_marker and a907_marker.resolved_position() in row_lookup:
                    r = call_named_marker(a907_marker, dict(row_lookup[a907_marker.resolved_position()]),
                                           nreads=row_lookup[a907_marker.resolved_position()].get("#Reads", 0))
                    has_907 = bool(r and r.get("fired")) and nreads_exon7_p422 >= self.panel.threshold("a2_min_reads")
                if a1032_marker and a1032_marker.resolved_position() in row_lookup:
                    r = call_named_marker(a1032_marker, dict(row_lookup[a1032_marker.resolved_position()]),
                                           nreads=row_lookup[a1032_marker.resolved_position()].get("#Reads", 0))
                    has_1032 = bool(r and r.get("fired")) and nreads_exon7_p422 >= self.panel.threshold("a2_min_reads")

                a297_marker = self.panel.role("a2_297")
                has_297 = False
                a297_phase_note = None
                if a297_marker and a297_marker.resolved_position() in row_lookup:
                    row = row_lookup[a297_marker.resolved_position()]
                    has_297 = float(row.get("G", 0) or 0) >= self.panel.threshold("c297_min_alt_pct")
                    if has_297:
                        # c.297A>G is dual-purpose (panel notes): also part of
                        # the shared B-lineage backbone. If phasing shows it
                        # riding on the SAME haplotype as the B primary
                        # marker, it's likely just B-backbone noise, not
                        # independently confirmatory for A2.01 -- advisory
                        # only, does not change has_297 or any call below.
                        b_marker = self.panel.role("b_796")
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
                    background_notes.append("c.907G>A seen (a pipeline background marker, not an ISBT-defined change; not used for allele naming)")
                if has_1032:
                    background_notes.append("c.1032G>A seen (a pipeline background marker, not an ISBT-defined change; not used for allele naming)")
                if has_297:
                    markers.append("c.297G")
                if a297_phase_note:
                    warns.append(a297_phase_note)

                o1v_fired = has_297 and not has_del and not has_907 and not has_1032
                if o1v_fired:
                    if a297_marker:
                        self._fired_markers.append((a297_marker, None))   # points at the O1v-type allele
                    for m_ in ("c.297G", "c.266T", "c.268C"):
                        if m_ in markers:
                            markers.remove(m_)
                    warns.append(
                        "Note: c.297G detected without c.1061del — possible A2-derived "
                        "O allele (O1v); manual review recommended"
                    )

                if not o1v_fired:
                    m266 = self.panel.role("a2_266")
                    m268 = self.panel.role("a2_268")
                    for m_marker, label in [(m266, "c.266T"), (m268, "c.268C")]:
                        if m_marker and m_marker.resolved_position() in row_lookup:
                            row = row_lookup[m_marker.resolved_position()]
                            row_values = {"A": row.get("A", 0), "G": row.get("G", 0),
                                          "C": row.get("C", 0), "T": row.get("T", 0)}
                            r = call_named_marker(m_marker, row_values, nreads=row.get("#Reads", 0))
                            if r and r.get("fired"):
                                markers.append(label)
                                self._fired_markers.append((m_marker, r.get("alt")))
                                self._fired_pct[m_marker.id] = r.get("pct")

                # an A2-panel position read with enough depth and not carrying its alt rules out the alleles it defines
                for m_ in self.panel.by_category("a2_panel"):
                    pos = m_.resolved_position()
                    if pos is None or pos not in row_lookup:
                        continue
                    row = row_lookup[pos]
                    r = call_named_marker(m_, {"A": row.get("A", 0), "G": row.get("G", 0),
                                               "C": row.get("C", 0), "T": row.get("T", 0)},
                                          nreads=row.get("#Reads", 0))
                    reads_here = row.get("#Reads", 0)
                    if (r and not r.get("fired") and not r.get("reason")
                            and reads_here and reads_here >= (m_.min_reads or 0)):
                        self._absent_markers.append(m_)

                if nreads_exon7_p422 >= self.panel.threshold("a2_panel_scan_min_reads") and not o1v_fired:
                    for m_ in self.panel.by_category("a2_panel"):
                        if m_.id in (self.panel.role_id("background_907"), self.panel.role_id("background_1032")):
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
                            self._fired_markers.append((m_, r.get("alt")))
                            self._fired_pct[m_.id] = r.get("pct")
                else:
                    # the A2 logic does not use these positions here (an O1v-type sample, or too few reads), but a
                    # marker that is clearly carried still decides the allele name, e.g. the c.771/c.829 cassette
                    for m_ in self.panel.by_category("a2_panel"):
                        if m_.id in (self.panel.role_id("background_907"), self.panel.role_id("background_1032"), self.panel.role_id("a2_297")):
                            continue
                        pos = m_.resolved_position()
                        if pos is None or pos not in row_lookup:
                            continue
                        row = row_lookup[pos]
                        if not row.get("#Reads", 0) or row.get("#Reads", 0) < max(self.panel.threshold("named_marker_min_reads"), m_.min_reads or 0):
                            continue
                        r = call_named_marker(m_, {"A": row.get("A", 0), "G": row.get("G", 0),
                                                   "C": row.get("C", 0), "T": row.get("T", 0)}, nreads=row.get("#Reads", 0))
                        if r and r.get("fired") and all(mk.id != m_.id for mk, _ in (_unpack_fired(x) for x in self._fired_markers)):
                            self._fired_markers.append((m_, r.get("alt")))
                            self._fired_pct[m_.id] = r.get("pct")

                return markers, warns

            a2_scanned = [False]

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
                a2_scanned[0] = True
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
                        # (c.297G) -- check ITS allele
                        # fraction directly via the same ref/alt/het
                        # classifier used for primary markers, since
                        # scan_a2_markers()/call_named_marker() only report
                        # whether the marker fired, not its zygosity.
                        driving_marker_id = {
                            "c.297G": self.panel.role_id("a2_297"),
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
                # ISBT: ABO*A1.02 = c.467C>T alone. Named only when c.1061delC is clearly absent on enough reads,
                # because a missed deletion looks the same and was the usual cause of this pattern in earlier audits.
                if (state_exon7_93 in ("alt", "het") and state_exon7_685 == "ref"
                        and nreads_exon7_p685 >= self.panel.threshold("a2_min_reads") and nreads_exon7_p422 >= self.panel.threshold("a2_min_reads")):
                    return "A1.02", warning_text, ("het" if state_exon7_93 == "het" else "hom")
                m467 = self.panel.role("a1a2_467")
                if m467 is not None and type_exon7_93 in (m467.alt_call_label, m467.ref_call_label):
                    return "A1", warning_text, "hom"
                return "", None, "hom"

            # ----- PART 1: PRIMARY PHENOTYPING LOGIC (GENOTYPE_TABLE) -----
            Phenotype, Genotype, ExtendedGenotype, a_subtype_warning = resolve_primary_genotype(
                (state_o1, state_796, state_802, state_803, state_804), determine_a_subtype)

            # ----- PART 1b: A2-panel markers for genotypes without an A allele (O and B samples) -----
            # They do not change the call, but a marker carried on the O or B allele (c.297, the c.646/681/771/829
            # cassette) decides which ISBT allele that allele is named.
            if not a2_scanned[0]:
                for m_ in list(self.panel.by_category("a2_panel")):
                    if m_.id in (self.panel.role_id("background_907"), self.panel.role_id("background_1032")):
                        continue
                    pos = m_.resolved_position()
                    if pos is None or pos not in row_lookup:
                        continue
                    row = row_lookup[pos]
                    reads_here = row.get("#Reads", 0) or 0
                    if reads_here < max(self.panel.threshold("named_marker_min_reads"), m_.min_reads or 0):
                        continue
                    r = call_named_marker(m_, {"A": row.get("A", 0), "G": row.get("G", 0),
                                               "C": row.get("C", 0), "T": row.get("T", 0)}, nreads=reads_here)
                    if r and r.get("fired"):
                        self._fired_markers.append((m_, r.get("alt")))
                        self._fired_pct[m_.id] = r.get("pct")
                    elif r and not r.get("reason"):
                        self._absent_markers.append(m_)

            # ----- PART 2: subtype-marker scan across ALL panel categories -----
            bw_markers, bw_warnings = self.scan_category_markers(
                "bw", row_lookup, require_b_allele=True,
                b_marker_states=(state_796, state_803),
            )
            ba_markers, ba_warnings = self.scan_category_markers(
                "ba", row_lookup,
                b_marker_states=(state_796, state_803),
            )
            cisab_markers, cisab_warnings = self.scan_category_markers(
                "cisab", row_lookup,
                b_marker_states=(state_796, state_803),
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
                if min_reads <= self.panel.threshold("reliability_very_low_reads"):
                    Reliability = self.reliability_label("very_low")
                elif min_reads <= self.panel.threshold("reliability_low_reads"):
                    Reliability = self.reliability_label("low")
                elif min_reads >= self.panel.threshold("reliability_robust_reads"):
                    Reliability = self.reliability_label("robust")
                else:
                    Reliability = "Normal"
                if max(valid_read_counts) / min(valid_read_counts) > 5:
                    Reliability += " (Variable coverage)"
            else:
                Reliability = "Unknown (no read data)"

            # ----- PART 5: QC / FAILURE CAPTURE (v2.1.0) -----
            # Annotates the call; never changes it. See the module docstring.
            qc_review: List[str] = []
            encoded_genotype, named_by_phase = self._encode_with_phase(ExtendedGenotype, resolved_phase_evidence)
            if named_by_phase:
                PhaseConfidence += "; alleles named per haplotype"
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
                return "hom" if frac >= self.panel.threshold("homozygous_min_pct") else ("het" if frac >= self.panel.threshold("heterozygous_min_pct") else "ref")

            s261 = vstate(alt_frac(self.panel.role_id("o1")))
            s796 = vstate(alt_frac(self.panel.role_id("b_796")))
            s802 = vstate(alt_frac(self.panel.role_id("o2_802")))
            s803 = vstate(alt_frac(self.panel.role_id("b_803")))
            s467 = vstate(alt_frac(self.panel.role_id("a1a2_467")))
            s1061 = vstate(alt_frac(self.panel.role_id("a2_1061del")))
            primary_states = [s for s in (s261, s796, s802, s803, s467, s1061) if s is not None]

            # evidence actually present at each primary position
            no_evidence, primary_reads, unscored = [], [], []
            for pid, short in ((self.panel.role_id("o1"), "c.261delG"), (self.panel.role_id("b_796"), "c.796C>A"),
                               (self.panel.role_id("o2_802"), "c.802G>A"), (self.panel.role_id("b_803"), "c.803G>C"),
                               (self.panel.role_id("a1a2_467"), "c.467C>T"), (self.panel.role_id("a2_1061del"), "c.1061delC")):
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

            qc_review.extend(background_notes)

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
                    if any(t == "A1.02" for t in tokens):
                        A2Trace = "2 A1.02 pattern: c.467C>T without c.1061delC (ABO*A1.02)"
                        if "ABO*AEL.08" not in encoded_genotype:    # an Ael call already explains the c.467C>T
                            qc_review.append(
                                f"A1.02 pattern: c.467C>T seen in {alt_frac('a1_a2_467') or 0:.0f}% of reads, c.1061delC in only "
                                f"{alt_frac('a1_a2_1061del') or 0:.0f}% "
                                f"({int(cell(self.panel.get('a1_a2_1061del').column_label(), '#Reads'))} reads); "
                                "ISBT defines A1.02 as c.467C>T alone. Rare allele: confirm that the deletion is truly absent. "
                                f"c.467C>T plus c.804dupG is ABO*AEL.08 (insertion signal here: "
                                f"{cell(self.panel.get('ael_804_indel').column_label(), 'Ins'):.0f}%); ONT can under-read an insertion in the 7-G run")
                    else:
                        A2Trace = "2 A2 partial: c.467 variant, c.1061delC not called"
                        qc_review.append("A2 partial: c.467 variant without c.1061delC (possible missed deletion or A1.02)")
                elif e1061:
                    if any(t == "A2.06" for t in tokens):
                        A2Trace = "3 A2.06 pattern: c.1061delC without c.467C>T (ABO*A2.06)"
                        qc_review.append(
                            f"A2.06 pattern: c.1061delC seen in {alt_frac('a1_a2_1061del') or 0:.0f}% of reads, c.467C>T in only "
                            f"{alt_frac('a1_a2_467') or 0:.0f}% ({int(cell(self.panel.get('a1_a2_467').column_label(), '#Reads'))} reads); "
                            "ISBT defines A2.06 as c.1061delC alone (A2.01 also needs c.467C>T). Rare allele: confirm that c.467C>T is truly absent")
                    else:
                        A2Trace = "3 A2 partial: c.1061delC variant, c.467 not called"
                        qc_review.append("A2 partial: c.1061delC variant without c.467")
                else:
                    A2Trace = "4 No A2 evidence (both primary markers reference)"
                a2_tokens = [t for t in tokens if t.startswith("A2")]
                a1_tokens = [t for t in tokens if t == "A1"]
                if any(t == "A" for t in tokens):
                    qc_review.append("A subtype unresolved (generic A)")
                if a1_tokens and not a2_tokens and not A2Trace.startswith("4") and not any(x == "A1.02" for x in tokens):
                    qc_review.append("reported A1 despite A2-lineage marker")
                if a2_tokens and A2Trace.startswith("4"):
                    qc_review.append("reported A2 without a primary A2 marker")
                if A2Trace.startswith("4"):
                    a2_reads = [cell(self.panel.get(p).column_label(), "#Reads")
                                for p in (self.panel.role_id("a1a2_467"), self.panel.role_id("a2_1061del"))
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

            try:
                _sid = str(df.at[0, ("", "Sequencing_ID")])
            except (KeyError, IndexError):
                _sid = ""

            # --- depth gate: reads at the primary positions. Annotates the sample; the call is never changed.
            #   normal sample     : Very low (<=20 reads) / Low (< --min-depth) -> REVIEW; otherwise OK
            #   negative control  : few reads is the expected, correct result; reads above --min-depth mean
            #                       possible contamination -> REVIEW
            #   positive control  : needs normal depth like a patient sample; thin reads mean a failed control -> REVIEW
            declared = self._current_control
            is_negative = declared == "negative" or bool(
                not declared and self.negative_control_re and self.negative_control_re.search(_sid))
            is_positive = declared == "positive" or bool(
                not declared and not is_negative and self.positive_control_re and self.positive_control_re.search(_sid))
            lowest = min(primary_reads) if primary_reads else 0
            highest = max(primary_reads) if primary_reads else 0
            if is_negative:
                if highest <= self.min_depth:
                    DepthCheck = f"Negative control: low depth as expected (max {int(highest)} reads)"
                    qc_review = []        # every other flag on a near-empty control is a consequence of the low depth
                else:
                    DepthCheck = f"Negative control: reads present (max {int(highest)} reads), check for contamination"
                    qc_review.append(f"negative control has reads (max {int(highest)} > {self.min_depth}): possible contamination")
            elif is_positive and (not primary_reads or highest <= 0 or lowest < self.min_depth):
                DepthCheck = f"Positive control: depth too low ({int(lowest)} reads), control failed"
                qc_review.append(f"positive control has too few reads (min {int(lowest)} < {self.min_depth})")
            elif is_positive:
                DepthCheck = f"Positive control: depth OK ({int(lowest)} reads)"
            elif not primary_reads or highest <= 0:
                DepthCheck = "No data"
                qc_review.append("no reads at the primary positions")
            elif lowest <= self.panel.threshold("reliability_very_low_reads"):
                DepthCheck = f"Very low ({int(lowest)} reads)"
                qc_review.append(f"low depth at primary positions (min {int(lowest)} reads < {self.min_depth})")
            elif lowest < self.min_depth:
                DepthCheck = f"Low ({int(lowest)} reads)"
                qc_review.append(f"low depth at primary positions (min {int(lowest)} reads < {self.min_depth})")
            else:
                DepthCheck = f"OK ({int(lowest)} reads)"

            if "ABO*AEL." in encoded_genotype:
                qc_review.append("Ael allele (c.804 homopolymer indel): very weak or silent A expression, so serology may read O; "
                                 "ONT reads are error-prone in this 7-G run, confirm the indel")

            QCStatus = "REVIEW" if qc_review else "PASS"
            QCFlags = "; ".join(qc_review)
            self.qc_records.append({
                "sample": _sid, "genotype": Genotype, "ext": ExtendedGenotype, "a2trace": A2Trace,
                "evidence": EvidenceBasis, "s261": s261, "reads_uniform": reads_uniform,
            })

            # A primary position with zero reads is missing evidence, not a position to skip.
            _zero = [r for r in primary_reads if r <= 0]
            if _zero:
                Reliability = self.reliability_label("very_low")

            ASubtype = a_subtype_warning or ""

            if any("cisAB" in w for w in all_subtype_warnings):
                Reliability += "; " + "; ".join(w for w in all_subtype_warnings if "cisAB" in w)

            df[("Result", "Phenotype")] = Phenotype
            df[("Result", "Genotype")] = Genotype
            df[("Result", "ExtendedGenotype")] = encoded_genotype
            df[("Notes", "ASubtype")] = ASubtype
            df[("Notes", "ReadReliability")] = Reliability
            df[("Notes", "DepthCheck")] = DepthCheck
            df[("Notes", "PhaseConfidence")] = PhaseConfidence
            df[("Notes", "BwSubtype")] = BwSubtype
            df[("Notes", "A2Trace")] = A2Trace
            df[("Notes", "EvidenceBasis")] = EvidenceBasis
            df[("Notes", "QCStatus")] = QCStatus
            df[("Notes", "QCFlags")] = QCFlags

            return df

        except Exception as e:
            # a bug or unexpected input in the typing logic: the sample is reported as failed (process_file records
            # it and the run exits non-zero) instead of carrying an "Error" call into the tables
            print(f"Error in assign_phenotype_genotype: {str(e)}")
            raise

    @staticmethod
    def _token_fits_haplotype(token: str, pos: set, ids: dict) -> bool:
        """Can this internal allele token (O1, B, A2.01, ...) sit on a haplotype that carries the ALT alleles at pos?"""
        has = lambda key: ids.get(key) in pos                     # noqa: E731
        o_or_b = has("o1") or has("o2") or has("o3") or has("b796") or has("b803")
        if token == "O1":
            return has("o1")
        if token == "O2":
            return has("o2")
        if token == "O3":
            return has("o3")
        if token == "B":
            return has("b796") and has("b803")
        if o_or_b:
            return False                                           # an A allele carries no O or B primary marker
        if token == "A2.01":
            return has("a467") and has("a1061")
        if token == "A2.06":
            return has("a1061") and not has("a467")
        if token == "A1.02":
            return has("a467") and not has("a1061")
        if token == "A1":
            return not has("a467") and not has("a1061")
        if token == "A2":
            return has("a467") or has("a1061")
        return True                                                # bare "A": no constraint

    def _name_alleles_by_haplotype(self, ext: str, evidence, fired, absent):
        """Name the two alleles of a genotype from the two phased haplotypes, or return None when the phasing cannot
        carry the decision (then the caller pools the markers of both alleles, as before).

        Each fired marker is given to the haplotype that carries it; a marker on the other haplotype counts as absent
        for this allele. This is what separates O.01.01 / O.01.02 in an O1/O1 sample and keeps the markers of the O
        allele from renaming the A allele."""
        if not evidence or not evidence.haplotype_counts or not isinstance(ext, str) or ext in ("", "Unknown", "Error"):
            return None
        tokens = ext.split("/")
        if len(tokens) != 2:
            return None
        ranked = sorted(evidence.haplotype_counts.items(), key=lambda x: x[1], reverse=True)
        if len(ranked) < 2 or evidence.total_reads <= 0:
            return None
        top2 = ranked[0][1] + ranked[1][1]
        if top2 / evidence.total_reads < 0.90 or ranked[1][1] / evidence.total_reads < 0.20:
            return None                                            # not two clean haplotypes (hom, or a broken phase set)
        haps = [set(int(p) for p in re.findall(r"p(\d+)=", h)) for h, _ in ranked[:2]]
        ids = {}
        for key, marker_id in (("o1", self.panel.role_id("o1")), ("o2", self.panel.role_id("o2_802")), ("o3", self.panel.role_id("o34_804")),
                               ("b796", self.panel.role_id("b_796")), ("b803", self.panel.role_id("b_803")),
                               ("a467", self.panel.role_id("a1a2_467")), ("a1061", self.panel.role_id("a2_1061del"))):
            m = self.panel.get(marker_id)
            if m is not None and m.resolved_position() is not None:
                ids[key] = m.resolved_position()
        # clair2haplotypes.py drops homozygous positions (the same on both sides), so a position both alleles need
        # (c.261delG in an O1/O1 sample) belongs to both haplotypes even though neither string lists it
        needs = [self._TOKEN_NEEDS.get(tok, set()) for tok in tokens]
        for key in needs[0] & needs[1]:
            if key in ids:
                haps[0].add(ids[key])
                haps[1].add(ids[key])
        fired_pos = [(mk, alt, mk.resolved_position()) for mk, alt in (_unpack_fired(x) for x in fired)]
        for mk, _, p in fired_pos:
            if p is not None and not (p in haps[0] or p in haps[1]):
                if (self._fired_pct.get(mk.id) or 0) >= self.panel.threshold("homozygous_min_pct"):        # homozygous marker: both alleles carry it
                    haps[0].add(p)
                    haps[1].add(p)
        if any(p is None or not (p in haps[0] or p in haps[1]) for _, _, p in fired_pos):
            return None                                            # a marker the phasing does not know about
        fits = [perm for perm in ((0, 1), (1, 0))
                if self._token_fits_haplotype(tokens[0], haps[perm[0]], ids)
                and self._token_fits_haplotype(tokens[1], haps[perm[1]], ids)]
        if not fits or (len(fits) == 2 and tokens[0] != tokens[1]):
            return None
        names = []
        for tok, h in zip(tokens, fits[0]):
            mine = [(mk, alt) for mk, alt, p in fired_pos if p in haps[h]]
            others = [mk for mk, alt, p in fired_pos if p not in haps[h]]
            names.append(encode_allele(tok, mine, list(absent) + others, primary_carriers(self.panel), pooled=False))
        return "/".join(sorted(names) if tokens[0] == tokens[1] else names)

    def _encode_with_phase(self, ext: str, evidence):
        """(ISBT-style name, True when the haplotypes decided it)."""
        try:
            named = self._name_alleles_by_haplotype(ext, evidence, self._fired_markers, self._absent_markers)
        except (KeyError, IndexError, ValueError, TypeError):
            named = None
        if named:
            return named, True
        return encode_extended_genotype(ext, self._fired_markers, self._absent_markers, primary_carriers(self.panel),
                                        fired_pct=self._fired_pct), False
