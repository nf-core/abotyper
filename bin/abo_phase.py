#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase evidence from the Haplotypes.tsv tables: phase confidence and cis/trans relationships.

Part of aggregate_abo_reports.py (nf-core/abotyper); split out of that script so each part can be read and
tested on its own. This file is part of nf-core/abotyper and is released under the MIT License (see LICENSE).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd




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
    if not haplo_file or not os.path.exists(haplo_file) or os.path.getsize(haplo_file) == 0:
        return None

    # a file that exists but cannot be parsed raises: the sample is reported as failed, not as unphased
    try:
        df = pd.read_csv(haplo_file, sep="\t", encoding="utf-8")
    except pd.errors.EmptyDataError:
        return None
    if df.empty:
        return None

    evidence = PhaseEvidence(total_reads=len(df))

    if "Haplotype" in df.columns:
        for haplo_str, count in df["Haplotype"].value_counts().items():
            evidence.haplotype_counts[str(haplo_str)] = int(count)

    return evidence


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
