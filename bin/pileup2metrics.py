#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pileup2metrics.py -- two lines of evidence for the ABO panel: Clair3 and the BAM pileup.

Reads a `samtools mpileup` text file taken over the panel positions, turns it into the
same per-position table clair2metrics.py writes from a Clair3 gVCF
(*.AlignmentStatistics.tsv), and compares the two at every panel position.

Outputs (only --merged-out is required):

  --merged-out        the table the pipeline scores downstream. Clair3 is the primary line:
                      every position it had evidence for keeps its Clair3 row. A position
                      Clair3 had no evidence for (no record, no reference block -- Depth 0),
                      or ALL positions when Clair3 produced nothing, is filled from the BAM
                      pileup, provided the pileup depth is >= --min-depth.
  --pileup-out        the BAM-derived table on its own
  --concordance-out   per panel position: Clair3 vs BAM mismatch %, their difference and a
                      flag (CONCORDANT / DISCORDANT / CLAIR3_ONLY / BAM_ONLY / NO_EVIDENCE)
  --status-out        one-line sample summary: evidence_source, positions covered by both /
                      Clair3 only / BAM only / neither, discordant positions (all and
                      primary-marker only, with names), primary positions filled from BAM and
                      the maximum pileup depth over the panel (0 = no reads over the ABO region)

evidence_source: CLAIR3+BAM (both lines cover at least one position) | CLAIR3_ONLY |
BAM_ONLY | NONE.

Pileup rows are noisier than Clair3's model-filtered calls, particularly for ONT
homopolymer indels, so disagreement is reported and, at the primary typing markers,
sent to review by aggregate_abo_reports.py rather than silently resolved. Indel-like
markers use a +-2bp window (alignment of a deletion inside a homopolymer is not unique);
SNP markers are exact.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from abo_panel import load_panel, VariantMarker

# Authorship and licence: provided as part of the nf-core/abotyper pipeline.
__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2026, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = [
    "Fredrick Mobegi",
    "Benedict Matern",
    "Mathijs Groeneweg",
    "Claude Sonnet 5",
]
__license__ = "GPL"
__version__ = "2.2.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"

TSV_HEADER = "\t".join([
    "Ref_Position_1based", "Ref_Base", "Match_Percent", "Mismatch_Percent",
    "Insertion_Percent", "Deletion_Percent", "A_Percent", "G_Percent",
    "C_Percent", "T_Percent", "Depth",
])
COLUMNS = TSV_HEADER.split("\t")
INDEL_LIKE_TYPES = {"deletion", "insertion", "indel", "snp_or_indel"}
INDEL_WINDOW = 2
INDEL_RE = re.compile(r"[+-](\d+)")


class PileupColumn:
    __slots__ = ("depth", "ref", "bases", "ins", "dele", "ins_seqs")

    def __init__(self, depth: int, ref: str):
        self.depth = depth
        self.ref = ref
        self.bases = {"A": 0, "C": 0, "G": 0, "T": 0}   # reads showing each base (ref reads included)
        self.ins = 0     # reads with an insertion starting right after this base
        self.dele = 0    # reads with a deletion placeholder ('*') at this base
        self.ins_seqs: List[str] = []   # the inserted sequence of every insertion starting right after this base


def parse_pileup_line(line: str) -> Optional[Tuple[str, int, PileupColumn]]:
    f = line.rstrip("\n").split("\t")
    if len(f) < 5:
        return None
    try:
        pos, depth = int(f[1]), int(f[3])
    except ValueError:
        return None
    ref = f[2].upper()
    col = PileupColumn(depth, ref)
    if depth == 0:
        return f[0], pos, col
    s, i = f[4], 0
    while i < len(s):
        c = s[i]
        if c == "^":
            i += 2                      # read start + mapping-quality character
            continue
        if c in "$<>":
            i += 1
            continue
        if c in "*#":
            col.dele += 1
            i += 1
            continue
        if c in ".,":
            if ref in col.bases:
                col.bases[ref] += 1
            i += 1
            continue
        if c.upper() in col.bases:
            col.bases[c.upper()] += 1
            i += 1
            continue
        if c in "+-":
            m = INDEL_RE.match(s, i)
            if not m:
                i += 1
                continue
            n = int(m.group(1))
            if c == "+":
                col.ins += 1
                col.ins_seqs.append(s[m.end():m.end() + n].upper())
            i = m.end() + n             # skip the inserted/deleted sequence
            continue
        i += 1
    return f[0], pos, col


def load_pileup(path: Path) -> Dict[int, PileupColumn]:
    cols: Dict[int, PileupColumn] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parsed = parse_pileup_line(line)
            if parsed:
                cols[parsed[1]] = parsed[2]
    return cols


def needs_window(marker: VariantMarker) -> bool:
    return marker.variant_type in INDEL_LIKE_TYPES or marker.ref_base == "REF"


def fmt(v: float) -> str:
    return str(int(v)) if float(v).is_integer() else f"{v:.2f}"


def homopolymer_span(pos: int, cols: Dict[int, PileupColumn]) -> Tuple[int, int, str]:
    """First and last reference position of the single-base run that contains pos, and the base."""
    here = cols.get(pos)
    base = here.ref if here else ""
    if base not in ("A", "C", "G", "T"):
        return pos, pos, ""
    first = last = pos
    while first - 1 in cols and cols[first - 1].ref == base:
        first -= 1
    while last + 1 in cols and cols[last + 1].ref == base:
        last += 1
    return first, last, base


def run_indel_percent(pos: int, cols: Dict[int, PileupColumn]) -> Tuple[float, float]:
    """Insertion and deletion % for a marker that sits in a homopolymer run longer than the +-INDEL_WINDOW.

    The aligner puts an indel anywhere in the run (usually at its left end: a c.804dupG insertion is anchored on the
    base before the 7-G run and a c.804delG deletion on its first G, six bases from c.804), so the whole run is searched:
    insertions of the run base anchored from one base before the run to its last base, and deletion placeholders on
    the run. Insertions are summed (one anchor per insertion); for deletions the busiest column is used."""
    first, last, base = homopolymer_span(pos, cols)
    here = cols.get(pos)
    if not base or here is None or not here.depth:
        return 0.0, 0.0
    ins = sum(1 for p in range(first - 1, last + 1) if p in cols for q in cols[p].ins_seqs if q and set(q) == {base})
    dele = max((cols[p].dele for p in range(first, last + 1) if p in cols), default=0)   # a 3-base deletion marks 3 columns: take the busiest, not the sum
    return min(100.0, 100.0 * ins / here.depth), min(100.0, 100.0 * dele / here.depth)


def pileup_row(marker: VariantMarker, cols: Dict[int, PileupColumn]) -> Tuple[List[str], int]:
    """One table row (COLUMNS order) for a panel marker; returns (row, depth)."""
    pos = marker.resolved_position()
    indel_like = needs_window(marker)
    window = INDEL_WINDOW if indel_like else 0
    lo, hi = marker.search_bounds(INDEL_WINDOW) if indel_like else (pos, pos)
    here = cols.get(pos)
    if here is None or here.depth == 0:
        return [str(pos), "N" if indel_like else (marker.ref_base or "N")] + ["0"] * 9, 0

    ins = dele = 0.0
    for p in range(lo, hi + 1):
        c = cols.get(p)
        if c and c.depth:
            ins = max(ins, 100.0 * c.ins / c.depth)
            dele = max(dele, 100.0 * c.dele / c.depth)
    if indel_like:
        first, last, _ = homopolymer_span(pos, cols)
        if first - 1 < pos - window or last + 1 > pos + window:     # the run is wider than the fixed window
            run_ins, run_del = run_indel_percent(pos, cols)
            ins, dele = max(ins, run_ins), max(dele, run_del)

    depth = here.depth
    if indel_like:
        match = max(0.0, 100.0 - ins - dele)
        pct = {b: 0.0 for b in "AGCT"}
        ref_base = "N"
        letter = (marker.ref_base or "").strip().upper()
        if letter in pct:               # e.g. c.1061delC (ref C): the reference fraction is scored on the C column
            pct[letter] = match
            ref_base = letter
        # the mismatch is the kind of indel the marker scores: a duplication row is not "mismatched" by the deletion
        # noise that homopolymer runs give (a 5-G run shows ~25% deletion placeholders in ONT reads)
        kinds = {t.lower() for t in marker.alt_bases()}
        scores_ins = bool(kinds & {"ins", "dup"})
        scores_del = "del" in kinds
        mismatch = (ins if scores_ins or not (scores_ins or scores_del) else 0.0) + (dele if scores_del or not (scores_ins or scores_del) else 0.0)
    else:
        ref = marker.ref_base.strip().upper() or here.ref
        pct = {b: 100.0 * here.bases[b] / depth for b in "AGCT"}
        ref_base = ref
        match = pct.get(ref, 0.0)
        mismatch = max(0.0, 100.0 - match)

    row = [str(pos), ref_base, fmt(round(match, 2)), fmt(round(mismatch, 2)),
           fmt(round(ins, 2)), fmt(round(dele, 2)),
           fmt(round(pct["A"], 2)), fmt(round(pct["G"], 2)),
           fmt(round(pct["C"], 2)), fmt(round(pct["T"], 2)), str(depth)]
    return row, depth


def read_table(path: Path) -> Dict[int, List[str]]:
    rows: Dict[int, List[str]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= len(COLUMNS):
                rows[int(f[0])] = f
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Combine Clair3 and BAM-pileup evidence for the ABO panel positions.")
    ap.add_argument("-p", "--pileup", required=True, help="samtools mpileup output over the panel positions")
    ap.add_argument("--panel", default="abo_variant_panel.yaml")
    ap.add_argument("--clair-metrics", default=None,
                     help="clair2metrics.py AlignmentStatistics.tsv; omit when Clair3 produced no gVCF")
    ap.add_argument("--merged-out", required=True)
    ap.add_argument("--pileup-out", default=None)
    ap.add_argument("--concordance-out", default=None)
    ap.add_argument("--status-out", default=None)
    ap.add_argument("--sample", default="")
    ap.add_argument("--min-depth", type=int, default=10,
                     help="minimum pileup depth for the BAM line to count as evidence (default: %(default)s)")
    ap.add_argument("--discordance-threshold", type=float, default=20.0,
                     help="Clair3-vs-BAM mismatch %% difference flagged DISCORDANT (default: %(default)s)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel)
    markers = sorted(panel.calibrated(), key=lambda m: m.resolved_position())
    cols = load_pileup(Path(args.pileup))
    clair_path = Path(args.clair_metrics) if args.clair_metrics else None
    clair = read_table(clair_path) if clair_path and clair_path.is_file() and clair_path.stat().st_size else {}

    merged, pile_rows, conc = [], [], []
    n_both = n_clair_only = n_bam_only = n_none = 0
    n_filled_primary = 0
    bam_max_depth = 0
    disc_all: List[str] = []
    disc_primary: List[str] = []
    for m in markers:
        pos = m.resolved_position()
        name = m.cdna_change.split()[0]
        primary = m.category == "primary"
        prow, pdepth = pileup_row(m, cols)
        pile_rows.append(prow)
        bam_max_depth = max(bam_max_depth, pdepth)
        crow = clair.get(pos)
        cdepth = int(crow[10]) if crow else 0
        has_c = bool(crow) and cdepth > 0
        has_b = pdepth >= args.min_depth

        cmis = float(crow[3]) if has_c else None
        bmis = float(prow[3]) if pdepth > 0 else None
        if has_c and has_b:
            n_both += 1
            merged.append(crow[:len(COLUMNS)])
            source = "CLAIR3+BAM"
            delta = bmis - cmis
            flag = "DISCORDANT" if abs(delta) >= args.discordance_threshold else "CONCORDANT"
            if flag == "DISCORDANT":
                disc_all.append(name)
                if primary:
                    disc_primary.append(name)
        elif has_c:
            n_clair_only += 1
            merged.append(crow[:len(COLUMNS)])
            source, delta, flag = "CLAIR3", None, "CLAIR3_ONLY"
        elif has_b:
            n_bam_only += 1
            n_filled_primary += primary
            merged.append(prow)
            source, delta, flag = "BAM_PILEUP", None, "BAM_ONLY"
        else:
            n_none += 1
            merged.append(crow[:len(COLUMNS)] if crow else prow)
            source, delta, flag = "NONE", None, "NO_EVIDENCE"
        conc.append([args.sample, name, str(pos), m.category,
                     "" if cmis is None else fmt(cmis), str(cdepth),
                     "" if bmis is None else fmt(bmis), str(pdepth),
                     "" if delta is None else fmt(round(delta, 2)), flag, source])

    def write_table(path: str, rows: List[List[str]]):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(TSV_HEADER + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")

    write_table(args.merged_out, merged)
    if args.pileup_out:
        write_table(args.pileup_out, pile_rows)
    if args.concordance_out:
        with open(args.concordance_out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(["sample", "variant", "position", "category", "clair3_mismatch_pct",
                                "clair3_depth", "bam_mismatch_pct", "bam_depth", "delta_pct",
                                "flag", "merged_source"]) + "\n")
            for r in conc:
                fh.write("\t".join(r) + "\n")

    if n_both:
        evidence = "CLAIR3+BAM"
    elif n_clair_only:
        evidence = "CLAIR3_ONLY"
    elif n_bam_only:
        evidence = "BAM_ONLY"
    else:
        evidence = "NONE"
    if args.status_out:
        with open(args.status_out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("sample\tevidence_source\tpositions\tboth\tclair3_only\tbam_only\tno_evidence\t"
                     "discordant\tdiscordant_primary\tdiscordant_primary_names\tfilled_primary_from_bam\tbam_max_depth\n")
            fh.write("\t".join([args.sample, evidence, str(len(markers)), str(n_both), str(n_clair_only),
                                str(n_bam_only), str(n_none), str(len(disc_all)), str(len(disc_primary)),
                                ";".join(disc_primary), str(int(n_filled_primary)),
                                str(bam_max_depth)]) + "\n")
    print(f"{args.sample or 'sample'}: source={evidence} both={n_both} clair3_only={n_clair_only} "
          f"bam_only={n_bam_only} none={n_none} discordant={len(disc_all)} "
          f"(primary {len(disc_primary)})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
