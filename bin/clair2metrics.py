#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
clair2metrics.py -- convert a Clair3 --gvcf output into the same
AlignmentStatistics.tsv format pysam_haploscan.py / stats_from_pileup.py
produce, so predict_abo_phenotype.py / aggregate_abo_reports.py can
consume Clair3 calls with zero changes to their scoring logic.

Only the panel's calibrated diagnostic positions are emitted (not a
full-reference scan): nothing downstream reads any other position, and
Clair3's gVCF only needs walking near those ~51 windows.

Works directly off Clair3's raw --gvcf output (no bcftools norm step
required) -- validated this way against a real 83-sample cohort:
Clair3's own indel-candidate anchor sits at the panel's declared
position minus 1 for every deletion tested, which the padded windows in
abo_variant_panel.bed / abo_variant_panel.yaml already account for.

gVCF non-variant blocks (ALT=<NON_REF>, INFO=END=...) give real
reference-depth coverage even where no variant is called -- this is
what recovers indel signal that Clair3's plain (non-gvcf) merge_output
silently drops (see c.1061delC / IMM-26-23099 in the panel notes).

No-evidence positions (v1.1.0)
------------------------------
Before v1.1.0, a position with neither a matching variant record nor a
reference block was written as 100% reference at the sample's MEDIAN block
depth. That made "Clair3 said nothing here" indistinguishable from a
confirmed reference call, and it gave every such position a depth it never
had. Now each panel position is classified as:

  variant   a record matched the panel alt allele      -> depth = record DP
  record    a record at the position, other allele     -> depth = record DP,
                                                         Match = real ref fraction
  block     covered by a gVCF reference block          -> depth = block MIN_DP, 100% ref
  none      nothing covers the position                -> Depth 0, Match 0 (no evidence)

Positions with no evidence are listed on stderr. The aggregator reports them
as "no evidence" rather than as reference. --legacy-fallback restores the old
behaviour for reproducing earlier runs.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from abo_panel import load_panel, Panel, VariantMarker

__version__ = "1.1.0"

TSV_HEADER = "\t".join([
    "Ref_Position_1based", "Ref_Base", "Match_Percent", "Mismatch_Percent",
    "Insertion_Percent", "Deletion_Percent", "A_Percent", "G_Percent",
    "C_Percent", "T_Percent", "Depth",
])

INDEL_LIKE_TYPES = {"deletion", "insertion", "indel", "snp_or_indel"}
# The +-15bp window is for FINDING an indel record that the caller anchored away from the
# nominal position. Evidence that a position was actually covered (a reference block, or a
# record of another allele) must be tight: a block 15bp away says nothing about this base.
EVIDENCE_WINDOW = 2
INFO_END_RE = re.compile(r"(?:^|;)END=(\d+)")
BUCKET_MAP = {"del": "Del", "dup": "Ins", "ins": "Ins"}


class VcfRecord:
    __slots__ = ("pos", "ref", "alts", "ad", "dp")

    def __init__(self, pos: int, ref: str, alts: List[str], ad: List[int], dp: int):
        self.pos = pos
        self.ref = ref
        self.alts = alts   # real ALT alleles, <NON_REF> placeholder stripped
        self.ad = ad        # [ref_depth, alt1_depth, alt2_depth, ...]
        self.dp = dp


RefBlocks = Tuple[List[int], List[int], List[int]]  # (starts, ends, min_dp), sorted by start


def parse_gvcf(path: Path) -> Tuple[Dict[int, List[VcfRecord]], RefBlocks]:
    """Parse a Clair3 --gvcf file into (variant records by POS, ref-block index)."""
    records: Dict[int, List[VcfRecord]] = {}
    block_starts, block_ends, block_dp = [], [], []

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]
            info = fields[7]
            fmt_keys = fields[8].split(":")
            sample_vals = fields[9].split(":")
            fmt = dict(zip(fmt_keys, sample_vals))

            raw_alts = alt_field.split(",") if alt_field not in (".", "") else []
            real_alts = [a for a in raw_alts if a != "<NON_REF>"]
            has_non_ref = "<NON_REF>" in raw_alts

            if not real_alts:
                m = INFO_END_RE.search(info)
                end = int(m.group(1)) if m else pos
                min_dp_raw = fmt.get("MIN_DP") or fmt.get("DP") or "0"
                try:
                    dp = int(min_dp_raw)
                except ValueError:
                    dp = 0
                block_starts.append(pos)
                block_ends.append(end)
                block_dp.append(dp)
                continue

            dp_raw = fmt.get("DP", "0")
            try:
                dp = int(dp_raw)
            except ValueError:
                dp = 0
            ad_raw = fmt.get("AD", "")
            try:
                ad_full = [int(x) for x in ad_raw.split(",")] if ad_raw not in ("", ".") else []
            except ValueError:
                ad_full = []
            if has_non_ref and len(ad_full) == len(real_alts) + 2:
                ad = ad_full[:-1]  # drop the trailing <NON_REF> allele depth
            else:
                ad = ad_full

            records.setdefault(pos, []).append(VcfRecord(pos, ref, real_alts, ad, dp))

    order = sorted(range(len(block_starts)), key=lambda i: block_starts[i])
    ref_blocks = (
        [block_starts[i] for i in order],
        [block_ends[i] for i in order],
        [block_dp[i] for i in order],
    )
    return records, ref_blocks


def dp_at(pos: int, ref_blocks: RefBlocks) -> int:
    starts, ends, dps = ref_blocks
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0 and starts[i] <= pos <= ends[i]:
        return dps[i]
    return 0


def median(values: List[int]) -> float:
    if not values:
        return 0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2


def needs_padding(marker: VariantMarker) -> bool:
    return marker.variant_type in INDEL_LIKE_TYPES or marker.ref_base == "REF"


def bucket_for_token(tok: str) -> str:
    return BUCKET_MAP.get(tok.lower(), tok.strip().upper())


def find_best_match(
    records_by_pos: Dict[int, List[VcfRecord]], center_pos: int, token: str, window: int,
) -> Tuple[Optional[VcfRecord], float]:
    """Search +-window around center_pos for a record matching `token`'s
    length profile (del/ins/dup -> indel by REF/ALT length comparison;
    a single letter -> SNP with that exact ALT base)."""
    tok_l = token.lower()
    is_indel_tok = tok_l in ("del", "ins", "dup")
    best_rec, best_pct = None, 0.0

    for p in range(center_pos - window, center_pos + window + 1):
        for rec in records_by_pos.get(p, []):
            if not rec.dp:
                continue
            for i, alt in enumerate(rec.alts):
                if is_indel_tok:
                    match = (len(rec.ref) > len(alt)) if tok_l == "del" else (len(rec.ref) < len(alt))
                else:
                    match = (len(rec.ref) == 1 and len(alt) == 1 and alt.upper() == tok_l.upper())
                if not match:
                    continue
                ad_idx = i + 1
                alt_depth = rec.ad[ad_idx] if ad_idx < len(rec.ad) else 0
                pct = 100.0 * alt_depth / rec.dp
                if pct > best_pct or best_rec is None:
                    best_rec, best_pct = rec, pct
    return best_rec, best_pct


def nearest_record(records_by_pos: Dict[int, List[VcfRecord]], center_pos: int,
                   window: int) -> Optional[VcfRecord]:
    """Closest VCF record (any allele, with depth) within +-window of center_pos."""
    for off in range(0, window + 1):
        for p in ((center_pos,) if off == 0 else (center_pos - off, center_pos + off)):
            for rec in records_by_pos.get(p, []):
                if rec.dp:
                    return rec
    return None


def block_dp_near(pos: int, window: int, ref_blocks: RefBlocks) -> int:
    """Reference-block depth at pos, or at the nearest covered position within +-window."""
    for off in range(0, window + 1):
        for p in ((pos,) if off == 0 else (pos - off, pos + off)):
            d = dp_at(p, ref_blocks)
            if d:
                return d
    return 0


def marker_row(
    marker: VariantMarker,
    records_by_pos: Dict[int, List[VcfRecord]],
    ref_blocks: RefBlocks,
    fallback_dp: float,
    legacy_fallback: bool = False,
) -> Tuple[str, str]:
    """Build one TSV row (matching TSV_HEADER) for a single panel marker.
    Returns (row, evidence); evidence is 'variant', 'record', 'block',
    'fallback' (legacy only) or 'none'."""
    row = {"A": 0.0, "G": 0.0, "C": 0.0, "T": 0.0, "Del": 0.0, "Ins": 0.0}
    pos = marker.resolved_position()
    window = 15 if needs_padding(marker) else 0

    matched_rec = None
    for alt_tok in marker.alt_bases():
        rec, pct = find_best_match(records_by_pos, pos, alt_tok, window)
        row[bucket_for_token(alt_tok)] = pct
        if rec is not None and matched_rec is None:
            matched_rec = rec

    ev_rec = matched_rec
    if matched_rec is not None:
        evidence, depth = "variant", matched_rec.dp
    else:
        ev_rec = nearest_record(records_by_pos, pos, min(window, EVIDENCE_WINDOW))
        if ev_rec is not None:
            evidence, depth = "record", ev_rec.dp
        else:
            bdp = block_dp_near(pos, min(window, EVIDENCE_WINDOW), ref_blocks)
            if bdp:
                evidence, depth = "block", bdp
            elif legacy_fallback:
                evidence, depth = "fallback", int(fallback_dp)
            else:
                evidence, depth = "none", 0

    ref_base = marker.ref_base if marker.ref_base not in ("REF", "") else "N"
    ref_letter = marker.ref_base.strip().upper() if marker.ref_base not in ("REF", "") else None

    def fmt(v: float) -> str:
        return str(int(v)) if float(v).is_integer() else f"{v:.2f}"

    if evidence == "none":
        # nothing covers this position: do not claim a reference call
        return "\t".join([str(pos), ref_base, "0", "0", "0", "0", "0", "0", "0", "0", "0"]), evidence

    alt_sum = sum(row.values())
    if ref_letter:
        if ev_rec is not None and ev_rec.dp and ev_rec.ad:
            ref_pct = 100.0 * ev_rec.ad[0] / ev_rec.dp
        else:
            ref_pct = 100.0  # gVCF reference block (or legacy fallback)
        row[ref_letter] = ref_pct
        mismatch = alt_sum
    else:
        ref_pct = None
        mismatch = alt_sum

    if evidence == "record":
        # a record sits here but carries an allele the panel does not score;
        # report the real reference fraction, and count the rest as mismatch
        match_pct = ref_pct if ref_pct is not None else max(0.0, 100.0 - alt_sum)
        mismatch = max(0.0, 100.0 - match_pct)
    else:
        match_pct = max(0.0, 100.0 - mismatch)
        mismatch = max(0.0, mismatch)

    return "\t".join([
        str(pos), ref_base,
        fmt(match_pct), fmt(mismatch),
        fmt(row["Ins"]), fmt(row["Del"]),
        fmt(row["A"]), fmt(row["G"]), fmt(row["C"]), fmt(row["T"]),
        str(depth),
    ]), evidence


def convert(gvcf_path: Path, panel: Panel, output_path: Path,
            legacy_fallback: bool = False,
            no_evidence_out: Optional[List[str]] = None) -> int:
    records_by_pos, ref_blocks = parse_gvcf(gvcf_path)
    fallback_dp = median(ref_blocks[2])

    markers = sorted(panel.calibrated(), key=lambda m: m.resolved_position())
    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(TSV_HEADER + "\n")
        for marker in markers:
            line, evidence = marker_row(marker, records_by_pos, ref_blocks, fallback_dp, legacy_fallback)
            if evidence == "none" and no_evidence_out is not None:
                no_evidence_out.append(marker.cdna_change.split()[0])
            f.write(line + "\n")
    return len(markers)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert a Clair3 --gvcf file into HAPLOSCAN-format "
                    "AlignmentStatistics.tsv metrics for the ABO panel positions.",
    )
    ap.add_argument("-i", "--gvcf", required=True, help="Clair3 merge_output.gvcf.gz")
    ap.add_argument("-o", "--output", required=True, help="Output *.AlignmentStatistics.tsv path")
    ap.add_argument("--panel", default="abo_variant_panel.yaml",
                     help="Path to the variant panel (.yaml/.json/.csv/.tsv). Default: %(default)s")
    ap.add_argument("--legacy-fallback", action="store_true",
                     help="Pre-1.1.0 behaviour: a position with no record and no reference block is written "
                          "as 100%% reference at the sample's median block depth instead of as 'no evidence'.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel)
    no_ev: List[str] = []
    n = convert(Path(args.gvcf), panel, Path(args.output),
                legacy_fallback=args.legacy_fallback, no_evidence_out=no_ev)
    print(f"Wrote {n} panel-position rows to {args.output}")
    if no_ev:
        print(f"WARNING: no gVCF evidence (no variant record and no reference block) at {len(no_ev)} "
              f"position(s): {', '.join(no_ev)}. Written as depth 0 / no match, not as reference.",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
