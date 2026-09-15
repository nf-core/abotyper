#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
clair2haplotypes.py -- convert a Clair3 phased VCF into the same
Haplotypes.tsv format pysam_haploscan.py produces, so aggregate_abo_reports.py's
existing PhaseConfidence scoring (parse_haplotypes_tsv / compute_phase_confidence)
works unchanged against Clair3 phasing.

Clair3's phased VCF (--enable_phasing) gives block-level phasing: GT uses
"|" (e.g. 0|1) and a PS FORMAT tag groups variants known to sit on the same
physical molecule, backed by its own read-based phasing algorithm. It does
NOT give per-read genotypes the way pysam_haploscan.py's BAM scan did, so
this is not a byte-for-byte reconstruction of the original per-read table --
it synthesizes one row per "virtual read" per haplotype side, in the same
proportions as each phase set's own allele depths (AD), which is the
information Clair3's phasing algorithm actually had available. Each phase
set becomes exactly 2 haplotype strings (one per side); a haplotype string
is built as pysam_haploscan.py's ReadHaplotype.to_string() does: "REF" if a
side carries no panel-position ALT, else "p<pos>=<base>;p<pos>=<base>..."
for whichever panel positions differ from reference on that side.

Only heterozygous, phased (GT contains "|") panel positions with the SAME
value on both sides ignored -- a homozygous position carries no phase
information and is dropped from the haplotype signature; it does not
distinguish the two sides.

If a sample's positions split across multiple phase sets (reads didn't
fully span the amplicon, or a low-confidence region broke phasing), each
phase set contributes its own 2-haplotype pair -- compute_phase_confidence
will then correctly report lower confidence, reflecting genuine phasing
fragmentation rather than a single clean phase block.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Tuple

from abo_panel import load_panel, Panel, VariantMarker

__version__ = "1.0.0"

INDEL_LIKE_TYPES = {"deletion", "insertion", "indel", "snp_or_indel"}


def needs_padding(marker: VariantMarker) -> bool:
    """Same rule as clair2metrics.py / panel_to_bed.py: only true indel,
    homopolymer, and multi-offset placeholder rows get a search window --
    plain SNP rows stay exact, since several panel positions sit only 2-5bp
    apart and a shared window would let one marker's record match another."""
    return marker.variant_type in INDEL_LIKE_TYPES or marker.ref_base == "REF"


class PhasedRecord:
    __slots__ = ("pos", "ref", "alts", "gt", "ps", "ad", "dp")

    def __init__(self, pos, ref, alts, gt, ps, ad, dp):
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.gt = gt   # tuple of two ints, e.g. (0, 1)
        self.ps = ps   # phase set id (string), or None if unphased/absent
        self.ad = ad   # [ref_depth, alt1_depth, ...]
        self.dp = dp


def parse_phased_vcf(path: Path) -> Dict[int, PhasedRecord]:
    """Return {POS: PhasedRecord} for every phased, heterozygous record."""
    records: Dict[int, PhasedRecord] = {}
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
            if alt_field in (".", ""):
                continue
            alts = alt_field.split(",")

            fmt_keys = fields[8].split(":")
            sample_vals = fields[9].split(":")
            fmt = dict(zip(fmt_keys, sample_vals))

            gt_raw = fmt.get("GT", "")
            if "|" not in gt_raw:
                continue  # not phased at this site
            try:
                gt = tuple(int(x) for x in gt_raw.split("|"))
            except ValueError:
                continue
            if len(gt) != 2 or gt[0] == gt[1]:
                continue  # homozygous -- carries no phase information

            ps = fmt.get("PS")

            ad_raw = fmt.get("AD", "")
            try:
                ad = [int(x) for x in ad_raw.split(",")] if ad_raw not in ("", ".") else []
            except ValueError:
                ad = []
            dp_raw = fmt.get("DP", "0")
            try:
                dp = int(dp_raw)
            except ValueError:
                dp = 0

            records[pos] = PhasedRecord(pos, ref, alts, gt, ps, ad, dp)
    return records


def find_phased_record(records: Dict[int, PhasedRecord], center_pos: int,
                        window: int) -> Optional[PhasedRecord]:
    """Closest phased record within +-window of center_pos (0 offset first)."""
    if center_pos in records:
        return records[center_pos]
    for offset in range(1, window + 1):
        if center_pos - offset in records:
            return records[center_pos - offset]
        if center_pos + offset in records:
            return records[center_pos + offset]
    return None


def allele_at(rec: PhasedRecord, side: int) -> str:
    """The base/event on one GT side: 'REF' or the ALT string (del/ins/base)."""
    idx = rec.gt[side]
    if idx == 0:
        return "REF"
    alt = rec.alts[idx - 1] if idx - 1 < len(rec.alts) else "?"
    if len(rec.ref) > len(alt):
        return "del"
    if len(rec.ref) < len(alt):
        return "ins"
    return alt.upper()


def haplotype_string(entries: List[Tuple[int, str]]) -> str:
    """Match ReadHaplotype.to_string(): 'REF' or 'p<pos>=<allele>;...'"""
    non_ref = [(pos, allele) for pos, allele in entries if allele != "REF"]
    if not non_ref:
        return "REF"
    return ";".join(f"p{pos}={allele}" for pos, allele in sorted(non_ref))


def build_rows(panel: Panel, records: Dict[int, PhasedRecord], window: int = 15):
    """Group calibrated panel markers' phased records by PS, and emit
    (Read_Name, Exon, Variant_Positions, Haplotype) rows -- 2 synthesized
    haplotype strings per phase set, repeated in proportion to that phase
    set's own median allele depth."""
    by_ps: Dict[str, List[Tuple[int, PhasedRecord, str]]] = {}

    for marker in panel.calibrated():
        pos = marker.resolved_position()
        marker_window = window if needs_padding(marker) else 0
        rec = find_phased_record(records, pos, marker_window)
        if rec is None:
            continue
        ps_key = rec.ps if rec.ps else f"unnamed_ps_at_{rec.pos}"
        by_ps.setdefault(ps_key, []).append((pos, rec, marker.exon))

    rows = []
    for ps_key, entries in by_ps.items():
        if len(entries) < 1:
            continue
        side0 = [(pos, allele_at(rec, 0)) for pos, rec, _ in entries]
        side1 = [(pos, allele_at(rec, 1)) for pos, rec, _ in entries]
        hap0 = haplotype_string(side0)
        hap1 = haplotype_string(side1)
        var_positions = ",".join(str(pos) for pos, _, _ in sorted(entries))
        exons = sorted({exon for _, _, exon in entries})
        exon_label = "+".join(exons) if len(exons) > 1 else (exons[0] if exons else "")

        ad0_vals = [rec.ad[0] for _, rec, _ in entries if rec.ad]
        ad1_vals = [rec.ad[rec.gt[1]] for _, rec, _ in entries
                    if rec.ad and rec.gt[1] < len(rec.ad)]
        n0 = int(median(ad0_vals)) if ad0_vals else 1
        n1 = int(median(ad1_vals)) if ad1_vals else 1

        for i in range(n0):
            rows.append((f"synthetic_{ps_key}_A{i}", exon_label, var_positions, hap0))
        for i in range(n1):
            rows.append((f"synthetic_{ps_key}_B{i}", exon_label, var_positions, hap1))

    return rows


def write_haplotypes_tsv(rows, output_path: Path) -> None:
    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("Read_Name\tExon\tVariant_Positions\tHaplotype\n")
        for read_name, exon, var_pos, haplo in rows:
            f.write(f"{read_name}\t{exon}\t{var_pos}\t{haplo}\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert a Clair3 phased VCF into a HAPLOSCAN-format "
                    "Haplotypes.tsv (synthesized from phase-set AD, not per-read).",
    )
    ap.add_argument("-i", "--phased-vcf", required=True, help="Clair3 phased_merge_output.vcf.gz")
    ap.add_argument("-o", "--output", required=True, help="Output *.Haplotypes.tsv path")
    ap.add_argument("--panel", default="abo_variant_panel.yaml",
                     help="Path to the variant panel (.yaml/.json/.csv/.tsv). Default: %(default)s")
    ap.add_argument("--window", type=int, default=15,
                     help="Search window (bp) around each panel position for a phased "
                          "record, matching panel BED padding. Default: %(default)s")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel)
    records = parse_phased_vcf(Path(args.phased_vcf))
    rows = build_rows(panel, records, window=args.window)
    write_haplotypes_tsv(rows, Path(args.output))

    n_ps = len({r[3] for r in rows}) // 2 if rows else 0
    print(f"Wrote {len(rows)} synthesized haplotype rows ({n_ps} phase set(s), "
          f"{len(records)} phased records matched) to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
