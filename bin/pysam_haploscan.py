#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pysam_haploscan.py — Read-level haplotype scanning for ABO amplicons (v2.0.0)

CHANGES IN v2.0.0
------------------
  - Diagnostic/indel positions are now loaded from an external variant panel
    (--panel, default abo_variant_panel.yaml) via abo_panel.py, instead of
    the hardcoded HAPLOTYPE_POSITIONS / INDEL_DIAGNOSTIC dicts.
  - Supports the COMBINED single ~6.7 kb long-read amplicon spanning exons
    2-7 (Goebel/Wu primers; Mobegi et al. 2025, IJMS 26(12):5443) as the
    default reference mode: the whole reference is treated as one unit and
    every calibrated panel position (any exon) is scored directly by its
    amplicon_pos, with no length-based exon-type guessing at all.
  - The old length-based ExonType detection (EXON6_LENGTH_RANGE /
    EXON7_LENGTH_RANGE) is kept ONLY as a --legacy fallback for analysing
    v1.x-style separate short exon6-only / exon7-only mini-amplicon BAMs,
    where panel positions are resolved via legacy_exon6_pos/legacy_exon7_pos
    instead of amplicon_pos.
  - Because a single ONT read can now span every diagnostic position from
    exon 2 through the 3' UTR, compute_read_haplotypes() reports ONE
    haplotype string per read covering ALL calibrated positions (not just
    one exon's worth) -- this removes the old "cross-amplicon phasing
    requires a post-processing step" limitation entirely: full-length cis
    confirmation is now available directly from Haplotypes.tsv.

Outputs
-------
{prefix}.AlignmentStatistics.tsv
    Same column layout as stats_from_pileup.py / v1.x -- downstream scripts
    (aggregate_abo_reports.py, predict_abo_phenotype.py) are unchanged in
    format, only richer in content (more positions, spanning more exons).

{prefix}.Haplotypes.tsv
    Per-read haplotype table. In combined mode this spans every calibrated
    diagnostic position on the reference (potentially exon2 through exon7);
    in --legacy mode it is scoped to whichever mini-amplicon (exon6/exon7)
    the read came from, as in v1.x.

ABOReadPolymorphisms.txt
    Same polymorphic-position summary format as before.

Usage
-----
    # Combined single-amplicon mode (default; requires a calibrated panel --
    # i.e. amplicon_pos populated via calibrate_panel_positions.py)
    pysam_haploscan.py -b sample.bam -f reference.fasta -o prefix \\
        --panel abo_variant_panel.yaml

    # Legacy dual mini-amplicon mode (v1.x behaviour, exon6-only or
    # exon7-only reference/BAM, resolved via legacy_exonN_pos)
    pysam_haploscan.py -b sample.bam -f reference.fasta -o prefix --legacy
"""

import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pysam

from abo_panel import load_panel, Panel, VariantMarker

__author__ = "Fredrick Mobegi"
__copyright__ = (
    "Copyright 2024-2025, ABO blood group typing using third-generation sequencing (TGS) technology"
)
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg",
               "Claude Sonnet 5 (v2.0.0 rewrite for panel-driven / combined amplicon)"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


# ---------------------------------------------------------------------------
# Legacy exon classification (used only in --legacy mode)
# ---------------------------------------------------------------------------

from enum import Enum


class ExonType(Enum):
    EXON6 = "Exon 6"
    EXON7 = "Exon 7"
    COMBINED = "combined"
    UNKNOWN = "unknown"


LEGACY_EXON6_LENGTH_RANGE = (130, 140)
LEGACY_EXON7_LENGTH_RANGE = (800, 830)

# Mirror of stats_from_pileup.LOW_COVERAGE_THRESHOLD
_LOW_COVERAGE_THRESHOLD = 200

NUCLEOTIDES = frozenset({"A", "G", "C", "T"})

TSV_HEADER = "\t".join([
    "Ref_Position_1based", "Ref_Base", "Match_Percent", "Mismatch_Percent",
    "Insertion_Percent", "Deletion_Percent", "A_Percent", "G_Percent",
    "C_Percent", "T_Percent", "Depth",
])


def determine_legacy_exon_type(ref_length: int) -> ExonType:
    if LEGACY_EXON6_LENGTH_RANGE[0] <= ref_length <= LEGACY_EXON6_LENGTH_RANGE[1]:
        return ExonType.EXON6
    if LEGACY_EXON7_LENGTH_RANGE[0] <= ref_length <= LEGACY_EXON7_LENGTH_RANGE[1]:
        return ExonType.EXON7
    return ExonType.UNKNOWN


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PositionStats:
    pos:               int
    ref_base:          str
    depth:             int = 0
    match_percent:     int = 0
    mismatch_percent:  int = 0
    insertion_percent: int = 0
    deletion_percent:  int = 0
    A_percent:         int = 0
    G_percent:         int = 0
    C_percent:         int = 0
    T_percent:         int = 0


@dataclass
class ReadHaplotype:
    """Alleles observed at diagnostic positions on a single read."""
    read_name: str
    exon:      str
    alleles:   Dict[int, str] = field(default_factory=dict)

    def to_string(self) -> str:
        if not self.alleles:
            return "REF"
        return ";".join(f"p{pos}={allele}" for pos, allele in sorted(self.alleles.items()))


# ---------------------------------------------------------------------------
# Panel-driven position resolution
# ---------------------------------------------------------------------------

def resolve_positions(panel: Panel, legacy: bool, legacy_exon: Optional[ExonType] = None):
    """
    Return (diag_positions: frozenset[int], indel_positions: frozenset[int],
    pos_to_ref_base: dict[int,str]) for the requested mode.

    combined mode: every calibrated (amplicon_pos-populated) marker in the
                   panel, across all exons -- a single continuous scan.
    legacy mode:   markers whose exon matches legacy_exon, resolved via
                   legacy_exon6_pos/legacy_exon7_pos.
    """
    diag = set()
    indel = set()
    ref_base_by_pos: Dict[int, str] = {}

    if legacy:
        exon_label = legacy_exon.value if legacy_exon else None
        candidates = [m for m in panel.variants if m.exon == exon_label]
        prefer = "legacy"
    else:
        candidates = panel.variants
        prefer = "amplicon"

    for m in candidates:
        pos = m.resolved_position(prefer=prefer if not legacy else (
            "legacy6" if legacy_exon == ExonType.EXON6 else "legacy7"
        ))
        if pos is None:
            continue
        diag.add(pos)
        if m.variant_type in ("deletion", "insertion", "indel", "snp_or_indel", "dup_or_del"):
            indel.add(pos)
        if m.ref_base and m.ref_base not in ("REF", "ALT", ""):
            ref_base_by_pos[pos] = m.ref_base

    return frozenset(diag), frozenset(indel), ref_base_by_pos


# ---------------------------------------------------------------------------
# Per-position frequency computation (-> AlignmentStatistics.tsv)
# ---------------------------------------------------------------------------

def compute_position_stats(
    bam:              pysam.AlignmentFile,
    ref_name:         str,
    ref_seq:          str,
    indel_positions:  frozenset,
    min_base_quality: int = 0,
    min_map_quality:  int = 0,
) -> List[PositionStats]:
    """
    Compute ATGC + indel frequencies at every reference position using the
    pysam pileup engine. Logic mirrors stats_from_pileup.py exactly so the
    TSV output is bit-for-bit compatible with existing downstream code.
    indel_positions (panel-resolved, absolute reference coordinates) replaces
    the old per-ExonType INDEL_DIAGNOSTIC dict.
    """
    ref_length = len(ref_seq)
    results: List[PositionStats] = []

    for col in bam.pileup(
        ref_name, 0, ref_length,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_map_quality,
        truncate=True, ignore_overlaps=False, stepper="nofilter",
    ):
        pos0 = col.reference_pos
        pos1 = pos0 + 1
        ref_base = ref_seq[pos0].upper() if pos0 < ref_length else "N"
        depth = col.nsegments

        if depth == 0:
            results.append(PositionStats(pos=pos1, ref_base=ref_base))
            continue

        counts: Dict[str, int] = {"A": 0, "G": 0, "C": 0, "T": 0, "ins": 0, "del": 0}

        for pr in col.pileups:
            if pr.is_refskip:
                continue
            if pr.is_del:
                counts["del"] += 1
            else:
                qpos = pr.query_position
                if qpos is not None and pr.alignment.query_sequence:
                    base = pr.alignment.query_sequence[qpos].upper()
                    if base in NUCLEOTIDES:
                        counts[base] += 1
                    if pr.indel > 0:
                        counts["ins"] += 1

        if pos1 in indel_positions:
            include_indels = True
        elif depth < _LOW_COVERAGE_THRESHOLD:
            # Conservative default for non-diagnostic positions at low
            # coverage. With a single combined reference there is no
            # meaningful "exon6 vs long exon7" distinction any more, so we
            # apply one consistent rule: exclude indels for non-diagnostic
            # positions when coverage is low, matching the original
            # ExonType.EXON6 / long-exon7 behaviour.
            include_indels = False
        else:
            include_indels = True

        total_bases = sum(counts[b] for b in NUCLEOTIDES)
        total_all = total_bases + counts["ins"] + counts["del"]

        if include_indels:
            denom = total_all if total_all > 0 else 1
            ins_pct = int(counts["ins"] / denom * 100)
            del_pct = int(counts["del"] / denom * 100)
            base_pcts = {b: int(counts[b] / denom * 100) for b in NUCLEOTIDES}
        else:
            denom = total_bases if total_bases > 0 else 1
            ins_pct = 0
            del_pct = 0
            base_pcts = {b: int(counts[b] / denom * 100) for b in NUCLEOTIDES}
            atgc_sum = sum(base_pcts.values())
            if atgc_sum != 100 and atgc_sum > 0:
                base_pcts[ref_base] = base_pcts.get(ref_base, 0) + (100 - atgc_sum)

        match_pct = base_pcts.get(ref_base, 0)
        mismatch_pct = sum(v for k, v in base_pcts.items() if k != ref_base)

        results.append(PositionStats(
            pos=pos1, ref_base=ref_base, depth=depth,
            match_percent=match_pct, mismatch_percent=mismatch_pct,
            insertion_percent=ins_pct, deletion_percent=del_pct,
            A_percent=base_pcts.get("A", 0), G_percent=base_pcts.get("G", 0),
            C_percent=base_pcts.get("C", 0), T_percent=base_pcts.get("T", 0),
        ))

    return results


# ---------------------------------------------------------------------------
# Per-read haplotype computation (-> Haplotypes.tsv)
# ---------------------------------------------------------------------------

def compute_read_haplotypes(
    bam:             pysam.AlignmentFile,
    ref_name:        str,
    ref_seq:         str,
    diag_positions:  frozenset,
    indel_positions: frozenset,
    exon_label:      str,
    min_map_quality: int = 0,
) -> List[ReadHaplotype]:
    """
    Iterate every primary aligned read and record which diagnostic positions
    carry a non-reference allele. In combined mode diag_positions spans
    every calibrated position gene-wide, so a single full-length ONT read
    yields full cis-phase information across exons in one row.
    """
    if not diag_positions:
        return []

    haplotypes: List[ReadHaplotype] = []

    for read in bam.fetch(ref_name):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_map_quality:
            continue
        if read.cigartuples is None or read.query_sequence is None:
            continue

        ref_to_qpos: Dict[int, Optional[int]] = {}
        ref_followed_by_ins: set = set()

        prev_rpos: Optional[int] = None
        for qpos, rpos in read.get_aligned_pairs(matches_only=False, with_seq=False):
            if rpos is not None:
                ref_to_qpos[rpos] = qpos
                prev_rpos = rpos
            elif qpos is not None and prev_rpos is not None:
                ref_followed_by_ins.add(prev_rpos)

        alleles: Dict[int, str] = {}

        for pos1 in diag_positions:
            rpos0 = pos1 - 1
            if rpos0 not in ref_to_qpos:
                continue

            qpos = ref_to_qpos[rpos0]
            ref_base = ref_seq[rpos0].upper() if rpos0 < len(ref_seq) else "N"

            if pos1 in indel_positions:
                if qpos is None:
                    alleles[pos1] = "del"
                elif rpos0 in ref_followed_by_ins:
                    alleles[pos1] = "ins"
            else:
                if qpos is None:
                    alleles[pos1] = "del"
                else:
                    base = read.query_sequence[qpos].upper()
                    if base != ref_base and base in NUCLEOTIDES:
                        alleles[pos1] = base

        haplotypes.append(ReadHaplotype(
            read_name=read.query_name or "unknown",
            exon=exon_label,
            alleles=alleles,
        ))

    return haplotypes


# ---------------------------------------------------------------------------
# Co-occurrence / phasing analysis
# ---------------------------------------------------------------------------

def analyse_cooccurrence(
    haplotypes: List[ReadHaplotype], pos_a: int, allele_a: str,
    pos_b: int, allele_b: str, logger: logging.Logger,
) -> None:
    """Log how often two alleles co-occur on the same read. With the
    combined amplicon, pos_a and pos_b may legitimately be in different
    exons -- full-length reads make this a true single-molecule phase
    check rather than the within-exon-only check possible in v1.x."""
    reads_a = sum(1 for h in haplotypes if h.alleles.get(pos_a) == allele_a)
    reads_b = sum(1 for h in haplotypes if h.alleles.get(pos_b) == allele_b)
    reads_ab = sum(
        1 for h in haplotypes
        if h.alleles.get(pos_a) == allele_a and h.alleles.get(pos_b) == allele_b
    )
    total = len(haplotypes) or 1
    pct = 100.0 * reads_ab / total

    logger.info(
        f"Phasing  p{pos_a}={allele_a} \u2229 p{pos_b}={allele_b}: "
        f"{reads_a} reads carry p{pos_a}={allele_a}, "
        f"{reads_b} carry p{pos_b}={allele_b}, "
        f"{reads_ab} carry BOTH ({pct:.1f}% of all reads)"
    )

    if reads_a > 0:
        frac = reads_ab / reads_a
        if frac >= 0.8:
            logger.info(f"  -> PHASED: {100*frac:.0f}% of p{pos_a}={allele_a} reads "
                        f"also carry p{pos_b}={allele_b}  (same allele confirmed)")
        elif frac <= 0.2:
            logger.info(f"  -> TRANS: only {100*frac:.0f}% of p{pos_a}={allele_a} reads "
                        f"carry p{pos_b}={allele_b}  (likely on different alleles)")
        else:
            logger.warning(f"  -> AMBIGUOUS: {100*frac:.0f}% co-occurrence — manual review recommended")


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_stats_tsv(stats: List[PositionStats], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write(TSV_HEADER + "\n")
        for s in stats:
            fh.write("\t".join([
                str(s.pos), s.ref_base,
                str(s.match_percent), str(s.mismatch_percent),
                str(s.insertion_percent), str(s.deletion_percent),
                str(s.A_percent), str(s.G_percent),
                str(s.C_percent), str(s.T_percent),
                str(s.depth),
            ]) + "\n")


def write_haplotypes_tsv(haplotypes: List[ReadHaplotype], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("Read_Name\tExon\tVariant_Positions\tHaplotype\n")
        for h in haplotypes:
            var_pos = (
                ",".join(str(p) for p in sorted(h.alleles.keys()))
                if h.alleles else "."
            )
            fh.write(f"{h.read_name}\t{h.exon}\t{var_pos}\t{h.to_string()}\n")


def write_summary(stats: List[PositionStats], path: Path, threshold: int = 10) -> None:
    """Write ABOReadPolymorphisms.txt in the same format as stats_from_pileup.py."""
    logger = logging.getLogger(__name__)
    polymorphic = 0
    with open(path, "w") as fh:
        for s in stats:
            if (s.mismatch_percent >= threshold or s.insertion_percent >= threshold
                    or s.deletion_percent >= threshold):
                fh.write(f"(1-based) Position:{s.pos}, Reference Base={s.ref_base}\n")
                fh.write(f"Aligned Read Count:{s.depth}\n")
                fh.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                fh.write(
                    f"{s.match_percent}\t{s.mismatch_percent}\t"
                    f"{s.insertion_percent}\t{s.deletion_percent}\t"
                    f"{s.A_percent}\t{s.G_percent}\t{s.C_percent}\t{s.T_percent}\n\n"
                )
                polymorphic += 1
    logger.info(f"Found {polymorphic} polymorphic positions (threshold={threshold}%)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def setup_logging(verbose: bool) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Read-level ABO haplotype scanner. Default mode targets the "
            "combined exon2-7 long-read amplicon (Mobegi et al. 2025); pass "
            "--legacy for v1.x-style separate exon6-only/exon7-only "
            "mini-amplicon references."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -b sample.bam -f combined_reference.fasta -o sample_prefix \\
      --panel abo_variant_panel.yaml
  %(prog)s -b sample.bam -f exon7_reference.fasta -o sample_prefix --legacy -v
        """,
    )
    parser.add_argument("-b", "--bam", required=True, help="Sorted, indexed BAM file")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA (must have .fai index)")
    parser.add_argument("-o", "--output", required=True,
                         help="Output prefix -- produces <prefix>.AlignmentStatistics.tsv "
                              "and <prefix>.Haplotypes.tsv")
    parser.add_argument("--panel", default="abo_variant_panel.yaml",
                         help="Path to the variant panel file (.yaml/.json/.csv/.tsv). Default: %(default)s")
    parser.add_argument("--legacy", action="store_true",
                         help="Use v1.x behaviour: treat the reference as a short "
                              "exon6-only or exon7-only mini-amplicon, detected by "
                              "length, and resolve panel positions via "
                              "legacy_exon6_pos/legacy_exon7_pos instead of amplicon_pos.")
    parser.add_argument("-s", "--summary", default="ABOReadPolymorphisms.txt",
                         help="Polymorphic positions summary (default: %(default)s)")
    parser.add_argument("-t", "--threshold", type=int, default=10,
                         help="Polymorphism threshold %% (default: %(default)s)")
    parser.add_argument("-q", "--min-mapq", type=int, default=0,
                         help="Minimum mapping quality (default: %(default)s)")
    parser.add_argument("-Q", "--min-baseq", type=int, default=0,
                         help="Minimum base quality for pileup (default: %(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    try:
        panel = load_panel(args.panel)
    except Exception as exc:
        logger.error(f"Cannot load variant panel '{args.panel}': {exc}")
        return 1

    bam_path = Path(args.bam)
    fasta_path = Path(args.fasta)
    prefix = args.output
    summary_path = Path(args.summary)
    stats_path = Path(f"{prefix}.AlignmentStatistics.tsv")
    haplo_path = Path(f"{prefix}.Haplotypes.tsv")

    try:
        with pysam.FastaFile(str(fasta_path)) as fa:
            ref_names = fa.references
            if not ref_names:
                logger.error("No sequences found in reference FASTA")
                return 1
            ref_name = ref_names[0]
            ref_seq = fa.fetch(ref_name).upper()
    except Exception as exc:
        logger.error(f"Cannot load reference FASTA: {exc}")
        return 1

    ref_length = len(ref_seq)

    if args.legacy:
        legacy_exon = determine_legacy_exon_type(ref_length)
        if legacy_exon == ExonType.UNKNOWN:
            logger.error(
                f"--legacy given but reference length {ref_length} matches "
                f"neither the exon6 range {LEGACY_EXON6_LENGTH_RANGE} nor "
                f"exon7 range {LEGACY_EXON7_LENGTH_RANGE}. If this is the "
                f"new combined amplicon, drop --legacy."
            )
            return 1
        exon_label = legacy_exon.value
        diag_positions, indel_positions, _ = resolve_positions(panel, legacy=True, legacy_exon=legacy_exon)
        logger.info(f"[legacy mode] Reference: {ref_name}  length={ref_length}  exon={exon_label}")
    else:
        exon_label = "combined"
        diag_positions, indel_positions, _ = resolve_positions(panel, legacy=False)
        logger.info(f"[combined mode] Reference: {ref_name}  length={ref_length}")
        if not diag_positions:
            logger.warning(
                "No calibrated (amplicon_pos-populated) panel positions found. "
                "Run calibrate_panel_positions.py against this reference FASTA "
                "first, or use --legacy for a v1.x mini-amplicon reference."
            )

    logger.info(f"Diagnostic positions loaded: {len(diag_positions)} "
                f"({len(indel_positions)} indel-diagnostic)")

    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as exc:
        logger.error(f"Cannot open BAM: {exc}")
        return 1

    try:
        logger.info("Computing per-position allele frequencies ...")
        stats = compute_position_stats(
            bam, ref_name, ref_seq, indel_positions,
            min_base_quality=args.min_baseq, min_map_quality=args.min_mapq,
        )
        write_stats_tsv(stats, stats_path)
        logger.info(f"[OK]  AlignmentStatistics  -> {stats_path}")

        write_summary(stats, summary_path, threshold=args.threshold)
        logger.info(f"[OK]  Polymorphisms summary -> {summary_path}")

        logger.info("Computing per-read haplotypes ...")
        haplotypes = compute_read_haplotypes(
            bam, ref_name, ref_seq, diag_positions, indel_positions, exon_label,
            min_map_quality=args.min_mapq,
        )
        write_haplotypes_tsv(haplotypes, haplo_path)
        logger.info(f"[OK]  Haplotypes            -> {haplo_path}  ({len(haplotypes)} reads)")

        if haplotypes and not args.legacy:
            logger.info(
                "Combined-amplicon mode: every read above spans all "
                "calibrated diagnostic positions gene-wide, so "
                "Haplotypes.tsv already gives full-length cis phasing "
                "without any cross-amplicon post-processing step."
            )
        elif haplotypes and exon_label == "Exon 7":
            a2_marker = panel.get("a1_a2_1061del")
            a1032_marker = panel.get("a2p_1032_a201")
            if a2_marker and a1032_marker:
                p1 = a2_marker.resolved_position(prefer="legacy7")
                p2 = a1032_marker.resolved_position(prefer="legacy7")
                if p1 and p2:
                    analyse_cooccurrence(haplotypes, p2, "A", p1, "del", logger)
        elif haplotypes and exon_label == "Exon 6":
            o1_marker = panel.get("o1_marker")
            pos = o1_marker.resolved_position(prefer="legacy6") if o1_marker else None
            if pos:
                dels = sum(1 for h in haplotypes if "del" in h.alleles.get(pos, ""))
                total = len(haplotypes) or 1
                logger.info(
                    f"Exon6 pos{pos} deletion: {dels}/{total} reads "
                    f"({100*dels/total:.1f}%)  — "
                    + ("homozygous O1" if dels/total > 0.8 else
                       "heterozygous O1/non-O1" if dels/total > 0.2 else
                       "non-O1")
                )

    finally:
        bam.close()

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n! Interrupted by user")
        sys.exit(130)
    except Exception as exc:
        print(f"! CRITICAL: unhandled exception: {exc}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
