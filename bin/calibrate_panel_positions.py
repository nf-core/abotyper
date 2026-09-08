#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calibrate_panel_positions.py — fill in amplicon_pos for the ABO variant panel
against your finalised combined exon2-7 reference FASTA (v2.0.0 utility).

WHY THIS EXISTS
---------------
The panel ships with amplicon_pos left blank for every new (exon2-6 / new
exon7 subtype-marker) position, because the exact coordinate depends on
YOUR reference FASTA -- which coordinate system it uses, exactly where its
first base sits relative to the ABO transcript, and exactly how long each
intron is in the sequence you built. Nobody outside your lab can respond
with correct numbers without that file, so rather than guess, this script
computes them for you, two ways:

MODE 1 -- anchor mode (no external tools required, always available)
    If you still have the OLD exon6-only and/or exon7-only mini-amplicon
    reference FASTAs from pipeline v1.x, this mode finds each one as a
    (possibly fuzzy) substring inside your new combined reference and uses
    that to shift every position that already has a legacy_exon6_pos /
    legacy_exon7_pos onto the new reference's coordinate system. This
    covers the 25 positions the assay already typed in v1.x immediately,
    with zero new dependencies.

MODE 2 -- spliced-alignment mode (requires a SAM/PAF alignment you supply)
    For genuinely new positions (exons 2-5, introns) that have no legacy
    coordinate, the only correct way to get their position in your new
    reference is to align the ABO mRNA/transcript (matching the cDNA
    numbering used by ISBT, e.g. NM_020469.3) against your reference with a
    SPLICED aligner (introns must be modelled as large gaps), e.g.:

        minimap2 -ax splice NM_020469.3.fasta your_combined_reference.fasta \\
            > abo_transcript_vs_reference.sam

    Then run:

        python3 calibrate_panel_positions.py \\
            --panel abo_variant_panel.yaml \\
            --from-sam abo_transcript_vs_reference.sam \\
            --out abo_variant_panel.calibrated.yaml

    This script parses the SAM CIGAR string itself (no aligner dependency
    at parse time) to build a cDNA-position -> reference-position map and
    fills in amplicon_pos for every panel row whose cdna_pos falls inside
    an aligned (matched/mismatched) block. Positions that fall inside a
    deletion relative to the reference, or outside the alignment entirely,
    are left blank with a warning -- that's a real "this variant isn't in
    your amplicon" signal worth investigating, not a bug to silence.

Both modes can be run in sequence (anchor mode first, then SAM mode to fill
in the rest) -- later modes only fill in blanks, never overwrite an
already-calibrated amplicon_pos, unless --force is given.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from abo_panel import load_panel, export_csv, FIELD_NAMES, Panel, VariantMarker

try:
    import yaml
except ImportError:
    yaml = None


# ---------------------------------------------------------------------------
# FASTA helpers (no external dependency)
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def find_best_substring_offset(haystack: str, needle: str, max_mismatch_frac: float = 0.05) -> Optional[int]:
    """
    Find the 0-based index in `haystack` where `needle` best matches,
    allowing up to max_mismatch_frac mismatches (handles the new reference
    having a slightly different consensus/basecall at a few positions).
    Returns None if no sufficiently good match is found.

    For a legacy mini-amplicon (~130-830 bp) against a ~6.7 kb haystack this
    naive O(n*m) scan is fast enough (a few hundred thousand comparisons).
    """
    n, m = len(haystack), len(needle)
    if m == 0 or m > n:
        return None
    max_mismatches = int(m * max_mismatch_frac)

    # Fast path: exact match first.
    exact = haystack.find(needle)
    if exact != -1:
        return exact

    best_offset, best_mismatches = None, max_mismatches + 1
    for offset in range(0, n - m + 1):
        window = haystack[offset:offset + m]
        mismatches = sum(1 for a, b in zip(window, needle) if a != b)
        if mismatches < best_mismatches:
            best_mismatches, best_offset = mismatches, offset
            if mismatches == 0:
                break
    if best_offset is not None and best_mismatches <= max_mismatches:
        return best_offset
    return None


def ng006669_calibrate(
    panel: Panel,
    ng_seq: str,
    new_ref_seq: str,
    anchor_len: int = 300,
    force: bool = False,
) -> Tuple[int, List[str]]:
    """
    Calibrate amplicon_pos directly from each row's ng006669_2_pos, given
    BOTH the full NG_006669.2 sequence and your new combined-amplicon
    reference sequence. This is the preferred calibration path whenever a
    row already has ng006669_2_pos populated (see the panel's meta.notes) --
    it needs no exon6/7 legacy references and, unlike SAM mode, works even
    for rows whose amplicon-mode position hasn't been filled by anything
    else yet.

    Method: locate the new reference as a (possibly reverse-complemented)
    substring of NG_006669.2 using its first `anchor_len` bases as an
    anchor (fast exact/near-exact search), establish the offset between the
    two coordinate systems, then apply that single offset to every row's
    ng006669_2_pos.
    """
    warnings: List[str] = []

    def revcomp(s: str) -> str:
        comp = str.maketrans("ACGTN", "TGCAN")
        return s.translate(comp)[::-1]

    anchor = new_ref_seq[:anchor_len]
    offset = find_best_substring_offset(ng_seq, anchor)
    strand = "+"
    if offset is None:
        anchor_rc = revcomp(anchor)
        offset = find_best_substring_offset(ng_seq, anchor_rc)
        strand = "-"

    if offset is None:
        warnings.append(
            "Could not locate your combined reference (forward or "
            "reverse-complement) inside NG_006669.2, even allowing 5% "
            "mismatch on the first "
            f"{anchor_len} bp. Your new reference may not be built from "
            "this RefSeqGene, or may need a longer/cleaner anchor region "
            "(try increasing --ng-anchor-len or trimming primer-only ends)."
        )
        return 0, warnings

    n_filled = 0
    if strand == "+":
        # ng_pos (1-based) = amplicon_pos (1-based) + offset
        for marker in panel.variants:
            if marker.amplicon_pos is not None and not force:
                continue
            if marker.ng006669_2_pos is None:
                continue
            amp_pos = marker.ng006669_2_pos - offset
            if 1 <= amp_pos <= len(new_ref_seq):
                marker.amplicon_pos = amp_pos
                n_filled += 1
    else:
        # New reference is the reverse complement of this NG_006669.2 span.
        # ng_pos = offset + (len(new_ref) - amplicon_pos + 1)  [1-based]
        for marker in panel.variants:
            if marker.amplicon_pos is not None and not force:
                continue
            if marker.ng006669_2_pos is None:
                continue
            amp_pos = len(new_ref_seq) - (marker.ng006669_2_pos - offset - 1)
            if 1 <= amp_pos <= len(new_ref_seq):
                marker.amplicon_pos = amp_pos
                n_filled += 1
        warnings.append(
            "Anchor matched on the REVERSE COMPLEMENT strand -- your "
            "combined reference appears to be oriented opposite to "
            "NG_006669.2. Positions were converted accordingly; spot-check "
            "a few ref_base values with --verify before trusting this."
        )

    return n_filled, warnings


def verify_ng006669_positions(panel: Panel, ng_seq: str) -> Tuple[int, int, List[str]]:
    """Sanity-check every row with ng006669_2_pos against the actual
    NG_006669.2 sequence: does the base there match ref_base? Returns
    (n_match, n_mismatch, mismatch_details)."""
    n_match, n_mismatch = 0, 0
    details = []
    for marker in panel.variants:
        ok = marker.ng_ref_base_matches(ng_seq)
        if ok is None:
            continue
        if ok:
            n_match += 1
        else:
            n_mismatch += 1
            idx = marker.ng006669_2_pos - 1
            actual = ng_seq[idx].upper() if 0 <= idx < len(ng_seq) else "?"
            details.append(
                f"{marker.id} ({marker.cdna_change}): expected ref_base="
                f"{marker.ref_base} at NG_006669.2:{marker.ng006669_2_pos}, "
                f"found '{actual}'"
            )
    return n_match, n_mismatch, details


def anchor_calibrate(
    panel: Panel,
    combined_ref_seq: str,
    legacy_exon6_seq: Optional[str],
    legacy_exon7_seq: Optional[str],
    force: bool = False,
) -> Tuple[int, List[str]]:
    """Shift legacy_exonN_pos coordinates onto the combined reference using
    substring anchoring. Returns (n_filled, warnings)."""
    warnings: List[str] = []
    n_filled = 0

    offsets = {}
    if legacy_exon6_seq:
        off = find_best_substring_offset(combined_ref_seq, legacy_exon6_seq)
        if off is None:
            warnings.append(
                "Could not locate the legacy exon6 reference sequence inside "
                "the combined reference (even allowing 5% mismatch). Exon6 "
                "positions were NOT calibrated by anchor mode."
            )
        else:
            offsets["Exon 6"] = off  # 0-based start of legacy seq in combined ref
    if legacy_exon7_seq:
        off = find_best_substring_offset(combined_ref_seq, legacy_exon7_seq)
        if off is None:
            warnings.append(
                "Could not locate the legacy exon7 reference sequence inside "
                "the combined reference (even allowing 5% mismatch). Exon7 "
                "positions were NOT calibrated by anchor mode."
            )
        else:
            offsets["Exon 7"] = off

    for marker in panel.variants:
        if marker.amplicon_pos is not None and not force:
            continue
        legacy_pos = None
        if marker.exon == "Exon 6" and marker.legacy_exon6_pos is not None:
            legacy_pos = marker.legacy_exon6_pos
        elif marker.exon == "Exon 7" and marker.legacy_exon7_pos is not None:
            legacy_pos = marker.legacy_exon7_pos
        if legacy_pos is None or marker.exon not in offsets:
            continue
        # legacy_pos is 1-based within the OLD mini-amplicon; offsets[exon]
        # is the 0-based start of that mini-amplicon sequence within the
        # new combined reference.
        marker.amplicon_pos = offsets[marker.exon] + legacy_pos
        n_filled += 1

    return n_filled, warnings


# ---------------------------------------------------------------------------
# SAM/CIGAR-based spliced-alignment calibration (no aligner dependency at
# parse time -- you supply the SAM produced by minimap2/your aligner of
# choice; this just walks the CIGAR).
# ---------------------------------------------------------------------------

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def build_cdna_to_ref_map(sam_path: Path) -> Dict[int, int]:
    """
    Parse a SAM file containing a spliced alignment of the ABO transcript
    (query, cDNA-numbered from its first base = c.1, matching ISBT
    numbering -- trim any 5' UTR from the transcript FASTA before aligning,
    or adjust with --cdna-offset) against the combined reference (subject).

    Returns {cdna_pos (1-based) : ref_pos (1-based)} for every transcript
    base that lands in an aligned (M/=/X) block. Bases in transcript
    insertions (I) or hard/soft clips have no reference position and are
    omitted. Bases the reference is missing entirely (large N gaps /
    deletions D) are also omitted -- if your positions of interest fall
    there, they are genuinely outside your amplicon/reference.
    """
    mapping: Dict[int, int] = {}
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            flag = int(fields[1])
            if flag & 0x4:  # unmapped
                continue
            if flag & 0x100 or flag & 0x800:  # secondary/supplementary -- skip, use primary only
                continue
            ref_start = int(fields[3])  # 1-based leftmost ref position
            cigar = fields[5]
            if cigar == "*":
                continue

            cdna_pos = 1  # walking position in the (unclipped) query/transcript
            ref_pos = ref_start

            for length_str, op in CIGAR_RE.findall(cigar):
                length = int(length_str)
                if op in ("M", "=", "X"):
                    for k in range(length):
                        mapping[cdna_pos + k] = ref_pos + k
                    cdna_pos += length
                    ref_pos += length
                elif op == "I":
                    cdna_pos += length  # consumes query, not reference
                elif op in ("D", "N"):
                    ref_pos += length  # consumes reference, not query
                elif op in ("S", "H"):
                    cdna_pos += length if op == "S" else 0
                # "P" (padding) consumes neither in practice here
            break  # only use the first (primary) alignment record found

    return mapping


def sam_calibrate(panel: Panel, sam_path: Path, cdna_offset: int = 0, force: bool = False) -> Tuple[int, List[str]]:
    """Fill amplicon_pos from a cDNA->reference map built from a SAM file.
    cdna_offset: add this to each panel cdna_pos before lookup, in case your
    transcript FASTA used for alignment didn't start exactly at c.1."""
    warnings: List[str] = []
    n_filled = 0
    cdna_map = build_cdna_to_ref_map(sam_path)
    if not cdna_map:
        warnings.append(f"No usable alignment records found in {sam_path}.")
        return 0, warnings

    for marker in panel.variants:
        if marker.amplicon_pos is not None and not force:
            continue
        if marker.cdna_pos is None:
            continue
        lookup_pos = marker.cdna_pos + cdna_offset
        ref_pos = cdna_map.get(lookup_pos)
        if ref_pos is None:
            warnings.append(
                f"{marker.id} (c.{marker.cdna_pos}, {marker.exon}) did not "
                f"land in an aligned block -- left uncalibrated. Verify it "
                f"is genuinely within your amplicon."
            )
            continue
        marker.amplicon_pos = ref_pos
        n_filled += 1

    return n_filled, warnings


# ---------------------------------------------------------------------------
# Writing the updated panel back out
# ---------------------------------------------------------------------------

def write_panel_yaml(panel: Panel, out_path: Path) -> None:
    if yaml is None:
        raise RuntimeError("pyyaml required to write a YAML panel")
    data = {
        "meta": panel.meta,
        "variants": [{k: getattr(v, k) for k in FIELD_NAMES} for v in panel.variants],
        "structural_variants": [
            {k: getattr(sv, k) for k in
             ["id", "cdna_change", "exon", "size_bp_approx", "associated_alleles", "notes", "source"]}
            for sv in panel.structural_variants
        ],
    }
    with open(out_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--panel", required=True, help="Panel file to calibrate (.yaml/.json/.csv/.tsv)")
    ap.add_argument("--out", required=True, help="Path to write the calibrated panel (.yaml or .csv)")
    ap.add_argument("--combined-reference", help="Your new combined exon2-7 reference FASTA (required for --legacy-exon6-reference/--legacy-exon7-reference anchor mode)")
    ap.add_argument("--legacy-exon6-reference", help="Old v1.x exon6-only mini-amplicon reference FASTA")
    ap.add_argument("--legacy-exon7-reference", help="Old v1.x exon7-only mini-amplicon reference FASTA")
    ap.add_argument("--from-sam", help="SAM file: ABO transcript (query) aligned with a SPLICED aligner "
                                        "(e.g. minimap2 -ax splice) against your combined reference (subject)")
    ap.add_argument("--ng006669-reference", help="Full NG_006669.2 FASTA (RefSeqGene). Combined with "
                                                   "--combined-reference, calibrates amplicon_pos directly "
                                                   "from each row's ng006669_2_pos -- no legacy exon6/7 "
                                                   "references or SAM alignment needed for any row that "
                                                   "already has ng006669_2_pos populated.")
    ap.add_argument("--ng-anchor-len", type=int, default=300,
                     help="Bases of --combined-reference's start used as the anchor when locating it "
                          "inside --ng006669-reference (default: %(default)s)")
    ap.add_argument("--verify", action="store_true",
                     help="Check every row's ng006669_2_pos against --ng006669-reference's actual "
                          "sequence (ref_base match) and report mismatches, without changing anything.")
    ap.add_argument("--cdna-offset", type=int, default=0,
                     help="Add this to every panel cdna_pos before looking it up in --from-sam's map "
                          "(use if your transcript FASTA's first base isn't cDNA position 1)")
    ap.add_argument("--force", action="store_true",
                     help="Recompute amplicon_pos even for rows that already have one")
    args = ap.parse_args()

    panel = load_panel(args.panel)
    total_warnings: List[str] = []
    total_filled = 0

    if args.legacy_exon6_reference or args.legacy_exon7_reference:
        if not args.combined_reference:
            print("! --combined-reference is required for anchor-mode calibration")
            sys.exit(1)
        combined_seqs = read_fasta(Path(args.combined_reference))
        combined_seq = next(iter(combined_seqs.values()))
        exon6_seq = None
        exon7_seq = None
        if args.legacy_exon6_reference:
            exon6_seq = next(iter(read_fasta(Path(args.legacy_exon6_reference)).values()))
        if args.legacy_exon7_reference:
            exon7_seq = next(iter(read_fasta(Path(args.legacy_exon7_reference)).values()))

        n, warns = anchor_calibrate(panel, combined_seq, exon6_seq, exon7_seq, force=args.force)
        print(f"[anchor mode] Calibrated {n} positions from legacy coordinates.")
        total_filled += n
        total_warnings.extend(warns)

    if args.from_sam:
        n, warns = sam_calibrate(panel, Path(args.from_sam), cdna_offset=args.cdna_offset, force=args.force)
        print(f"[SAM mode] Calibrated {n} positions from spliced alignment.")
        total_filled += n
        total_warnings.extend(warns)

    if args.ng006669_reference and args.verify:
        ng_seq = next(iter(read_fasta(Path(args.ng006669_reference)).values()))
        n_match, n_mismatch, details = verify_ng006669_positions(panel, ng_seq)
        print(f"\n[verify] ng006669_2_pos ref_base check: {n_match} match, {n_mismatch} mismatch")
        for d in details:
            print(f"  MISMATCH: {d}")
        if not (args.legacy_exon6_reference or args.legacy_exon7_reference or args.from_sam
                or (args.ng006669_reference and args.combined_reference)):
            out_path = Path(args.out)
            if out_path.suffix.lower() in (".csv", ".tsv"):
                export_csv(panel, out_path)
            else:
                write_panel_yaml(panel, out_path)
            print(f"\nWrote unmodified panel (verify-only run) -> {out_path}")
            return

    if args.ng006669_reference and args.combined_reference:
        ng_seq = next(iter(read_fasta(Path(args.ng006669_reference)).values()))
        combined_seq = next(iter(read_fasta(Path(args.combined_reference)).values()))
        n, warns = ng006669_calibrate(panel, ng_seq, combined_seq,
                                       anchor_len=args.ng_anchor_len, force=args.force)
        print(f"[NG_006669.2 mode] Calibrated {n} positions via direct genomic-coordinate offset.")
        total_filled += n
        total_warnings.extend(warns)

    ran_anything = (
        args.legacy_exon6_reference or args.legacy_exon7_reference or args.from_sam
        or (args.ng006669_reference and args.combined_reference)
    )
    if not ran_anything:
        print("Nothing to do: pass --legacy-exon6-reference/--legacy-exon7-reference "
              "(with --combined-reference), --from-sam, and/or "
              "--ng006669-reference (with --combined-reference).")
        sys.exit(1)

    out_path = Path(args.out)
    if out_path.suffix.lower() in (".csv", ".tsv"):
        export_csv(panel, out_path)
    else:
        write_panel_yaml(panel, out_path)

    print(f"\nTotal newly calibrated: {total_filled}")
    print(f"Still uncalibrated: {len(panel.uncalibrated())}")
    if total_warnings:
        print("\nWarnings:")
        for w in total_warnings:
            print(f"  - {w}")
    print(f"\nWrote calibrated panel -> {out_path}")


if __name__ == "__main__":
    main()
