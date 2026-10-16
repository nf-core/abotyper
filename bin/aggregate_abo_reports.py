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
__version__ = "2.4.0"
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

Typing sheet and ISBT allele names (v2.2.0)
-------------------------------------------
ABO_result.xlsx now has two sheets: "Typing" (one NGSengine-style row per
sample: reads, read length, insert size, mappability, median depth, median
analysed depth, number of unresolved ("?") primary positions, typing result,
phenotype, genotype, QC) and "Variants" (the per-position table, formerly
"ABO_Result"). ExtendedGenotype is written with ISBT-style names (A1 ->
ABO*A1.01, O1 -> ABO*O.01.01, B -> ABO*B.01 are only reference defaults; when
non-primary panel markers fire, the allele is narrowed to the one they point
to). The decision tree and QC logic still work on the internal tokens and only the written value is encoded.
Alignment metrics are read from samtools stats / flagstat files found under
each sample directory; a missing file leaves the cell blank.

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
# ISBT allele names and alignment metrics (v2.2.0)
# ===========================================================================

# Reference (default) ISBT allele per internal token: used only when no panel marker
# points to a more specific allele of the same family. Checked against the ISBT ABO allele
# table v9 (blooddatabase.isbtweb.org/system/ABO): ABO*A1.01 (reference, no changes),
# ABO*B.01 (B backbone), ABO*O.01.01 (c.261delG only), ABO*O.02.01 (c.53G>T, c.220C>T,
# c.297A>G, c.526C>G, c.802G>A) and ABO*O.03 (no numeric suffix in the table).
RESEARCH_USE_STATEMENT = (
    "RESEARCH USE ONLY. Not validated for, and not to be used for, clinical diagnosis or "
    "transfusion decisions; confirm results with an accredited method before any clinical use."
)
ISBT_ALLELE_NAMES = {
    "A1": "ABO*A1.01",
    "A2.01": "ABO*A2.01",
    "A2.06": "ABO*A2.06",
    "A2": "ABO*A2",
    "A": "ABO*A",
    "B": "ABO*B.01",
    "O1": "ABO*O.01.01",
    "O2": "ABO*O.02.01",
    "O3": "ABO*O.03",
}

# Family prefix an allele name must start with to be a candidate for a token.
_ISBT_FAMILY = {
    "A1": "ABO*A1.", "A2.01": "ABO*A2.", "A2.06": "ABO*A2.", "A2": "ABO*A2.",
    "B": "ABO*B.", "O1": "ABO*O.01.", "O2": "ABO*O.02.", "O3": "ABO*O.03",
}
_ISBT_NAME_RE = re.compile(r"(?:ABO\*)?([A-Za-z]+[0-9]*\.\d+(?:\.\d+)?)(?![\d-])")


def _marker_candidate_alleles(marker) -> set:
    """Exact (non-wildcard) ISBT allele names listed in a panel marker's associated_alleles."""
    out = set()
    for m in _ISBT_NAME_RE.finditer(getattr(marker, "associated_alleles", "") or ""):
        out.add("ABO*" + m.group(1).upper())
    return out


def encode_allele(token: str, fired_markers=()) -> str:
    """Internal allele token (A1, O2, ...) -> ISBT-style name.

    Starts from the reference allele of the token's family, then narrows using the
    non-primary panel markers that actually fired in this sample: the allele names
    those markers are associated with, restricted to the token's family. One surviving
    candidate -> that allele; several -> the family name only (not resolved further);
    none -> the reference allele. Unknown tokens pass through.
    """
    t = str(token).strip()
    default = ISBT_ALLELE_NAMES.get(t, t)
    fam = _ISBT_FAMILY.get(t)
    if not fam:
        return default
    sets = []
    for mk in fired_markers:
        cand = {a for a in _marker_candidate_alleles(mk) if a.startswith(fam) or a == fam.rstrip(".")}
        if cand:
            sets.append(cand)
    if not sets:
        return default
    common = set.intersection(*sets)
    cands = common if common else set.union(*sets)
    if len(cands) == 1:
        return next(iter(cands))
    return fam.rstrip(".")   # ambiguous: family only, e.g. ABO*O.01


def encode_extended_genotype(ext: str, fired_markers=()) -> str:
    """'A2.01/O1' -> 'ABO*A2.01/ABO*O.01.01'. 'Unknown', 'Error' and empty values are left untouched."""
    if not isinstance(ext, str) or ext in ("", "Unknown", "Error"):
        return ext
    tokens = ext.split("/")
    markers = list(fired_markers)
    if "B" in tokens:
        # B-lineage backbone SNPs (associated alleles include ABO*B.*) come from the B
        # allele of an AB/BO sample; they must not rename the non-B allele (e.g. to A2.04)
        backbone = [m for m in markers if "ABO*B.*" in (getattr(m, "associated_alleles", "") or "")]
        non_b = [m for m in markers if m not in backbone]
        return "/".join(encode_allele(t, markers if t == "B" else non_b) for t in tokens)
    return "/".join(encode_allele(t, markers) for t in tokens)


def _first_file(sample_dir: str, patterns) -> Optional[str]:
    """First non-empty file under sample_dir matching any glob pattern (recursive).
    Prefers a combined/ location, then the shallowest path, then alphabetical."""
    found = []
    for pat in patterns:
        found.extend(glob.glob(os.path.join(sample_dir, "**", pat), recursive=True))
    found = [f for f in found if os.path.isfile(f) and os.path.getsize(f) > 0]
    if not found:
        return None

    def rank(path):
        rel = os.path.relpath(path, sample_dir)
        parts = rel.replace("\\", "/").split("/")
        return (0 if "combined" in parts[:-1] else 1, len(parts), rel)

    return sorted(found, key=rank)[0]


def _weighted_median(pairs) -> Optional[float]:
    """Median of values given as (value, count) pairs."""
    pairs = sorted((v, c) for v, c in pairs if c > 0)
    total = sum(c for _, c in pairs)
    if total <= 0:
        return None
    half = total / 2.0
    running = 0
    for i, (v, c) in enumerate(pairs):
        running += c
        if running > half:
            return float(v)
        if running == half:
            nxt = pairs[i + 1][0] if i + 1 < len(pairs) else v
            return (v + nxt) / 2.0
    return float(pairs[-1][0])


def parse_alignment_metrics(sample_dir: str) -> Dict[str, Optional[float]]:
    """Per-sample alignment metrics from samtools stats (preferred) and flagstat (fallback).

    Keys: reads, read_length, insert_size, mappability, median_depth (None when unavailable).
    Never raises. Note samtools stats caps the COV histogram at its last bin
    (default 1000+), so a median depth above the cap is reported as the cap.
    """
    m: Dict[str, Optional[float]] = {
        "reads": None, "read_length": None, "insert_size": None,
        "mappability": None, "median_depth": None,
    }
    mapped = None
    try:
        stats_file = _first_file(sample_dir, ["*.stats"])
        if stats_file:
            cov = []
            with open(stats_file, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith("SN\t"):
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3:
                            continue
                        key, val = parts[1].rstrip(":").strip(), parts[2].strip()
                        try:
                            if key == "raw total sequences":
                                m["reads"] = float(val)
                            elif key == "reads mapped":
                                mapped = float(val)
                            elif key == "average length":
                                m["read_length"] = float(val)
                            elif key == "insert size average":
                                m["insert_size"] = float(val)
                        except ValueError:
                            pass
                    elif line.startswith("COV\t"):
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 4:
                            try:
                                cov.append((int(parts[2]), int(parts[3])))
                            except ValueError:
                                pass
            m["median_depth"] = _weighted_median(cov)
    except Exception:
        pass

    try:
        if m["reads"] is None or mapped is None:
            flag_file = _first_file(sample_dir, ["*.flagstat"])
            if flag_file:
                with open(flag_file, "r", encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        t = re.match(r"^(\d+) \+ \d+ in total", line)
                        if t and m["reads"] is None:
                            m["reads"] = float(t.group(1))
                        g = re.match(r"^(\d+) \+ \d+ mapped \(", line)
                        if g and mapped is None:
                            mapped = float(g.group(1))
    except Exception:
        pass

    m["mapped"] = mapped
    if m["reads"] and mapped is not None:
        m["mappability"] = 100.0 * mapped / m["reads"]
    if m["insert_size"] is not None and m["insert_size"] <= 0:
        m["insert_size"] = None   # 0 = unpaired/long reads: not meaningful
    return m


def parse_coverage_file(sample_dir: str) -> Dict[str, Optional[float]]:
    """ABO-specific figures from a samtools coverage table (*.coverage.txt, one row per contig).

    abo_reads     reads mapped to the ABO reference (numreads)
    covered_depth mean depth over the bases that have any coverage (meandepth * length / covbases)
    Both None when the file is absent. In an ABO-HLA duplex run the HLA reads do not map to the
    ABO reference, so these describe the ABO amplicon rather than the whole library.
    """
    out: Dict[str, Optional[float]] = {"abo_reads": None, "covered_depth": None}
    try:
        f = _first_file(sample_dir, ["*.coverage.txt"])
        if not f:
            return out
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]
        reads = depth_x_len = cov = 0.0
        for r in rows:
            if len(r) < 7:
                continue
            start, end, numreads, covbases, meandepth = (float(r[1]), float(r[2]), float(r[3]),
                                                         float(r[4]), float(r[6]))
            reads += numreads
            cov += covbases
            depth_x_len += meandepth * (end - start + 1)
        out["abo_reads"] = reads
        if cov > 0:
            out["covered_depth"] = depth_x_len / cov
    except (OSError, ValueError):
        pass
    return out


def parse_evidence_status(sample_dir: str) -> Dict[str, object]:
    """Two-line evidence summary written by the pipeline's ABO_EVIDENCE step (*.evidence_status.tsv).

    evidence_source  CLAIR3+BAM | CLAIR3_ONLY | BAM_ONLY | NONE ('' when no file was found)
    filled_primary   primary typing positions taken from the BAM pileup (Clair3 had no evidence)
    discordant       panel positions where Clair3 and the BAM pileup disagree
    primary_names    the primary typing positions among them
    """
    out: Dict[str, object] = {"evidence_source": "", "filled_primary": None, "bam_only": None,
                              "discordant": None, "discordant_primary": None, "primary_names": "",
                              "bam_max_depth": None}
    try:
        f = _first_file(sample_dir, ["*.evidence_status.tsv"])
        if not f:
            return out
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
        if len(lines) >= 2:
            row = dict(zip(lines[0], lines[1]))
            out["evidence_source"] = row.get("evidence_source", "")
            out["filled_primary"] = int(row.get("filled_primary_from_bam", "0"))
            out["bam_only"] = int(row.get("bam_only", "0"))
            if row.get("bam_max_depth", "") != "":
                out["bam_max_depth"] = int(row["bam_max_depth"])
            out["discordant"] = int(row.get("discordant", "0"))
            out["discordant_primary"] = int(row.get("discordant_primary", "0"))
            out["primary_names"] = row.get("discordant_primary_names", "").replace(";", ", ")
    except (OSError, ValueError):
        pass
    return out


def parse_nanostats(path: str) -> Dict[str, Optional[float]]:
    """Read count and mean read length from a NanoPlot NanoStats.txt (None when absent)."""
    out: Dict[str, Optional[float]] = {"reads": None, "read_length": None}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                key, _, val = line.partition(":")
                val = val.strip().replace(",", "")
                try:
                    if key.strip() == "Number of reads":
                        out["reads"] = float(val)
                    elif key.strip() == "Mean read length":
                        out["read_length"] = float(val)
                except ValueError:
                    pass
    except OSError:
        pass
    return out


def find_nanostats(sample_dir: str, nanoplot_dir: Optional[str], sample_name: str,
                   barcode: str) -> Optional[str]:
    """NanoStats.txt inside the sample folder, else under nanoplot_dir in a path that
    contains the sample name (preferred) or the barcode token."""
    found = _first_file(sample_dir, ["*NanoStats.txt"])
    if found or not nanoplot_dir or not os.path.isdir(nanoplot_dir):
        return found
    cands = [f for f in glob.glob(os.path.join(nanoplot_dir, "**", "*NanoStats.txt"), recursive=True)
             if os.path.getsize(f) > 0]
    for token in (sample_name, barcode):
        hits = sorted(f for f in cands if token and token.lower() in
                      os.path.relpath(f, nanoplot_dir).lower())
        if token == barcode and len(hits) > 1:
            continue   # barcode alone is ambiguous across runs
        if hits:
            return hits[0]
    return None


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
        df = pd.read_csv(haplo_file, sep="\t", encoding="utf-8")
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


_WORKER_PARSER = None


def _init_worker(parser):
    """Pool initializer: each worker process keeps its own copy of the configured parser."""
    global _WORKER_PARSER
    _WORKER_PARSER = parser


def _process_one(filename):
    """Run process_file for one sample in a worker; return everything it recorded plus its log."""
    import contextlib
    import io
    p = _WORKER_PARSER
    p.results, p.qc_records, p.sample_metrics, p.failed_samples = [], [], {}, []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            p.process_file(filename)
        except Exception as e:    # process_file handles its own errors; this is a backstop
            p.failed_samples.append({"sample": filename, "reason": f"Processing error: {e}"})
    return filename, buf.getvalue(), p.results, p.qc_records, p.sample_metrics, p.failed_samples


class ABOReportParser:
    """
    Collates all ABO phenotype results from each sample into a general table
    and generates an Excel worksheet and a CSV file for export to LIS soft
    or other general purpose lab management systems. Column layout, position
    parsing, and subtype marker scanning are driven by an external variant
    panel (see abo_panel.py).
    """

    def __init__(self, input_dir, panel: Panel, default_barcode="barcode00",
                 min_depth=30, qc_report="ABO_run_qc.txt", nanoplot_dir=None, threads=1):
        self.input_dir = input_dir
        self.panel = panel
        self.default_barcode = default_barcode
        self.results = []
        self.failed_samples = []
        self.min_depth = min_depth            # reads below which an A1 call is flagged
        self.qc_report_path = qc_report       # run-level QC text report ("" disables)
        self.qc_records = []                  # one dict per sample, for the run summary
        self.threads = max(1, int(threads or 1))   # worker processes for per-sample parsing
        self.nanoplot_dir = nanoplot_dir      # optional folder of NanoPlot output (NanoStats.txt)
        self._fired_markers = []
        self.sample_metrics = {}             # (sample_name, int barcode) -> Typing-sheet metrics
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
                              "A2Trace", "EvidenceBasis", "QCStatus", "QCFlags", "ReplicateNote"]
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
                self._fired_markers.append(marker)
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
            self._fired_markers = []   # non-primary markers that fired; refine ISBT allele names

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
                    if a297_marker:
                        self._fired_markers.append(a297_marker)   # points at the O1v-type allele
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
                                self._fired_markers.append(m_marker)

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
                            self._fired_markers.append(m_)

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
            df[("Result", "ExtendedGenotype")] = encode_extended_genotype(ExtendedGenotype, self._fired_markers)
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

            self._record_sample_metrics(result_df, sample_name, barcode, sample_dir)

            self.results.append(result_df)
            print(f"Successfully processed {filename}")

        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            import traceback
            traceback.print_exc()

    def _record_sample_metrics(self, result_df, sample_name, barcode, sample_dir):
        """Collect the Typing-sheet metrics for one sample (never raises)."""
        try:
            metrics = parse_alignment_metrics(sample_dir)
        except Exception:
            metrics = {}
        # Duplex (ABO-HLA) runs: HLA reads never map to the ABO reference, so whole-library
        # mappability is misleading. Report reads in the library (NanoPlot, else samtools),
        # reads on the ABO reference (samtools coverage) and on-target % between them.
        try:
            metrics.update(parse_coverage_file(sample_dir))
        except Exception:
            pass
        try:
            status = parse_evidence_status(sample_dir)
            metrics.update({"evidence_source": status["evidence_source"], "evidence_label": None,
                            "bam_only": status["bam_only"], "discordant": status["discordant"],
                            "discordant_primary": status["discordant_primary"]})
            msgs = []
            label = "OK"
            if status["evidence_source"] == "NONE":
                mx = status["bam_max_depth"]
                if mx == 0:
                    label = "No reads over ABO region (likely amplicon dropout)"
                elif mx is not None:
                    label = f"ABO region barely covered (max pileup depth {mx})"
                else:
                    label = "No evidence from either Clair3 or BAM pileup"
                # replace the generic per-position message with the specific cause
                old = str(result_df.at[0, ("Notes", "QCFlags")] or "")
                kept = [f for f in old.split("; ") if f and not f.startswith("no evidence at primary position")]
                result_df.at[0, ("Notes", "QCFlags")] = "; ".join([label] + kept)
                result_df.at[0, ("Notes", "QCStatus")] = "REVIEW"
            elif status["evidence_source"] == "BAM_ONLY":
                label = "Typed from BAM pileup only (Clair3 produced no usable calls)"
            if status["evidence_source"] == "BAM_ONLY":
                msgs.append("Clair3 produced no usable calls; typed from BAM pileup only")
            elif status["filled_primary"]:
                msgs.append(f"{status['filled_primary']} primary position(s) rest on BAM pileup "
                            f"alone (no Clair3 evidence)")
            if status["discordant_primary"]:
                msgs.append(f"Clair3 and BAM pileup disagree at primary position(s): "
                            f"{status['primary_names']}")
            metrics["evidence_label"] = label
            if msgs:
                # the two lines of evidence do not support the call equally: always review
                old = str(result_df.at[0, ("Notes", "QCFlags")] or "")
                result_df.at[0, ("Notes", "QCFlags")] = (old + "; " if old else "") + "; ".join(msgs)
                result_df.at[0, ("Notes", "QCStatus")] = "REVIEW"
        except Exception:
            pass
        try:
            # NanoPlot runs on the raw reads, so it is the better source for read count/length
            ns_file = find_nanostats(sample_dir, self.nanoplot_dir, str(sample_name), str(barcode))
            if ns_file:
                ns = parse_nanostats(ns_file)
                for k in ("reads", "read_length"):
                    if ns.get(k) is not None:
                        metrics[k] = ns[k]
        except Exception:
            pass
        try:
            # the coverage table counts secondary/supplementary alignments; samtools stats'
            # primary mapped-read count is the ceiling when we have it
            if metrics.get("abo_reads") is not None and metrics.get("mapped") is not None:
                metrics["abo_reads"] = min(metrics["abo_reads"], metrics["mapped"])
            if metrics.get("abo_reads") is not None and metrics.get("reads"):
                metrics["mappability"] = min(100.0, 100.0 * metrics["abo_reads"] / metrics["reads"])
        except Exception:
            pass
        try:
            depths = []
            for col_label in self.position_columns:
                v = pd.to_numeric(result_df.at[0, (col_label, "#Reads")], errors="coerce")
                if pd.notna(v) and v > 0:
                    depths.append(float(v))
            metrics["median_analyzed_depth"] = float(pd.Series(depths).median()) if depths else None

            qmarks = 0
            for marker in self.primary_markers:
                t = result_df.at[0, (marker.column_label(), "Type")]
                if pd.isna(t) or not str(t).strip():
                    qmarks += 1
            metrics["question_marks"] = qmarks
        except Exception:
            metrics.setdefault("median_analyzed_depth", None)
            metrics.setdefault("question_marks", None)
        try:
            key = (str(sample_name), int(str(barcode).replace("barcode", "")))
        except ValueError:
            return
        self.sample_metrics[key] = metrics

    def _build_typing_df(self, final_df) -> pd.DataFrame:
        """One NGSengine-style row per sample, same order as final_df (single-level columns)."""
        def num(v, nd=None):
            if v is None or pd.isna(v):
                return None
            return round(float(v), nd) if nd is not None else int(round(float(v)))

        rows = []
        for _, r in final_df.iterrows():
            sid = str(r["Sequencing_ID"])
            try:
                bc = int(r["Barcode"])
            except (TypeError, ValueError):
                bc = None
            m = self.sample_metrics.get((sid, bc), {})
            rows.append({
                "Sample name": sid,
                "Barcode": f"barcode{bc:02d}" if bc is not None else "",
                "Number of reads": num(m.get("reads")),
                "Read Length": num(m.get("read_length")),
                "Insert Size": num(m.get("insert_size")),
                "ABO reads": num(m.get("abo_reads")),
                "Mappability %": num(m.get("mappability"), 1),
                "Mean depth (covered bases)": num(m.get("covered_depth")),
                "Median Read depth": num(m.get("median_depth")),
                "Median analyzed Read depth": num(m.get("median_analyzed_depth")),
                "Number of question mark positions": m.get("question_marks"),
                "Evidence source": m.get("evidence_source") or "",
                "BAM-only positions": m.get("bam_only"),
                "Clair3/BAM discordant": m.get("discordant"),
                "Discordant primary positions": m.get("discordant_primary"),
                "Typing Result": r["ExtendedGenotype"],
                "Phenotype": r["Phenotype"],
                "Genotype": r["Genotype"],
                "QCStatus": r["QCStatus"],
                "QCFlags": r["QCFlags"],
                "Replicate note": r.get("ReplicateNote"),
            })
        return pd.DataFrame(rows, columns=[
            "Sample name", "Barcode", "Number of reads", "Read Length", "Insert Size",
            "ABO reads", "Mappability %", "Mean depth (covered bases)", "Median Read depth",
            "Median analyzed Read depth",
            "Number of question mark positions", "Evidence source", "BAM-only positions",
            "Clair3/BAM discordant", "Discordant primary positions", "Typing Result", "Phenotype", "Genotype",
            "QCStatus", "QCFlags", "Replicate note"]).astype(object)

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

        parallel = self.threads > 1 and len(valid_directories) > 1
        if parallel:
            self._process_files_parallel(valid_directories)

        for filename in ([] if parallel else valid_directories):
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

    def _process_files_parallel(self, valid_directories):
        """Parse samples in a process pool (--threads N). Output order does not matter: the
        merged table is sorted by sample and barcode, and each sample's log is printed whole,
        in directory order, once its worker finishes."""
        from concurrent.futures import ProcessPoolExecutor

        for filename in valid_directories:        # cheap bookkeeping stays serial
            self.processing_stats["pattern_matched"] += 1
            try:
                _, _, pattern_type = self.extract_sample_info(filename)
                self.processing_stats["with_barcode" if pattern_type == "explicit" else "without_barcode"] += 1
            except Exception:
                pass

        workers = min(self.threads, len(valid_directories))
        print(f"Parsing {len(valid_directories)} samples with {workers} worker processes")
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(self,)) as ex:
            for filename, log, results, qc_records, sample_metrics, failed in ex.map(
                    _process_one, valid_directories):
                print(f"\nProcessing file: {filename}")
                print(log, end="")
                print(f"Finished processing file: {filename}")
                self.results.extend(results)
                self.qc_records.extend(qc_records)
                self.sample_metrics.update(sample_metrics)
                self.failed_samples.extend(failed)

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
            final_df.to_csv("./ABO_result.txt", sep="\t", index=False, encoding="utf-8")
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

            workbook = writer.book
            data_format = workbook.add_format({"bg_color": "white", "font_color": "black", "border": 1})
            header_format = workbook.add_format({"bold": True, "fg_color": "#007399", "border": 1, "font_color": "white"})
            red_bg_format = workbook.add_format({"bg_color": "#e2725b", "font_color": "black"})
            orange_bg_format = workbook.add_format({"bg_color": "#ff9a00", "font_color": "black"})
            header_format.set_align("center")
            header_format.set_align("vcenter")

            # Sheet 1: Typing (one NGSengine-style row per sample). Sheet order = creation order.
            try:
                typing_df = self._build_typing_df(final_df)
                typing_df.to_excel(writer, sheet_name="Typing", index=False)
                ws_t = writer.sheets["Typing"]
                for c, name in enumerate(typing_df.columns):
                    ws_t.write(0, c, name, header_format)
                int_fmt = workbook.add_format({"border": 1, "num_format": "#,##0"})
                dec_fmt = workbook.add_format({"border": 1, "num_format": "0.0"})
                for r in range(len(typing_df)):
                    for c, name in enumerate(typing_df.columns):
                        v = typing_df.iat[r, c]
                        if v is None or (not isinstance(v, str) and pd.isna(v)):
                            ws_t.write_blank(r + 1, c, None, data_format)
                        elif name == "Mappability %":
                            ws_t.write_number(r + 1, c, float(v), dec_fmt)
                        elif name in ("Number of reads", "Read Length", "Insert Size", "Median Read depth",
                                      "ABO reads", "Mean depth (covered bases)",
                                      "Median analyzed Read depth", "Number of question mark positions",
                                      "BAM-only positions", "Clair3/BAM discordant",
                                      "Discordant primary positions"):
                            ws_t.write_number(r + 1, c, float(v), int_fmt)
                        else:
                            ws_t.write(r + 1, c, v, data_format)
                for c, name in enumerate(typing_df.columns):
                    longest = max([len(str(name))] + [len(str(x)) for x in typing_df[name] if x is not None])
                    ws_t.set_column(c, c, min(max(longest + 2, 10), 80 if name == "Replicate note" else 60 if name == "QCFlags" else 40))
                ws_t.set_row(0, 45, header_format)
                ws_t.freeze_panes(1, 0)
                last = len(typing_df)
                if last:
                    qc_t = xl_col_to_name(list(typing_df.columns).index("QCStatus"))
                    ws_t.conditional_format(
                        f"{qc_t}2:{qc_t}{last + 1}",
                        {"type": "cell", "criteria": "==", "value": '"REVIEW"', "format": orange_bg_format},
                    )
                    rn_t = xl_col_to_name(list(typing_df.columns).index("Replicate note"))
                    ws_t.conditional_format(
                        f"A2:{xl_col_to_name(len(typing_df.columns) - 1)}{last + 1}",
                        {"type": "formula", "criteria": f"=LEN(${rn_t}2)>0", "format": red_bg_format},
                    )
                    rd_t = xl_col_to_name(list(typing_df.columns).index("Number of reads"))
                    ws_t.conditional_format(
                        f"{rd_t}2:{rd_t}{last + 1}",
                        {"type": "cell", "criteria": "<=", "value": self.min_depth, "format": orange_bg_format},
                    )
            except Exception as typing_err:
                print(f"Warning: could not write Typing sheet: {typing_err}")

            # Sheet 2: Variants (per-position table)
            final_df.to_excel(writer, sheet_name="Variants", header=True, index=False, startrow=1)
            worksheet = writer.sheets["Variants"]

            about = workbook.add_worksheet("About")
            bold = workbook.add_format({"bold": True})
            about.set_column(0, 0, 34)
            about.set_column(1, 1, 110)
            about.write(0, 0, "ABO typing results", bold)
            about.write(1, 0, RESEARCH_USE_STATEMENT, bold)
            for r_, (k_, v_) in enumerate(self._about_items(), start=3):
                about.write(r_, 0, k_, bold)
                about.write(r_, 1, v_)

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

            if "ReplicateNote" in _cols and num_rows:
                rn_letter = xl_col_to_name(_cols.index("ReplicateNote"))
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",
                    {"type": "formula", "criteria": f"=LEN(${rn_letter}3)>0", "format": red_bg_format},
                )

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
                           ("A2Trace", "A2 Trace"), ("EvidenceBasis", "Evidence Basis"),
                           ("ReplicateNote", "Replicate Note")):
            if _src in final_df.columns:
                self.df_for_lis_soft[_dst] = final_df[_src]

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")
        print(f"LIS export file created successfully with {len(self.df_for_lis_soft)} samples")

    # -----------------------------------------------------------------
    # QC: replicate cross-check and run-level summary (v2.1.0)
    # -----------------------------------------------------------------
    @staticmethod
    def _zygosity(genotype) -> str:
        g = str(genotype)
        if len(g) != 2 or g in ("Unknown", "nan"):
            return "unresolved"
        return "homozygous" if g[0] == g[1] else "heterozygous"

    def _apply_replicate_checks(self, final_df):
        """Samples sharing a Sequencing_ID (the same sample on different barcodes/runs) must
        agree on genotype, A2 evidence and allele-level typing. Disagreement marks EVERY member
        of the group: REVIEW + a QC flag, and a ReplicateNote that lists each replicate's
        genotype, zygosity and typing (the sheets colour these rows red). Nothing is changed.
        Rows are already sorted by sample then barcode, so replicates sit together."""
        sid = final_df[("", "Sequencing_ID")].astype(str)
        stats = {"groups": 0, "discordant_genotype": 0, "discordant_a2": 0, "discordant_allele": 0}
        for _, idx in sid.groupby(sid).groups.items():
            idx = list(idx)
            if len(idx) < 2:
                continue
            stats["groups"] += 1
            grp = final_df.loc[idx]
            genos = set(grp[("Result", "Genotype")].astype(str))
            exts = set(grp[("Result", "ExtendedGenotype")].astype(str))
            a2 = grp[("Notes", "A2Trace")].astype(str)
            a2_has = set(not x.startswith("4") for x in a2 if "n/a" not in x and not x.startswith("5"))
            rows = list(zip(grp[("", "Barcode")], grp[("Result", "Genotype")].astype(str),
                            grp[("Result", "ExtendedGenotype")].astype(str)))
            members = ", ".join(f"barcode{int(bc):02d}={e}" for bc, _, e in rows)
            msg = None
            if len(genos) > 1:
                stats["discordant_genotype"] += 1
                msg = f"replicates disagree on genotype ({members})"
            elif len(a2_has) > 1:
                stats["discordant_a2"] += 1
                msg = f"replicates disagree on A2 evidence ({members})"
            elif len(exts) > 1:
                stats["discordant_allele"] += 1
                msg = f"replicates agree on genotype but differ at allele level ({members})"
            if not msg:
                continue

            # ReplicateNote: every replicate with its zygosity and typing, plus what the
            # mismatch most likely means
            parts = [f"barcode{int(bc):02d}: {g} ({self._zygosity(g)}), {e}" for bc, g, e in rows]
            note = f"Replicate mismatch ({len(rows)} runs of this sample): " + " | ".join(parts) + "."
            hets = [(bc, g) for bc, g, _ in rows if self._zygosity(g) == "heterozygous"]
            homs = [(bc, g) for bc, g, _ in rows if self._zygosity(g) == "homozygous"]
            dropout = [f"barcode{int(bc):02d}" for bc, g in homs
                       if any(set(g) <= set(hg) for _, hg in hets)]
            if dropout:
                note += (f" Homozygous call in {', '.join(dropout)} is contained in a heterozygous replicate: "
                         "possible allele dropout; confirm with another method.")
            unresolved = [f"barcode{int(bc):02d}" for bc, g, _ in rows if self._zygosity(g) == "unresolved"]
            if unresolved:
                note += f" Unresolved (no usable typing): {', '.join(unresolved)}."
            for i in idx:
                old = str(final_df.at[i, ("Notes", "QCFlags")] or "")
                final_df.at[i, ("Notes", "QCFlags")] = (old + "; " if old else "") + msg
                final_df.at[i, ("Notes", "QCStatus")] = "REVIEW"
                final_df.at[i, ("Notes", "ReplicateNote")] = note
        self.replicate_stats = stats

    def _about_items(self):
        """Version / provenance lines for the About sheet and the QC report."""
        meta = getattr(self.panel, "meta", {}) or {}
        items = [("aggregate_abo_reports.py", __version__),
                 ("Panel version", str(meta.get("panel_version", "unknown"))),
                 ("Panel source", str(meta.get("source", "unknown")))]
        if meta.get("isbt_allele_table_version"):
            items.append(("ISBT ABO allele table", str(meta["isbt_allele_table_version"])))
        if meta.get("isbt_allele_table_checked"):
            items.append(("ISBT table checked", str(meta["isbt_allele_table_checked"])))
        if meta.get("isbt_allele_table_release_date"):
            items.append(("ISBT table release date", str(meta["isbt_allele_table_release_date"])))
        return items

    def _evidence_summary(self, final_df):
        """Per-sample evidence rows (None when no sample carries an evidence status)."""
        rows = []
        for _, r in final_df.iterrows():
            try:
                key = (str(r["Sequencing_ID"]), int(r["Barcode"]))
            except (TypeError, ValueError):
                continue
            m = self.sample_metrics.get(key, {})
            if not m.get("evidence_source"):
                continue
            rows.append({
                "Sample": key[0],
                "Barcode": f"barcode{key[1]:02d}",
                "Evidence source": m["evidence_source"],
                "Evidence status": m.get("evidence_label") or "OK",
                "BAM-only positions": m.get("bam_only"),
                "Clair3/BAM discordant": m.get("discordant"),
                "Discordant primary positions": m.get("discordant_primary"),
                "QC status": r.get("QCStatus", ""),
                "Typing result": r.get("ExtendedGenotype", ""),
            })
        return rows or None

    def write_evidence_summary(self, final_df, path="ABO_evidence_summary_mqc.tsv"):
        """MultiQC custom-content table: where each sample's evidence came from, and which
        samples had failed or empty Clair3 output, no reads over the ABO region, or disagreement
        between Clair3 and the BAM pileup."""
        rows = self._evidence_summary(final_df)
        if not rows:
            return
        cols = list(rows[0].keys())
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("# id: 'abo_evidence_summary'\n")
                fh.write("# section_name: 'ABO evidence summary'\n")
                fh.write("# description: 'Agreement between the two lines of evidence (Clair3 and BAM pileup) "
                         "per sample. Samples with no reads over the ABO region, Clair3-only or BAM-only "
                         "evidence, or discordant primary positions need review. "
                         + RESEARCH_USE_STATEMENT.replace("'", "") + "'\n")
                fh.write("# plot_type: 'table'\n")
                fh.write("\t".join(cols) + "\n")
                for row in rows:
                    fh.write("\t".join("" if row[c] is None else str(row[c]) for c in cols) + "\n")
            print(f"Evidence summary written to {path}")
        except OSError as e:
            print(f"Could not write evidence summary: {e}")

    def report_qc(self, final_df):
        """Print (and optionally write) a run-level QC summary. Called after
        the results are saved, so final_df has single-level column names."""
        lines = ["=== QC SUMMARY ===", RESEARCH_USE_STATEMENT]
        lines += [f"{k}: {v}" for k, v in self._about_items()]
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

        ev_rows = self._evidence_summary(final_df)
        if ev_rows:
            lines.append("Evidence (Clair3 + BAM pileup):")
            for k, v in Counter(r["Evidence source"] for r in ev_rows).most_common():
                lines.append(f"  {v:6d}  {k}")
            problems = [r for r in ev_rows if r["Evidence status"] != "OK" or r["Discordant primary positions"]]
            if problems:
                lines.append("Samples needing attention:")
                for r in problems:
                    extra = (f"; discordant at {r['Discordant primary positions']} primary position(s)"
                             if r["Discordant primary positions"] else "")
                    lines.append(f"  {r['Sample']} ({r['Barcode']}): {r['Evidence status']}{extra}")

        if "EvidenceBasis" in final_df.columns:
            lines.append("Evidence basis: " + ", ".join(
                f"{k}={v}" for k, v in final_df["EvidenceBasis"].value_counts().items()))
        if "A2Trace" in final_df.columns:
            lines.append("A2 trace:")
            for k, v in final_df["A2Trace"].value_counts().sort_index().items():
                lines.append(f"  {v:6d}  {k}")
        rs = self.replicate_stats
        if rs["groups"]:
            bad = rs["discordant_genotype"] + rs["discordant_a2"] + rs.get("discordant_allele", 0)
            lines.append(f"Replicate groups: {rs['groups']}   disagreeing: {bad} "
                         f"(genotype {rs['discordant_genotype']}, A2 evidence {rs['discordant_a2']}, "
                         f"allele level {rs.get('discordant_allele', 0)})")

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
        self.write_evidence_summary(final_df)

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
                        # tokens are ISBT names (ABO*A1.01, ABO*B.01, ...); "ABO*" itself contains a "B"
                        toks = str(genotype).split("/")
                        has_a = any(t.startswith("ABO*A") for t in toks)
                        has_b = any(t.startswith("ABO*B") for t in toks)
                        if has_a and has_b:
                            ab_subtypes[genotype] = count
                        elif has_a:
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
  - ABO_result.xlsx has two sheets: "Typing" (one row per sample with read /
    mappability / depth metrics from samtools stats + flagstat, the ISBT-style
    typing result and QC) and "Variants" (per-position table)
  - ExtendedGenotype uses ISBT-style names, e.g. ABO*A2.01/ABO*O.01

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
    parser.add_argument("--nanoplot-dir", default=None, metavar="DIR",
                         help="Folder holding NanoPlot output (*NanoStats.txt); used for Number of reads / Read Length "
                              "in the Typing sheet when a sample folder has no NanoStats.txt of its own")
    parser.add_argument("--threads", "-t", type=int, default=1, metavar="N",
                         help="Worker processes for per-sample parsing (default: %(default)s, serial). "
                              "Useful for large runs; results are identical.")
    parser.add_argument("--min-depth", type=int, default=30, metavar="N",
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
                                      min_depth=args.min_depth, qc_report=args.qc_report,
                                      nanoplot_dir=args.nanoplot_dir, threads=args.threads)
    parser_instance.run()
    print("All done!\n")
