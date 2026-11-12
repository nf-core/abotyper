#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import re
import sys
import glob
from typing import Dict, List, Optional

import pandas as pd

from abo_panel import load_panel, Panel, VariantMarker
from abo_calling import CallingMixin, missing_required_markers  # noqa: F401
from abo_calling import GENOTYPE_TABLE, resolve_primary_genotype  # noqa: F401  (re-exported for tests)
from abo_inputs import (find_nanostats, parse_alignment_metrics, parse_coverage_file,
                        parse_evidence_status, parse_nanostats, read_controls)
from abo_phase import PhaseEvidence, parse_haplotypes_tsv
from abo_report import ReportMixin

# Authorship and licence: provided as part of the nf-core/abotyper pipeline.
__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2026, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = [
    "Fredrick Mobegi",
    "Benedict Matern",
    "Mathijs Groeneweg",
    "Claude Sonnet 5",
]
__license__ = "MIT"
__version__ = "2.5.0"
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
phenotype, genotype, the variants detected, QC) and "Variants" (one column per
panel position holding the call, e.g. "A 50% (215x)"; the full per-base table
is kept in ABO_result.txt, formerly "ABO_Result"). ExtendedGenotype is written with ISBT-style names (A1 ->
ABO*A1.01, O1 -> ABO*O.01.01, B -> ABO*B.01 are only reference defaults; when
non-primary panel markers fire, the allele is narrowed to the one they point
to). The decision tree and QC logic still work on the internal tokens and only
 the written value is encoded.
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


# ===========================================================================
# Phase Confidence Scoring — reads Haplotypes.tsv from pysam_haploscan.py
# ===========================================================================


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


# ===========================================================================
# Primary genotype rules
# ===========================================================================


class ABOReportParser(CallingMixin, ReportMixin):
    """
    Collates all ABO phenotype results from each sample into a general table
    and generates an Excel worksheet and a CSV file for export to LIS soft
    or other general purpose lab management systems. Column layout, position
    parsing, and subtype marker scanning are driven by an external variant
    panel (see abo_panel.py).
    """

    def __init__(self, input_dir, panel: Panel, default_barcode="barcode00",
                 min_depth=30, qc_report="ABO_run_qc.txt", nanoplot_dir=None, threads=1,
                 negative_control_pattern=None, positive_control_pattern=None, controls=None):
        self.input_dir = input_dir
        self.script_version = __version__      # shown in the About sheet
        self.panel = panel
        self.default_barcode = default_barcode
        self.results = []
        self.failed_samples = []
        self.output_errors = []               # outputs that could not be written: the run must fail
        self.controls = dict(controls or {})  # {sample folder name: 'negative' | 'positive'}, from --controls
        self._current_control = ""            # control type of the sample being processed
        self.min_depth = min_depth            # reads below which an A1 call is flagged and a sample fails the depth check
        # optional, off by default: sample names matching this pattern are also negative controls
        self.negative_control_re = re.compile(negative_control_pattern, re.IGNORECASE) if negative_control_pattern else None
        # positive controls must show normal depth: a thin positive control means the run or the library failed
        self.positive_control_re = re.compile(positive_control_pattern, re.IGNORECASE) if positive_control_pattern else None
        self.qc_report_path = qc_report       # run-level QC text report ("" disables)
        self.qc_records = []                  # one dict per sample, for the run summary
        self.threads = max(1, int(threads or 1))   # worker processes for per-sample parsing
        self.nanoplot_dir = nanoplot_dir      # optional folder of NanoPlot output (NanoStats.txt)
        self._fired_markers = []
        self._absent_markers = []
        self._fired_pct = {}                  # marker id -> alt % of the call that fired it
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

        # A sample folder is named by the samplesheet ID. Any name is accepted; a trailing
        # "_barcode<digits>" is reported as the barcode, anything else is part of the sample name.
        self.barcode_suffix = re.compile(r"^(?P<name>.+)_barcode(?P<num>\d+)$", re.IGNORECASE)

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

        self.notes_columns = ["ASubtype", "ReadReliability", "DepthCheck", "PhaseConfidence", "BwSubtype",
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
    def reliability_label(self, band: str) -> str:
        """Read-reliability label, worded from the panel's thresholds ("Very Low(<=20 reads)" by default)."""
        t = self.panel.threshold
        return {"very_low": f"Very Low(<={t('reliability_very_low_reads')} reads)",
                "low": f"Low (<={t('reliability_low_reads')} reads)",
                "robust": f"Robust(>={t('reliability_robust_reads')} reads)"}[band]

    def extract_sample_info(self, filename):
        """Sample name and barcode from a sample folder name.

        "S1_barcode07" gives ("S1", "barcode07"). A name without a trailing _barcode<digits> keeps the whole
        name and gets the default barcode, so "patient-01", "S1.r1" and "X_barcode06_rep2" are all valid.
        """
        m = self.barcode_suffix.match(filename)
        if m:
            num = m.group("num")
            barcode = f"barcode{num.zfill(2)}" if int(num) <= 99 else f"barcode{num}"
            return m.group("name"), barcode, "explicit"
        print(f"No barcode found in '{filename}', using default barcode: {self.default_barcode}")
        return filename, self.default_barcode, "default"

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

        except OSError as e:
            print(f"Error reading {exon_label} file {filename}: {str(e)}")
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


    # -----------------------------------------------------------------
    # Subtype marker scanner
    # -----------------------------------------------------------------


    # -----------------------------------------------------------------
    # Phenotype/genotype assignment
    # -----------------------------------------------------------------


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
            self._current_control = self.controls.get(filename, "")

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
            self.failed_samples.append({"sample": filename, "reason": f"Processing error: {e}"})

    def _record_sample_metrics(self, result_df, sample_name, barcode, sample_dir):
        """Collect the Typing-sheet metrics for one sample (never raises)."""
        try:
            metrics = parse_alignment_metrics(sample_dir)
        except (OSError, ValueError):
            metrics = {}
        # Duplex (ABO-HLA) runs: HLA reads never map to the ABO reference, so whole-library
        # mappability is misleading. Report reads in the library (NanoPlot, else samtools),
        # reads on the ABO reference (samtools coverage) and on-target % between them.
        try:
            metrics.update(parse_coverage_file(sample_dir))
        except (OSError, ValueError):
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
        except (OSError, ValueError, KeyError, IndexError, TypeError) as status_err:
            # without the evidence summary the Clair3/BAM checks cannot be made: never report PASS
            print(f"Warning: could not read the evidence status for {sample_name}: {status_err}")
            old = str(result_df.at[0, ("Notes", "QCFlags")] or "")
            result_df.at[0, ("Notes", "QCFlags")] = ((old + "; " if old else "")
                                                     + f"evidence status unreadable ({status_err})")
            result_df.at[0, ("Notes", "QCStatus")] = "REVIEW"
        try:
            # NanoPlot runs on the raw reads, so it is the better source for read count/length
            ns_file = find_nanostats(sample_dir, self.nanoplot_dir, str(sample_name), str(barcode))
            if ns_file:
                ns = parse_nanostats(ns_file)
                for k in ("reads", "read_length"):
                    if ns.get(k) is not None:
                        metrics[k] = ns[k]
        except (OSError, ValueError):
            pass
        try:
            # the coverage table counts secondary/supplementary alignments; samtools stats'
            # primary mapped-read count is the ceiling when we have it
            if metrics.get("abo_reads") is not None and metrics.get("mapped") is not None:
                metrics["abo_reads"] = min(metrics["abo_reads"], metrics["mapped"])
            if metrics.get("abo_reads") is not None and metrics.get("reads"):
                metrics["mappability"] = min(100.0, 100.0 * metrics["abo_reads"] / metrics["reads"])
        except (TypeError, ValueError, ZeroDivisionError):
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
        except (KeyError, IndexError, TypeError, ValueError):
            metrics.setdefault("median_analyzed_depth", None)
            metrics.setdefault("question_marks", None)
        try:
            key = (str(sample_name), int(str(barcode).replace("barcode", "")))
        except ValueError:
            return
        self.sample_metrics[key] = metrics


    # -----------------------------------------------------------------
    # Allele naming by haplotype (Clair3 phasing)
    # -----------------------------------------------------------------
    _TOKEN_NEEDS = {"O1": {"o1"}, "O2": {"o2"}, "O3": {"o3"}, "B": {"b796", "b803"},
                    "A2.01": {"a467", "a1061"}, "A2.06": {"a1061"}, "A1.02": {"a467"}}


    def process_files(self):
        """Process all files in the input directory that match expected patterns."""
        print(f"Scanning directory: {self.input_dir}")

        valid_directories = []
        for filename in os.listdir(self.input_dir):
            if os.path.isdir(os.path.join(self.input_dir, filename)):
                valid_directories.append(filename)

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

        print("\n--- Processing Statistics ---")
        print(f"Files with explicit barcode: {self.processing_stats['with_barcode']}")
        print(f"Files using default barcode: {self.processing_stats['without_barcode']}")
        print(f"Total pattern matches: {self.processing_stats['pattern_matched']}")
        print(f"Pattern match failures: {self.processing_stats['pattern_failed']}")
        print("----------------------------")

    def _process_files_parallel(self, valid_directories):
        """Parse samples in a process pool (--threads N). Output order does not matter: the
        merged table is sorted by sample and barcode, and each sample's log is printed whole,
        in directory order, once its worker finishes."""
        from concurrent.futures import ProcessPoolExecutor

        for filename in valid_directories:        # cheap bookkeeping stays serial
            self.processing_stats["pattern_matched"] += 1
            _, _, pattern_type = self.extract_sample_info(filename)
            self.processing_stats["with_barcode" if pattern_type == "explicit" else "without_barcode"] += 1

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


    # -----------------------------------------------------------------
    # QC: replicate cross-check and run-level summary (v2.1.0)
    # -----------------------------------------------------------------
    @staticmethod
    def _alleles_compatible(x: str, y: str) -> bool:
        """Two ISBT-style allele names do not conflict when they are equal or one is a less specific form of
        the other (ABO*A2 / ABO*A2.01, ABO*O.01 / ABO*O.01.01, or a bare ABO*A / ABO*A1.01)."""
        a, b = x.replace("ABO*", "").strip(), y.replace("ABO*", "").strip()
        if a == b:
            return True
        for s, t in ((a, b), (b, a)):
            if t.startswith(s + "."):
                return True
            if len(s) == 1 and s.isalpha() and t[:1] == s and (len(t) == 1 or not t[1].isalpha()):
                return True
        return False

    @classmethod
    def _typings_compatible(cls, e1: str, e2: str) -> bool:
        """Do two ExtendedGenotype strings describe the same pair of alleles, allowing for one being less specific?"""
        t1, t2 = str(e1).split("/"), str(e2).split("/")
        if len(t1) != len(t2):
            return str(e1) == str(e2)
        if len(t1) == 2:
            return ((cls._alleles_compatible(t1[0], t2[0]) and cls._alleles_compatible(t1[1], t2[1]))
                    or (cls._alleles_compatible(t1[0], t2[1]) and cls._alleles_compatible(t1[1], t2[0])))
        return all(cls._alleles_compatible(a, c) for a, c in zip(t1, t2))

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
            elif len(exts) > 1 and not all(self._typings_compatible(p, q) for p in exts for q in exts):
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


    def run(self) -> int:
        """Run the ABOReportParser. Returns the exit status: 1 when any sample failed or an output
        could not be written (the outputs that could be written are still written), else 0."""
        self.process_files()
        final_df = self.merge_dataframes()
        if final_df.empty:
            print("No results to report. Exiting.")
            self._print_failures()
            return 1
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

        print("\n=== PROCESSING SUMMARY ===")
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
                print("\nPhenotype Distribution:")
                for phenotype, count in phenotype_counts.items():
                    print(f"  {phenotype}: {count}")

                try:
                    extended_genotype_counts = final_df["ExtendedGenotype"].value_counts()
                    print("\nExtended Genotype Distribution (with A subtypes):")
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
                        print("  A Subtypes:")
                        for genotype, count in sorted(a_subtypes.items()):
                            print(f"    {genotype}: {count}")
                    if ab_subtypes:
                        print("  AB Subtypes:")
                        for genotype, count in sorted(ab_subtypes.items()):
                            print(f"    {genotype}: {count}")
                    if other_genotypes:
                        print("  Other Genotypes:")
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
                        print("\nPhase Confidence Distribution:")
                        for val, count in pc_col.value_counts().head(10).items():
                            print(f"    {val}: {count}")
                except Exception as e:
                    print(f"Could not analyze phase confidence: {e}")

            except Exception as e:
                print(f"Could not analyze phenotype distribution: {e}")

        print("============================")
        self.report_qc(final_df)
        return self._print_failures()

    def _print_failures(self) -> int:
        """Summarise failed samples and output errors; return the exit status."""
        if not self.failed_samples and not self.output_errors:
            return 0
        print("\nERROR: the run is incomplete.", file=sys.stderr)
        for sample in self.failed_samples:
            print(f"  failed sample {sample['sample']}: {sample['reason']}", file=sys.stderr)
        for err in self.output_errors:
            print(f"  output not written: {err}", file=sys.stderr)
        return 1


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
    typing result, variants detected and QC) and "Variants" (one column per panel
    position)
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
                         help="Flag A1 calls made on fewer than N reads at c.467/c.1061delC, and samples with fewer than "
                              "N reads at a primary position (depth check) (default: %(default)s)")
    parser.add_argument("--controls", default=None, metavar="FILE",
                         help="Tab-separated file of sample folder name and control type ('negative', 'positive' or "
                              "empty), from the samplesheet `control` column. Negative controls are expected to have "
                              "few reads, and reads above --min-depth are flagged as possible contamination; positive "
                              "controls need normal depth.")
    parser.add_argument("--positive-control-pattern", default=None, metavar="REGEX",
                         help="Optional: also treat sample names matching this regular expression (case-insensitive) "
                              "as positive controls when --controls does not list them (default: off)")
    parser.add_argument("--negative-control-pattern", default=None, metavar="REGEX",
                         help="Optional: also treat sample names matching this regular expression (case-insensitive) "
                              "as negative controls when --controls does not list them (default: off)")
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
        print("Verbose mode enabled")
        print(f"Input directory: {args.input_directory}")
        print(f"Default barcode: {args.default_barcode}")
        print(f"Panel: {args.panel}")
        print(f"Script version: {__version__}")

    try:
        panel = load_panel(args.panel)
    except Exception as exc:
        print(f"! CRITICAL: could not load variant panel '{args.panel}': {exc}")
        sys.exit(1)

    missing = missing_required_markers(panel)
    if missing:
        print(f"! CRITICAL: the variant panel '{args.panel}' does not define marker(s) this script needs: "
              f"{', '.join(missing)}. Restore them (they may be recalibrated, not renamed).")
        sys.exit(1)

    try:
        controls = read_controls(args.controls)
    except (OSError, ValueError) as exc:
        print(f"! CRITICAL: could not read the controls file '{args.controls}': {exc}")
        sys.exit(1)

    print(f"Using default barcode: {args.default_barcode}")
    print(f"Loaded panel v{panel.meta.get('panel_version', '?')} with {len(panel.variants)} variants "
          f"({len(panel.calibrated())} usable now, {len(panel.uncalibrated())} awaiting calibration)")

    parser_instance = ABOReportParser(args.input_directory, panel, args.default_barcode,
                                      min_depth=args.min_depth, qc_report=args.qc_report,
                                      nanoplot_dir=args.nanoplot_dir, threads=args.threads,
                                      negative_control_pattern=args.negative_control_pattern,
                                      positive_control_pattern=args.positive_control_pattern,
                                      controls=controls)
    status = parser_instance.run()
    print("All done!\n" if status == 0 else "Finished with errors.\n")
    sys.exit(status)
