#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Readers for the per-sample files the aggregation uses besides the phenotype report: samtools stats/flagstat/coverage, the ABO_EVIDENCE status file, NanoPlot NanoStats and the controls table.

Part of aggregate_abo_reports.py (nf-core/abotyper); split out of that script so each part can be read and
tested on its own. This file is part of nf-core/abotyper and is released under the MIT License (see LICENSE).
"""
from __future__ import annotations

import glob
import os
import re
from typing import Dict, Optional





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
    except OSError:
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
    except OSError:
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
            rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]
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
            lines = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
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


def read_controls(path: Optional[str]) -> Dict[str, str]:
    """Controls declared in the samplesheet: {sample id: 'negative' | 'positive'}.

    Two tab-separated columns (sample id, control type); an empty type marks a sample.
    """
    controls: Dict[str, str] = {}
    if not path:
        return controls
    with open(path, encoding="utf-8") as fh:
        for n, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            sample = parts[0].strip()
            kind = parts[1].strip().lower() if len(parts) > 1 else ""
            if kind not in ("", "negative", "positive"):
                raise ValueError(f"{path} line {n}: control must be empty, 'negative' or 'positive', got '{kind}'")
            if kind:
                controls[sample] = kind
    return controls
