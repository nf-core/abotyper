"""Shared fixtures for the bin/ script tests.

The golden tests run the pipeline's Python steps, with the same command lines as the modules, on the 10 samples
of the pipeline test (tests/python/data/<sample>/) and compare every output with tests/python/golden/.
Set UPDATE_GOLDEN=1 to rewrite the golden files after an intended change.
"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"
DATA = Path(__file__).resolve().parent / "data"
GOLDEN = Path(__file__).resolve().parent / "golden"
PANEL = ROOT / "assets" / "refs" / "abo_variant_panel.yaml"
SAMPLES = sorted(p.name for p in DATA.iterdir() if p.is_dir())
UPDATE = os.environ.get("UPDATE_GOLDEN") == "1"

sys.path.insert(0, str(BIN))   # the scripts import abo_panel from their own folder, as in the pipeline


def run(script, *args, cwd):
    """Run a bin/ script like a pipeline task; fail with its output when it exits non-zero."""
    r = subprocess.run([sys.executable, str(BIN / script), *map(str, args)], cwd=cwd, capture_output=True, text=True,
                       env={**os.environ, "PYTHONHASHSEED": "0"})
    if r.returncode != 0:
        raise AssertionError(f"{script} exited {r.returncode}\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    return r


def dump_xlsx(path):
    """Every cell of every sheet, as tab-separated text (one block per sheet)."""
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True)
    out = []
    for ws in wb.worksheets:
        out.append(f"## sheet {ws.title}")
        for row in ws.iter_rows(values_only=True):
            out.append("\t".join("" if v is None else str(v) for v in row).rstrip("\t"))
    return "\n".join(out) + "\n"


def normalise(text, tmp):
    """Mask what legitimately differs between runs: temporary paths, run dates and script versions."""
    text = text.replace(str(tmp), "<TMP>")
    text = re.sub(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}(:\d{2}(\.\d+)?)?", "<DATETIME>", text)
    text = re.sub(r"(\.py\W{0,3}(v|version)?\s*)\d+\.\d+\.\d+", r"\1<VERSION>", text)
    return text


@pytest.fixture(scope="session")
def chain(tmp_path_factory):
    """Run every Python step of the pipeline on the 10 fixture samples; return {relative name: text}."""
    tmp = tmp_path_factory.mktemp("chain")
    out = {}
    run("panel_to_bed.py", "--panel", PANEL, "--output", "abo_variant_panel.bed", cwd=tmp)
    out["abo_variant_panel.bed"] = (tmp / "abo_variant_panel.bed").read_text(encoding="utf-8")
    agg = tmp / "aggregate"
    for s in SAMPLES:
        d = DATA / s
        w = tmp / s
        w.mkdir()
        run("clair2metrics.py", "-i", d / f"{s}_clair3.merge_output.gvcf.gz", "-o", f"{s}.clair.AlignmentStatistics.tsv",
            "--panel", PANEL, cwd=w)
        run("clair2haplotypes.py", "-i", d / f"{s}_clair3.phased_merge_output.vcf.gz", "-o", f"{s}.Haplotypes.tsv",
            "--panel", PANEL, cwd=w)
        run("pileup2metrics.py", "--pileup", d / f"{s}.pileup.txt", "--clair-metrics", f"{s}.clair.AlignmentStatistics.tsv",
            "--panel", PANEL, "--sample", s, "--merged-out", f"{s}.AlignmentStatistics.tsv",
            "--pileup-out", f"{s}.pileup_metrics.tsv", "--concordance-out", f"{s}.concordance.tsv",
            "--status-out", f"{s}.evidence_status.tsv", cwd=w)
        run("predict_abo_phenotype.py", "-i", f"{s}.AlignmentStatistics.tsv", "-o", f"{s}.ABOPhenotype.txt",
            "-c", d / f"{s}.coverage.txt", "--panel", PANEL, cwd=w)
        for f in ("clair.AlignmentStatistics.tsv", "Haplotypes.tsv", "AlignmentStatistics.tsv", "pileup_metrics.tsv",
                  "concordance.tsv", "evidence_status.tsv", "ABOPhenotype.txt"):
            out[f"{s}/{s}.{f}"] = (w / f"{s}.{f}").read_text(encoding="utf-8")
        # the per_sample_processing/<id>/combined layout ABO_SNPS2PHENO builds in its task
        c = agg / "per_sample_processing" / s / "combined"
        c.mkdir(parents=True)
        for f in ("ABOPhenotype.txt", "Haplotypes.tsv", "evidence_status.tsv", "concordance.tsv"):
            shutil.copy(w / f"{s}.{f}", c)
        for f in ("coverage.txt", "stats", "flagstat"):
            shutil.copy(d / f"{s}.{f}", c)
    run("aggregate_abo_reports.py", "per_sample_processing", "--panel", PANEL, "--threads", "1", cwd=agg)
    for f in ("ABO_result.txt", "final_export.csv", "ABO_run_qc.txt", "ABO_evidence_summary_mqc.tsv"):
        out[f"aggregate/{f}"] = (agg / f).read_text(encoding="utf-8")
    out["aggregate/ABO_result.xlsx.tsv"] = dump_xlsx(agg / "ABO_result.xlsx")
    return {k: normalise(v, tmp) for k, v in out.items()}
