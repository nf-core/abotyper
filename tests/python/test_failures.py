"""A sample whose files cannot be parsed fails the run (non-zero exit) but does not stop the other samples' results."""
import shutil
import subprocess
import sys

from conftest import BIN, PANEL, STATE


def test_corrupt_haplotype_table_fails_the_run_but_reports_the_rest(tmp_path, chain):
    src = STATE["aggregate_dir"] / "per_sample_processing"
    dst = tmp_path / "per_sample_processing"
    shutil.copytree(src, dst)
    bad = dst / "ABO_BO1_barcode06" / "combined" / "ABO_BO1_barcode06.Haplotypes.tsv"
    bad.write_text('Read_Name\tHaplotype\n"unterminated\tx\n', encoding="utf-8")
    r = subprocess.run([sys.executable, str(BIN / "aggregate_abo_reports.py"), str(dst), "--panel", str(PANEL),
                        "--qc-report", ""], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 1, r.stdout[-2000:]
    assert "failed sample ABO_BO1_barcode06" in r.stderr
    rows = (tmp_path / "final_export.csv").read_text(encoding="utf-8").splitlines()[1:]
    assert len(rows) == 9 and not any(row.startswith("ABO_BO1,") for row in rows)
