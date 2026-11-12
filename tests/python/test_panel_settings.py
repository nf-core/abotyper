"""Marker roles and thresholds come from the panel (meta.roles / meta.thresholds)."""
import copy

import pytest
import yaml

from conftest import PANEL, STATE, run

import abo_panel


@pytest.fixture(scope="module")
def panel_yaml():
    return yaml.safe_load(PANEL.read_text(encoding="utf-8"))


def aggregate_with(panel_data, tmp_path, chain):
    """Rerun the aggregation of the 10 samples with another panel; return final_export.csv and ABO_result.txt."""
    p = tmp_path / "panel.yaml"
    p.write_text(yaml.safe_dump(panel_data, sort_keys=False), encoding="utf-8")
    run("aggregate_abo_reports.py", STATE["aggregate_dir"] / "per_sample_processing", "--panel", p, "--threads", "1",
        "--qc-report", "", cwd=tmp_path)
    return (tmp_path / "final_export.csv").read_text(encoding="utf-8"), (tmp_path / "ABO_result.txt").read_text(encoding="utf-8")


def test_renamed_marker_with_its_role_gives_the_same_results(panel_yaml, tmp_path, chain):
    data = copy.deepcopy(panel_yaml)
    for v in data["variants"]:
        if v["id"] == "b_vs_ao_796":
            v["id"] = "c796_renamed"
    data["meta"]["roles"] = {"b_796": "c796_renamed"}
    csv_text, _ = aggregate_with(data, tmp_path, chain)
    assert csv_text == (STATE["aggregate_dir"] / "final_export.csv").read_text(encoding="utf-8")


def test_unknown_role_or_threshold_key_fails(panel_yaml, tmp_path):
    for key, value in (("roles", {"b796": "b_vs_ao_796"}), ("thresholds", {"a2_minreads": 30}), ("thresholds", {"a2_min_reads": "x"})):
        data = copy.deepcopy(panel_yaml)
        data["meta"][key] = value
        p = tmp_path / "panel.yaml"
        p.write_text(yaml.safe_dump(data), encoding="utf-8")
        with pytest.raises(abo_panel.PanelError):
            abo_panel.load_panel(p)


def test_threshold_override_is_used(panel_yaml, tmp_path, chain):
    data = copy.deepcopy(panel_yaml)
    data["meta"]["thresholds"] = {"reliability_very_low_reads": 100000}
    csv_text, _ = aggregate_with(data, tmp_path, chain)
    rows = csv_text.splitlines()[1:]
    assert rows and all("Very Low(<=100000 reads)" in r for r in rows)
