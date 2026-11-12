"""Every output of the pipeline's Python steps on the 10 test samples matches the stored golden copy.

A refactor must leave these unchanged. After an intended change, rerun with UPDATE_GOLDEN=1 and review the diff.
"""
import difflib

import pytest

from conftest import GOLDEN, UPDATE


def test_golden_outputs(chain):
    missing, changed = [], []
    for name, text in sorted(chain.items()):
        g = GOLDEN / name
        if UPDATE:
            g.parent.mkdir(parents=True, exist_ok=True)
            g.write_text(text, encoding="utf-8", newline="\n")
            continue
        if not g.exists():
            missing.append(name)
            continue
        old = g.read_text(encoding="utf-8")
        if old != text:
            diff = "".join(difflib.unified_diff(old.splitlines(True), text.splitlines(True), f"golden/{name}", name, n=1))
            changed.append(diff[:4000])
    if UPDATE:
        pytest.skip(f"golden files rewritten ({len(chain)})")
    assert not missing, f"no golden file for: {missing}"
    assert not changed, "outputs differ from the golden copies:\n" + "\n".join(changed)


def test_golden_set_is_complete(chain):
    """No stale golden files: every stored file is still produced."""
    stored = {str(p.relative_to(GOLDEN)).replace("\\", "/") for p in GOLDEN.rglob("*") if p.is_file()}
    assert stored <= set(chain), f"golden files no longer produced: {sorted(stored - set(chain))}"
