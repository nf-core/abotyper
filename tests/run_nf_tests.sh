#!/usr/bin/env bash
# Run the abotyper nf-test suite on a Linux host (not from PowerShell/Windows, where nf-test and the
# container runtimes are not supported).
#
#   tests/run_nf_tests.sh [docker|singularity|conda] [--update-snapshot]
#
# First run on a new machine: pass --update-snapshot once so nf-test writes the .snap files for the
# new tests (commit them); later runs compare against those snapshots.
set -euo pipefail

PROFILE="${1:-docker}"
shift || true
EXTRA=("$@")

cd "$(git rev-parse --show-toplevel)"

echo "== environment =="
command -v nextflow >/dev/null || { echo "nextflow not found"; exit 2; }
command -v nf-test  >/dev/null || { echo "nf-test not found";  exit 2; }
nextflow -version | head -3
nf-test version

echo "== sanity: every script in bin/ must be executable (exit status 126 otherwise) =="
bad=0
for f in bin/*.py; do
    if [ ! -x "$f" ]; then echo "NOT EXECUTABLE: $f"; bad=1; fi
done
[ "$bad" -eq 0 ] || { echo "fix with: git update-index --chmod=+x <file> && chmod +x <file>"; exit 3; }

echo "== local modules and subworkflows, stub then real (profile: test,${PROFILE}) =="
# nf-test combines -profile test (from nf-test.config) with the profile given here
nf-test test --tag abotyper --profile "+${PROFILE}" --verbose "${EXTRA[@]}"
