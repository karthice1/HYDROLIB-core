#!/bin/bash
# Install this repository in a new Python environment and run the release checks:
# the GSA unit tests, the published-number checks and (if DATA_DIR is given) the
# source-data checksums. Writes a report to research/periyar/test_report/.
#
# Usage: bash research/periyar/release_checks.sh [DATA_DIR]
# Run from the repository root. Needs python3.13 and network access to PyPI.
set -u
ROOT=$(pwd)
OUT=$ROOT/research/periyar/test_report
VENV=$(mktemp -d)/venv
mkdir -p "$OUT"
python3.13 -m venv "$VENV" && . "$VENV/bin/activate"
pip install -q --upgrade pip
pip install -q -e . && pip install -q -r research/periyar/requirements.txt pytest
{
  echo "# Release check report"
  echo
  echo "- Date (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "- Commit: $(git rev-parse HEAD 2>/dev/null || cat COMMIT.txt)"
  echo "- Python: $(python --version 2>&1)"
  echo "- Platform: $(uname -sm)"
  echo "- Environment: new virtual environment, installed with pip install -e . and requirements.txt"
} > "$OUT/REPORT.md"
python -m pytest research/gsa/tests -v --junitxml="$OUT/gsa_tests.xml" > "$OUT/gsa_tests.log" 2>&1
TESTS=$?
python -m research.periyar.check_results > "$OUT/check_results.log" 2>&1
CHECKS=$?
DATA="not run (no DATA_DIR given)"
if [ $# -ge 1 ]; then
  python -m research.periyar.verify_data "$1" > "$OUT/verify_data.log" 2>&1
  DATA="exit code $? ($(tail -1 "$OUT/verify_data.log"))"
fi
pip freeze --exclude-editable > "$OUT/pip_freeze.txt"
{
  echo
  echo "| Check | Result |"
  echo "|---|---|"
  echo "| GSA unit tests | exit code $TESTS ($(grep -Eo '[0-9]+ passed[^=]*' "$OUT/gsa_tests.log" | tail -1)) |"
  echo "| Published numbers (check_results) | exit code $CHECKS ($(grep -c '^PASS' "$OUT/check_results.log") PASS, $(grep -c '^FAIL' "$OUT/check_results.log") FAIL) |"
  echo "| Source data (verify_data) | $DATA |"
  echo
  echo "These checks do not run D-Flow FM. Per-run results, figures, structural tests and"
  echo "the drainage diagnosis are reproduced only by rerunning the model (REPRODUCE.md)."
} >> "$OUT/REPORT.md"
cat "$OUT/REPORT.md"
exit $((TESTS + CHECKS))
