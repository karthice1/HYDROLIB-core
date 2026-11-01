#!/bin/bash
# Build the Zenodo archive from one exact commit, so that the archive, COMMIT.txt and the
# clean-environment test report all refer to the same source snapshot.
#
# Usage: bash research/periyar/zenodo/build_archive.sh OUT_DIR [DATA_DIR]
# Run from the repository root on a clean working tree. The commit archived is HEAD; tag it
# (e.g. v1.0-periyar) and push the tag so that GitHub and Zenodo point at the same snapshot.
set -euo pipefail
OUT=$(realpath "$1"); DATA=${2:-}
[ -z "$(git status --porcelain)" ] || { echo "working tree not clean"; exit 1; }
COMMIT=$(git rev-parse HEAD)
NAME=periyar-flood-workflow-v1.0
WORK=$(mktemp -d)
git clone -q "$(pwd)" "$WORK/src" && git -C "$WORK/src" checkout -q "$COMMIT"
(cd "$WORK/src" && bash research/periyar/release_checks.sh $DATA)
mkdir -p "$OUT"
echo "$COMMIT" > "$WORK/COMMIT.txt"
cp research/periyar/zenodo/README.md "$WORK/README_FIRST.md"
git archive --format=zip --prefix="$NAME/" --add-file="$WORK/COMMIT.txt" \
    --add-file="$WORK/README_FIRST.md" -o "$OUT/$NAME.zip" "$COMMIT"
python3 - "$OUT/$NAME.zip" "$WORK/src/research/periyar/test_report" "$NAME" <<'PY'
import sys, zipfile, pathlib
archive, report, name = sys.argv[1], pathlib.Path(sys.argv[2]), sys.argv[3]
with zipfile.ZipFile(archive, "a", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(report.iterdir()):
        z.write(f, f"{name}/research/periyar/test_report/{f.name}")
PY
# Verify the archive as a user would receive it.
mkdir "$WORK/check" && (cd "$WORK/check" && unzip -q "$OUT/$NAME.zip" && cd "$NAME" \
  && [ "$(cat COMMIT.txt)" = "$COMMIT" ] \
  && grep -q "$COMMIT" research/periyar/test_report/REPORT.md \
  && python3 -m research.periyar.check_results | tail -1)
sha256sum "$OUT/$NAME.zip" | tee "$OUT/$NAME.zip.sha256"
echo "archived commit $COMMIT"
rm -rf "$WORK"
