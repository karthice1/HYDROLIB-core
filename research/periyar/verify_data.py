"""Check downloaded source data against data_manifest.csv (file size and SHA-256).

Usage: python -m research.periyar.verify_data DATA_DIR

Files listed under research/periyar/data are checked in the repository. The exit code is
the number of missing or mismatched files.
"""

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parent.parent


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(data_dir):
    problems = 0
    with open(HERE / "data_manifest.csv", newline="") as f:
        for row in csv.DictReader(f):
            name = row["file (under DATA_DIR unless stated)"]
            path = REPO / name if name.startswith("research/") else Path(data_dir) / name
            if not path.exists():
                print(f"MISSING  {name}")
                problems += 1
            elif path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
                print(f"MISMATCH {name}")
                problems += 1
            else:
                print(f"OK       {name}")
    print(f"{problems} problem(s)")
    return problems


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
