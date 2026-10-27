from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
INVENTORY_EXCLUSIONS = {"MANIFEST.csv", "SHA256SUMS.txt"}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

failures = []
checked = 0
covered = set()
for line in (ROOT / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    expected, rel = line.split("  ", 1)
    covered.add(rel)
    path = ROOT / rel
    if not path.is_file():
        failures.append({"file": rel, "issue": "missing"})
        continue
    observed = sha256(path)
    if observed != expected:
        failures.append({"file": rel, "issue": "SHA256 mismatch", "observed": observed})
    checked += 1
for path in ROOT.rglob("*"):
    rel = path.relative_to(ROOT).as_posix()
    if not rel.isascii():
        failures.append({"file": rel, "issue": "non-ASCII archive path"})
    if path.is_dir() and not any(path.iterdir()):
        failures.append({"file": rel, "issue": "empty directory"})
    if path.is_file() and rel not in covered and rel not in INVENTORY_EXCLUSIONS:
        failures.append({"file": rel, "issue": "file absent from SHA256SUMS.txt"})
report = {"status": "PASS" if not failures else "FAIL", "files_checked": checked, "failures": failures}
print(json.dumps(report, indent=2))
sys.exit(0 if not failures else 1)
