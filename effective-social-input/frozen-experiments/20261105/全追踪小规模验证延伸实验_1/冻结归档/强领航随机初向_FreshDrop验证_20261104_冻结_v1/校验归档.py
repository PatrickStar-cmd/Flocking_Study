"""Verify frozen archive without changing its contents. Python standard library only."""
from pathlib import Path
import csv, hashlib, json, sys, zipfile

ROOT = Path(__file__).resolve().parent

def digest_file(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

try:
    for line in (ROOT/"SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ",1)
        assert digest_file(ROOT/relative) == expected, f"Outer file mismatch: {relative}"
    metadata = json.loads((ROOT/"FREEZE_METADATA.json").read_text(encoding="utf-8"))
    with (ROOT/"FROZEN_INVENTORY.csv").open(encoding="utf-8-sig",newline="") as stream:
        rows = list(csv.DictReader(stream))
    archive = ROOT/metadata["archive_file"]
    assert archive.stat().st_size == metadata["archive_bytes"]
    assert digest_file(archive) == metadata["archive_sha256"]
    assert len(rows) == metadata["payload_file_count"]
    with zipfile.ZipFile(archive) as package:
        assert len(package.infolist()) == len(rows)
        assert set(package.namelist()) == {row["archive_path"] for row in rows}
        for row in rows:
            h = hashlib.sha256()
            size = 0
            with package.open(row["archive_path"]) as stream:
                for block in iter(lambda:stream.read(1024*1024),b""):
                    size += len(block)
                    h.update(block)
            assert size == int(row["bytes"]), row["archive_path"]
            assert h.hexdigest() == row["sha256"], row["archive_path"]
    print(f"PASS: {len(rows)} payload files verified; outer hashes and ZIP contents match.")
    print(f"Archive: {archive}")
except Exception as error:
    print(f"FAIL: {type(error).__name__}: {error}",file=sys.stderr)
    sys.exit(1)
