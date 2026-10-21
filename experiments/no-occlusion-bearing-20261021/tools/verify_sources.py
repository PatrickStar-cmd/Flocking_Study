"""Verify exact exported frozen source bytes. Python standard library only."""
from pathlib import Path
import json, hashlib
root=Path(__file__).resolve().parents[1]
rows=json.loads((root/'provenance/frozen-source-manifest.json').read_text(encoding='utf-8'))
for r in rows:
    assert hashlib.sha256((root/r['file']).read_bytes()).hexdigest()==r['sha256'],r['file']
print(f'PASS: {len(rows)} frozen source files match archived hashes')
