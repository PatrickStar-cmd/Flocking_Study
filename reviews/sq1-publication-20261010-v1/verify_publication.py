"""Check only the latest C publication; never run simulations or rewrite evidence."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREFIXES = ['reviews/sq1-b-local-smoke-audit-v1/', 'reviews/sq1-b-r2-audit-v1/',
            'reviews/sq1-b-r3-audit-v1/', 'reviews/sq1-c-batch-freeze-v1/',
            'reviews/sq1-publication-20261010-v1/']
ENTRIES = ['reviews/README.md', 'reviews/B_HANDOFF.md']


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--revision')
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    rev = args.revision
    if rev:
        paths = git('diff-tree', '--no-commit-id', '--name-only', '-r', '-z', rev)
        ref = rev+':'
    else:
        paths = git('diff', '--cached', '--name-only', '-z')
        ref = ':'
    names = [n for n in paths.decode('utf-8').split('\0') if n]
    assert names, 'No publication files selected.'
    assert all(n in ENTRIES or any(n.startswith(p) for p in PREFIXES) for n in names)
    assert not any('/source/' in n or '/downloads/' in n or Path(n).suffix.lower() in ['.mat', '.zip', '.7z', '.exe', '.dll'] for n in names)
    blobs = {n: git('show', ref+n) for n in names}
    assert all(blobs[n] == (ROOT/n).read_bytes() for n in names), 'Selected bytes differ from working files.'
    assert max(map(len, blobs.values())) < 5_000_000
    freeze = 'reviews/sq1-c-batch-freeze-v1/'
    for r in json.loads(blobs[freeze+'SHA256SUMS.json'])['files']:
        raw = blobs[freeze+r['path']]
        assert hashlib.sha256(raw).hexdigest() == r['sha256'] and len(raw) == r['bytes']
    cfg = json.loads(blobs[freeze+'execution_config.json'])
    runs = json.loads(blobs[freeze+'RUN_BINDINGS.json'])['runs']
    t0 = json.loads(blobs[freeze+'T0_BINDINGS.json'])['runs']
    assert len(runs) == len({r['run_id'] for r in runs}) == 560
    assert len(t0) == len({r['run_id'] for r in t0}) == 80
    assert not cfg['authorization']['full_batch_allowed']
    audit = 'reviews/sq1-b-r3-audit-v1/'
    checks = json.loads(blobs[audit+'source_audit.json'])['delivery_hash_checks']
    assert len(checks) == 96 and all(r['hash_match'] and r['size_match'] for r in checks)
    tests = json.loads(blobs[audit+'C_tests.json'])
    assert len(tests) == 15 and all(t['status'] == 'PASS' for t in tests)
    readers = ENTRIES + ['reviews/sq1-publication-20261010-v1/README.md',
                          audit+'REVIEW.md', freeze+'README.md', freeze+'B_HANDOFF.md', freeze+'BATCH_ACCEPTANCE.md']
    count = 0
    for name in readers:
        for target in re.findall(r'\]\(([^)]+)\)', blobs[name].decode('utf-8')):
            if '://' in target or target.startswith('#'):
                continue
            target = target.split('#')[0]
            dest = ((ROOT/name).parent/target).resolve()
            relative = dest.relative_to(ROOT).as_posix()
            if args.record and relative == 'reviews/sq1-publication-20261010-v1/verification.json' and not dest.exists():
                continue  # The only forward reference is this newly generated audit record.
            assert dest.is_file(), (name, target)
            git('cat-file', '-e', ref+relative)
            count += 1
    result = {'status': 'PASS', 'date': '2026-10-10', 'evidence_type': 'publication_validation',
              'published_files': len(names) + int(args.record and 'reviews/sq1-publication-20261010-v1/verification.json' not in names),
              'selected_bytes_excluding_record': sum(len(b) for n,b in blobs.items() if n != 'reviews/sq1-publication-20261010-v1/verification.json'),
              'reader_links_checked': count, 'freeze_hashes_match': True,
              'main_keys': len(runs), 'T0_keys': len(t0), 'r3_delivery_hash_matches': 96,
              'r3_independent_tests': '15 PASS', 'raw_binaries_staged': 0,
              'full_batch_allowed': False,
              'inventory': [{'path': n, 'sha256': hashlib.sha256(b).hexdigest(), 'bytes': len(b)}
                            for n, b in sorted(blobs.items()) if n != 'reviews/sq1-publication-20261010-v1/verification.json']}
    if args.record:
        with (HERE/'verification.json').open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'inventory'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
