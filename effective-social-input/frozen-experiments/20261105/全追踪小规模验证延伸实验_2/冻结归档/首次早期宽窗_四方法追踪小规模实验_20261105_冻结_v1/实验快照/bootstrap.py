"""Locate optional, already-installed dependencies without changing the old experiments."""
from pathlib import Path
import os
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

ROOT = Path(__file__).resolve().parent
for candidate in (ROOT / '.runtime_deps', ROOT.parent / '.runtime_deps',
                  ROOT.parents[1] / '追踪验证' / '.runtime_deps'):
    if candidate.is_dir():
        sys.path.insert(0, str(candidate))
        break
