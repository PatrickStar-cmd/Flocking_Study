"""Freeze reviewed executable sources before the scientific replay batch."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def run(root):
    here = Path(__file__).resolve().parent
    destination = root/'sources/p2'
    assert not destination.exists()
    destination.mkdir()
    for path in sorted(here.iterdir()):
        if path.suffix in {'.py','.m'} or path.name == 'requirements.txt':
            shutil.copyfile(path, destination/path.name)
    sources = []
    for path in sorted((root/'sources').rglob('*')):
        if path.is_file():
            sources.append({'path':path.relative_to(root).as_posix(),
                'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest(), 'bytes':path.stat().st_size})
    record = {'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
              'branch':subprocess.check_output(['git','branch','--show-current'],text=True).strip(),
              'status_before_replay':subprocess.check_output(['git','status','--short'],text=True),
              'sources':sources}
    path = root/'source_hashes.json'; assert not path.exists()
    path.write_text(json.dumps(record,indent=2),encoding='utf-8')
    print('Frozen source files:',len(sources))


if __name__ == '__main__':
    run(Path(sys.argv[1]))
