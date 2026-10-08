import json
from pathlib import Path
from hashlib import sha256
import platform
import time
from .experiment import run
from .audit import audit

ROOT = Path(__file__).resolve().parent


def main():
    start = time.perf_counter()
    data = run()
    files = sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md', ROOT/'PRIOR_ART.md']
    data['source_manifest'] = {p.name: sha256(p.read_bytes()).hexdigest() for p in files}
    data['runtime'] = dict(python=platform.python_version(), platform=platform.platform())
    data['audit_at_generation'] = audit(data)
    data['elapsed_seconds'] = time.perf_counter()-start
    directory = ROOT/'results'; directory.mkdir(exist_ok=True)
    destination = directory/'verification.json'
    destination.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(audit=data['audit_at_generation'], elapsed_seconds=data['elapsed_seconds'],
                          sha256=sha256(destination.read_bytes()).hexdigest()), indent=2), flush=True)


if __name__ == '__main__':
    main()
