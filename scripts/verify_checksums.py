"""Verify files against a SHA-256 manifest on any operating system (replaces `sha256sum -c`).
Usage: python scripts/verify_checksums.py <folder containing MANIFEST.sha256 or SHA256SUMS>
"""
import hashlib, os, sys
folder = sys.argv[1] if len(sys.argv) > 1 else '.'
man = next((os.path.join(folder, n) for n in ('MANIFEST.sha256', 'SHA256SUMS') if os.path.exists(os.path.join(folder, n))), None)
if man is None: sys.exit(f'no MANIFEST.sha256 or SHA256SUMS in {folder}')
ok = bad = 0
for line in open(man, encoding='utf-8'):
    if not line.strip(): continue
    digest, name = line.split(maxsplit=1); name = name.strip().lstrip('*')
    path = os.path.join(folder, name.replace('/', os.sep))
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 20), b''): h.update(chunk)
        good = h.hexdigest() == digest
    except FileNotFoundError:
        good = False
    ok += good; bad += not good
    if not good: print(f'FAILED  {name}')
print(f'{ok}/{ok + bad} OK'); sys.exit(1 if bad else 0)
