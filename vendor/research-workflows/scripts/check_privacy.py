#!/usr/bin/env python3
"""Check staged/tracked portable files before commits and pushes."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

PATTERNS=[rb'-----BEGIN (?:OPENSSH |RSA |EC |DSA )?PRIVATE KEY-----',
          rb'\bgh[pousr]_[A-Za-z0-9]{20,}\b',rb'\bgithub_pat_[A-Za-z0-9_]{30,}\b',
          rb'\bAKIA[A-Z0-9]{16}\b',rb'/Users/[A-Za-z][A-Za-z0-9_-]*/',
          rb'/(?:scratch|data)/[A-Za-z][A-Za-z0-9_-]*/']

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--tracked',action='store_true');args=p.parse_args()
    root=args.root.resolve();private=root/'local/fleet.json';values=[]
    if private.exists():
        home=Path(json.loads(private.read_text())['fleet_home']);registry=home/'local/machines.json'
        if registry.exists():
            for machine in json.loads(registry.read_text())['machines'].values():
                for key in ['hostname','identity_file']:
                    value=machine.get(key)
                    if isinstance(value,str) and len(value)>5:values.append(value.encode())
                values.extend(s.encode() for s in machine.get('research_roots',[]))
    command=['git','ls-files','-z'] if args.tracked else ['git','diff','--cached','--name-only','--diff-filter=ACMR','-z']
    names=subprocess.check_output(command,cwd=root).decode().split('\0');bad=[]
    for name in filter(None,names):
        path=Path(name)
        mode=subprocess.check_output(['git','ls-files','-s','--',name],cwd=root).decode().split()[0]
        if any(part in {'local','reference','.ssh','.aws','__pycache__','drafts'} or part.startswith('.env') for part in path.parts) or mode=='120000':
            bad.append(name);continue
        # A gitlink records a commit in another repository, not a blob in this
        # index. Its path is checked above; .gitmodules is scanned normally.
        if mode=='160000':continue
        body=subprocess.check_output(['git','show',':'+name],cwd=root)
        if any(re.search(pattern,body) for pattern in PATTERNS) or any(value in body for value in values):bad.append(name)
    if bad:print('Privacy check refused: '+', '.join(bad),file=sys.stderr);return 1
    print('Privacy check passed.');return 0

if __name__=='__main__':sys.exit(main())
