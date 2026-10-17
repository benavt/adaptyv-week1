"""Resolve project-pinned workflow code before host-level legacy bindings."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import subprocess


SUBMODULE_PATH = 'external/research-workflows'


def _adopted(project: Path) -> bool:
    modules = project / '.gitmodules'
    if not modules.is_file() or modules.is_symlink():
        return False
    result = subprocess.run(
        ['git', 'config', '--file', str(modules), '--get-regexp', r'^submodule\..*\.path$'],
        capture_output=True, text=True, check=False)
    if result.returncode not in (0, 1):
        raise ValueError('Could not parse project .gitmodules')
    for row in result.stdout.splitlines():
        fields = row.split(None, 1)
        if len(fields) == 2 and fields[1] == SUBMODULE_PATH:
            return True
    return False


def resolve_workflow_home(project, private, environ=None):
    project = Path(project).resolve()
    if _adopted(project):
        checkout = project / SUBMODULE_PATH
        script = checkout / 'scripts/workflow.py'
        if checkout.is_symlink() or not script.is_file() or script.is_symlink():
            raise ValueError('Initialize the pinned workflow: git submodule update --init external/research-workflows')
        lock_path = project / 'workflow.lock.json'
        try:
            lock = json.loads(lock_path.read_text())
            expected_commit = lock['workflow']['commit']
        except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
            raise ValueError('Project workflow.lock.json has no valid workflow commit pin') from error
        if not isinstance(expected_commit, str) or not re.fullmatch(r'[0-9a-f]{40}', expected_commit):
            raise ValueError('Project workflow.lock.json has an invalid workflow commit pin')
        result = subprocess.run(['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
                                capture_output=True, text=True, check=False)
        if result.returncode or result.stdout.strip() != expected_commit:
            raise ValueError('Initialized research-workflows commit differs from workflow.lock.json')
        return checkout.resolve()

    environment = os.environ if environ is None else environ
    configured = environment.get('WORKFLOWS_HOME') or private.get('workflow_home')
    if not configured:
        return None
    home = Path(configured).expanduser()
    if not home.is_absolute() or home.is_symlink() or not (home / 'scripts/workflow.py').is_file():
        raise ValueError('Shared workflow unavailable at configured host binding.')
    return home.resolve()
