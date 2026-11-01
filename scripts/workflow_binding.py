"""Resolve project-pinned workflow code before host-level legacy bindings."""
from __future__ import annotations

import json
import os
from pathlib import Path


VENDOR_PATH = 'vendor/research-workflows'


def resolve_workflow_home(project, private, environ=None):
    """
    Resolve workflow home directory.
    
    Priority order:
    1. Vendored workflow (for public distribution)
    2. WORKFLOWS_HOME environment variable
    3. private['workflow_home'] from local configuration
    """
    project = Path(project).resolve()
    
    # Check for vendored workflow first
    vendor = project / VENDOR_PATH
    script = vendor / 'scripts/workflow.py'
    if vendor.is_dir() and script.is_file() and not script.is_symlink():
        return vendor.resolve()
    
    # Fall back to configured workflow home
    environment = os.environ if environ is None else environ
    configured = environment.get('WORKFLOWS_HOME') or private.get('workflow_home')
    if not configured:
        return None
    home = Path(configured).expanduser()
    if not home.is_absolute() or home.is_symlink() or not (home / 'scripts/workflow.py').is_file():
        raise ValueError('Shared workflow unavailable at configured host binding.')
    return home.resolve()
