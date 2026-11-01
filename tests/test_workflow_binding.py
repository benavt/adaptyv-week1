import json
from pathlib import Path

import pytest

from scripts.workflow_binding import resolve_workflow_home


def init_project_with_vendored_workflow(tmp_path):
    """Initialize project with vendored workflow (new public distribution model)."""
    project = tmp_path / 'project'
    vendor = project / 'vendor/research-workflows'
    (vendor / 'scripts').mkdir(parents=True)
    (vendor / 'scripts/workflow.py').write_text('# vendored workflow\n')
    return project, vendor


def test_vendored_workflow_takes_priority_over_host_config(tmp_path):
    """Vendored workflow is preferred for public distribution."""
    project, vendor = init_project_with_vendored_workflow(tmp_path)
    resolved = resolve_workflow_home(project, {'workflow_home': '/legacy/workflow'},
                                     {'WORKFLOWS_HOME': '/other/workflow'})
    assert resolved == vendor.resolve()


def test_fallback_to_environment_and_private_workflow_home(tmp_path):
    """When no vendored workflow exists, fall back to configured workflow home."""
    env_home = tmp_path / 'environment'
    private_home = tmp_path / 'private'
    for home in (env_home, private_home):
        (home / 'scripts').mkdir(parents=True)
        (home / 'scripts/workflow.py').write_text('# host workflow\n')
    project = tmp_path / 'project'
    project.mkdir()
    
    # Environment variable takes priority
    assert resolve_workflow_home(project, {'workflow_home': str(private_home)},
                                 {'WORKFLOWS_HOME': str(env_home)}) == env_home.resolve()
    # Private config used when no environment variable
    assert resolve_workflow_home(project, {'workflow_home': str(private_home)}, {}) == private_home.resolve()
    # None when neither configured
    assert resolve_workflow_home(project, {}, {}) is None
