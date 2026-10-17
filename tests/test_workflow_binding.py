import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from scripts.workflow_binding import resolve_workflow_home


def init_project_with_workflow(tmp_path):
    project = tmp_path / 'project'
    checkout = project / 'external/research-workflows'
    (checkout / 'scripts').mkdir(parents=True)
    (checkout / 'scripts/workflow.py').write_text('# pinned fixture\n')
    subprocess.run(['git', 'init', '-q', str(checkout)], check=True)
    subprocess.run(['git', '-C', str(checkout), 'config', 'user.name', 'Test'], check=True)
    subprocess.run(['git', '-C', str(checkout), 'config', 'user.email', 'test@example.invalid'], check=True)
    subprocess.run(['git', '-C', str(checkout), 'add', 'scripts/workflow.py'], check=True)
    subprocess.run(['git', '-C', str(checkout), 'commit', '-qm', 'fixture'], check=True)
    commit = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    (project / '.gitmodules').write_text(
        '[submodule "external/research-workflows"]\n'
        '\tpath = external/research-workflows\n'
        '\turl = https://github.com/benavt/research-workflows.git\n')
    (project / 'workflow.lock.json').write_text(json.dumps({'workflow': {'commit': commit}}))
    return project, checkout, commit


def test_adopted_project_uses_exact_project_pin_and_ignores_host_override(tmp_path):
    project, checkout, _ = init_project_with_workflow(tmp_path)
    resolved = resolve_workflow_home(project, {'workflow_home': '/legacy/workflow'},
                                     {'WORKFLOWS_HOME': '/other/workflow'})
    assert resolved == checkout.resolve()


def test_adopted_project_fails_closed_when_submodule_is_missing(tmp_path):
    project = tmp_path / 'project'
    project.mkdir()
    (project / '.gitmodules').write_text(
        '[submodule "external/research-workflows"]\npath = external/research-workflows\n')
    with pytest.raises(ValueError, match='submodule update --init'):
        resolve_workflow_home(project, {'workflow_home': '/legacy/workflow'},
                              {'WORKFLOWS_HOME': '/other/workflow'})


def test_adopted_project_rejects_commit_lock_mismatch(tmp_path):
    project, _, _ = init_project_with_workflow(tmp_path)
    lock = json.loads((project / 'workflow.lock.json').read_text())
    lock['workflow']['commit'] = '0' * 40
    (project / 'workflow.lock.json').write_text(json.dumps(lock))
    with pytest.raises(ValueError, match='differs from workflow.lock'):
        resolve_workflow_home(project, {})


def test_unadopted_project_preserves_environment_and_private_fallbacks(tmp_path):
    env_home = tmp_path / 'environment'
    private_home = tmp_path / 'private'
    for home in (env_home, private_home):
        (home / 'scripts').mkdir(parents=True)
        (home / 'scripts/workflow.py').write_text('# legacy\n')
    project = tmp_path / 'unadopted'
    project.mkdir()
    assert resolve_workflow_home(project, {'workflow_home': str(private_home)},
                                 {'WORKFLOWS_HOME': str(env_home)}) == env_home.resolve()
    assert resolve_workflow_home(project, {'workflow_home': str(private_home)}, {}) == private_home.resolve()
    assert resolve_workflow_home(project, {}, {}) is None
