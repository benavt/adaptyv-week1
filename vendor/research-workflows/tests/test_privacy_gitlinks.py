"""A submodule is index metadata; its URL and surrounding files remain scanned."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CHECKER = Path(__file__).resolve().parents[1]/'scripts/check_privacy.py'


class GitlinkPrivacyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.git('init', '-q')
        (self.root/'.gitmodules').write_text('[submodule "external/research-workflows"]\n'
            ' path = external/research-workflows\n url = https://github.com/example/research-workflows.git\n')
        self.git('add', '.gitmodules')
        self.git('update-index', '--add', '--cacheinfo', '160000,'+'a'*40+',external/research-workflows')

    def git(self, *args):
        return subprocess.run(['git', *args], cwd=self.root, capture_output=True, check=True)

    def check(self, *args):
        return subprocess.run([sys.executable, str(CHECKER), '--root', str(self.root), *args], capture_output=True)

    def test_gitlink_with_no_local_commit_object_passes_staged_and_tracked_checks(self):
        for args in ((), ('--tracked',)):
            result = self.check(*args)
            self.assertEqual(result.returncode, 0, result.stderr.decode())

    def test_secret_in_submodule_url_still_fails(self):
        file = self.root/'.gitmodules'
        file.write_text(file.read_text().replace('https://github.com/', 'https://github_pat_'+'x'*35+'@github.com/'))
        self.git('add', '.gitmodules')
        result = self.check()
        self.assertEqual(result.returncode, 1)
        self.assertIn(b'.gitmodules', result.stderr)

    def test_private_gitlink_path_still_fails(self):
        self.git('update-index', '--add', '--cacheinfo', '160000,'+'a'*40+',local/private-workflow')
        result = self.check()
        self.assertEqual(result.returncode, 1)
        self.assertIn(b'local/private-workflow', result.stderr)


if __name__ == '__main__':
    unittest.main()
