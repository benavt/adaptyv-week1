"""Transport contract tests; all transfers stay in temporary local fixtures."""

import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import shlex
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("remote_sync", SCRIPTS / "remote_sync.py")
sync = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sync)


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.index = Path(self.tmp.name) / "results.jsonl"
        self.root = PurePosixPath("/fixture/campaign")

    def write(self, rows):
        self.index.write_text("\n".join(json.dumps(row) for row in rows) + "\n")

    def test_preserves_shards_sidecars_and_filters(self):
        artifact = str(self.root / "shards/007/folding/native/model.cif")
        row = {"status": "complete", "rank": 1, "engine": "boltz2", "target_species": "human",
               "structure_path": artifact, "pae_path": artifact + ".npz"}
        self.write([row, row, {**row, "rank": 2}, {**row, "status": "missing_output"},
                    {**row, "target_species": "mouse"}, {"status": "invalid_input"}])
        self.assertEqual(sync.selected_files(self.index, self.root, engine="boltz2", species="human"),
                         ["shards/007/folding/native/model.cif", "shards/007/folding/native/model.cif.npz"])

    def test_rejects_unsafe_and_external_paths(self):
        for path in ("/fixture/other/model.cif", "/fixture/campaign/../other.cif",
                     "/fixture/campaign/a\nb.cif", "/fixture/campaign", "relative.cif"):
            with self.subTest(path=path):
                self.write([{"status": "complete", "rank": 1, "structure_path": path}])
                with self.assertRaises(ValueError):
                    sync.selected_files(self.index, self.root)

    def test_missing_rank_does_not_match(self):
        self.write([{"status": "complete", "structure_path": "/fixture/campaign/model.cif"}])
        self.assertEqual(sync.selected_files(self.index, self.root), [])

    def test_private_file_list_rejects_absolute_and_traversal(self):
        for value in (b"../outside.cif\n", b"/absolute.cif\0", b".\n"):
            self.index.write_bytes(value)
            with self.assertRaises(ValueError):
                sync.read_file_list(self.index)

    def test_config_rejects_host_injection_and_path_traversal(self):
        for alias in ("-oProxyCommand=bad", "user@host", "some.host", "host;bad"):
            with self.subTest(alias=alias), patch.dict(os.environ, {"REMOTE": alias}):
                with self.assertRaises(ValueError):
                    sync.Transport()
        for path in ("/", "relative", "/fixture/../archive", "/fixture\narchive"):
            with self.subTest(path=path), patch.dict(os.environ, {"REMOTE_INPUT_ROOT": path}):
                with self.assertRaises(ValueError):
                    sync.remote_path("REMOTE_INPUT_ROOT")

    def test_modern_rsync_gets_raw_protected_path(self):
        with patch.dict(os.environ, {"REMOTE": "fixture"}):
            transport = sync.Transport()
        transport.protected_args = True
        path = "/fixture/space's (2)"
        self.assertEqual(transport.operand(path), "fixture:" + path)

    def test_project_staging_parents_preserve_run_identity(self):
        with tempfile.TemporaryDirectory() as project:
            for parent in ("01_Staging", "01_Staging/runs"):
                with self.subTest(parent=parent), patch.dict(os.environ, {"PROJECT_ROOT": project, "LOCAL_INPUT_ROOT": parent}):
                    run_id = sync.identifier("2026-10-04_exact_run_001")
                    staged = sync.local_root("LOCAL_INPUT_ROOT", "inputs") / run_id
                    self.assertEqual(staged, Path(project).resolve() / parent / run_id)


@unittest.skipUnless(shutil.which("rsync"), "rsync is required for local transport fixtures")
class LocalRsyncTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.bin = self.root / "bin"
        self.bin.mkdir()
        # Only the SSH network hop is replaced. Both rsync processes are real.
        ssh = self.bin / "ssh"
        ssh.write_text(
            "#!/usr/bin/env python3\n"
            "import os, sys\n"
            "args = sys.argv[1:]\n"
            "while args and args[0] == '-o': args = args[2:]\n"
            "args = args[1:]\n"
            "os.execv('/bin/sh', ['/bin/sh', '-c', ' '.join(args)])\n"
        )
        ssh.chmod(0o755)
        self.remote = self.root / "remote space's $(false) (2)"
        self.remote.mkdir()
        self.inputs = self.root / "local inputs"
        self.results = self.root / "local results"
        self.logs = self.root / "local logs"
        values = {
            "REMOTE": "fixture", "REMOTE_INPUT_ROOT": str(self.remote / "incoming"),
            "REMOTE_FOLDING_ROOT": str(self.remote / "runs"),
            "REMOTE_CAMPAIGN_ROOT": str(self.remote / "runs"),
            "REMOTE_FOUNDRY_ROOT": str(self.remote / "foundry"),
            "REMOTE_LOG_ROOT": str(self.remote / "logs"),
            "REMOTE_RFD3_ROOT": "", "REMOTE_BUNDLE": str(self.remote / "bundle (2)"),
            "REMOTE_CMD": "exit 7", "LOCAL_INPUT_ROOT": str(self.inputs),
            "LOCAL_RESULT_ROOT": str(self.results), "LOCAL_LOG_ROOT": str(self.logs),
        }
        self.config = self.root / "test.env.local"
        self.config.write_text("\n".join(f"{key}={shlex.quote(value)}" for key, value in values.items()) + "\n")
        self.env = {**os.environ, "PATH": str(self.bin) + os.pathsep + os.environ["PATH"],
                    "SYNC_ENV_FILE": str(self.config)}

    def run_sync(self, *args, success=True):
        result = subprocess.run([str(SCRIPTS / "sync.sh"), *args], env=self.env,
                                cwd=self.root, text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_upload_preview_real_upload_and_symlink_rejection(self):
        staged = self.inputs / "experiment/folding/molecules"
        staged.mkdir(parents=True)
        (staged / "input.fasta").write_text(">A\nMKTAYI\n")
        (self.remote / "incoming").mkdir()
        self.run_sync("push", "experiment", "--dry-run")
        self.assertFalse((self.remote / "incoming/experiment").exists())
        self.run_sync("push", "experiment")
        self.assertEqual((self.remote / "incoming/experiment/folding/molecules/input.fasta").read_text(), ">A\nMKTAYI\n")
        (staged / "linked.fasta").symlink_to(staged / "input.fasta")
        result = self.run_sync("push", "experiment", success=False)
        self.assertIn("symlinks", result.stderr)

    def test_upload_rejects_linked_remote_ancestor(self):
        staged = self.inputs / "experiment"
        staged.mkdir(parents=True)
        (staged / "input.fasta").write_text("test")
        archive = self.remote / "archive"
        archive.mkdir()
        (self.remote / "incoming").symlink_to(archive, target_is_directory=True)
        self.run_sync("push", "experiment", success=False)
        self.assertFalse((archive / "experiment").exists())

    def test_prepare_only_creates_staging_destination(self):
        self.assertFalse((self.remote / "incoming").exists())
        self.run_sync("prepare", "experiment")
        destination = self.remote / "incoming/experiment"
        self.assertTrue(destination.is_dir())
        self.assertEqual(list(destination.iterdir()), [])
        self.assertFalse((self.remote / "runs/experiment").exists())

    def test_private_file_list_limits_foundry_download(self):
        output = self.remote / "foundry/experiment/mpnn/outputs"
        output.mkdir(parents=True)
        (output / "selected.cif").write_text("selected")
        (output / "other.cif").write_text("other")
        manifest = self.root / "files.txt"
        manifest.write_bytes(b"mpnn/outputs/selected.cif\0")
        self.run_sync("pull", "foundry", "experiment", "--files-from", str(manifest))
        local = self.results / "foundry/experiment/mpnn/outputs"
        self.assertEqual((local / "selected.cif").read_text(), "selected")
        self.assertFalse((local / "other.cif").exists())

    def test_sharded_metadata_and_selected_sidecars(self):
        campaign = self.remote / "runs/campaign"
        artifact = campaign / "shards/007/folding/native/#model space.cif"
        artifact.parent.mkdir(parents=True)
        artifact.write_text("coordinates")
        sidecar = artifact.with_suffix(".npz")
        sidecar.write_bytes(b"pae")
        (artifact.parent / "unselected.cif").write_text("unused")
        index = campaign / "folding/results.jsonl"
        index.parent.mkdir()
        index.write_text(json.dumps({"status": "complete", "rank": 1,
                                     "structure_path": str(artifact), "pae_path": str(sidecar)}) + "\n")
        (artifact.parent / "state.json").write_text('{"status":"complete"}')
        self.run_sync("pull", "campaign", "campaign", "--metadata")
        local = self.results / "campaigns/campaign"
        self.assertTrue((local / "folding/results.jsonl").exists())
        self.assertTrue((local / "shards/007/folding/native/state.json").exists())
        self.assertFalse((local / artifact.relative_to(campaign)).exists())
        self.run_sync("pull", "campaign", "campaign", "--selected", "--dry-run")
        self.assertFalse((local / artifact.relative_to(campaign)).exists())
        self.run_sync("pull", "campaign", "campaign", "--selected")
        self.assertEqual((local / artifact.relative_to(campaign)).read_text(), "coordinates")
        self.assertEqual((local / sidecar.relative_to(campaign)).read_bytes(), b"pae")
        self.assertFalse((local / "shards/007/folding/native/unselected.cif").exists())

    def test_foundry_link_materialization_and_compressed_export(self):
        run = self.remote / "foundry/experiment"
        output = run / "rfd3/outputs"
        output.mkdir(parents=True)
        (output / "design.cif.gz").write_bytes(b"compressed")
        (output / "design.cif").write_text("expanded")
        (output / "design.json").write_text("{}")
        log = self.remote / "external.out"
        log.write_text("run log")
        (run / "slurm.out").symlink_to(log)
        self.run_sync("pull", "rfd3", "experiment")
        compressed = self.results / "foundry/experiment/rfd3_compressed"
        self.assertTrue((compressed / "design.cif.gz").exists())
        self.assertTrue((compressed / "design.json").exists())
        self.assertFalse((compressed / "design.cif").exists())
        self.run_sync("pull", "foundry", "experiment")
        local_log = self.results / "foundry/experiment/slurm.out"
        self.assertFalse(local_log.is_symlink())
        self.assertEqual(local_log.read_text(), "run log")

    def test_checksum_finds_same_size_same_mtime_change(self):
        source = self.remote / "foundry/experiment/model.cif"
        source.parent.mkdir(parents=True)
        source.write_text("original")
        self.run_sync("pull", "foundry", "experiment")
        local = self.results / "foundry/experiment/model.cif"
        timestamp = source.stat().st_mtime_ns
        source.write_text("modified")
        os.utime(source, ns=(timestamp, timestamp))
        self.run_sync("pull", "foundry", "experiment", "--checksum")
        self.assertEqual(local.read_text(), "modified")

    def test_bundle_and_log_keep_exact_paths(self):
        bundle = self.remote / "bundle (2)"
        (bundle / "provenance").mkdir(parents=True)
        (bundle / "provenance/receipt.json").write_text("{}")
        self.run_sync("pull", "bundle", "archive_label")
        self.assertTrue((self.results / "bundles/archive_label/provenance/receipt.json").exists())
        (self.remote / "logs").mkdir()
        (self.remote / "logs/folding-12345.out").write_text("job log")
        self.run_sync("log", "12345")
        self.assertEqual((self.logs / "folding-12345.out").read_text(), "job log")

    def test_failure_status_propagates(self):
        self.assertEqual(self.run_sync("run", success=False).returncode, 7)
        result = self.run_sync("pull", "foundry", "missing", success=False)
        self.assertIn("Sync failed: rsync exited", result.stderr)


if __name__ == "__main__":
    unittest.main()
