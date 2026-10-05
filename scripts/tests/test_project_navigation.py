"""Verify navigation freshness, status evidence and safe source paths."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_project_navigation import build, load_campaigns, render


class ProjectNavigationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = self.root / "03_Filtering/campaigns/example/manifest.json"
        self.source.parent.mkdir(parents=True)
        self.stage = self.root / "01_Staging/exact_run_001"
        self.stage.mkdir(parents=True)
        (self.root / "01_Staging/README.md").write_text("# Staging\n")
        (self.root / "registry").mkdir()
        (self.root / "registry/README.md").write_text("# Registry\n")
        (self.stage / "submission_status.json").write_text(json.dumps({"status": "running", "completed_models": 0}))
        (self.stage / "RUN_RECORD.md").write_text("# Exact run\n\nRecorded receipt, not live status.\n")
        self.artifact = self.root / "results/folding/exact_run_001/native.cif"
        self.artifact.parent.mkdir(parents=True)
        self.artifact.write_text("data_original\n")
        self.manifest = {"schema_version": 1, "campaign_id": "example", "title": "Example campaign", "summary": "Source-backed navigation.", "locations": [{"role": "Native result", "label": "Original", "path": "results/folding/exact_run_001/native.cif", "evidence_path": "01_Staging/exact_run_001/RUN_RECORD.md"}], "current_outputs": [], "notes": ["No new prediction."], "runs": [{"run_id": "exact_run_001", "phase": "folding_submission", "path": "01_Staging/exact_run_001", "status_evidence": "01_Staging/exact_run_001/submission_status.json", "status_locator": "/status"}]}
        self.save()

    def tearDown(self):
        self.tmp.cleanup()

    def save(self):
        self.source.write_text(json.dumps(self.manifest))

    def test_build_and_read_only_check_preserve_scientific_files(self):
        original = self.artifact.read_bytes()
        build(self.root)
        before = {p: p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
        self.assertEqual(build(self.root, check=True)["status"], "passed")
        self.assertEqual(before, {p: p.read_bytes() for p in self.root.rglob("*") if p.is_file()})
        self.assertEqual(original, self.artifact.read_bytes())
        with (self.root / "registry/campaign_runs.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["run_id"], "exact_run_001")
        self.assertEqual(row["recorded_status"], "running")
        self.assertEqual(row["completed_models"], "0")
        self.assertEqual(row["submitted_requests"], "")

    def test_changed_receipt_is_stale_until_regenerated(self):
        build(self.root)
        (self.stage / "submission_status.json").write_text(json.dumps({"status": "complete", "completed_models": 5}))
        with self.assertRaisesRegex(ValueError, "stale"):
            build(self.root, check=True)
        build(self.root)
        self.assertEqual(build(self.root, check=True)["status"], "passed")
        self.assertIn("complete", (self.root / "registry/campaign_runs.csv").read_text())

    def test_mixed_outcomes_do_not_promote_first_job_to_run_status(self):
        self.manifest["runs"][0]["status_locator"] = "/outcomes"
        self.save()
        (self.stage / "submission_status.json").write_text(json.dumps({"outcomes": [{"status": "complete"}, {"status": "failed_before_prediction"}]}))
        build(self.root)
        with (self.root / "registry/campaign_runs.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["recorded_status"], "mixed_outcomes")

    def test_missing_source_rejected_before_partial_writes(self):
        self.artifact.unlink()
        with self.assertRaisesRegex(ValueError, "Missing navigation reference"):
            build(self.root)
        self.assertFalse((self.root / "registry/campaigns.csv").exists())
        self.assertFalse((self.source.parent / "README.md").exists())

    def test_unavailable_report_is_not_a_current_output(self):
        self.manifest["documented_unavailable"] = [{"expected_path": "results/folding/exact_run_001/promised.pptx", "evidence_path": "01_Staging/exact_run_001/RUN_RECORD.md"}]
        self.save()
        build(self.root)
        self.assertIn("unavailable locally", (self.source.parent / "README.md").read_text())
        (self.artifact.parent / "promised.pptx").write_bytes(b"new report")
        with self.assertRaisesRegex(ValueError, "review authority"):
            build(self.root)

    def test_traversal_and_external_symlink_rejected(self):
        self.manifest["locations"][0]["path"] = "../outside"
        self.save()
        with self.assertRaisesRegex(ValueError, "project-relative"):
            load_campaigns(self.root)
        with tempfile.TemporaryDirectory() as external:
            (self.root / "outside_link").symlink_to(external, target_is_directory=True)
            self.manifest["locations"][0]["path"] = "outside_link"
            self.save()
            with self.assertRaisesRegex(ValueError, "escapes project"):
                load_campaigns(self.root)

    def test_supersession_requires_named_evidence_and_no_cycles(self):
        first = self.manifest["runs"][0]
        first.update(superseded_by="exact_run_002", supersession_evidence="01_Staging/exact_run_001/RUN_RECORD.md")
        second = dict(first, run_id="exact_run_002", superseded_by="exact_run_001")
        self.manifest["runs"].append(second)
        self.save()
        with self.assertRaisesRegex(ValueError, "not named in evidence"):
            render(self.root)
        (self.stage / "RUN_RECORD.md").write_text("exact_run_001 and exact_run_002\n")
        with self.assertRaisesRegex(ValueError, "Supersession cycle"):
            render(self.root)

    def test_missing_local_document_link_rejected(self):
        (self.root / "README.md").write_text("[Missing](01_Staging/nonexistent.fasta)\n")
        with self.assertRaisesRegex(ValueError, "Broken local link"):
            build(self.root)

    def test_generated_symlink_cannot_overwrite_scientific_source(self):
        (self.source.parent / "README.md").symlink_to(self.artifact)
        before = self.artifact.read_bytes()
        with self.assertRaisesRegex(ValueError, "symlink"):
            build(self.root)
        self.assertEqual(before, self.artifact.read_bytes())


if __name__ == "__main__":
    unittest.main()
