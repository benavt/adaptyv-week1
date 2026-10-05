"""Exercise identity collisions, swapped species, mutation checks and refreshes."""
import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_model_registry import Registry, parse_structure
from validate_model_registry import validate

THREE = {"A": "ALA", "D": "ASP", "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "R": "ARG"}


def cif(chains):
    columns = "group_PDB label_atom_id label_comp_id label_asym_id label_seq_id pdbx_PDB_model_num".split()
    text = "data_test\n_entry.date 2026-09-30\nloop_\n" + "\n".join("_atom_site." + k for k in columns) + "\n"
    for chain, sequence in chains.items():
        for i, aa in enumerate(sequence, 1):
            text += f"ATOM CA {THREE[aa]} {chain} {i} 1\n"
    return text + "#\n"


class ModelRegistryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.output = self.root / "registry"
        self.source_rel = "02_Design/foundry/runs/DD1_Human_EGFR_hotspots_09_29/rfd3/outputs/shared_model.cif.gz"
        self.refold_rel = "03_Filtering/Refolding/Putative_DD1_EGFR_D3/DD1_001_Human_shared_model/DD1_001_Human_shared_model_model_0.cif"
        self.put(self.source_rel, gzip.compress(cif({"A": "AAA", "B": "GHIK"}).encode()))
        self.put(self.source_rel.replace(".cif.gz", ".json"), json.dumps({"metrics": {"score": 0.3}, "seed": None, "ckpt_path": "/recorded/model.ckpt", "specification": {"contig": "3,/0,A1-4"}}))
        self.put("01_Staging/Putative_DD1/FASTA/DD1_001_Human_shared_model.fasta", ">binder\nAAA\n>Human EGFR\nGHIK\n")
        self.put("01_Staging/Putative_DD1/FASTA/DD1_001_Mouse_shared_model.fasta", ">binder\nAAA\n>Mouse EGFR\nGHIR\n")
        self.put(self.refold_rel, cif({"A": "AAA", "B": "GHIR"}))
        self.csv("03_Filtering/RFD3_filtering/out/dd1_asp_near_both_his_sequences.csv", [{"run": "human", "design_id": "shared_model", "binder_sequence": "AAA", "structure_path": str(self.root / self.source_rel)}])
        self.csv("03_Filtering/Boltz2_metrics/out/boltz2_dd1_metrics_long.csv", [{"design_id": "shared_model", "source_run": "human", "target": "Mouse", "boltz_dir": str((self.root / self.refold_rel).parent), "iptm": "0.8"}])
        self.put("02_Design/foundry/runs/DD1_Human_protonated_EGFR_hotspots_09_29/rfd3/outputs/shared_model.cif.gz", gzip.compress(cif({"A": "DDD", "B": "GHIK"}).encode()))
        self.put("03_Filtering/JustHIpKA/empty.pdb", "")

    def tearDown(self):
        self.temp.cleanup()

    def put(self, rel, value):
        p = self.root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(value) if isinstance(value, bytes) else p.write_text(value)

    def csv(self, rel, rows):
        p = self.root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    def build(self):
        r = Registry(self.root, self.output)
        with contextlib.redirect_stdout(io.StringIO()):
            r.build()
        return r

    def test_species_identity_statistics_and_run_name_collision(self):
        before = (self.root / self.source_rel).read_bytes()
        r = self.build()
        meta = r.models[r.path_model[self.refold_rel]]
        self.assertEqual(meta["design_id"], "DD1_001")
        self.assertEqual(meta["target_species"], "Mouse")
        self.assertTrue(any("Filename says Human" in warning for warning in meta["warnings"]))
        self.assertEqual(meta["statistics"][0]["values"]["iptm"], "0.8")
        self.assertEqual(r.resolve_design("shared_model"), "DD1_001")
        source = r.models[r.path_model[self.source_rel]]
        self.assertEqual(source["provenance"]["checkpoint"], "/recorded/model.ckpt")
        self.assertIsNone(source["provenance"]["seed"])
        self.assertIn("executed_command", source["provenance"]["unknowns"])
        self.assertEqual((self.root / self.source_rel).read_bytes(), before)
        self.assertEqual(validate(self.root, self.output)["status"], "passed")

    def test_new_tool_environments_excluded_without_dropping_campaign_sources(self):
        self.put("tools/environments/new_runtime/fixture.pdb", "fixture")
        self.put("tools/cache/new_runtime/fixture.cif", "fixture")
        self.put("03_Filtering/campaigns/example/analyses/metric/run_001/source.cif", cif({"A": "AAA", "B": "GHIK"}))
        paths = {p.relative_to(self.root.resolve()).as_posix() for p in Registry(self.root, self.output).discover()}
        self.assertNotIn("tools/environments/new_runtime/fixture.pdb", paths)
        self.assertNotIn("tools/cache/new_runtime/fixture.cif", paths)
        self.assertIn("03_Filtering/campaigns/example/analyses/metric/run_001/source.cif", paths)
        self.assertIn(self.source_rel, paths)
        self.assertIn(self.refold_rel, paths)

    def test_refresh_stability_and_annotations(self):
        first = self.build()
        self.put("registry/annotations.json", json.dumps({"models": {}, "designs": {"DD1_001": {"operator_note": "verified by operator"}}}))
        second = self.build()
        self.assertEqual(set(first.models), set(second.models))
        meta = json.loads((self.output / "designs/DD1_001/metadata.json").read_text())
        self.assertEqual(meta["annotations"]["operator_note"], "verified by operator")
        self.assertEqual(validate(self.root, self.output)["status"], "passed")

    def test_empty_structure_and_changed_source_are_explicit(self):
        r = self.build()
        meta = r.models[r.path_model["03_Filtering/JustHIpKA/empty.pdb"]]
        self.assertEqual(meta["structure_status"], "empty_file")
        self.assertEqual(meta["chains"], [])
        self.put(self.source_rel, gzip.compress(cif({"A": "DDD", "B": "GHIK"}).encode()))
        result = validate(self.root, self.output)
        self.assertEqual(result["status"], "failed")
        self.assertTrue(any("Source content changed" in e for e in result["errors"]))

    def test_mutation_labels_checked_against_sequence(self):
        self.put("01_Staging/DD1_001_A3D_EGFR_species_intersection/DD1_001_A3D_Human_EGFR.fasta", ">binder\nAAD\n>target\nGHIK\n")
        r = self.build()
        self.assertTrue(r.designs["DD1_001_A3D"]["provenance"]["mutations_verified"])
        self.assertEqual(r.designs["DD1_001_A3D"]["provenance"]["observed_mutations"], ["A3D"])

    def test_wrapped_cif_rows_use_full_parser(self):
        self.put("wrapped.cif", "data_test\nloop_\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\nCA ALA\nA 1\n#\n")
        self.assertEqual(parse_structure(self.root / "wrapped.cif")["chains"][0]["sequence"], "A")

    def test_coordinate_cache_checks_current_source_hash(self):
        self.build()
        self.put(self.source_rel, gzip.compress(cif({"A": "DDD", "B": "GHIK"}).encode()))
        registry = Registry(self.root, self.output, reuse_coordinates=True)
        registry.load_designs()
        registry.build_model(self.root / self.source_rel)
        model = registry.models[registry.path_model[self.source_rel]]
        binder = next(c for c in model["chains"] if c["role"] == "binder")
        self.assertEqual(binder["sequence"], "DDD")
        self.assertNotEqual(model["design_id"], "DD1_001")

    def test_target_only_pka_provenance_and_result(self):
        rel = "03_Filtering/JustHIpKA/putative_DD1_EGFR_D3/shared_model_Mouse/egfr.cif"
        self.put(rel, cif({"A": "GHIR"}))
        self.put("03_Filtering/JustHIpKA/putative_DD1_EGFR_D3/shared_model_Mouse/predictions_A_HIS37_HIS100/HIS37.txt", "PKA=6.9 SID=37 RES=HIS\n")
        r = self.build()
        prep = r.models[r.path_model[rel]]
        self.assertEqual(prep["design_id"], "DD1_001")
        self.assertEqual(prep["chains"][0]["role"], "EGFR_D3")
        refold = r.models[r.path_model[self.refold_rel]]
        pka = next(s for s in refold["statistics"] if "JustHIpKA" in s["source_path"])
        self.assertEqual(pka["values"][0]["pKa"], 6.9)
        self.assertIn("target-only", pka["association"])
        self.assertEqual(validate(self.root, self.output)["status"], "passed")


if __name__ == "__main__":
    unittest.main()
