"""Species must follow receptor sequence, independently of source-design naming."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

from Bio.PDB import Atom, Chain, Model, Residue, Structure, MMCIFIO
import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "DD1_ASP_EGFR_D3_HIS_proximity_filter.py"
spec = importlib.util.spec_from_file_location("proximity", SCRIPT)
proximity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(proximity)
TARGETS = {"human": "AHH", "mouse": "VHH"}


def write_complex(path, target_residue):
    structure = Structure.Structure("test")
    model = Model.Model(0)
    structure.add(model)
    for chain_id, residues in (("A", [(1, "ASP")]), ("B", [(1, target_residue), (37, "HIS"), (100, "HIS")])):
        chain = Chain.Chain(chain_id)
        model.add(chain)
        for number, name in residues:
            residue = Residue.Residue((" ", number, " "), name, " ")
            chain.add(residue)
            residue.add(Atom.Atom("CA", np.array([1.0 if chain_id == "B" else 0.0, 0.0, 0.0]), 1.0, 1.0, " ", " CA ", number, element="C"))
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(path))


class RefoldSpeciesTest(unittest.TestCase):
    def test_human_target_from_mouse_design(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "DD1_001_Human_Mouse_EGFR_Mouse_EGFR_1_model_0" / "model_0.cif"
            write_complex(path, "ALA")
            row = proximity.analyze_refolded(path, 5.0, TARGETS)
            self.assertEqual(row["species"], "human")
            self.assertEqual(row["design_id"], "DD1_001")
            self.assertEqual(row["asp_near_his37"], "1")
            self.assertEqual(row["asp_near_his100"], "1")
            self.assertTrue(row["has_asp_near_both"])

    def test_mouse_target_from_human_design(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "DD1_001_Mouse_Human_EGFR_Human_EGFR_1_model_0" / "model_0.cif"
            write_complex(path, "VAL")
            self.assertEqual(proximity.analyze_refolded(path, 5.0, TARGETS)["species"], "mouse")

    def test_species_and_geometry_do_not_depend_on_folder_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = []
            for folder in ("DD1_001_Human_Mouse_EGFR", "DD1_001_Mouse_Human_EGFR", "DD1_001_no_species_label"):
                path = Path(tmp) / folder / "model_0.cif"
                write_complex(path, "ALA")
                rows.append(proximity.analyze_refolded(path, 5.0, TARGETS))
            for row in rows[1:]:
                self.assertEqual({k: v for k, v in row.items() if k not in {"model_id", "structure_path"}}, {k: v for k, v in rows[0].items() if k not in {"model_id", "structure_path"}})

    def test_unknown_or_ambiguous_target_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "DD1_001_Human_label" / "model_0.cif"
            write_complex(path, "ALA")
            for targets in ({"human": "VHH", "mouse": "GHH"}, {"human": "AHH", "mouse": "AHH"}):
                with self.assertRaisesRegex(ValueError, "cannot assign species"):
                    proximity.analyze_refolded(path, 5.0, targets)


if __name__ == "__main__":
    unittest.main()
