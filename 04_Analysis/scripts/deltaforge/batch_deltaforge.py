import argparse
import json
import sys
from pathlib import Path

LIGAND_AI = Path(__file__).resolve().parents[3] / "design" / "ligand_ai"
if str(LIGAND_AI) not in sys.path:
    sys.path.insert(0, str(LIGAND_AI))

from ligandai import LigandAI

from src.ligandai_local.config import api_key_debug_info, get_api_key


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch score PDB files with LigandAI DeltaForge and save JSON outputs."
    )
    parser.add_argument("input_dir", help="Directory containing input PDB files")
    parser.add_argument("output_dir", help="Directory where output JSON files will be saved")
    parser.add_argument(
        "--glob",
        default="*.pdb",
        help="Glob for input files inside input_dir (default: *.pdb)",
    )
    parser.add_argument(
        "--receptor-chain",
        dest="receptor_chains",
        action="append",
        required=True,
        help="Repeat for multiple receptor chains, e.g. --receptor-chain A",
    )
    parser.add_argument(
        "--peptide-chain",
        required=True,
        help="Peptide chain identifier, e.g. B",
    )
    parser.add_argument(
        "--scorer",
        default="auto",
        help="DeltaForge scorer to use (default: auto)",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()

    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    pdb_files = sorted(path for path in input_dir.glob(args.glob) if path.is_file())
    if not pdb_files:
        raise SystemExit(f"No files matched {args.glob!r} in {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    key_info = api_key_debug_info()
    print(
        f"Using LigandAI API key from {key_info['source']} "
        f"(prefix={key_info['prefix']}, length={key_info['length']})"
    )
    print(f"Scoring {len(pdb_files)} file(s) from {input_dir}")

    client = LigandAI(api_key=get_api_key())

    for pdb_file in pdb_files:
        print(f"Scoring {pdb_file.name}...")
        score = client.peptides.score_pdb(
            pdb_file=str(pdb_file),
            receptor_chains=args.receptor_chains,
            peptide_chain=args.peptide_chain,
            scorer=args.scorer,
        )
        output_path = output_dir / f"{pdb_file.stem}.json"
        output_path.write_text(json.dumps(score.model_dump(), indent=2), encoding="utf-8")
        print(f"Saved {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
