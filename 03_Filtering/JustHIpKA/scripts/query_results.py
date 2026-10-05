#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--csv", type=Path, default=Path(__file__).parents[1] / "hispka_results.csv")
p.add_argument("--model", required=True); p.add_argument("--chain", required=True)
p.add_argument("--res-name", required=True); p.add_argument("--res-num", required=True)
a = p.parse_args()
with a.csv.open(newline="") as f:
    rows = csv.DictReader(f)
    for row in rows:
        if (row["model"], row["chain_id"], row["res_name"], row["res_num"]) == (a.model, a.chain, a.res_name.upper(), a.res_num):
            print(row)
