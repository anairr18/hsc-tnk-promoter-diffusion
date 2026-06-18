#!/usr/bin/env python3
"""Run or prepare sequence predictor ensemble validation.

The script always writes a standardized FASTA and manifest. If prediction
tables already exist, it merges them. Heavy models are optional because Enformer,
Borzoi and ChromBPNet require GPU/runtime-specific installation.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from project_utils import read_sequence_file, validate_seq


def write_fasta(seqs: list[str], out: Path) -> pd.DataFrame:
    rows = []
    with out.open("w") as f:
        for i, seq in enumerate(seqs, 1):
            sid = f"seq_{i:06d}"
            rows.append({"sequence_id": sid, "sequence": seq})
            f.write(f">{sid}\n{seq}\n")
    return pd.DataFrame(rows)


def merge_predictions(index: pd.DataFrame, prediction_files: list[str]) -> pd.DataFrame:
    merged = index.copy()
    for path in prediction_files:
        pred = pd.read_csv(path, sep="\t")
        if "sequence_id" not in pred.columns and "sequence" in pred.columns:
            pred = pred.merge(index, on="sequence", how="left")
        if "sequence_id" not in pred.columns:
            raise SystemExit(f"Prediction file lacks sequence_id or sequence column: {path}")
        prefix = Path(path).stem
        rename = {c: f"{prefix}__{c}" for c in pred.columns if c not in {"sequence_id", "sequence"}}
        pred = pred.rename(columns=rename)
        merged = merged.merge(pred.drop(columns=["sequence"], errors="ignore"), on="sequence_id", how="left")
    return merged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True, help="TXT/TSV/CSV containing 200bp sequences.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--prediction-files", nargs="*", default=[])
    parser.add_argument("--enformer-command", default=None, help="Optional command template with {fasta} and {out}.")
    parser.add_argument("--borzoi-command", default=None, help="Optional command template with {fasta} and {out}.")
    parser.add_argument("--chrombpnet-command", default=None, help="Optional command template with {fasta} and {out}.")
    parser.add_argument("--run-commands", action="store_true")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    seqs = [s for s in read_sequence_file(args.sequences) if validate_seq(s)]
    index = write_fasta(seqs, out / "candidate_sequences.fa")
    index.to_csv(out / "candidate_sequence_index.tsv", sep="\t", index=False)

    commands = {
        "enformer": args.enformer_command,
        "borzoi": args.borzoi_command,
        "chrombpnet": args.chrombpnet_command,
    }
    manifest = {
        "sequence_count": len(index),
        "fasta": str(out / "candidate_sequences.fa"),
        "expected_prediction_schema": "TSV with sequence_id or sequence plus numeric score columns",
        "commands": commands,
    }
    (out / "predictor_ensemble_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    generated_prediction_files = []
    if args.run_commands:
        for name, template in commands.items():
            if not template:
                continue
            pred_out = out / f"{name}_predictions.tsv"
            cmd = template.format(fasta=(out / "candidate_sequences.fa"), out=pred_out)
            subprocess.run(cmd, shell=True, check=True)
            generated_prediction_files.append(str(pred_out))

    all_prediction_files = args.prediction_files + generated_prediction_files
    if all_prediction_files:
        merged = merge_predictions(index, all_prediction_files)
        merged.to_csv(out / "merged_predictor_scores.tsv", sep="\t", index=False)
        print(f"Wrote merged scores to {out / 'merged_predictor_scores.tsv'}")
    else:
        print(f"Wrote predictor input FASTA and manifest to {out}")


if __name__ == "__main__":
    main()
