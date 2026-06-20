#!/usr/bin/env python3
"""Scan candidate sequences with simple JASPAR/HOCOMOCO-style PFM files."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import pandas as pd

from project_utils import read_sequence_file, revcomp, validate_seq


BASES = "ACGT"


def parse_numbers(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", text)]


def parse_pfm(path: Path) -> dict[str, dict[str, list[float]]]:
    motifs: dict[str, dict[str, list[float]]] = {}
    current = None
    rows: dict[str, list[float]] = {}
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current and set(rows) == set(BASES):
                motifs[current] = rows
            current = line[1:].strip().split()[0]
            rows = {}
            continue
        if current and line[0].upper() in BASES:
            rows[line[0].upper()] = parse_numbers(line)
    if current and set(rows) == set(BASES):
        motifs[current] = rows
    return motifs


def pfm_to_pwm(pfm: dict[str, list[float]], pseudocount: float = 0.1) -> list[dict[str, float]]:
    length = min(len(v) for v in pfm.values())
    pwm = []
    for i in range(length):
        col = {b: pfm[b][i] + pseudocount for b in BASES}
        total = sum(col.values())
        pwm.append({b: math.log2((col[b] / total) / 0.25) for b in BASES})
    return pwm


def best_pwm_score(seq: str, pwm: list[dict[str, float]]) -> float:
    k = len(pwm)
    best = float("-inf")
    for scan_seq in [seq, revcomp(seq)]:
        for i in range(0, len(scan_seq) - k + 1):
            window = scan_seq[i : i + k]
            if set(window) <= set(BASES):
                score = sum(pwm[j][base] for j, base in enumerate(window))
                best = max(best, score)
    return best


def load_sequences(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".csv"}:
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(path, sep=sep)
        if "sequence" in df.columns:
            keep = [c for c in ["candidate_id", "sequence", "rank_score"] if c in df.columns]
            return df[keep].copy()
    return pd.DataFrame({"sequence": read_sequence_file(path)})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--motifs", required=True, nargs="+")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--score-threshold", type=float, default=6.0)
    args = parser.parse_args()

    seqs = load_sequences(Path(args.sequences))
    seqs["sequence"] = seqs["sequence"].astype(str).str.upper().str.strip()
    seqs = seqs[seqs["sequence"].map(validate_seq)].drop_duplicates("sequence").reset_index(drop=True)
    if "candidate_id" not in seqs.columns:
        seqs.insert(0, "candidate_id", [f"seq_{i+1:06d}" for i in range(len(seqs))])

    motif_pwms = []
    for item in args.motifs:
        for motif_id, pfm in parse_pfm(Path(item)).items():
            motif_pwms.append((motif_id, pfm_to_pwm(pfm)))
    if not motif_pwms:
        raise SystemExit("No motifs parsed. Expected JASPAR-style PFM records.")

    rows = []
    summary_rows = []
    for motif_id, pwm in motif_pwms:
        scores = [best_pwm_score(seq, pwm) for seq in seqs["sequence"]]
        hit_count = sum(score >= args.score_threshold for score in scores)
        summary_rows.append({"motif_id": motif_id, "hit_count": hit_count, "sequence_count": len(seqs), "hit_rate": hit_count / len(seqs) if len(seqs) else 0})
        for candidate_id, seq, score in zip(seqs["candidate_id"], seqs["sequence"], scores):
            if score >= args.score_threshold:
                rows.append({"candidate_id": candidate_id, "sequence": seq, "motif_id": motif_id, "best_score": score})

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_prefix.with_suffix(".motif_hits.tsv"), sep="\t", index=False)
    pd.DataFrame(summary_rows).sort_values("hit_rate", ascending=False).to_csv(out_prefix.with_suffix(".motif_summary.tsv"), sep="\t", index=False)
    print(f"Wrote motif scan outputs with prefix {out_prefix}")


if __name__ == "__main__":
    main()
