"""Shared helpers for the DNA-Diffusion HSC/T-NK project scripts."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

DNA_RE = re.compile(r"^[ACGTacgt]{200}$")
BASES = "ACGT"
RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_json(path: str | Path, data) -> None:
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def validate_seq(seq: object, length: int = 200) -> bool:
    seq = str(seq).strip()
    return len(seq) == length and bool(DNA_RE.match(seq))


def normalize_seq(seq: object) -> str:
    return str(seq).strip().upper()


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def gc_content(seq: str) -> float:
    seq = normalize_seq(seq)
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else math.nan


def cpg_count(seq: str) -> int:
    return normalize_seq(seq).count("CG")


def kmer_counts(seqs: Iterable[str], k: int = 3) -> Counter[str]:
    counts: Counter[str] = Counter()
    for seq in seqs:
        seq = normalize_seq(seq)
        for i in range(0, len(seq) - k + 1):
            mer = seq[i : i + k]
            if set(mer) <= set(BASES):
                counts[mer] += 1
    return counts


def jensen_shannon_distance(a: Counter[str], b: Counter[str]) -> float:
    keys = sorted(set(a) | set(b))
    if not keys:
        return math.nan
    total_a = sum(a.values())
    total_b = sum(b.values())
    if total_a == 0 or total_b == 0:
        return math.nan
    pa = [a[k] / total_a for k in keys]
    pb = [b[k] / total_b for k in keys]
    pm = [(x + y) / 2 for x, y in zip(pa, pb)]

    def kl(p, q):
        return sum(x * math.log2(x / y) for x, y in zip(p, q) if x > 0 and y > 0)

    return math.sqrt((kl(pa, pm) + kl(pb, pm)) / 2)


def read_sequence_file(path: str | Path, seq_col: str | None = None) -> list[str]:
    """Read 200bp DNA sequences from plain text, TSV, or CSV-like files."""
    path = Path(path)
    text = path.read_text(errors="replace").splitlines()
    if not text:
        return []

    # Fast path: one sequence per line.
    one_per_line = [normalize_seq(line) for line in text if validate_seq(line)]
    if len(one_per_line) >= max(1, len(text) // 2):
        return one_per_line

    delimiter = "\t" if "\t" in text[0] else ","
    rows = list(csv.reader(text, delimiter=delimiter))
    if not rows:
        return []

    header = rows[0]
    start = 1
    col_idx = None
    if seq_col and seq_col in header:
        col_idx = header.index(seq_col)
    elif any(h.lower() == "sequence" for h in header):
        col_idx = [h.lower() for h in header].index("sequence")
    else:
        start = 0

    seqs: list[str] = []
    if col_idx is not None:
        for row in rows[start:]:
            if col_idx < len(row) and validate_seq(row[col_idx]):
                seqs.append(normalize_seq(row[col_idx]))
    else:
        for row in rows[start:]:
            for value in row:
                if validate_seq(value):
                    seqs.append(normalize_seq(value))
                    break
    return seqs


def infer_cell_type(name: str) -> str | None:
    lowered = name.lower()
    if "k562" in lowered:
        return "K562"
    if "hepg2" in lowered:
        return "HepG2"
    if "gm12878" in lowered:
        return "GM12878"
    if "hesc" in lowered:
        return "hESCT0"
    if "nk" in lowered:
        return "NK"
    if re.search(r"\bt[\W_]?cell|\bcd[48]\b|\btn\b", lowered):
        return "T"
    return None


def sequence_summary(seqs: Iterable[str]) -> dict[str, float | int]:
    seq_list = [normalize_seq(s) for s in seqs if validate_seq(s)]
    n = len(seq_list)
    unique = len(set(seq_list))
    if n == 0:
        return {
            "n": 0,
            "unique": 0,
            "duplicate_rate": math.nan,
            "gc_mean": math.nan,
            "gc_sd": math.nan,
            "cpg_mean": math.nan,
        }
    gcs = [gc_content(s) for s in seq_list]
    cpgs = [cpg_count(s) for s in seq_list]
    gc_mean = sum(gcs) / n
    gc_sd = math.sqrt(sum((x - gc_mean) ** 2 for x in gcs) / n)
    return {
        "n": n,
        "unique": unique,
        "duplicate_rate": 1 - unique / n,
        "gc_mean": gc_mean,
        "gc_sd": gc_sd,
        "cpg_mean": sum(cpgs) / n,
    }
