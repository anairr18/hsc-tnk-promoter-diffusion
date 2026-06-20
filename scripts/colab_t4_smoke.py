#!/usr/bin/env python3
"""Fast Colab/T4 smoke test for the DNA-Diffusion training/sampling stack.

This intentionally does not build the biological ENCODE dataset. It creates a
tiny synthetic 200bp three-label dataset, verifies uv/DNA-Diffusion setup,
patches pretrained fine-tuning if requested, trains for a few epochs, samples a
few sequences, and writes a smoke-test report. Use this before spending A100
compute units on the full pipeline.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import random
import shutil
import subprocess
from pathlib import Path

import pandas as pd


CELL_TYPES = ["K562", "HepG2", "GM12878"]
MOTIFS = {"K562": "GATA", "HepG2": "TGTTT", "GM12878": "TTTGCAT"}
BASES = "ACGT"


def run(cmd: list[str], cwd: Path | None = None, env: dict | None = None) -> None:
    print(">", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def make_seq(label: str, rng: random.Random) -> str:
    seq = [rng.choice(BASES) for _ in range(200)]
    motif = MOTIFS[label]
    for pos in [35, 90, 145]:
        seq[pos : pos + len(motif)] = list(motif)
    return "".join(seq)


def ensure_uv() -> str:
    uv = shutil.which("uv")
    if uv:
        return uv
    run(["python", "-m", "pip", "install", "-q", "uv"])
    uv = shutil.which("uv")
    if not uv:
        raise RuntimeError("uv install did not put uv on PATH")
    return uv


def ensure_repo(repo: Path) -> None:
    if not (repo / ".git").exists():
        if repo.exists():
            shutil.rmtree(repo)
        run(["git", "clone", "https://github.com/pinellolab/DNA-Diffusion.git", str(repo)])


def create_smoke_data(repo: Path, n_per_label: int, seed: int) -> Path:
    rng = random.Random(seed)
    rows = []
    for label in CELL_TYPES:
        for i in range(n_per_label):
            rows.append({"chr": f"chr{(i % 6) + 3}", "sequence": make_seq(label, rng), "TAG": label})
    df = pd.DataFrame(rows).sample(frac=1, random_state=seed).reset_index(drop=True)
    data_dir = repo / "data" / "hsc_tnk_smoke"
    data_dir.mkdir(parents=True, exist_ok=True)
    train = df.sample(frac=0.8, random_state=seed)
    rest = df.drop(train.index)
    val = rest.sample(frac=0.5, random_state=seed)
    test = rest.drop(val.index)
    train[["chr", "sequence", "TAG"]].to_csv(data_dir / "train.txt", sep="\t", index=False)
    val[["chr", "sequence", "TAG"]].to_csv(data_dir / "val.txt", sep="\t", index=False)
    test[["chr", "sequence", "TAG"]].to_csv(data_dir / "test.txt", sep="\t", index=False)
    with (data_dir / "encode_data.pkl").open("wb") as f:
        pickle.dump(
            {
                "train_df": train[["chr", "sequence", "TAG"]].reset_index(drop=True),
                "validation_df": val[["chr", "sequence", "TAG"]].reset_index(drop=True),
                "test_df": test[["chr", "sequence", "TAG"]].reset_index(drop=True),
            },
            f,
        )
    return data_dir


def patch_pretrained(project_root: Path, repo: Path) -> None:
    patcher = project_root / "scripts" / "patch_dnadiffusion_pretrained_finetune.py"
    run(["python", str(patcher), "--dnadiffusion-repo", str(repo)])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", default=Path(__file__).resolve().parents[1])
    p.add_argument("--dnadiffusion-repo", default=Path.cwd().parent / "DNA-Diffusion")
    p.add_argument("--epochs", type=int, default=3)
    p.add_argument("--min-epochs", type=int, default=1)
    p.add_argument("--n-per-label", type=int, default=24)
    p.add_argument("--samples-per-label", type=int, default=10)
    p.add_argument("--patch-pretrained", action="store_true")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    project_root = Path(args.project_root)
    repo = Path(args.dnadiffusion_repo)
    uv = ensure_uv()
    ensure_repo(repo)
    if args.patch_pretrained:
        patch_pretrained(project_root, repo)
    data_dir = create_smoke_data(repo, args.n_per_label, args.seed)

    run([uv, "sync"], cwd=repo)
    env = os.environ.copy()
    env["WANDB_DISABLED"] = "true"
    env["PYTHONUNBUFFERED"] = "1"
    model_override = ["model=unet_pretrained"] if args.patch_pretrained else []
    run(
        [
            uv,
            "run",
            "python",
            "-u",
            "train.py",
            *model_override,
            f"data.data_path={data_dir / 'train.txt'}",
            "data.load_saved_data=True",
            f"data.saved_data_path={data_dir / 'encode_data.pkl'}",
            "training.use_wandb=False",
            f"training.num_epochs={args.epochs}",
            f"training.min_epochs={args.min_epochs}",
            "training.patience=2",
            "training.batch_size=12",
            "training.sample_epoch=999999",
        ],
        cwd=repo,
        env=env,
    )

    ckpts = sorted((repo / "checkpoints").glob("*.pt"), key=lambda x: x.stat().st_mtime)
    if not ckpts:
        raise RuntimeError("Smoke training finished but produced no checkpoint")
    latest = ckpts[-1]
    run(
        [
            uv,
            "run",
            "python",
            "-u",
            "sample.py",
            f"sampling.checkpoint_path={latest}",
            "data.load_saved_data=True",
            f"data.saved_data_path={data_dir / 'encode_data.pkl'}",
            "data.cell_types=K562,HepG2,GM12878",
            f"sampling.number_of_samples={args.samples_per_label}",
            "sampling.sample_batch_size=5",
        ],
        cwd=repo,
        env=env,
    )

    reports = project_root / "reports"
    reports.mkdir(exist_ok=True)
    output_files = sorted((repo / "data" / "outputs").glob("*.txt"))
    report = {
        "status": "pass",
        "repo": str(repo),
        "checkpoint": str(latest),
        "outputs": [str(x) for x in output_files],
        "epochs": args.epochs,
        "samples_per_label": args.samples_per_label,
    }
    (reports / "t4_smoke_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
