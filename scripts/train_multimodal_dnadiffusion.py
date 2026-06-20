#!/usr/bin/env python3
"""Prepare and optionally launch multimodal DNA-Diffusion training.

This script exports continuous promoter activity profiles for each sequence and
creates a compatible proxy TAG for the current DNA-Diffusion code. The proxy
mode lets the project progress immediately while preserving continuous profiles
for the follow-up model patch that projects profile vectors into the UNet time
embedding.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-set", required=True)
    parser.add_argument("--dnadiffusion-repo", default=Path.home() / "DNA-Diffusion")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--profile-prefixes", default="accessibility__,initiation__,expression__")
    parser.add_argument("--num-profile-bins", type=int, default=4)
    parser.add_argument("--run", action="store_true", help="Actually launch uv run train.py after exporting data.")
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--min-epochs", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    args = parser.parse_args()

    repo = Path(args.dnadiffusion_repo)
    out_dir = Path(args.output_dir) if args.output_dir else repo / "data" / "hsc_tnk_multimodal"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.training_set, sep="\t")
    prefixes = tuple(x.strip() for x in args.profile_prefixes.split(",") if x.strip())
    profile_cols = [c for c in df.columns if c.startswith(prefixes)]
    if not profile_cols:
        raise SystemExit(f"No profile columns found with prefixes: {prefixes}")
    required = {"chr", "sequence", "TAG"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"training set missing required columns: {sorted(missing)}")

    profile = df[profile_cols].astype(float).to_numpy(dtype=np.float32)
    np.save(out_dir / "condition_profiles.npy", profile)
    (out_dir / "condition_columns.json").write_text(json.dumps(profile_cols, indent=2) + "\n")

    # Compatibility proxy for current DNA-Diffusion, which only supports class labels.
    score = df.get("tnk_specificity_score", pd.Series(np.linspace(0, 1, len(df))))
    bins = pd.qcut(score.rank(method="first"), q=args.num_profile_bins, labels=False, duplicates="drop")
    proxy = df[["chr", "sequence", "TAG", *profile_cols]].copy()
    proxy["TAG"] = proxy["TAG"].astype(str) + "_PROFILE_Q" + bins.astype(int).astype(str)
    train_path = out_dir / "train_multimodal_proxy.txt"
    proxy[["chr", "sequence", "TAG"]].to_csv(train_path, sep="\t", index=False)

    # Preserve explicit split for DNA-Diffusion dataloader.
    train_df = proxy.sample(frac=0.8, random_state=42)
    rest = proxy.drop(train_df.index)
    val_df = rest.sample(frac=0.5, random_state=42) if len(rest) else rest
    test_df = rest.drop(val_df.index)
    split_pkl = out_dir / "hsc_tnk_multimodal_proxy.pkl"
    with split_pkl.open("wb") as f:
        pickle.dump(
            {
                "train_df": train_df.reset_index(drop=True),
                "validation_df": val_df.reset_index(drop=True),
                "test_df": test_df.reset_index(drop=True),
            },
            f,
        )

    manifest = {
        "mode": "profile-binned proxy for current DNA-Diffusion",
        "continuous_profile_files": ["condition_profiles.npy", "condition_columns.json"],
        "proxy_train_path": str(train_path),
        "proxy_split_pickle": str(split_pkl),
        "next_model_patch": "replace class-only label embedding with label embedding + profile MLP in UNet time embedding",
    }
    (out_dir / "training_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    if args.run:
        uv = shutil.which("uv")
        if not uv:
            raise SystemExit("uv not found on PATH")
        cmd = [
            uv,
            "run",
            "python",
            "-u",
            "train.py",
            f"data.data_path={train_path.as_posix()}",
            "data.load_saved_data=True",
            f"data.saved_data_path={split_pkl.as_posix()}",
            "training.use_wandb=False",
        ]
        if args.epochs is not None:
            cmd.append(f"training.num_epochs={args.epochs}")
        if args.min_epochs is not None:
            cmd.append(f"training.min_epochs={args.min_epochs}")
        if args.batch_size is not None:
            cmd.append(f"training.batch_size={args.batch_size}")
        cmd.append("training.sample_epoch=999999")
        env = os.environ.copy()
        env["WANDB_DISABLED"] = "true"
        env["PYTHONUNBUFFERED"] = "1"
        subprocess.run(cmd, cwd=repo, env=env, check=True)

    print(f"Wrote multimodal training export to {out_dir}")


if __name__ == "__main__":
    main()
