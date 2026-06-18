#!/usr/bin/env python3
"""Generate T/NK promoter candidates from a trained DNA-Diffusion checkpoint."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnadiffusion-repo", default=Path.home() / "DNA-Diffusion")
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--saved-data-path", required=True)
    parser.add_argument("--cell-types", default="TNK_HIGH,TNK_HIGH_PROFILE_Q3")
    parser.add_argument("--number-of-samples", type=int, default=10000)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    desired_profile = {
        "target": "high T/NK promoter activity",
        "offtarget": "low HSC/HSPC, B, myeloid, erythroid activity",
        "cell_types": args.cell_types.split(","),
        "number_of_samples": args.number_of_samples,
    }
    (out / "desired_generation_profile.json").write_text(json.dumps(desired_profile, indent=2) + "\n")

    cmd = [
        "uv",
        "run",
        "python",
        "-u",
        "sample.py",
        f"sampling.checkpoint_path={Path(args.checkpoint).as_posix()}",
        "data.load_saved_data=True",
        f"data.saved_data_path={Path(args.saved_data_path).as_posix()}",
        f"data.cell_types={args.cell_types}",
        f"sampling.number_of_samples={args.number_of_samples}",
    ]
    (out / "sample_command.txt").write_text(" ".join(cmd) + "\n")
    if args.run:
        uv = shutil.which("uv")
        if not uv:
            raise SystemExit("uv not found on PATH")
        cmd[0] = uv
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        env["WANDB_DISABLED"] = "true"
        subprocess.run(cmd, cwd=args.dnadiffusion_repo, env=env, check=True)
        for dname in ["samples", "generated_sequences"]:
            d = Path(args.dnadiffusion_repo) / dname
            if d.exists():
                for f in d.glob("*"):
                    if f.is_file():
                        shutil.copy(f, out / f.name)
    print(f"Wrote generation profile and command to {out}")


if __name__ == "__main__":
    main()
