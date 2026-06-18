#!/usr/bin/env python3
"""Check whether the local DNA-Diffusion repo supports true pretrained training."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def read(path: Path) -> str:
    return path.read_text(errors="replace") if path.exists() else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnadiffusion-repo", default=Path.home() / "DNA-Diffusion")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    repo = Path(args.dnadiffusion_repo)
    train_py = read(repo / "train.py")
    sample_hf = read(repo / "sample_hf.py")
    unet_pretrained = read(repo / "configs" / "model" / "unet_pretrained.yaml")
    pretrained_model = read(repo / "src" / "dnadiffusion" / "models" / "pretrained_unet.py")

    evidence = {
        "repo": str(repo),
        "train_py_exists": bool(train_py),
        "sample_hf_exists": bool(sample_hf),
        "unet_pretrained_config_exists": bool(unet_pretrained),
        "pretrained_unet_wrapper_exists": bool(pretrained_model),
        "train_mentions_from_pretrained": "from_pretrained" in train_py,
        "train_mentions_load_state_dict": "load_state_dict" in train_py,
        "train_mentions_checkpoint_path": "checkpoint_path" in train_py,
        "sample_hf_mentions_from_pretrained": "from_pretrained" in sample_hf,
        "unet_pretrained_mentions_from_pretrained": "from_pretrained" in unet_pretrained,
        "unet_pretrained_target": re.findall(r"_target_:\s*(.+)", unet_pretrained),
    }
    supported = (
        evidence["train_mentions_from_pretrained"]
        or evidence["train_mentions_checkpoint_path"]
        or evidence["train_mentions_load_state_dict"]
    )
    result = {
        "true_pretrained_training_supported_by_current_train_py": bool(supported),
        "recommended_poc_baseline": "true_pretrained_finetune" if supported else "training_from_scratch_or_apply_model_patch",
        "evidence": evidence,
        "notes": [
            "sample_hf.py can instantiate PretrainedUNet.from_pretrained for sampling.",
            "train.py must load pretrained weights before optimizer construction to perform true fine-tuning.",
            "If unsupported, record the first POC as a from-scratch baseline and run a second patched fine-tune later.",
        ],
    }
    output = Path(args.output) if args.output else repo / "pretrained_initialization_check.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {output}")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
