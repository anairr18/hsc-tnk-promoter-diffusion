#!/usr/bin/env python3
"""Patch DNA-Diffusion train.py so model=unet_pretrained performs true fine-tuning.

Upstream sample_hf.py unwraps PretrainedUNet.from_pretrained(...).model before
building the Diffusion wrapper. The training entrypoint does not. This patch
adds the same unwrapping in train.py so this command becomes valid:

uv run python -u train.py model=unet_pretrained data.load_saved_data=True ...
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


OLD = """    model = hydra.utils.instantiate(cfg.model)
    data = hydra.utils.instantiate(cfg.data)
    optimizer = hydra.utils.instantiate(cfg.optimizer, model.parameters())
    diffusion = hydra.utils.instantiate(cfg.diffusion, model=model)
"""

NEW = """    model = hydra.utils.instantiate(cfg.model)
    if hasattr(model, "model") and model.__class__.__name__ == "PretrainedUNet":
        print("Using pretrained UNet weights for fine-tuning")
        model = model.model
    data = hydra.utils.instantiate(cfg.data)
    optimizer = hydra.utils.instantiate(cfg.optimizer, model.parameters())
    diffusion = hydra.utils.instantiate(cfg.diffusion, model=model)
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnadiffusion-repo", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    train_py = Path(args.dnadiffusion_repo) / "train.py"
    if not train_py.exists():
        raise SystemExit(f"Missing train.py: {train_py}")
    text = train_py.read_text()
    if NEW in text:
        print("Pretrained fine-tune patch already applied.")
        return
    if OLD not in text:
        raise SystemExit("Could not find expected train.py block; inspect manually before patching.")
    if args.dry_run:
        print(f"Dry-run OK: would patch {train_py}")
        return
    backup = train_py.with_suffix(".py.bak_pretrained_finetune")
    if not backup.exists():
        shutil.copy(train_py, backup)
    train_py.write_text(text.replace(OLD, NEW, 1))
    print(f"Patched {train_py}; backup at {backup}")


if __name__ == "__main__":
    main()
