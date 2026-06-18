#!/usr/bin/env python3
"""Patch a DNA-Diffusion checkout for continuous promoter-profile conditioning.

The patch keeps backward compatibility with class-only datasets. When a split
pickle contains numeric columns whose names include "__" or start with
"profile_", the dataloader returns `(x, y, profile)`, training passes profiles
through the Diffusion wrapper, and UNet adds a learned profile embedding to its
time/class conditioning.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def backup(path: Path) -> None:
    backup_path = path.with_suffix(path.suffix + ".bak_hsc_tnk")
    if not backup_path.exists():
        shutil.copy(path, backup_path)


def replace_once(text: str, old: str, new: str, path: Path) -> str:
    if old not in text:
        raise RuntimeError(f"Expected text not found in {path}: {old[:120]!r}")
    return text.replace(old, new, 1)


def patch_unet(path: Path) -> None:
    backup(path)
    text = path.read_text()
    text = replace_once(
        text,
        "        num_classes: int = 10,\n        output_attention: bool = False,\n",
        "        num_classes: int = 10,\n        profile_dim: int = 0,\n        output_attention: bool = False,\n",
        path,
    )
    text = replace_once(
        text,
        "        if num_classes is not None:\n            self.label_emb = nn.Embedding(num_classes, time_dim)\n\n        # layers\n",
        "        if num_classes is not None:\n            self.label_emb = nn.Embedding(num_classes, time_dim)\n        self.profile_dim = profile_dim\n        self.profile_mlp = None\n        if profile_dim and profile_dim > 0:\n            self.profile_mlp = nn.Sequential(\n                nn.Linear(profile_dim, time_dim),\n                nn.GELU(),\n                nn.Linear(time_dim, time_dim),\n            )\n\n        # layers\n",
        path,
    )
    text = replace_once(
        text,
        "    def forward(self, x: torch.Tensor, time: torch.Tensor, classes: torch.Tensor):\n",
        "    def forward(self, x: torch.Tensor, time: torch.Tensor, classes: torch.Tensor, profiles: torch.Tensor | None = None):\n",
        path,
    )
    text = replace_once(
        text,
        "        if classes is not None:\n            t_start += self.label_emb(classes)\n            t_mid += self.label_emb(classes)\n            t_end += self.label_emb(classes)\n            t_cross += self.label_emb(classes)\n\n        h = []\n",
        "        if classes is not None:\n            t_start += self.label_emb(classes)\n            t_mid += self.label_emb(classes)\n            t_end += self.label_emb(classes)\n            t_cross += self.label_emb(classes)\n        if profiles is not None and self.profile_mlp is not None:\n            profiles = profiles.to(t_start.device, dtype=t_start.dtype)\n            p = self.profile_mlp(profiles)\n            t_start += p\n            t_mid += p\n            t_end += p\n            t_cross += p\n\n        h = []\n",
        path,
    )
    path.write_text(text)


def patch_diffusion(path: Path) -> None:
    backup(path)
    text = path.read_text()
    text = replace_once(
        text,
        "    def p_losses(self, x_start, t, classes, noise=None, loss_type=\"huber\", p_uncond=0.1):\n",
        "    def p_losses(self, x_start, t, classes, profiles=None, noise=None, loss_type=\"huber\", p_uncond=0.1):\n",
        path,
    )
    text = replace_once(
        text,
        "        predicted_noise = self.model(x_noisy, t, classes)\n",
        "        predicted_noise = self.model(x_noisy, t, classes, profiles=profiles)\n",
        path,
    )
    text = replace_once(
        text,
        "    def forward(self, x, classes):\n",
        "    def forward(self, x, classes, profiles=None):\n",
        path,
    )
    text = replace_once(
        text,
        "        return self.p_losses(x, t, classes)\n",
        "        return self.p_losses(x, t, classes, profiles=profiles)\n",
        path,
    )
    path.write_text(text)


def patch_train_util(path: Path) -> None:
    backup(path)
    text = path.read_text()
    text = replace_once(
        text,
        "    x: torch.Tensor,\n    y: torch.Tensor,\n    model: torch.nn.Module,\n",
        "    batch,\n    model: torch.nn.Module,\n",
        path,
    )
    text = replace_once(
        text,
        "    x = x.to(device, dtype=torch.float32)\n    y = y.to(device)\n    with torch.autocast(device_type=device, dtype=dtype):\n        loss = model(x, y)\n",
        "    if len(batch) == 3:\n        x, y, profiles = batch\n        profiles = profiles.to(device, dtype=torch.float32)\n    else:\n        x, y = batch\n        profiles = None\n    x = x.to(device, dtype=torch.float32)\n    y = y.to(device)\n    with torch.autocast(device_type=device, dtype=dtype):\n        loss = model(x, y, profiles=profiles)\n",
        path,
    )
    text = replace_once(
        text,
        "    x: torch.Tensor,\n    y: torch.Tensor,\n    model: torch.nn.Module,\n",
        "    batch,\n    model: torch.nn.Module,\n",
        path,
    )
    text = replace_once(
        text,
        "    x = x.to(device, dtype=dtype)\n    y = y.to(device)\n    loss = model(x, y)\n",
        "    if len(batch) == 3:\n        x, y, profiles = batch\n        profiles = profiles.to(device, dtype=torch.float32)\n    else:\n        x, y = batch\n        profiles = None\n    x = x.to(device, dtype=dtype)\n    y = y.to(device)\n    loss = model(x, y, profiles=profiles)\n",
        path,
    )
    path.write_text(text)


def patch_train(path: Path) -> None:
    backup(path)
    text = path.read_text()
    text = text.replace(
        "        for x, y in train_dl:\n            loss = train_step(x, y, model, optimizer, device, precision)\n",
        "        for batch in train_dl:\n            loss = train_step(batch, model, optimizer, device, precision)\n",
    )
    text = text.replace(
        "        for x, y in val_dl:\n            val_loss = val_step(x, y, model, device, precision)\n",
        "        for batch in val_dl:\n            val_loss = val_step(batch, model, device, precision)\n",
    )
    path.write_text(text)


def patch_dataloader(path: Path) -> None:
    backup(path)
    text = path.read_text()
    text = replace_once(
        text,
        "    train_data = SequenceDataset(x_data, y_data)\n    val_data = SequenceDataset(x_val_data, y_val_data)\n",
        "    train_profiles = encode_data.get(\"train_profiles\")\n    val_profiles = encode_data.get(\"val_profiles\")\n    train_data = SequenceDataset(x_data, y_data, profiles=train_profiles)\n    val_data = SequenceDataset(x_val_data, y_val_data, profiles=val_profiles)\n",
        path,
    )
    text = replace_once(
        text,
        "    # Collecting variables into a dict\n    encode_data_dict = {\n",
        "    profile_cols = [c for c in df.columns if \"__\" in c or c.startswith(\"profile_\")]\n    train_profiles = None\n    val_profiles = None\n    if profile_cols:\n        train_profiles = torch.tensor(df[profile_cols].astype(float).values, dtype=torch.float32)\n        val_profiles = torch.tensor(val_df[profile_cols].astype(float).values, dtype=torch.float32)\n\n    # Collecting variables into a dict\n    encode_data_dict = {\n",
        path,
    )
    text = replace_once(
        text,
        "        \"x_val_cell_type\": x_val_cell_type,\n",
        "        \"x_val_cell_type\": x_val_cell_type,\n        \"profile_cols\": profile_cols,\n        \"train_profiles\": train_profiles,\n        \"val_profiles\": val_profiles,\n",
        path,
    )
    text = replace_once(
        text,
        "        transform: T.Compose = T.Compose([T.ToTensor()]),\n    ):\n",
        "        transform: T.Compose = T.Compose([T.ToTensor()]),\n        profiles: torch.Tensor | None = None,\n    ):\n",
        path,
    )
    text = replace_once(
        text,
        "        self.transform = transform\n",
        "        self.transform = transform\n        self.profiles = profiles\n",
        path,
    )
    text = replace_once(
        text,
        "        return x, y\n",
        "        if self.profiles is not None:\n            return x, y, self.profiles[index]\n        return x, y\n",
        path,
    )
    path.write_text(text)


def patch_unet_config(path: Path, profile_dim: int) -> None:
    backup(path)
    text = path.read_text()
    if "profile_dim:" not in text:
        text += f"\nprofile_dim: {profile_dim}\n"
    else:
        lines = [f"profile_dim: {profile_dim}" if line.startswith("profile_dim:") else line for line in text.splitlines()]
        text = "\n".join(lines) + "\n"
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnadiffusion-repo", required=True)
    parser.add_argument("--profile-dim", type=int, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    repo = Path(args.dnadiffusion_repo)
    targets = [
        repo / "src" / "dnadiffusion" / "models" / "unet.py",
        repo / "src" / "dnadiffusion" / "models" / "diffusion.py",
        repo / "src" / "dnadiffusion" / "utils" / "train_util.py",
        repo / "src" / "dnadiffusion" / "data" / "dataloader.py",
        repo / "train.py",
        repo / "configs" / "model" / "unet.yaml",
    ]
    missing = [str(p) for p in targets if not p.exists()]
    if missing:
        raise SystemExit(f"Missing expected DNA-Diffusion files: {missing}")
    if args.dry_run:
        print("Patch dry-run OK. Files that would be patched:")
        for p in targets:
            print(f"  {p}")
        return

    patch_unet(targets[0])
    patch_diffusion(targets[1])
    patch_train_util(targets[2])
    patch_dataloader(targets[3])
    patch_train(targets[4])
    patch_unet_config(targets[5], args.profile_dim)
    print("Applied multimodal conditioning patch. Backups use suffix .bak_hsc_tnk")


if __name__ == "__main__":
    main()
