"""
Cell 2: Generate sequences from fine-tuned checkpoint + backup to Drive

Run this AFTER cell1_dataset_and_training.py completes.
Re-establishes paths/env independently so it works even in a fresh runtime
(as long as Drive has the checkpoints from cell 1).
"""

import os, shutil, subprocess
from pathlib import Path

try:
    from google.colab import drive
    IN_COLAB = True
except ImportError:
    drive = None
    IN_COLAB = False

if IN_COLAB and os.environ.get("DNA_DIFFUSION_CACHE", "").startswith("/content/drive"):
    drive.mount('/content/drive', force_remount=False)

DEFAULT_CACHE = (
    Path("/content/drive/MyDrive/dna_diffusion_cache")
    if IN_COLAB else Path.home() / "dna_diffusion_cache"
)
CACHE = Path(os.environ.get("DNA_DIFFUSION_CACHE", DEFAULT_CACHE)).expanduser()
OUTPUT = CACHE / "training_final"
DNA_REPO = Path(os.environ.get(
    "DNA_DIFFUSION_REPO",
    "/content/DNA-Diffusion" if IN_COLAB else str(Path.home() / "DNA-Diffusion"),
)).expanduser()

os.chdir(DNA_REPO)
uv_bin = shutil.which("uv")
if uv_bin is None:
    raise RuntimeError("uv not found on PATH. Install uv, then rerun this script.")

env = os.environ.copy()
env["WANDB_DISABLED"] = "true"
env["PYTHONUNBUFFERED"] = "1"

# -----------------------------------------------------------------------
# Locate checkpoint
# -----------------------------------------------------------------------

ckpt_dir = Path("checkpoints")
ckpts = sorted(ckpt_dir.glob("*.pt"), key=lambda p: p.stat().st_mtime)
if not ckpts:
    raise FileNotFoundError("No checkpoint found in checkpoints/ - did training complete?")

best_ckpt = ckpts[-1]
print(f"Using checkpoint: {best_ckpt}")

# -----------------------------------------------------------------------
# Generate sequences
# -----------------------------------------------------------------------

gen_dir = OUTPUT / "generated"
gen_dir.mkdir(exist_ok=True, parents=True)

print("Generating 1000 sequences per cell type (K562, HepG2, GM12878)...")

subprocess.run([
    uv_bin, "run", "python", "-u", "sample.py",
    f"sampling.checkpoint_path={best_ckpt.as_posix()}",
    "data.load_saved_data=True",
    "data.saved_data_path=data/hsc_encode_data.pkl",
    "data.cell_types=K562,HepG2,GM12878",
    "sampling.number_of_samples=1000",
], check=True, env=env)

# Copy generated outputs to Drive
copied = []
for f in Path(".").glob("*generated*"):
    shutil.copy(f, gen_dir / f.name)
    copied.append(f.name)
if Path("samples").exists():
    for f in Path("samples").glob("*"):
        shutil.copy(f, gen_dir / f.name)
        copied.append(f.name)

print(f"  copied: {copied}")
print(f"  saved to {gen_dir}")

# -----------------------------------------------------------------------
# Backup checkpoints and logs
# -----------------------------------------------------------------------

print("\nBacking up checkpoints and logs to Drive...")

ckpt_backup = CACHE / "checkpoints"
ckpt_backup.mkdir(exist_ok=True, parents=True)
for f in ckpt_dir.glob("*.pt"):
    shutil.copy(f, ckpt_backup / f.name)
    print(f"  {f.name} -> {ckpt_backup}")

if Path("outputs").exists():
    logs_backup = CACHE / "training_logs"
    if logs_backup.exists():
        shutil.rmtree(logs_backup)
    shutil.copytree("outputs", logs_backup)
    print(f"  logs -> {logs_backup}")

# -----------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------

print()
print("="*60)
print("COMPLETE")
print("="*60)
print(f"  Generated sequences: {gen_dir}")
print(f"  Checkpoint backup:   {ckpt_backup}")
print(f"  Training logs:       {CACHE / 'training_logs'}")
print()
print("Next steps:")
print("  1. Run Enformer validation on generated sequences")
print("     (compare against endogenous + pretrained-model outputs)")
print("  2. Review training loss curve in training_logs/")
print("  3. Share results with Giacomo")
print("  4. If results look good, select top candidates for MPRA")
