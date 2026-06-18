"""
HSC Promoter Training Dataset Pipeline
Pinello Lab - Harvard Medical School

Generates a training dataset for DNA-Diffusion fine-tuning
using cell-type-specific regulatory sequences from ENCODE.

Cell types:  K562, HepG2, GM12878
Output:      200bp sequences linked to differentially expressed genes
Format:      DNA-Diffusion compatible (chr, sequence, TAG)

Parameters:
    MAX_PEAKS_PER_CELL  = 50000  (top peaks by DNase score)
    DIFF_EXPR_THRESHOLD = 2.0    (TPM ratio vs other cell types)
    MIN_TPM             = 1.0
    PEAK_GENE_WINDOW    = 10000  (bp from TSS)
    STRIDE              = 50     (bp, for peaks > 200bp)
    SPLIT               = 80/10/10 random stratified
"""

import subprocess, sys, os, time, tempfile, pickle
import gzip, shutil
from pathlib import Path

try:
    from google.colab import drive
    IN_COLAB = True
except ImportError:
    drive = None
    IN_COLAB = False

if IN_COLAB:
    drive.mount('/content/drive', force_remount=False)
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q',
        'pyfaidx', 'tqdm', 'pandas', 'numpy', 'scikit-learn', 'requests'])

try:
    import pandas as pd
    import numpy as np
    import requests
    from tqdm.auto import tqdm
    from sklearn.model_selection import train_test_split
except ImportError as exc:
    raise SystemExit(
        "Missing Python dependencies. Install with: "
        "python -m pip install pyfaidx tqdm pandas numpy scikit-learn requests"
    ) from exc

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

MAX_PEAKS_PER_CELL   = 50000
DIFF_EXPR_THRESHOLD  = 2.0
MIN_TPM              = 1.0
PEAK_GENE_WINDOW     = 10000
STRIDE               = 50

DEFAULT_CACHE = (
    Path("/content/drive/MyDrive/dna_diffusion_cache")
    if IN_COLAB else Path.home() / "dna_diffusion_cache"
)
CACHE    = Path(os.environ.get("DNA_DIFFUSION_CACHE", DEFAULT_CACHE)).expanduser()
DATA_DIR = CACHE / "data"
OUTPUT   = CACHE / "training_final"
TMP_DIR  = Path(os.environ.get("DNA_DIFFUSION_TMP", tempfile.gettempdir())).expanduser()
DNA_REPO = Path(os.environ.get(
    "DNA_DIFFUSION_REPO",
    "/content/DNA-Diffusion" if IN_COLAB else str(Path.home() / "DNA-Diffusion"),
)).expanduser()
for d in [CACHE, DATA_DIR, OUTPUT, TMP_DIR]:
    d.mkdir(exist_ok=True, parents=True)

EXPERIMENTS = {
    "K562":    {"dnase": "ENCSR000EOT", "rnaseq": ["ENCSR000AEQ", "ENCSR000COK"]},
    "HepG2":   {"dnase": "ENCSR000ENQ", "rnaseq": ["ENCSR000EYR", "ENCSR931WGT"]},
    "GM12878": {"dnase": "ENCSR000EMT", "rnaseq": ["ENCSR000AEF", "ENCSR000AEG"]},
}

# -----------------------------------------------------------------------
# Download utilities
# -----------------------------------------------------------------------

def download(url, cache_path, label=""):
    """Download to local /content, verify, copy to Drive cache."""
    if cache_path.exists():
        try:
            if str(cache_path).endswith('.gz'):
                with gzip.open(cache_path, 'rt') as f:
                    f.read(1000)
            return
        except:
            cache_path.unlink()

    local = TMP_DIR / cache_path.name
    for attempt in range(3):
        try:
            r = requests.get(url, stream=True, timeout=120)
            r.raise_for_status()
            total = int(r.headers.get('content-length', 0))
            with open(local, 'wb') as f:
                with tqdm(total=total, unit='B', unit_scale=True,
                          desc=label or cache_path.name) as pbar:
                    for chunk in r.iter_content(8192):
                        f.write(chunk)
                        pbar.update(len(chunk))
            if str(local).endswith('.gz'):
                with gzip.open(local, 'rt') as f:
                    f.read(1000)
            shutil.copy(local, cache_path)
            local.unlink()
            return
        except Exception as e:
            if local.exists():
                local.unlink()
            if attempt == 2:
                raise RuntimeError(f"Download failed after 3 attempts: {url}") from e
            time.sleep(5)

def encode_tsv_files(exp_id):
    r = requests.get(
        f"https://www.encodeproject.org/experiments/{exp_id}/?format=json",
        headers={'accept': 'application/json'}, timeout=30)
    return [f['accession'] for f in r.json().get('files', [])
            if 'tsv' in f.get('file_format', '')
            and 'gene quantifications' in f.get('output_type', '')]

def encode_bed_files(exp_id):
    r = requests.get(
        f"https://www.encodeproject.org/experiments/{exp_id}/?format=json",
        headers={'accept': 'application/json'}, timeout=30)
    return [f['accession'] for f in r.json().get('files', [])
            if 'bed' in f.get('file_format', '')]

# -----------------------------------------------------------------------
# Step 1: Gene annotations (GENCODE v38, protein-coding only)
# -----------------------------------------------------------------------

print("Step 1: Gene annotations")

genes_pkl = DATA_DIR / "genes.pkl"
if genes_pkl.exists():
    genes = pd.read_pickle(genes_pkl)
else:
    gtf = DATA_DIR / "gencode.gtf.gz"
    download(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_38/gencode.v38.annotation.gtf.gz",
        gtf, "GENCODE v38")

    rows = []
    with gzip.open(gtf, 'rt') as f:
        for line in tqdm(f, desc="  parsing GTF"):
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            p = line.split("\t")
            gid = gnm = gtp = ""
            for item in p[8].split(";"):
                if "gene_id"   in item: gid = item.split('"')[1].split(".")[0]
                elif "gene_name" in item: gnm = item.split('"')[1]
                elif "gene_type" in item: gtp = item.split('"')[1]
            if gtp != "protein_coding":
                continue
            ch  = p[0] if p[0].startswith("chr") else f"chr{p[0]}"
            tss = int(p[3]) if p[6] == "+" else int(p[4])
            rows.append({"chr": ch, "start": int(p[3]), "end": int(p[4]),
                         "gene_id": gid, "gene_name": gnm, "TSS": tss})
    genes = pd.DataFrame(rows)
    genes.to_pickle(genes_pkl)

print(f"  {len(genes):,} protein-coding genes\n")

# -----------------------------------------------------------------------
# Step 2: RNA-seq expression (2 datasets per cell type, median TPM)
# -----------------------------------------------------------------------

print("Step 2: RNA-seq expression")

expr_pkl = DATA_DIR / "expression_merged_median.pkl"
if expr_pkl.exists():
    expr = pd.read_pickle(expr_pkl)
else:
    cell_expr = {}
    for ct, cfg in EXPERIMENTS.items():
        dfs = []
        for exp_id in cfg["rnaseq"]:
            fids = encode_tsv_files(exp_id)
            if not fids:
                continue
            fid = fids[0]
            cp  = DATA_DIR / "rnaseq" / f"{ct}_{fid}.tsv"
            cp.parent.mkdir(exist_ok=True, parents=True)
            if not cp.exists():
                r = requests.get(
                    f"https://www.encodeproject.org/files/{fid}/@@download/{fid}.tsv",
                    timeout=120)
                cp.write_bytes(r.content)
            df = pd.read_csv(cp, sep="\t", comment="#")
            if "TPM" in df.columns and "gene_id" in df.columns:
                df["gene_id"] = df["gene_id"].astype(str).str.split(".").str[0]
                dfs.append(df[["gene_id", "TPM"]])

        if not dfs:
            continue
        merged = dfs[0].rename(columns={"TPM": "TPM_0"})
        for i, d in enumerate(dfs[1:], 1):
            merged = merged.merge(d.rename(columns={"TPM": f"TPM_{i}"}),
                                  on="gene_id", how="outer")
        tpm_cols = [c for c in merged.columns if c.startswith("TPM")]
        merged[f"{ct}_TPM"] = merged[tpm_cols].median(axis=1)
        cell_expr[ct] = merged[["gene_id", f"{ct}_TPM"]].fillna(0)
        print(f"  {ct}: {len(cell_expr[ct]):,} genes, "
              f"{(cell_expr[ct][f'{ct}_TPM'] > MIN_TPM).sum():,} expressed")

    expr = cell_expr["K562"]
    for ct in ["HepG2", "GM12878"]:
        expr = expr.merge(cell_expr[ct], on="gene_id", how="outer")
    expr = expr.fillna(0)
    expr.to_pickle(expr_pkl)

print(f"  merged matrix: {len(expr):,} genes\n")

# -----------------------------------------------------------------------
# Step 3: Differential expression
#   CELLTYPE: ratio = TPM(target) / max(TPM(others)) > 2.0 AND TPM > 1.0
# -----------------------------------------------------------------------

print("Step 3: Differential expression")

diff_genes = {}
for ct in ["K562", "HepG2", "GM12878"]:
    others    = [c for c in ["K562", "HepG2", "GM12878"] if c != ct]
    df        = expr.copy()
    df["ratio"] = df[f"{ct}_TPM"] / (df[[f"{c}_TPM" for c in others]].max(axis=1) + 0.01)
    diff      = df[(df["ratio"] > DIFF_EXPR_THRESHOLD) & (df[f"{ct}_TPM"] > MIN_TPM)]
    diff_genes[ct] = diff
    print(f"  {ct}: {len(diff):,} genes  "
          f"(median ratio {diff['ratio'].median():.1f}x, "
          f"mean TPM {diff[f'{ct}_TPM'].mean():.1f})")

print()

# -----------------------------------------------------------------------
# Step 4: DNase-seq peaks (top 50k by score)
# -----------------------------------------------------------------------

print("Step 4: DNase-seq peaks")

peaks       = {}
valid_chrs  = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

for ct, cfg in EXPERIMENTS.items():
    pk = DATA_DIR / f"{ct}_peaks_top50k.pkl"
    if pk.exists():
        peaks[ct] = pd.read_pickle(pk)
        print(f"  {ct}: {len(peaks[ct]):,} peaks (cached)")
        continue

    exp_id    = cfg["dnase"]
    fids      = encode_bed_files(exp_id)
    fid       = fids[0]
    local_gz  = TMP_DIR / f"{ct}.bed.gz"
    local_bed = TMP_DIR / f"{ct}.bed"

    r = requests.get(
        f"https://www.encodeproject.org/files/{fid}/@@download/{fid}.bed.gz",
        stream=True)
    with open(local_gz, 'wb') as f:
        for chunk in tqdm(r.iter_content(8192), desc=f"  {ct} DNase"):
            f.write(chunk)
    with gzip.open(local_gz, 'rb') as fi, open(local_bed, 'wb') as fo:
        shutil.copyfileobj(fi, fo)
    local_gz.unlink()

    df = pd.read_csv(local_bed, sep="\t", header=None,
                     names=["chr","start","end","name","score",
                            "strand","sv","pv","qv","peak"])
    local_bed.unlink()

    df["chr"] = df["chr"].astype(str).apply(
        lambda c: c if c.startswith("chr") else f"chr{c}")
    df = df[df["chr"].isin(valid_chrs)]

    if len(df) > MAX_PEAKS_PER_CELL:
        df = df.nlargest(MAX_PEAKS_PER_CELL, "score")

    df.to_pickle(pk)
    peaks[ct] = df
    print(f"  {ct}: {len(df):,} peaks")

print()

# -----------------------------------------------------------------------
# Step 5: Link peaks to differentially expressed genes (10kb window)
# -----------------------------------------------------------------------

print("Step 5: Peak-gene linking")

linked_all = []
for ct in ["K562", "HepG2", "GM12878"]:
    pk      = peaks[ct]
    tgt_ids = set(diff_genes[ct]["gene_id"].values)

    genes_by_chr = {
        ch: genes[(genes["chr"] == ch) & (genes["gene_id"].isin(tgt_ids))]
        for ch in pk["chr"].unique()
    }

    linked = []
    for _, peak in tqdm(pk.iterrows(), total=len(pk), desc=f"  {ct}"):
        cg = genes_by_chr.get(peak["chr"], pd.DataFrame())
        if cg.empty:
            continue
        center = (peak["start"] + peak["end"]) / 2
        nearby = cg[abs(cg["TSS"] - center) < PEAK_GENE_WINDOW]
        if not nearby.empty:
            gene = nearby.iloc[(nearby["TSS"] - center).abs().argmin()]
            linked.append({
                "chr":       peak["chr"],
                "start":     peak["start"],
                "end":       peak["end"],
                "label":     ct,
                "gene_name": gene["gene_name"],
                "gene_id":   gene["gene_id"],
            })

    df = pd.DataFrame(linked)
    linked_all.append(df)
    print(f"  {ct}: {len(df):,} peaks linked")

all_peaks = pd.concat(linked_all)
print(f"  total: {len(all_peaks):,} peaks\n")

# -----------------------------------------------------------------------
# Step 6: hg38 reference genome
# -----------------------------------------------------------------------

print("Step 6: hg38 genome")

fa = CACHE / "hg38.fa"
if not fa.exists():
    local_gz = TMP_DIR / "hg38.fa.gz"
    local_fa = TMP_DIR / "hg38.fa"
    download(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        local_gz, "hg38")
    with gzip.open(local_gz, 'rb') as fi, open(local_fa, 'wb') as fo:
        shutil.copyfileobj(fi, fo)
    local_gz.unlink()
    shutil.copy(local_fa, fa)
    local_fa.unlink()

from pyfaidx import Fasta
genome = Fasta(str(fa))
print(f"  loaded\n")

# -----------------------------------------------------------------------
# Step 7: Extract 200bp sequences (stride = 50bp for long peaks)
# -----------------------------------------------------------------------

print("Step 7: Sequence extraction")

seqs = []
for _, peak in tqdm(all_peaks.iterrows(), total=len(all_peaks)):
    plen = peak["end"] - peak["start"]
    if plen < 200:
        center = (peak["start"] + peak["end"]) // 2
        s, e   = center - 100, center + 100
        try:
            seq = str(genome[peak["chr"]][s:e]).upper()
            if len(seq) == 200 and "N" not in seq:
                seqs.append({"chr": peak["chr"], "sequence": seq,
                             "TAG": peak["label"], "gene": peak["gene_name"]})
        except:
            pass
    else:
        for i in range((plen - 200) // STRIDE + 1):
            s, e = peak["start"] + i * STRIDE, peak["start"] + i * STRIDE + 200
            if e > peak["end"]:
                break
            try:
                seq = str(genome[peak["chr"]][s:e]).upper()
                if len(seq) == 200 and "N" not in seq:
                    seqs.append({"chr": peak["chr"], "sequence": seq,
                                 "TAG": peak["label"], "gene": peak["gene_name"]})
            except:
                pass

sequences = pd.DataFrame(seqs)
print(f"  {len(sequences):,} sequences ({len(sequences)/len(all_peaks):.1f} per peak)")
print(sequences["TAG"].value_counts().to_string())
print()

# -----------------------------------------------------------------------
# Step 8: Train/val/test split (random stratified 80/10/10)
# -----------------------------------------------------------------------

print("Step 8: Train/val/test split")

train_val, test = train_test_split(
    sequences, test_size=0.1, random_state=42, stratify=sequences["TAG"])
train, val = train_test_split(
    train_val, test_size=0.111, random_state=42, stratify=train_val["TAG"])

print(f"  train: {len(train):,}  val: {len(val):,}  test: {len(test):,}")
print(train["TAG"].value_counts().to_string())

for name, df in [("train", train), ("val", val), ("test", test)]:
    df[["chr", "sequence", "TAG"]].to_csv(
        OUTPUT / f"{name}.txt", sep="\t", index=False, header=False)
train.to_csv(OUTPUT / "train_with_genes.csv", index=False)
print(f"  saved to {OUTPUT}\n")

# -----------------------------------------------------------------------
# DNA-Diffusion setup
# -----------------------------------------------------------------------

print("Setting up DNA-Diffusion")

DNA_REPO.parent.mkdir(exist_ok=True, parents=True)
os.chdir(DNA_REPO.parent)
if not DNA_REPO.exists():
    subprocess.run(["git", "clone",
        "https://github.com/pinellolab/DNA-Diffusion.git", str(DNA_REPO)], check=True)

os.chdir(DNA_REPO)

uv_bin = shutil.which("uv")
if uv_bin is None:
    raise RuntimeError("uv not found on PATH. Install uv, then rerun this script.")

subprocess.run([uv_bin, "sync"], check=True)
print("  installed\n")

# -----------------------------------------------------------------------
# Format data for DNA-Diffusion (requires header row)
# -----------------------------------------------------------------------

print("Formatting data")

os.makedirs("data", exist_ok=True)
for split in ["train", "val", "test"]:
    df = pd.read_csv(OUTPUT / f"{split}.txt",
                     sep="\t", header=None, names=["chr", "sequence", "TAG"])
    df.to_csv(f"data/{split}.txt", sep="\t", index=False)
    print(f"  data/{split}.txt  ({len(df):,} sequences)")

with open("data/hsc_encode_data.pkl", "wb") as f:
    pickle.dump({
        "train_df": train[["chr", "sequence", "TAG"]].reset_index(drop=True),
        "validation_df": val[["chr", "sequence", "TAG"]].reset_index(drop=True),
        "test_df": test[["chr", "sequence", "TAG"]].reset_index(drop=True),
    }, f)
print("  data/hsc_encode_data.pkl  (explicit random stratified split)")

print()
print("Pipeline complete.")
print(f"  sequences: {len(sequences):,}")
print(f"  labels:    {sorted(sequences['TAG'].unique())}")
print()
print("To train:")
print("  !uv run train.py \\")
print("      data.data_path=data/train.txt \\")
print("      data.load_saved_data=True \\")
print("      data.saved_data_path=data/hsc_encode_data.pkl")

# -----------------------------------------------------------------------
# Launch training
# -----------------------------------------------------------------------

import os
env = os.environ.copy()
env["WANDB_DISABLED"] = "true"
env["PYTHONUNBUFFERED"] = "1"

print()
print("Launching training (this runs for many hours)...")
print()

subprocess.run([
    uv_bin, "run", "python", "-u", "train.py",
    "data.data_path=data/train.txt",
    "data.load_saved_data=True",
    "data.saved_data_path=data/hsc_encode_data.pkl",
    "training.use_wandb=False",
] + (
    [f"training.num_epochs={os.environ['DNA_DIFFUSION_EPOCHS']}"]
    if os.environ.get("DNA_DIFFUSION_EPOCHS") else []
) + (
    [f"training.min_epochs={os.environ['DNA_DIFFUSION_MIN_EPOCHS']}"]
    if os.environ.get("DNA_DIFFUSION_MIN_EPOCHS") else []
), check=True, env=env)

print()
print("Training finished. Checkpoints are in checkpoints/")
