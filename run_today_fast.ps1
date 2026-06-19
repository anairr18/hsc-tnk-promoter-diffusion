param(
  [string]$ProjectRoot = "$PSScriptRoot",
  [string]$DNADiffusionRepo = "$HOME\DNA-Diffusion",
  [string]$CacheDir = "$HOME\dna_diffusion_cache",
  [switch]$AllowCpu,
  [switch]$PatchPretrained,
  [int]$Epochs = 5000,
  [int]$MinEpochs = 2000
)

$ErrorActionPreference = "Stop"
$env:DNA_DIFFUSION_REPO = $DNADiffusionRepo
$env:DNA_DIFFUSION_CACHE = $CacheDir
$env:DNA_DIFFUSION_EPOCHS = "$Epochs"
$env:DNA_DIFFUSION_MIN_EPOCHS = "$MinEpochs"

Write-Host "== Checking GPU =="
$gpuOk = $false
try {
  python -c "import torch; print('cuda', torch.cuda.is_available()); print('device_count', torch.cuda.device_count()); raise SystemExit(0 if torch.cuda.is_available() else 2)"
  $gpuOk = $true
} catch {
  $gpuOk = $false
}

if (-not $gpuOk -and -not $AllowCpu) {
  Write-Error "No CUDA GPU detected. Re-run on a GPU machine, or pass -AllowCpu for a slow CPU smoke run."
}

Write-Host "== Installing lightweight requirements =="
python -m pip install -r "$ProjectRoot\requirements.txt"

Write-Host "== Ensuring uv and DNA-Diffusion repo =="
$uv = Get-Command uv -ErrorAction SilentlyContinue
if (-not $uv) {
  python -m pip install -q uv
}
if (-not (Test-Path (Join-Path $DNADiffusionRepo ".git"))) {
  if (Test-Path $DNADiffusionRepo) {
    Remove-Item -Recurse -Force $DNADiffusionRepo
  }
  git clone https://github.com/pinellolab/DNA-Diffusion.git $DNADiffusionRepo
}

Write-Host "== Downloading metadata =="
python "$ProjectRoot\scripts\download_project_assets.py" --output-dir "$ProjectRoot\downloads"

Write-Host "== Checking pretrained fine-tune support =="
python "$ProjectRoot\scripts\check_pretrained_initialization.py" `
  --dnadiffusion-repo "$DNADiffusionRepo" `
  --output "$ProjectRoot\reports\pretrained_initialization_check.json"

if ($PatchPretrained) {
  Write-Host "== Patching DNA-Diffusion train.py for pretrained fine-tuning =="
  python "$ProjectRoot\scripts\patch_dnadiffusion_pretrained_finetune.py" `
    --dnadiffusion-repo "$DNADiffusionRepo"
}

Write-Host "== Running cell-line POC dataset + training =="
python "$ProjectRoot\pipeline\cell1_dataset_and_training.py"

Write-Host "== Generating sequences =="
python "$ProjectRoot\pipeline\cell2_generate_and_backup.py"

Write-Host "== Validating POC =="
python "$ProjectRoot\scripts\validate_cellline_poc.py" --project-root "$ProjectRoot"

Write-Host "DONE. Open reports\cellline_poc_report.md"
