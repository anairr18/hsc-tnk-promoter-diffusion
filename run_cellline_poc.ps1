param(
  [string]$ProjectRoot = "$PSScriptRoot",
  [string]$DNADiffusionRepo = "$HOME\DNA-Diffusion",
  [string]$CacheDir = "$HOME\dna_diffusion_cache",
  [switch]$PatchPretrained,
  [switch]$RunTraining
)

$ErrorActionPreference = "Stop"
$env:DNA_DIFFUSION_REPO = $DNADiffusionRepo
$env:DNA_DIFFUSION_CACHE = $CacheDir

python "$ProjectRoot\scripts\download_project_assets.py" --output-dir "$ProjectRoot\downloads"
python "$ProjectRoot\scripts\check_pretrained_initialization.py" --dnadiffusion-repo "$DNADiffusionRepo" --output "$ProjectRoot\reports\pretrained_initialization_check.json"

if ($PatchPretrained) {
  python "$ProjectRoot\scripts\patch_dnadiffusion_pretrained_finetune.py" --dnadiffusion-repo "$DNADiffusionRepo"
}

if ($RunTraining) {
  python "$ProjectRoot\pipeline\cell1_dataset_and_training.py"
  python "$ProjectRoot\pipeline\cell2_generate_and_backup.py"
}

python "$ProjectRoot\scripts\validate_cellline_poc.py" --project-root "$ProjectRoot"

$GeneratedDir = Join-Path $CacheDir "training_final\generated"
$SplitDir = Join-Path $CacheDir "training_final"
if ((Test-Path "$SplitDir\train.txt") -and (Test-Path "$SplitDir\test.txt")) {
  $GeneratedFiles = @()
  if (Test-Path $GeneratedDir) {
    $GeneratedFiles = Get-ChildItem $GeneratedDir -File | ForEach-Object { $_.FullName }
  }
  python "$ProjectRoot\scripts\leakage_and_novelty_audit.py" `
    --train "$SplitDir\train.txt" `
    --val "$SplitDir\val.txt" `
    --test "$SplitDir\test.txt" `
    --generated $GeneratedFiles `
    --output "$ProjectRoot\reports\cellline_poc_leakage_audit.tsv"
}

Write-Host "Cell-line POC workflow complete. Read reports\cellline_poc_report.md"
