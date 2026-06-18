param(
  [string]$ProjectRoot = "$PSScriptRoot",
  [string]$DNADiffusionRepo = "$HOME\DNA-Diffusion",
  [string]$Hg38Fasta,
  [string]$GencodeGtf,
  [string]$FantomTssBed,
  [string]$ActivityLongTsv
)

$ErrorActionPreference = "Stop"
$DataDir = Join-Path $ProjectRoot "data\hsc_tnk"
New-Item -ItemType Directory -Force -Path $DataDir | Out-Null

python "$ProjectRoot\scripts\catalog_public_hematopoietic_data.py"

if ($Hg38Fasta -and ($GencodeGtf -or $FantomTssBed)) {
  $cmd = @("$ProjectRoot\scripts\build_promoter_windows.py", "--fasta", $Hg38Fasta, "--output", "$DataDir\promoter_windows.tsv")
  if ($GencodeGtf) { $cmd += @("--gencode-gtf", $GencodeGtf) }
  if ($FantomTssBed) { $cmd += @("--fantom-tss-bed", $FantomTssBed) }
  python @cmd
}

if ($ActivityLongTsv) {
  python "$ProjectRoot\scripts\merge_promoter_activity_matrix.py" `
    --promoters "$DataDir\promoter_windows.tsv" `
    --activity-long "$ActivityLongTsv" `
    --output "$DataDir\promoter_activity_matrix.tsv"

  python "$ProjectRoot\scripts\select_tnk_specific_training_set.py" `
    --activity-matrix "$DataDir\promoter_activity_matrix.tsv" `
    --output "$DataDir\tnk_specific_training_set.tsv"

  python "$ProjectRoot\scripts\train_multimodal_dnadiffusion.py" `
    --training-set "$DataDir\tnk_specific_training_set.tsv" `
    --dnadiffusion-repo "$DNADiffusionRepo"
}

Write-Host "HSC/T-NK backbone prepared under $DataDir"
