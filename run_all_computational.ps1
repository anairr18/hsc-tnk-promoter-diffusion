param(
  [string]$ProjectRoot = "$PSScriptRoot",
  [string]$BaseDir = "$PSScriptRoot\..",
  [string]$DNADiffusionRepo = "",
  [string]$CacheDir = "",
  [string]$Stage2Mode = "demo",
  [switch]$BuildRealStage2Inputs,
  [switch]$DownloadReferences,
  [string]$Hg38Fasta = "",
  [string]$GencodeGtf = "",
  [string]$FantomTssBed = "",
  [string]$SignalManifest = "",
  [string]$ExpressionLongTsv = "",
  [string]$TssActivityLongTsv = "",
  [int]$Stage1Epochs = 3000,
  [int]$Stage1MinEpochs = 1000,
  [int]$Stage2Epochs = 500,
  [int]$Stage2MinEpochs = 100,
  [int]$Stage2Samples = 1000
)

$ErrorActionPreference = "Stop"
if (-not $DNADiffusionRepo) { $DNADiffusionRepo = Join-Path $BaseDir "DNA-Diffusion" }
if (-not $CacheDir) { $CacheDir = Join-Path $BaseDir "dna_diffusion_cache" }

$env:PROJECT_ROOT = $ProjectRoot
$env:BASE_DIR = $BaseDir
$env:DNA_DIFFUSION_REPO = $DNADiffusionRepo
$env:DNA_DIFFUSION_CACHE = $CacheDir
$env:STAGE2_MODE = $Stage2Mode
$env:BUILD_REAL_STAGE2_INPUTS = if ($BuildRealStage2Inputs) { "1" } else { "0" }
$env:DOWNLOAD_REFERENCES = if ($DownloadReferences) { "1" } else { "0" }
if ($Hg38Fasta) { $env:HG38_FASTA = $Hg38Fasta }
if ($GencodeGtf) { $env:GENCODE_GTF = $GencodeGtf }
if ($FantomTssBed) { $env:FANTOM_TSS_BED = $FantomTssBed }
if ($SignalManifest) { $env:SIGNAL_MANIFEST = $SignalManifest }
if ($ExpressionLongTsv) { $env:EXPRESSION_LONG_TSV = $ExpressionLongTsv }
if ($TssActivityLongTsv) { $env:TSS_ACTIVITY_LONG_TSV = $TssActivityLongTsv }
$env:PATCH_PRETRAINED = "1"
$env:DNA_DIFFUSION_EPOCHS = "$Stage1Epochs"
$env:DNA_DIFFUSION_MIN_EPOCHS = "$Stage1MinEpochs"
$env:STAGE2_EPOCHS = "$Stage2Epochs"
$env:STAGE2_MIN_EPOCHS = "$Stage2MinEpochs"
$env:STAGE2_SAMPLES = "$Stage2Samples"

bash "$ProjectRoot/run_all_computational.sh"
