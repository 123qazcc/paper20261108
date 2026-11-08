param(
    [string]$Rscript = (Get-Command Rscript -ErrorAction Stop).Source,
    [string]$StageRoot = 'private\stage'
)

$ErrorActionPreference = 'Stop'

$projectRoot = Split-Path -Parent $PSScriptRoot
$cleanInput = Join-Path $StageRoot 'output\data_clean'
$analysisOutput = Join-Path $StageRoot 'analysis_output'

foreach ($path in @($Rscript, $cleanInput)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required path not found: $path"
    }
}

New-Item -ItemType Directory -Force -Path $analysisOutput | Out-Null

$env:GA_STAGE_ROOT = $StageRoot.Replace('\', '/')
$env:GA_CLEAN_INPUT_ROOT = $cleanInput.Replace('\', '/')
$env:GA_ANALYSIS_OUTPUT_ROOT = $analysisOutput.Replace('\', '/')
$env:LC_ALL = 'C.UTF-8'
$env:LANG = 'C.UTF-8'
$env:LC_CTYPE = 'C.UTF-8'

$sourceScript = Join-Path $projectRoot 'R\03_build_harmonized_exposures_baseline.R'
$stageScript = Join-Path $StageRoot '03_build_harmonized_exposures_baseline.R'
Copy-Item -LiteralPath $sourceScript -Destination $stageScript -Force

& $Rscript $stageScript 2>&1 | Tee-Object -FilePath (Join-Path $StageRoot 'exposure_baseline_console.log')
if ($LASTEXITCODE -ne 0) {
    throw "Exposure/baseline script failed with exit code $LASTEXITCODE"
}

$verifySource = Join-Path $projectRoot 'R\04_verify_exposures_baseline.R'
$verifyStage = Join-Path $StageRoot '04_verify_exposures_baseline.R'
Copy-Item -LiteralPath $verifySource -Destination $verifyStage -Force
& $Rscript $verifyStage 2>&1 | Tee-Object -FilePath (Join-Path $StageRoot 'exposure_baseline_verification_console.log')
if ($LASTEXITCODE -ne 0) {
    throw "Exposure/baseline verification failed with exit code $LASTEXITCODE"
}

foreach ($folder in @('data_derived', 'qc')) {
    $sourceFolder = Join-Path $analysisOutput $folder
    $destinationFolder = Join-Path $projectRoot $folder
    New-Item -ItemType Directory -Force -Path $destinationFolder | Out-Null
    Get-ChildItem -LiteralPath $sourceFolder -File | ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $destinationFolder -Force
    }
}

$logFolder = Join-Path $projectRoot 'logs'
New-Item -ItemType Directory -Force -Path $logFolder | Out-Null
Copy-Item -LiteralPath (Join-Path $StageRoot 'exposure_baseline_console.log') -Destination (Join-Path $logFolder 'exposure_baseline_console.log') -Force
Copy-Item -LiteralPath (Join-Path $StageRoot 'exposure_baseline_verification_console.log') -Destination (Join-Path $logFolder 'exposure_baseline_verification_console.log') -Force

Write-Host "Exposure construction and baseline outputs copied to: $projectRoot"
