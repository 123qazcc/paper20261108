param(
    [string]$SourceRoot = 'private\source',
    [string]$Rscript = (Get-Command Rscript -ErrorAction Stop).Source,
    [string]$StageRoot = 'private\stage'
)

$ErrorActionPreference = 'Stop'

$projectRoot = Split-Path -Parent $PSScriptRoot
$rawRoot = Join-Path $SourceRoot 'data_raw'
$teachingRoot = Join-Path $SourceRoot 'reference_scripts'
$elsaRawRoot = Join-Path $rawRoot '2_ELSA_UK'
$elsaStage = Join-Path $StageRoot 'elsa'
$elsaEolStage = Join-Path $StageRoot 'elsa_eol'
$elsaHarmonizedStage = Join-Path $StageRoot 'elsa_harmonized'
$stageOutput = Join-Path $StageRoot 'output'

foreach ($path in @($SourceRoot, $rawRoot, $teachingRoot, $elsaRawRoot, $Rscript)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required path not found: $path"
    }
}

New-Item -ItemType Directory -Force -Path $elsaStage, $elsaEolStage, $elsaHarmonizedStage, $stageOutput | Out-Null

$elsaFiles = @(
    @{ Wave = 4; Name = 'wave_4_elsa_data_v3.dta' },
    @{ Wave = 5; Name = 'wave_5_elsa_data_v4.dta' },
    @{ Wave = 6; Name = 'wave_6_elsa_data_v2.dta' },
    @{ Wave = 7; Name = 'wave_7_elsa_data.dta' },
    @{ Wave = 8; Name = 'wave_8_elsa_data_eul_v2.dta' },
    @{ Wave = 9; Name = 'wave_9_elsa_data_eul_v1.dta' },
    @{ Wave = 10; Name = 'wave_10_elsa_data_eul_v1.dta' }
)

foreach ($item in $elsaFiles) {
    $source = Join-Path (Join-Path $elsaRawRoot ("wave" + $item.Wave)) $item.Name
    $destination = Join-Path $elsaStage ("w" + $item.Wave + '.dta')
    if (-not (Test-Path -LiteralPath $source)) {
        throw "ELSA source file not found: $source"
    }
    $copyNeeded = -not (Test-Path -LiteralPath $destination)
    if (-not $copyNeeded) {
        $copyNeeded = (Get-Item -LiteralPath $source).Length -ne (Get-Item -LiteralPath $destination).Length
    }
    if ($copyNeeded) {
        Copy-Item -LiteralPath $source -Destination $destination -Force
    }
}

$elsaEolSource = Join-Path $SourceRoot 'harmonized\ELSA\Harmonized_ELSA_EOL_A2\h_elsa_eol_a2.dta'
$elsaEolDestination = Join-Path $elsaEolStage 'h_elsa_eol_a2.dta'
if (Test-Path -LiteralPath $elsaEolSource) {
    Copy-Item -LiteralPath $elsaEolSource -Destination $elsaEolDestination -Force
}

# GlobalAgingR::elsa is the longitudinal spine, but its packaged Wave-5
# `angine` column is empty.  Stage the exact Harmonized ELSA source used by the
# package so the cleaning script can restore that component and independently
# verify the broad heart-problem phenotype against all component variables.
$elsaHarmonizedSource = Join-Path $SourceRoot 'harmonized\ELSA\Harmonized_ELSA_G3\h_elsa_g3.dta'
$elsaHarmonizedDestination = Join-Path $elsaHarmonizedStage 'h_elsa_g3.dta'
if (-not (Test-Path -LiteralPath $elsaHarmonizedSource)) {
    throw "Harmonized ELSA source file not found: $elsaHarmonizedSource"
}
$harmonizedCopyNeeded = -not (Test-Path -LiteralPath $elsaHarmonizedDestination)
if (-not $harmonizedCopyNeeded) {
    $harmonizedCopyNeeded = (Get-Item -LiteralPath $elsaHarmonizedSource).Length -ne (Get-Item -LiteralPath $elsaHarmonizedDestination).Length
}
if ($harmonizedCopyNeeded) {
    Copy-Item -LiteralPath $elsaHarmonizedSource -Destination $elsaHarmonizedDestination -Force
}

$env:GA_STAGE_ROOT = $StageRoot.Replace('\', '/')
$env:GA_OUTPUT_ROOT = $stageOutput.Replace('\', '/')
$env:LC_ALL = 'C.UTF-8'
$env:LANG = 'C.UTF-8'
$env:LC_CTYPE = 'C.UTF-8'

$rScriptSource = Join-Path $projectRoot 'R\01_clean_three_cohorts.R'
$rScriptStage = Join-Path $StageRoot '01_clean_three_cohorts.R'
Copy-Item -LiteralPath $rScriptSource -Destination $rScriptStage -Force

& $Rscript $rScriptStage 2>&1 | Tee-Object -FilePath (Join-Path $StageRoot 'cleaning_console.log')
if ($LASTEXITCODE -ne 0) {
    throw "R cleaning script failed with exit code $LASTEXITCODE"
}

$verifyScriptSource = Join-Path $projectRoot 'R\02_verify_clean_outputs.R'
$verifyScriptStage = Join-Path $StageRoot '02_verify_clean_outputs.R'
Copy-Item -LiteralPath $verifyScriptSource -Destination $verifyScriptStage -Force
& $Rscript $verifyScriptStage 2>&1 | Tee-Object -FilePath (Join-Path $StageRoot 'verification_console.log')
if ($LASTEXITCODE -ne 0) {
    throw "R output verification failed with exit code $LASTEXITCODE"
}

foreach ($folder in @('data_clean', 'qc', 'logs')) {
    $sourceFolder = Join-Path $stageOutput $folder
    $destinationFolder = Join-Path $projectRoot $folder
    New-Item -ItemType Directory -Force -Path $destinationFolder | Out-Null
    Get-ChildItem -LiteralPath $sourceFolder -File | ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $destinationFolder -Force
    }
}
Copy-Item -LiteralPath (Join-Path $StageRoot 'cleaning_console.log') -Destination (Join-Path $projectRoot 'logs\cleaning_console.log') -Force
Copy-Item -LiteralPath (Join-Path $StageRoot 'verification_console.log') -Destination (Join-Path $projectRoot 'logs\verification_console.log') -Force

Write-Host "Cleaning completed. Outputs copied to: $projectRoot"
