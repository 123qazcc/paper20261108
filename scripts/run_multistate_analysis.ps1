[CmdletBinding()]
param(
    [ValidateSet('Core', 'Full')]
    [string]$Profile = 'Full',
    [string]$Rscript = (Get-Command Rscript -ErrorAction Stop).Source,
    [string]$StageRoot = 'private\stage'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$projectRoot = Split-Path -Parent $PSScriptRoot
$exposureInput = Join-Path $StageRoot 'analysis_output\data_derived'
$multistateOutput = Join-Path $StageRoot 'multistate_output'
$multistateData = Join-Path $multistateOutput 'data_multistate'
$coxMiceCache = Join-Path $multistateData 'cox_mice_cache'
$markovResultRoot = Join-Path $multistateOutput 'results\markov'
$markovModelRoot = Join-Path $multistateOutput 'models\markov'
$markovLogRoot = Join-Path $multistateOutput 'logs\markov'
$stagedScriptRoot = Join-Path $StageRoot 'multistate_scripts'
$profileKey = $Profile.ToLowerInvariant()
$runStamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$pipelineLogRoot = Join-Path $multistateOutput (Join-Path 'logs\pipeline' (Join-Path $profileKey $runStamp))

foreach ($path in @($Rscript, $exposureInput)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required path not found: $path"
    }
}

$steps = @(
    '05_prepare_multistate_data.R',
    '05b_verify_disease_states.R',
    '05c_audit_elsa_w10_heart_mapping.R',
    '06_fit_transition_cox_models.R',
    '07_fit_markov_models.R',
    '07b_protocol_sensitivities.R',
    '08_collate_multistate_outputs.R',
    '08_make_multistate_figures.R',
    '09_verify_multistate_outputs.R'
)

# Validate the complete pipeline before starting a long model run.
foreach ($step in $steps) {
    $source = Join-Path $projectRoot (Join-Path 'R' $step)
    if (-not (Test-Path -LiteralPath $source)) {
        throw "Analysis step not found: $source"
    }
}

foreach ($path in @($multistateOutput, $multistateData, $stagedScriptRoot, $pipelineLogRoot)) {
    New-Item -ItemType Directory -Force -Path $path | Out-Null
}

function Convert-ToRPath {
    param([Parameter(Mandatory = $true)][string]$Path)
    return $Path.Replace('\', '/')
}

function Copy-OutputTree {
    param(
        [Parameter(Mandatory = $true)][string]$SourceRoot,
        [Parameter(Mandatory = $true)][string]$DestinationRoot
    )

    if (-not (Test-Path -LiteralPath $SourceRoot)) {
        return
    }

    New-Item -ItemType Directory -Force -Path $DestinationRoot | Out-Null
    Get-ChildItem -LiteralPath $SourceRoot -File -Recurse | ForEach-Object {
        $relative = $_.FullName.Substring($SourceRoot.Length).TrimStart([char[]]'\/')
        $destination = Join-Path $DestinationRoot $relative
        $destinationParent = Split-Path -Parent $destination
        if (-not (Test-Path -LiteralPath $destinationParent)) {
            New-Item -ItemType Directory -Force -Path $destinationParent | Out-Null
        }
        Copy-Item -LiteralPath $_.FullName -Destination $destination -Force
    }
}

function Copy-AllOutputs {
    foreach ($folder in @('data_multistate', 'results', 'figures', 'reports', 'qc', 'models', 'logs')) {
        Copy-OutputTree `
            -SourceRoot (Join-Path $multistateOutput $folder) `
            -DestinationRoot (Join-Path $projectRoot $folder)
    }
}

# All R scripts run from an ASCII-only staging path because this Windows R
# installation does not reliably translate the Chinese project/source paths.
$env:GA_STAGE_ROOT = Convert-ToRPath $StageRoot
$env:GA_EXPOSURE_INPUT_ROOT = Convert-ToRPath $exposureInput
$env:GA_MULTISTATE_OUTPUT_ROOT = Convert-ToRPath $multistateOutput
$env:GA_MULTISTATE_DATA_ROOT = Convert-ToRPath $multistateData
$env:GA_COX_OUTPUT_ROOT = Convert-ToRPath $multistateOutput
$env:GA_COX_MICE_CACHE_ROOT = Convert-ToRPath $coxMiceCache
$env:GA_MARKOV_RESULT_ROOT = Convert-ToRPath $markovResultRoot
$env:GA_MARKOV_MODEL_ROOT = Convert-ToRPath $markovModelRoot
$env:GA_MARKOV_LOG_ROOT = Convert-ToRPath $markovLogRoot
$env:GA_PROJECT_ROOT = ''
$env:GA_COX_COHORTS = 'HRS,ELSA'
$env:GA_COX_TRANSITIONS = 'ALL'
$env:GA_COX_WEIGHTED_ONLY = 'false'
$env:GA_RUN_ZPH = 'true'
$env:GA_RUN_GLOBALAGINGR_BASELINE = 'true'
$env:GA_MARKOV_ANALYSIS_REGEX = '.'
$env:GA_MARKOV_FORCE_REFIT = 'false'
$env:GA_MARKOV_MAX_PANEL_ROWS = '0'
$env:GA_MARKOV_MAX_MODELS = '0'
$env:GA_MARKOV_MAX_SECONDS_PER_MODEL = '0'
$env:GA_MARKOV_MAX_ITER = '1000'
$env:GA_MARKOV_NUISANCE_CONSTRAINT = 'common'
$env:GA_COX_SEED = '20260908'
$env:GA_MARKOV_MI_SEED = '20260908'
$env:GA_ELSA_DEATH_CUTOFF = '2012.999'
$env:GA_PROTOCOL_SENS_COHORTS = 'HRS,ELSA'
$env:GA_PROTOCOL_SENS_SCENARIOS = 'rcs,alternate_death,attained_age_period,ipcw,proxy_exclusion,expanded_hypertension,semi_markov'
$env:GA_PROTOCOL_SENS_MAX_MODELS = '0'
$env:GA_PROTOCOL_SENS_SAVE_MODELS = 'true'
$env:GA_FIGURE_RUN_PDF_QA = 'true'
$env:GA_MULTISTATE_PROFILE = $profileKey
$env:NATURE_FIGURE_SKILL_ROOT = 'private/render_tools/nature-figure'
$env:GA_FIGURE_AUDIT_PYTHON = 'python'
$env:LC_ALL = 'C.UTF-8'
$env:LANG = 'C.UTF-8'
$env:LC_CTYPE = 'C.UTF-8'

if ($Profile -eq 'Full') {
    # Complete requested analysis: every Cox sensitivity module, unweighted and
    # survey-weighted fits, MICE m=20, and the complete Markov/MI model registry.
    $env:GA_COX_TIER = 'full'
    $env:GA_COX_MODULES = 'all'
    $env:GA_COX_RUN_LABEL = 'full_all_weighted_mi20'
    $env:GA_RUN_MICE = 'true'
    $env:GA_RUN_WEIGHTED = 'true'
    $env:GA_MICE_M = '20'
    $env:GA_MICE_MAXIT = '10'
    $env:GA_SAVE_COX_OBJECTS = 'true'
    $env:GA_SAVE_IMPUTATION_ESTIMATES = 'true'
    $env:GA_MARKOV_RUN_MODE = 'full'
    $env:GA_MARKOV_MODEL_REGEX = '.'
    $env:GA_MARKOV_MI_MODE = 'run'
    $env:GA_MARKOV_MI_M = '20'
    $env:GA_MARKOV_MI_MAXIT = '10'
    $env:GA_PROTOCOL_SENS_RUN_MODE = 'full'
    $env:GA_FIGURE_STRICT_COMPLETENESS = 'true'
} else {
    # Development profile: primary M1-M4 complete-case models only. It omits
    # MI, weighting, secondary exposures, dynamic patterns, and interactions.
    $env:GA_COX_TIER = 'core'
    $env:GA_COX_MODULES = 'main'
    $env:GA_COX_RUN_LABEL = 'core_main_unweighted_complete_case'
    $env:GA_RUN_MICE = 'false'
    $env:GA_RUN_WEIGHTED = 'false'
    $env:GA_MICE_M = '20'
    $env:GA_MICE_MAXIT = '10'
    $env:GA_SAVE_COX_OBJECTS = 'false'
    $env:GA_SAVE_IMPUTATION_ESTIMATES = 'false'
    $env:GA_MARKOV_RUN_MODE = 'full'
    $env:GA_MARKOV_MODEL_REGEX = '^(GA_M0|M1|M2|M3|M4)$'
    $env:GA_MARKOV_MI_MODE = 'off'
    $env:GA_MARKOV_MI_M = '20'
    $env:GA_MARKOV_MI_MAXIT = '10'
    $env:GA_PROTOCOL_SENS_RUN_MODE = 'audit'
    $env:GA_FIGURE_STRICT_COMPLETENESS = 'false'
}

$pipelineFailed = $false
try {
    foreach ($step in $steps) {
        $source = Join-Path $projectRoot (Join-Path 'R' $step)
        $staged = Join-Path $stagedScriptRoot $step
        Copy-Item -LiteralPath $source -Destination $staged -Force

        $logName = [System.IO.Path]::GetFileNameWithoutExtension($step) + '.log'
        $logPath = Join-Path $pipelineLogRoot $logName
        Write-Host "[$Profile] Running $step"
        & $Rscript $staged 2>&1 | Tee-Object -FilePath $logPath
        $stepExitCode = $LASTEXITCODE
        if ($stepExitCode -ne 0) {
            throw "Analysis step failed: $step (exit $stepExitCode). Log: $logPath"
        }
    }
} catch {
    $pipelineFailed = $true
    throw
} finally {
    # Preserve successful and partial artefacts plus diagnostic logs even when a
    # later step fails; existing project files are updated but never deleted.
    Copy-AllOutputs
    if ($pipelineFailed) {
        Write-Warning "Pipeline stopped before completion. Partial outputs and logs were copied to: $projectRoot"
    }
}

Write-Host "Multistate $Profile profile completed. Outputs copied to: $projectRoot"
