param([string]$OutputName='sq1-fourbranch-protocol12-seed1-r3-results-v4')
$packageRoot=$PSScriptRoot
$resultRoot=Join-Path (Join-Path $packageRoot 'results') $OutputName
$started=Get-Date
$modes=@('fresh','drop','drop-normalized','hold')
$completed=Import-Csv -LiteralPath (Join-Path $resultRoot 'T0_TEST_MANIFEST.csv')
$pending=foreach($seed in 2..20) {
    foreach($mode in $modes) {
        [pscustomobject]@{
            run_id=$mode+'_seed'+$seed+'_T0'; seed=$seed; mode=$mode
            mask_rule='random'; duration_steps=0; reference_id='fresh_seed'+$seed
            checkpoint_id='seed'+$seed+'_step6000'; file_path=''
            status='NOT_RUN'; error_reason='Formal seed batch postponed until C independent acceptance'
            is_shared_reference=($mode -eq 'fresh')
        }
    }
}
@($completed)+@($pending) | Export-Csv -LiteralPath (Join-Path $resultRoot 'T0_TEST_MANIFEST.csv') -NoTypeInformation -Encoding utf8
$environment=[ordered]@{
    peak_memory_measured=$false
    peak_memory_bytes=$null
    reason='Peak process memory was not instrumented during the MATLAB run'
    proposed_worker_count=0
    scheduling='serial; formal-batch performance not validated'
    t0_completed=4
    t0_not_run=76
}
$environment | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $resultRoot 'resource_and_batch_status.json') -Encoding utf8
foreach($name in @('R3_STATUS.md','README.md','DATA_DICTIONARY.md','FIXES.md','load_pilot_example.m')) {
    Copy-Item -LiteralPath (Join-Path $packageRoot $name) -Destination $resultRoot -Force
}
$sourceHash=(Get-FileHash -LiteralPath (Join-Path $resultRoot 'source_manifest.json') -Algorithm SHA256).Hash.ToLowerInvariant()
$audit=[ordered]@{
    command='./complete_delivery.ps1 -OutputName '+$OutputName
    script_sha256=(Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash.ToLowerInvariant()
    source_manifest_sha256=$sourceHash
    implementation_commit=$null
    implementation_identity='Frozen source file SHA256 manifest; local package not committed or pushed'
    start_time=$started.ToString('o')
    end_time=(Get-Date).ToString('o')
    scope='Add pending T0 statuses, resource limitations, document copies and delivery audit; MAT and test evidence unchanged'
}
$audit | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $resultRoot 'delivery_audit.json') -Encoding utf8
Write-Output 'Delivery metadata completed. Regenerate and verify SHA256SUMS.csv after this step.'
