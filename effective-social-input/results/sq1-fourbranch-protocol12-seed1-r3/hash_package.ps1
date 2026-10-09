param([switch]$Verify)
$packageRoot = $PSScriptRoot
$manifestPath = Join-Path $packageRoot 'SHA256SUMS.csv'
if ($Verify) {
    $hashRows = Import-Csv -LiteralPath $manifestPath
    $failures = @()
    foreach ($row in $hashRows) {
        $target = Join-Path $packageRoot $row.relative_path
        if (-not (Test-Path -LiteralPath $target -PathType Leaf)) {
            $failures += $row.relative_path
            continue
        }
        $info = Get-Item -LiteralPath $target
        $actual = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($actual -ne $row.sha256 -or $info.Length -ne [long]$row.bytes) { $failures += $row.relative_path }
    }
    if ($failures.Count -gt 0) { throw ('Hash validation failed: ' + ($failures -join ', ')) }
    Write-Output ('Verified ' + $hashRows.Count + ' files, all hashes and sizes match.')
} else {
    $hashRows = Get-ChildItem -LiteralPath $packageRoot -Recurse -File |
        Where-Object { $_.FullName -ne $manifestPath } |
        Sort-Object FullName |
        ForEach-Object {
            [pscustomobject]@{
                relative_path = $_.FullName.Substring($packageRoot.Length + 1).Replace('\','/')
                bytes = $_.Length
                sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
            }
        }
    $hashRows | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding utf8
    Write-Output ('Wrote ' + $hashRows.Count + ' file hashes. The manifest excludes itself.')
}
