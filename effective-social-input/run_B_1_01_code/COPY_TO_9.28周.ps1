$src = 'C:\Users\WHQ\Documents\Codex\2026-09-16\new-chat\run_B_1'
$dst = 'C:\Users\WHQ\Studying\大三上学期\无人机群\9.28周\run_B_1'
New-Item -ItemType Directory -Force -Path $dst | Out-Null
Copy-Item -LiteralPath (Join-Path $src '*') -Destination $dst -Recurse -Force
Write-Host "已复制到 $dst"
