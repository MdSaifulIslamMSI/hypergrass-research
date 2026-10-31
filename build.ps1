$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
New-Item -ItemType Directory -Force -Path bin | Out-Null
javac --release 21 -Xlint:all -d bin src/hypergrass/*.java
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Built $((Get-ChildItem bin/hypergrass/*.class | Measure-Object).Count) classes into bin/"
