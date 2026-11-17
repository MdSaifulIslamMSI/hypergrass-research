$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
New-Item -ItemType Directory -Force -Path bin | Out-Null
# Default to --release 25 (latest Java 25 LTS baseline, matching GraalVM 25 & Kaggle E2E).
# Set $env:JAVA_RELEASE = '27' to target Java 27 bytecode explicitly.
$release = if ($env:JAVA_RELEASE) { $env:JAVA_RELEASE } else { '25' }
$jdk = @(
    (Join-Path $env:USERPROFILE 'tools\jdk-27'),
    (Join-Path $env:USERPROFILE 'tools\jdk-25'),
    $env:JAVA_HOME
) | Where-Object { $_ -and (Test-Path (Join-Path $_ 'bin\javac.exe')) } | Select-Object -First 1
$javac = if ($jdk) { Join-Path $jdk 'bin\javac.exe' } else { 'javac' }
& $javac --release $release -Xlint:all -d bin src/hypergrass/*.java
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Built $((Get-ChildItem bin/hypergrass/*.class | Measure-Object).Count) classes (release $release) into bin/"
