$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
New-Item -ItemType Directory -Force -Path bin | Out-Null
# Default to --release 21 (minimum supported LTS runtime, also matching Kaggle E2E & release JAR).
# Set $env:JAVA_RELEASE = '27' to target a newer bytecode release explicitly.
$release = if ($env:JAVA_RELEASE) { $env:JAVA_RELEASE } else { '21' }
$jdk = if (Test-Path (Join-Path $env:USERPROFILE 'tools\jdk-27\bin\javac.exe')) {
    Join-Path $env:USERPROFILE 'tools\jdk-27'
} elseif ($env:JAVA_HOME) {
    $env:JAVA_HOME
} else {
    $null
}
$javac = if ($jdk) { Join-Path $jdk 'bin\javac.exe' } else { 'javac' }
& $javac --release $release -Xlint:all -d bin src/hypergrass/*.java
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Built $((Get-ChildItem bin/hypergrass/*.class | Measure-Object).Count) classes (release $release) into bin/"
