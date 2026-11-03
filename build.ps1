$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
New-Item -ItemType Directory -Force -Path bin | Out-Null
# Build with JDK 27; --release 27 targets that feature release regardless of the
# host JDK, so a newer or older javac on PATH cannot silently change the output.
$jdk = if (Test-Path (Join-Path $env:USERPROFILE 'tools\jdk-27\bin\javac.exe')) {
    Join-Path $env:USERPROFILE 'tools\jdk-27'
} elseif ($env:JAVA_HOME) {
    $env:JAVA_HOME
} else {
    $null
}
$javac = if ($jdk) { Join-Path $jdk 'bin\javac.exe' } else { 'javac' }
& $javac --release 27 -Xlint:all -d bin src/hypergrass/*.java
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Built $((Get-ChildItem bin/hypergrass/*.class | Measure-Object).Count) classes into bin/"
