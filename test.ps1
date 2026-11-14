$ErrorActionPreference = 'Stop'
if (-not ($env:CI -eq 'true' -and $env:GITHUB_ACTIONS -eq 'true')) {
    Write-Host 'MANDATORY POLICY: Never run tests locally.' -ForegroundColor Red
    Write-Host 'Run all verification suites and world-compressor tests end-to-end on Kaggle:' -ForegroundColor Yellow
    Write-Host '  python benchmarks/kaggle/run_kaggle_e2e.py' -ForegroundColor Cyan
    exit 1
}
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
& "$root/build.ps1"

$jdk = if (Test-Path (Join-Path $env:USERPROFILE 'tools\jdk-27\bin\java.exe')) {
    Join-Path $env:USERPROFILE 'tools\jdk-27'
} elseif ($env:JAVA_HOME) {
    $env:JAVA_HOME
} else {
    $null
}
$java = if ($jdk) { Join-Path $jdk 'bin\java.exe' } else { 'java' }

Write-Host '== VerifyAllEdgeCases =='
& $java -cp bin hypergrass.VerifyAllEdgeCases
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyDonorIndex =='
& $java -cp bin hypergrass.VerifyDonorIndex
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyAnchorCodecs =='
& $java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyArchiveAccounting =='
& $java -Xmx512m -cp bin hypergrass.VerifyArchiveAccounting
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyDonorCosts =='
& $java -Xmx512m -cp bin hypergrass.VerifyDonorCosts
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyCompressionLaws =='
& $java -Xmx512m -cp bin hypergrass.VerifyCompressionLaws
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyFusionEngine =='
& $java -Xmx512m -cp bin hypergrass.VerifyFusionEngine
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host 'All verification suites passed.'
