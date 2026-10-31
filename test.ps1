$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
& "$root/build.ps1"

Write-Host '== VerifyAllEdgeCases =='
& java -cp bin hypergrass.VerifyAllEdgeCases
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyDonorIndex =='
& java -cp bin hypergrass.VerifyDonorIndex
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyAnchorCodecs =='
& java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyArchiveAccounting =='
& java -Xmx512m -cp bin hypergrass.VerifyArchiveAccounting
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host '== VerifyDonorCosts =='
& java -Xmx512m -cp bin hypergrass.VerifyDonorCosts
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host 'All verification suites passed.'
