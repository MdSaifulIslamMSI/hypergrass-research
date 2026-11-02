# Fetch and build the third-party codecs used by the specialized harness.
# The checkouts and binaries are gitignored; this script recreates them.
$ErrorActionPreference = "Stop"

$vendor = $PSScriptRoot
$project = (Resolve-Path (Join-Path $vendor "..\..\..")).Path
$grassCommit = "6d93700f5a2de8841f7cf2c20064413c80df2c4a"
$libbscCommit = "1dc676540cb40a41562af6f2db04bd1adce97f90"

function Sync-Repo([string]$name, [string]$url, [string]$commit) {
    $path = Join-Path $vendor $name
    if (-not (Test-Path $path)) {
        git clone $url $path
    }
    git -C $path fetch --depth 1 origin $commit
    git -C $path checkout --detach $commit
}

Sync-Repo "GraSS" "https://github.com/AnshuOishik/GraSS" $grassCommit
Sync-Repo "libbsc" "https://github.com/IlyaGrebnov/libbsc" $libbscCommit

$image = (Get-Content (Join-Path $project "benchmarks\corpus\image.lock.json") |
    ConvertFrom-Json).image_id
docker run --rm --network=none --cpus 2 --memory 2g `
    --mount "type=bind,source=$project,target=/work" `
    --workdir /work --entrypoint sh $image `
    /work/benchmarks/specialized/vendor/build_vendor.sh
