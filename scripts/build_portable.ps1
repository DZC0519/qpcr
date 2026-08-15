param([switch]$SkipTests)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$projectRoot = Split-Path -Parent $PSScriptRoot
$workspaceRoot = Split-Path -Parent (Split-Path -Parent $projectRoot)
$outputs = Join-Path $workspaceRoot 'outputs'
$buildRoot = Join-Path $projectRoot 'build'
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$product = 'qPCR' + (-join ([char[]](0x5206, 0x6790, 0x52A9, 0x624B)))
$portableFolder = Join-Path $outputs ($product + '-portable-win64')
$portableZip = Join-Path $outputs ($product + '-portable-win64.zip')
$sourceZip = Join-Path $outputs ($product + '-source.zip')
$rootManifest = Join-Path $outputs 'SHA256SUMS.txt'

function Remove-SafeBuildDirectory([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return }
    $resolvedBuild = [IO.Path]::GetFullPath($buildRoot).TrimEnd('\') + '\'
    $resolvedTarget = [IO.Path]::GetFullPath($Path).TrimEnd('\') + '\'
    if (-not $resolvedTarget.StartsWith($resolvedBuild, [StringComparison]::OrdinalIgnoreCase)) {
        throw ('Refusing to remove a directory outside build: ' + $resolvedTarget)
    }
    Remove-Item -LiteralPath $Path -Recurse -Force
}

function Remove-SafeOutputDirectory([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return }
    $resolvedOutputs = [IO.Path]::GetFullPath($outputs).TrimEnd('\') + '\'
    $resolvedTarget = [IO.Path]::GetFullPath($Path).TrimEnd('\') + '\'
    if (-not $resolvedTarget.StartsWith($resolvedOutputs, [StringComparison]::OrdinalIgnoreCase)) {
        throw ('Refusing to remove a directory outside outputs: ' + $resolvedTarget)
    }
    Remove-Item -LiteralPath $Path -Recurse -Force
}

function Write-DirectoryManifest([string]$Root, [string]$ManifestPath) {
    $resolvedRoot = [IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    $lines = Get-ChildItem -LiteralPath $Root -File -Recurse |
        Where-Object { $_.FullName -ne $ManifestPath } |
        Sort-Object FullName |
        ForEach-Object {
            $relative = $_.FullName.Substring($resolvedRoot.Length).Replace('\', '/')
            $hash = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
            $hash + '  ' + $relative
        }
    Set-Content -LiteralPath $ManifestPath -Value $lines -Encoding utf8
}

if (-not (Test-Path -LiteralPath $python)) {
    throw 'Virtual environment not found. Create .venv and install requirements first.'
}

New-Item -ItemType Directory -Force -Path $outputs, $buildRoot | Out-Null
Push-Location $projectRoot
try {
    & $python 'scripts\generate_docs.py'
    if ($LASTEXITCODE -ne 0) { throw 'Documentation generation failed.' }
    $env:PYTHONPATH = Join-Path $projectRoot 'src'
    & $python -m qpcr_analyzer.demo
    if ($LASTEXITCODE -ne 0) { throw 'Synthetic example generation failed.' }
    if (-not $SkipTests) {
        & $python -m pytest -q
        if ($LASTEXITCODE -ne 0) { throw 'Tests failed.' }
    }

    $pyiDist = Join-Path $buildRoot 'pyinstaller-dist'
    $pyiWork = Join-Path $buildRoot 'pyinstaller-work'
    Remove-SafeBuildDirectory $pyiDist
    Remove-SafeBuildDirectory $pyiWork
    & $python -m PyInstaller `
        --noconfirm --clean --windowed --onedir `
        --name $product `
        --paths 'src' `
        --distpath $pyiDist `
        --workpath $pyiWork `
        --specpath $buildRoot `
        --hidden-import matplotlib.backends.backend_agg `
        --hidden-import matplotlib.backends.backend_svg `
        --hidden-import matplotlib.backends.backend_pdf `
        'src\qpcr_analyzer\__main__.py'
    if ($LASTEXITCODE -ne 0) { throw 'PyInstaller build failed.' }

    $packageStage = Join-Path $buildRoot 'portable-staging'
    Remove-SafeBuildDirectory $packageStage
    New-Item -ItemType Directory -Force -Path $packageStage | Out-Null
    $appSource = Join-Path $pyiDist $product
    $appTarget = Join-Path $packageStage $product
    Copy-Item -LiteralPath $appSource -Destination $appTarget -Recurse
    Copy-Item -LiteralPath 'assets' -Destination (Join-Path $appTarget 'assets') -Recurse
    Copy-Item -LiteralPath 'docs\QUICKSTART_zh-CN.md' -Destination $appTarget
    Copy-Item -LiteralPath 'LICENSE' -Destination $appTarget
    Copy-Item -LiteralPath 'THIRD_PARTY_LICENSES.md' -Destination $appTarget
    Write-DirectoryManifest $packageStage (Join-Path $packageStage 'SHA256SUMS.txt')
    Remove-SafeOutputDirectory $portableFolder
    Copy-Item -LiteralPath $packageStage -Destination $portableFolder -Recurse
    if (Test-Path -LiteralPath $portableZip) { Remove-Item -LiteralPath $portableZip -Force }
    Compress-Archive -Path (Join-Path $packageStage '*') -DestinationPath $portableZip -CompressionLevel Optimal

    $sourceStage = Join-Path $buildRoot 'source-staging'
    Remove-SafeBuildDirectory $sourceStage
    $sourceTarget = Join-Path $sourceStage 'qpcr_analyzer'
    New-Item -ItemType Directory -Force -Path $sourceTarget | Out-Null
    foreach ($directory in @('src', 'tests', 'scripts', 'assets', 'docs')) {
        Copy-Item -LiteralPath $directory -Destination (Join-Path $sourceTarget $directory) -Recurse
    }
    foreach ($file in @('.gitignore', 'LICENSE', 'README.md', 'THIRD_PARTY_LICENSES.md',
            'pyproject.toml', 'requirements.txt', 'requirements-dev.txt',
            'requirements-lock.txt', 'run_qpcr_analyzer.bat')) {
        Copy-Item -LiteralPath $file -Destination $sourceTarget
    }
    $cacheDirectories = @(Get-ChildItem -LiteralPath $sourceTarget -Directory -Recurse -Force |
        Where-Object { $_.Name -in @('__pycache__', '.pytest_cache') })
    foreach ($cacheDirectory in $cacheDirectories) {
        Remove-SafeBuildDirectory $cacheDirectory.FullName
    }
    if (Test-Path -LiteralPath $sourceZip) { Remove-Item -LiteralPath $sourceZip -Force }
    Compress-Archive -Path (Join-Path $sourceStage '*') -DestinationPath $sourceZip -CompressionLevel Optimal

    Copy-Item -LiteralPath 'docs\QUICKSTART_zh-CN.md' -Destination $outputs -Force
    Copy-Item -LiteralPath 'assets\synthetic_qpcr_example.xlsx' -Destination $outputs -Force
    Copy-Item -LiteralPath 'assets\synthetic_qpcr_calculation.xlsx' -Destination $outputs -Force
    Copy-Item -LiteralPath 'assets\synthetic_plate_template.json' -Destination $outputs -Force
    Copy-Item -LiteralPath 'THIRD_PARTY_LICENSES.md' -Destination $outputs -Force
    $rootHashFiles = @(
        $portableZip,
        $sourceZip,
        (Join-Path $outputs 'QUICKSTART_zh-CN.md'),
        (Join-Path $outputs 'synthetic_qpcr_example.xlsx'),
        (Join-Path $outputs 'synthetic_qpcr_calculation.xlsx'),
        (Join-Path $outputs 'synthetic_plate_template.json'),
        (Join-Path $outputs 'THIRD_PARTY_LICENSES.md')
    )
    $rootHashLines = $rootHashFiles | ForEach-Object {
        $hash = (Get-FileHash -LiteralPath $_ -Algorithm SHA256).Hash.ToLowerInvariant()
        $hash + '  ' + (Split-Path -Leaf $_)
    }
    Set-Content -LiteralPath $rootManifest -Value $rootHashLines -Encoding utf8
    Write-Host ('Portable folder:  ' + $portableFolder)
    Write-Host ('Portable package: ' + $portableZip)
    Write-Host ('Source package:   ' + $sourceZip)
}
finally {
    Pop-Location
}
