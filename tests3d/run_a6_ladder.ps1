param([string]$Only = '')
$ErrorActionPreference = 'Stop'
$root = 'D:\CodeProjects\cfd'
Set-Location $root
$env:WING3D_ABLATE_RAM_GB = '1.5'
$env:WING3D_TSD_SPAN = '1.6'
$env:WING3D_CTEST_EP = '8'

function Invoke-Probe([string]$Tag, [string]$Mach, [hashtable]$Knobs) {
    $out = "tests3d\validation\results\ablate_mp_${Tag}_M$Mach.json"
    if (Test-Path $out) {
        Write-Host "SKIP $Tag M$Mach (result exists: $out)"
        return
    }
    foreach ($key in @('WING3D_FREE_K', 'WING3D_BC_ROWS', 'WING3D_TE_GAP')) {
        Remove-Item "Env:$key" -ErrorAction SilentlyContinue
    }
    foreach ($key in $Knobs.Keys) { Set-Item "Env:$key" $Knobs[$key] }
    Write-Host "Starting $Tag M$Mach at $(Get-Date -Format o)"
    & python tests3d/tsd_ablate_mp.py $Tag $Mach --span 1.6
    if ($LASTEXITCODE -ne 0) { throw "$Tag M$Mach failed (exit $LASTEXITCODE); stopping ladder" }
}

$steps = @(
    @{Tag='u1_freek'; Knobs=@{WING3D_FREE_K='1'}},
    @{Tag='u2_bcrows'; Knobs=@{WING3D_BC_ROWS='1'}},
    @{Tag='u4_tegap0'; Knobs=@{WING3D_TE_GAP='0'}}
)
foreach ($step in $steps) {
    if ($Only -and $Only -ne $step.Tag) { continue }
    foreach ($mach in @('0.850', '0.925', '1.000')) {
        Invoke-Probe $step.Tag $mach $step.Knobs
    }
}
if ($Only -and $Only -ne 'u3') { exit 0 }
Remove-Item Env:WING3D_RAW_SECTION -ErrorAction SilentlyContinue
& python tests3d/tsd_ablate_u3.py --mach 0.85 --epochs 8
if ($LASTEXITCODE -ne 0) { throw "U3 failed (exit $LASTEXITCODE)" }
