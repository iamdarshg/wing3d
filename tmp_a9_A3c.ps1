$repo = 'D:\CodeProjects\cfd'
Set-Location $repo
$env:WING3D_IMPRINT = 'll'
$env:WING3D_ENTROPY = '0'
$env:WING3D_BL = '0'
$env:A9_GEO = 'dfvlr_r4'
$env:WING3D_KUTTA = '3'
"=== START $(Get-Date) KUTTA=$env:WING3D_KUTTA IMPRINT=$env:WING3D_IMPRINT ===" | Out-File -Encoding utf8 tmp_a9_A3c.log
python tests3d/tsd_a9_ladder.py A *>> tmp_a9_A3c.log
"=== EXITCODE=$LASTEXITCODE $(Get-Date) ===" | Out-File -Encoding utf8 -Append tmp_a9_A3c.log
