$repo = 'D:\CodeProjects\cfd'
Set-Location $repo

# 1) Gamma-scan: map f(Gamma) = dCp(TE) with phi converged at each Gamma.
$env:WING3D_IMPRINT = 'll'
$env:WING3D_KUTTA = '3'
"=== START $(Get-Date) gscan mode3 M0.3 ===" | Out-File -Encoding utf8 tmp_a9_gscan.log
python tests3d/tsd_kutta_scan.py 0.3 3 -1.0 -0.5 0.0 0.5 1.0 2.0 *>> tmp_a9_gscan.log
"=== EXITCODE=$LASTEXITCODE $(Get-Date) ===" | Out-File -Encoding utf8 -Append tmp_a9_gscan.log

# 2) B ladder for the formulation as specified: mode 1 + entropy, ep8.
$env:WING3D_KUTTA = '1'
$env:WING3D_ENTROPY = '1'
$env:WING3D_BL = '0'
$env:A9_GEO = 'dfvlr_r4'
$env:A9_EP = '8'
"=== START $(Get-Date) B mode1 entropy ep8 ===" | Out-File -Encoding utf8 tmp_a9_B1.log
python tests3d/tsd_a9_ladder.py B *>> tmp_a9_B1.log
"=== EXITCODE=$LASTEXITCODE $(Get-Date) ===" | Out-File -Encoding utf8 -Append tmp_a9_B1.log

# 3) C: BL on top at M0.85 (mode 1), ep8.
$env:WING3D_BL = '1'
"=== START $(Get-Date) C mode1 entropy+BL ===" | Out-File -Encoding utf8 tmp_a9_C1.log
python tests3d/tsd_a9_ladder.py C *>> tmp_a9_C1.log
"=== EXITCODE=$LASTEXITCODE $(Get-Date) ===" | Out-File -Encoding utf8 -Append tmp_a9_C1.log

# 4) D: legacy free_k at M0.85/M1.0, ep4 (mode 0).
$env:WING3D_BL = '0'
$env:WING3D_KUTTA = '0'
$env:A9_EP = '4'
"=== START $(Get-Date) D free_k ===" | Out-File -Encoding utf8 tmp_a9_D.log
python tests3d/tsd_a9_ladder.py D *>> tmp_a9_D.log
"=== EXITCODE=$LASTEXITCODE $(Get-Date) ===" | Out-File -Encoding utf8 -Append tmp_a9_D.log
