@echo off
rem Instrumented conv2 run wrapper: %1=geo %2=M [%3 %4 %5=grid %6=nz? %7=epochs]
rem Usage: tmp_run_conv2.cmd <geo> <M> [nx ny nz [ep]]
cd /d D:\CodeProjects\cfd
set SUF=
if not "%3"=="" set SUF=_g%3x%4x%5
set LOG=tmp_conv2_%1_M%2%SUF%.log
set ERR=tmp_conv2_%1_M%2%SUF%.err
echo === START %DATE% %TIME% === > "%LOG%"
python tests3d\tsd_conv2.py %1 %2 %3 %4 %5 %6 %7 >> "%LOG%" 2> "%ERR%"
echo === EXITCODE=%ERRORLEVEL% %DATE% %TIME% === >> "%LOG%"
