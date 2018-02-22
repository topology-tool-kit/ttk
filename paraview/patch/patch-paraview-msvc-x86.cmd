@echo off
echo Patching ParaView sources for compiling on MSVC x86
set cdir=%cd%
set paraviewsrc=%1
:loop
if "%paraviewsrc%"=="" goto error
if not exist %paraviewsrc% goto error

set paraviewpath=%paraviewsrc%/ThirdParty/cgns/vtkcgns/src/adf/ADF_internals.c
if exist %paraviewpath% goto patch
echo File %paraviewpath% not exist
goto exit
:patch
%~d1
cd %paraviewsrc%
git apply --ignore-whitespace %cdir%/ADF_internals.c.patch
cd %cdir%
echo Finished patching.

goto exit
:error
echo Error: Missing ParaView source tree path.
echo Usage: patch-paraview-msvc-x86.cmd "<ParaView source tree path>"
set /P "query=Continue (y/n)?: "
if %query%==n goto exit
set /P "paraviewsrc=ParaView source tree path: "
goto loop
:exit
