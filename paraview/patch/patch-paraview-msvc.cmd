@echo off
echo Patching ParaView sources
set cdir=%cd%
set paraviewsrc=%1
:loop
if "%paraviewsrc%"=="" goto error
if not exist %paraviewsrc% goto error

%~d1
cd %paraviewsrc%

git apply --ignore-whitespace %cdir%/paraview-msvc-CMakeLists.txt.patch

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
