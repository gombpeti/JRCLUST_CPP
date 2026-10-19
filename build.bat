@echo off
REM ---------------------------------------------------------------------------
REM Configures and builds jrclust with MSVC + CUDA + Ninja.
REM
REM Pins the v142 toolset (MSVC 14.29). The 14.16 that also ships with this
REM VS2019 install is a v141 toolset and is not what we want for CUDA 12.8.
REM ---------------------------------------------------------------------------

setlocal

set "VSDEV=C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"
set "CMAKE=C:\Qt\Tools\CMake_64\bin\cmake.exe"
set "NINJA=C:\Qt\Tools\Ninja\ninja.exe"
set "ROOT=%~dp0"

if not exist "%VSDEV%" (
    echo ERROR: cannot find vcvars64.bat at "%VSDEV%"
    exit /b 1
)

call "%VSDEV%" -vcvars_ver=14.29 >nul
if errorlevel 1 (
    echo ERROR: failed to initialise the MSVC environment
    exit /b 1
)

echo Using compiler:
cl 2>&1 | findstr /C:"Version"

"%CMAKE%" -S "%ROOT%." -B "%ROOT%build" -G Ninja -DCMAKE_MAKE_PROGRAM="%NINJA%" -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 exit /b 1

"%CMAKE%" --build "%ROOT%build"
if errorlevel 1 exit /b 1

echo.
echo Built: %ROOT%build\bin\jrclust.exe
endlocal
