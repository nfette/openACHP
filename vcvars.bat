@REM Call or it won't come back here
@Call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86 10.0.15063.0
@REM @Call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86 8.1

@if "%WindowsSdkVersion:~0,3%" == "10." (
  @echo You are using Windows SDK matching version 10.*
  @echo so let's add the correct SDK\bin to your path,
  @echo and then you can call rc.exe, etc.
  call :setWindowsSdkActualBin
) else (
  @echo You are using Windows SDK that does not match version 10.*
  @echo No change required to your path.
)
@REM @endlocal
exit /B 0

:setWindowsSdkActualBin
@set "WindowsSdkActualBin="
@REM eg. C:\Program Files (x86)\Windows Kits\10\bin\10.0.15063.0\x86
@set "WindowsSdkActualBin=%WindowsSdkDir%bin\%WindowsSDKVersion%x86"
@REM @echo "%WindowsSdkActualBin%"
@set "PATH=%WindowsSdkActualBin%;%PATH%"
exit /B 0

:end