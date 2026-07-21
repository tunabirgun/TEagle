; TEagle Windows installer (Inno Setup 6).
; Build:  ISCC.exe /DMyAppVersion=1.0.0 installer\teagle.iss   (version injected by build_installer.ps1)
; Per-user install (no admin/UAC), self-contained bundle from dist\TEagle.
; Detects a previous version and upgrades in place; offers an opt-in clean install
; that also wipes previously generated environments (app data + WSL annotation backend).

#ifndef MyAppVersion
  #define MyAppVersion "1.0.0"
#endif
#define MyAppName "TEagle"
#define MyAppPublisher "Tuna Birgun"
#define MyAppExeName "TEagle.exe"
#define MyAppId "{A7F3C2E9-5D41-4B8A-9E2F-1C6D8B3A0F51}"
; must stay byte-identical to SetCurrentProcessExplicitAppUserModelID(...) in app/native/main.py, so Windows
; binds the running window's taskbar button to this shortcut's icon. NEVER change the value (invalidates pins).
#define MyAppAumid "TEagle.desktop.2"

[Setup]
AppId={{#MyAppId}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
VersionInfoVersion={#MyAppVersion}
VersionInfoProductVersion={#MyAppVersion}
DefaultDirName={autopf}\{#MyAppName}
UsePreviousAppDir=yes
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
DisableReadyPage=no
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
OutputDir=..\dist
OutputBaseFilename=TEagle-Setup-{#MyAppVersion}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
CloseApplications=yes
RestartApplications=no
AppMutex=Global\TEagle_native_single_instance
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
LicenseFile=..\LICENSE
UninstallDisplayName={#MyAppName} {#MyAppVersion}
UninstallDisplayIcon={app}\{#MyAppExeName}
SetupIconFile=teagle.ico
; eagle brand logo in the wizard: top-right small image (interior pages) + left image (Welcome/Finish),
; each a DPI ladder so Inno picks the crisp size for the current scaling.
WizardSmallImageFile=wizard-small.bmp,wizard-small-83.bmp,wizard-small-110.bmp,wizard-small-138.bmp
WizardImageFile=wizard-large.bmp,wizard-large-246.bmp,wizard-large-328.bmp,wizard-large-410.bmp
WizardImageStretch=no
WizardImageBackColor=$FFFFFF

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Shortcuts:"
Name: "cleaninstall"; Description: "Clean install — remove previous TEagle settings, cache, and any downloaded databases / environments (including the WSL annotation backend)"; GroupDescription: "Installation type:"; Flags: unchecked

[Files]
Source: "..\dist\TEagle\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; AppUserModelID: "{#MyAppAumid}"
Name: "{userdesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon; AppUserModelID: "{#MyAppAumid}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch TEagle now"; Flags: nowait postinstall skipifsilent

[Code]
var
  PrevVersion: String;

function GetInstalledVersion(): String;
var
  v: String;
begin
  Result := '';
  if RegQueryStringValue(HKCU, 'Software\Microsoft\Windows\CurrentVersion\Uninstall\{#MyAppId}_is1', 'DisplayVersion', v) then
    Result := v
  else if RegQueryStringValue(HKLM, 'Software\Microsoft\Windows\CurrentVersion\Uninstall\{#MyAppId}_is1', 'DisplayVersion', v) then
    Result := v;
end;

function InitializeSetup(): Boolean;
begin
  PrevVersion := GetInstalledVersion();
  Result := True;
end;

{ resolve wsl.exe from a 32-bit setup process on 64-bit Windows (Sysnative alias), fall back to System32 }
function WslExePath(): String;
begin
  Result := ExpandConstant('{win}\Sysnative\wsl.exe');
  if not FileExists(Result) then
    Result := ExpandConstant('{sys}\wsl.exe');
end;

{ remove the generated WSL backend (env + Dfam) — targets the distro the app recorded, never unregisters it }
procedure CleanWslBackend();
var
  rc: Integer;
  wsl, distro, distroArg, marker: String;
  raw: AnsiString;
begin
  wsl := WslExePath();
  if not FileExists(wsl) then
    exit;
  distroArg := '';                                  { default distro if the app never recorded one }
  marker := ExpandConstant('{localappdata}\TEagle\wsl_distro.txt');
  if FileExists(marker) and LoadStringFromFile(marker, raw) then
  begin
    distro := Trim(String(raw));
    if distro <> '' then
      distroArg := '-d "' + distro + '" ';
  end;
  Exec(wsl, distroArg + '-- bash -lc "rm -rf ~/micromamba ~/bin/micromamba ~/teagle_wsl_install.log ~/.teagle_install.lock"',
       '', SW_HIDE, ewWaitUntilTerminated, rc);
end;

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if CurStep = ssInstall then
  begin
    if WizardIsTaskSelected('cleaninstall') then
    begin
      { wipe Windows-side app data: cache, timings, environment signature }
      DelTree(ExpandConstant('{localappdata}\TEagle'), True, True, True);
      { wipe the generated WSL annotation backend }
      CleanWslBackend();
    end;
  end;
end;

function UpdateReadyMemo(Space, NewLine, MemoUserInfoInfo, MemoDirInfo, MemoTypeInfo, MemoComponentsInfo, MemoGroupInfo, MemoTasksInfo: String): String;
var
  S: String;
begin
  if PrevVersion <> '' then
    S := 'Existing TEagle detected: version ' + PrevVersion + NewLine
       + 'Action: upgrade to {#MyAppVersion} (settings and downloaded databases are kept unless clean install is selected).' + NewLine + NewLine
  else
    S := 'Fresh install of TEagle {#MyAppVersion}.' + NewLine + NewLine;
  if WizardIsTaskSelected('cleaninstall') then
    S := S + 'CLEAN INSTALL selected: previous settings, cache, and downloaded environments (including the WSL backend) will be removed.' + NewLine + NewLine;
  if MemoDirInfo <> '' then S := S + MemoDirInfo + NewLine + NewLine;
  if MemoTasksInfo <> '' then S := S + MemoTasksInfo;
  Result := S;
end;

{ on uninstall, offer to also remove downloaded databases + settings + WSL backend }
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
begin
  { only prompt in an interactive uninstall; a silent uninstall keeps data (never auto-wipes the WSL backend).
    Default button is No, so an accidental Enter keeps the user's downloaded databases. }
  if (CurUninstallStep = usUninstall) and (not UninstallSilent()) then
  begin
    if MsgBox('Also remove TEagle settings, cache, and downloaded databases / environments (including the WSL annotation backend)?'
              + #13#10 + '(Choose No to keep them for a future reinstall.)',
              mbConfirmation, MB_YESNO or MB_DEFBUTTON2) = IDYES then
    begin
      DelTree(ExpandConstant('{localappdata}\TEagle'), True, True, True);
      CleanWslBackend();
    end;
  end;
end;
