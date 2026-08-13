; TEagle Windows installer (Inno Setup 6).
; Build:  ISCC.exe /DMyAppVersion=<version> installer\teagle.iss   (injected by build_installer.ps1)
; Per-user install (no admin/UAC), self-contained bundle from dist\TEagle.
; Detects a previous version and upgrades in place; offers an opt-in clean install
; that also removes TEagle-generated environments and data (app data + WSL annotation backend).

#ifndef MyAppVersion
  #error "MyAppVersion must be injected by installer/build_installer.ps1"
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
Name: "cleaninstall"; Description: "Clean install — remove previous TEagle settings and cache; remove TEagle-owned WSL data only from a valid recorded distribution"; GroupDescription: "Installation type:"; Flags: unchecked

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

{ A distro name is read from a writable per-user marker and later used as a wsl.exe argument. Accept only
  the conservative character set used by registered TEagle distro names; an invalid marker must make
  cleanup stop rather than alter the command line or select an unrelated distribution. }
function IsSafeDistroName(const Distro: String): Boolean;
var
  I: Integer;
begin
  Result := False;
  if (Length(Distro) < 1) or (Length(Distro) > 128) then
    exit;
  for I := 1 to Length(Distro) do
  begin
    if Pos(Copy(Distro, I, 1), 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-') = 0 then
      exit;
  end;
  Result := True;
end;

function TryReadRecordedWslDistro(var Distro: String): Boolean;
var
  marker: String;
  raw: AnsiString;
begin
  Result := False;
  Distro := '';
  marker := ExpandConstant('{localappdata}\TEagle\wsl_distro.txt');
  if not FileExists(marker) then
  begin
    Log('Skipping TEagle WSL cleanup: no recorded distro marker exists.');
    exit;
  end;
  if not LoadStringFromFile(marker, raw) then
  begin
    Log('Skipping TEagle WSL cleanup: the recorded distro marker could not be read.');
    exit;
  end;
  Distro := Trim(String(raw));
  if not IsSafeDistroName(Distro) then
  begin
    Log('Skipping TEagle WSL cleanup: the recorded distro name is empty or invalid.');
    Distro := '';
    exit;
  end;
  Result := True;
end;

{ Remove only TEagle-owned WSL environments, caches, and logs. The micromamba root and binary may be
  shared with another application, so a TEagle clean install must never delete either one. The caller
  captures and validates the app-recorded distro before deleting LocalAppData. No default-distro fallback
  is permitted for this destructive operation, and this procedure never unregisters a distribution. }
function CleanWslBackend(const Distro: String): Boolean;
var
  rc: Integer;
  wsl: String;
begin
  Result := False;
  if not IsSafeDistroName(Distro) then
  begin
    Log('Skipping TEagle WSL cleanup: no valid recorded distro was supplied.');
    exit;
  end;
  wsl := WslExePath();
  if not FileExists(wsl) then
  begin
    Log('TEagle WSL cleanup failed: wsl.exe was not found.');
    exit;
  end;
  if not Exec(wsl, '-d "' + Distro + '" -- bash -lc "rm -rf ~/micromamba/envs/te ~/micromamba/envs/teagle-vrna ~/teagle_genomes ~/teagle_libs ~/teagle_wsl_install.sh ~/teagle_wsl_install.log ~/teagle_genome_prepare.log ~/teagle_genome_annotate.log ~/.teagle_install.lock"',
              '', SW_HIDE, ewWaitUntilTerminated, rc) then
  begin
    Log('TEagle WSL cleanup failed: wsl.exe could not be launched.');
    exit;
  end;
  if rc <> 0 then
  begin
    Log(Format('TEagle WSL cleanup failed in distribution %s (exit code %d).', [Distro, rc]));
    exit;
  end;
  Result := True;
end;

function CleanTeagleData(): Boolean;
var
  Distro: String;
  HasRecordedDistro: Boolean;
  WslCleanupOk: Boolean;
  WindowsCleanupOk: Boolean;
  MarkerRestored: Boolean;
  marker: String;
  FailureText: String;
begin
  Result := False;
  FailureText := '';
  { Capture the only authorized WSL cleanup target before its marker is removed. }
  HasRecordedDistro := TryReadRecordedWslDistro(Distro);
  WslCleanupOk := not HasRecordedDistro;
  if HasRecordedDistro then
    WslCleanupOk := CleanWslBackend(Distro);
  WindowsCleanupOk := DelTree(ExpandConstant('{localappdata}\TEagle'), True, True, True);
  if not WindowsCleanupOk then
    FailureText := 'Windows-side TEagle settings or cache could not be removed. Close TEagle and any program using its files, then retry.';
  if HasRecordedDistro and (not WslCleanupOk) then
  begin
    { Preserve the validated target for a later repair instead of erasing the only safe retry record. }
    marker := ExpandConstant('{localappdata}\TEagle\wsl_distro.txt');
    MarkerRestored := ForceDirectories(ExtractFileDir(marker));
    if MarkerRestored then
      MarkerRestored := SaveStringToFile(marker, Distro, False);
    if FailureText <> '' then
      FailureText := FailureText + #13#10 + #13#10;
    if MarkerRestored then
      FailureText := FailureText + 'The optional WSL backend could not be cleaned. No other distribution was touched. '
        + 'The validated target (' + Distro + ') was retained for a later retry.'
    else
    begin
      Log('WARNING: could not restore the TEagle WSL cleanup marker after cleanup failed.');
      FailureText := FailureText + 'The optional WSL backend could not be cleaned. No other distribution was touched, '
        + 'but the retry marker could not be restored. The validated target was ' + Distro + '; inspect the install log '
        + 'before removing TEagle data manually from that distribution.';
    end;
  end;
  Result := WindowsCleanupOk and WslCleanupOk;
  if FailureText <> '' then
    SuppressibleMsgBox(FailureText, mbError, MB_OK, IDOK);
end;

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if CurStep = ssInstall then
  begin
    if WizardIsTaskSelected('cleaninstall') then
    begin
      if not CleanTeagleData() then
        RaiseException('Clean install stopped because the selected TEagle data could not all be removed. Correct the reported problem and retry.');
    end;
  end;
end;

function UpdateReadyMemo(Space, NewLine, MemoUserInfoInfo, MemoDirInfo, MemoTypeInfo, MemoComponentsInfo, MemoGroupInfo, MemoTasksInfo: String): String;
var
  S: String;
begin
  if PrevVersion <> '' then
    S := 'Existing TEagle detected: version ' + PrevVersion + NewLine
       + 'Action: upgrade to {#MyAppVersion} (settings and downloaded TEagle databases are kept unless clean install is selected).' + NewLine + NewLine
  else
    S := 'Fresh install of TEagle {#MyAppVersion}.' + NewLine + NewLine;
  if WizardIsTaskSelected('cleaninstall') then
    S := S + 'CLEAN INSTALL selected: previous Windows settings and cache will be removed. TEagle-owned WSL data is removed only when cleanup succeeds in a valid recorded distribution; Linux cleanup is skipped without a valid record, and a failed cleanup retains that record for retry.' + NewLine + NewLine;
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
    if MsgBox('Also remove TEagle settings and cache, and attempt to remove TEagle-owned WSL data from its valid recorded distribution?'
              + #13#10 + '(Linux cleanup is skipped if no valid distro record exists. A failed cleanup retains its target for retry. Choose No to keep all data for a future reinstall.)',
              mbConfirmation, MB_YESNO or MB_DEFBUTTON2) = IDYES then
    begin
      if not CleanTeagleData() then
        Log('TEagle data cleanup was incomplete; uninstall continues after reporting the retained data.');
    end;
  end;
end;
