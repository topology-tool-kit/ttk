function append_to_env_var($var, $val) {
    # target the current user environment variables registry
    $dest = [EnvironmentVariableTarget]::User
    # existing value of the current environment variable
    $old_val = [Environment]::GetEnvironmentVariable($var, $dest)
    If ($null -eq $old_val) {
        # set environment variable to value
        $new_val = $val
    } else {
        # append new value to existing ones
        $new_val = $old_val + ";" + $val
    }
    [Environment]::SetEnvironmentVariable($var, $new_val, $dest)
}

function remove_from_env_var($var, $val) {
    # target the current user environment variables registry
    $dest = [EnvironmentVariableTarget]::User
    # existing value of the current environment variable
    $old_val = [Environment]::GetEnvironmentVariable($var, $dest)
    # remove $val from $old_val
    If ($null -ne $old_val) {
        $new_val = $old_val -replace [Regex]::Escape($val), ""
        # remove trailing semi-colon
        $new_val = $new_val -replace ";$", ""
        [Environment]::SetEnvironmentVariable($var, $new_val, $dest)
    }
}

If ($args[0] -eq "install") {
    $main = ${function:append_to_env_var}
} else {
    $main = ${function:remove_from_env_var}
}

$ttk_dir = $MyInvocation.MyCommand.Definition | Split-Path | Split-Path
$prefix = $ttk_dir | Split-Path
# assumption: TTK and TTK-ParaView installed in the same prefix
$pv_dir = "${prefix}\TTK-ParaView"

function find_conda_dir() {
    # assumption: Anaconda environment variables set
    $python_exe = Get-Command python
    If ($null -eq $python_exe) {
        Write-Output "Error: Python executable not accessible in PATH"
        exit
    }
    $conda_dir = $python_exe.Definition | Split-Path
    If (-not ($conda_dir -Match "conda")) {
        Write-Output "Error: Python not part of an Anaconda distribution"
        exit
    }
    return $conda_dir
}

$conda_dir = find_conda_dir

& $main "PATH" "${ttk_dir}\bin;${pv_dir}\bin"
& $main "PV_PLUGIN_PATH" "${ttk_dir}\bin\plugins"
& $main "QT_PLUGIN_PATH" "${conda_dir}\Library\plugins"
& $main "PYTHONPATH" "${conda_dir}\Lib;${ttk_dir}\bin\Lib\site-packages;${pv_dir}\bin\Lib\site-packages"
& $main "CMAKE_PREFIX_PATH" "${ttk_dir}\lib\cmake;${pv_dir}\lib\cmake"
