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

& $main "PATH" "C:\Program Files\TTK\bin;C:\Program Files\TTK-ParaView\bin"
& $main "PV_PLUGIN_PATH" "C:\Program Files\TTK\bin\plugins"
& $main "QT_PLUGIN_PATH" "C:\ProgramData\Anaconda3\Library\plugins"
& $main "PYTHONPATH ""C:\ProgramData\Anaconda3\Lib;C:\Program Files\TTK\bin\Lib\site-packages;C:\Program Files\TTK-ParaView\bin\Lib\site-packages"
& $main "CMAKE_PREFIX_PATH" "C:\Program Files\TTK\lib\cmake;C:\Program Files\TTK-ParaView\lib\cmake"
