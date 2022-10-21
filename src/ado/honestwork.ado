capture program drop honestwork
program honestwork
    if "${HONEST_CALLER}" != "honestdid" {
        disp as err "honestwork is an internal honestdid command"
        exit 198
    }
    tempname results
    mata {
        `results' = _honestPLLLoad(st_sdata(1, "resfile"))
        (void) HonestDiDPLL(`results', st_data(., "mindex"))
        _honestPLLSave(st_sdata(., "parfile")[1], `results')
    }
end

if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) {
    local c_os_ macosx

    cap program drop honestosqp_plugin
    cap program drop honestecos_plugin

    local rc = 0
    cap program honestosqp_plugin, plugin using("honestosqp_`c_os_'.plugin")
    local rc = _rc | `rc'
    cap program honestecos_plugin, plugin using("honestecos_`c_os_'.plugin")
    local rc = _rc | `rc'

    if `rc' {
        local c_os_ macosxarm64
        cap program honestosqp_plugin, plugin using("honestosqp_`c_os_'.plugin")
        cap program honestecos_plugin, plugin using("honestecos_`c_os_'.plugin")
    }
}
else {
    local c_os_: di lower("`c(os)'")

    cap program drop honestosqp_plugin
    cap program honestosqp_plugin, plugin using("honestosqp_`c_os_'.plugin")

    cap program drop honestecos_plugin
    cap program honestecos_plugin, plugin using("honestecos_`c_os_'.plugin")
}
