capture program drop honestdid
program honestdid
    version 14.1
    plugin call honestosqp_plugin, _plugin_check
    plugin call honestecos_plugin, _plugin_check
end


if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
else local c_os_: di lower("`c(os)'")

cap program drop honestosqp_plugin
cap program honestosqp_plugin, plugin using("honestosqp_`c_os_'.plugin")

cap program drop honestecos_plugin
cap program honestecos_plugin, plugin using("honestecos_`c_os_'.plugin")
