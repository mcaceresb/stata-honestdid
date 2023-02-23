if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) {
    local c_os_ macosx
}
else {
    local c_os_: di lower("`c(os)'")
}

* stata14-mp
cap program drop honesttest
program honesttest, plugin using("honesttest_`c_os_'.plugin")
plugin call honesttest
* exit, clear
