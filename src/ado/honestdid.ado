*! version 0.5.2 02Oct2022 Mauricio Caceres Bravo, mauricio.caceres.bravo@gmail.com
*! HonestDiD R to Stata translation

capture program drop honestdid
program honestdid, sclass
    version 14.1
    cap plugin call honestosqp_plugin, _plugin_check
    if ( _rc ) {
        disp as err "Failed to load ECOS plugin"
        exit _rc
    }
    cap plugin call honestecos_plugin, _plugin_check
    if ( _rc ) {
        disp as err "Failed to load OSQP plugin"
        exit _rc
    }

    syntax,                                ///
    [                                      ///
        b(str)                             /// name of coefficient matrix; default is e(b)
        vcov(str)                          /// name of vcov matrix; default is e(V)
        l_vec(str)                         /// Vector with parameters of interest (default is first period post event)
        mvec(str)                          /// Vector or list with with M-values
        grid_lb(str)                       /// Lower bound for grid
        grid_ub(str)                       /// Upper bound for grid
        gridPoints(str)                    /// Number of grid points
        alpha(passthru)                    /// 1 - level of the CI
                                           ///
        omit                               /// Omit levels parsing b vector column names
        delta(str)                         /// Delta to use (rm for relative magnitudes or sd)
        NUMPREperiods(int 0)               /// number of pre-treatment periods
        PREperiodindices(numlist)          /// pre-period indices
        POSTperiodindices(numlist)         /// post-period indices
        method(str)                        /// FLCI, Conditional, C-F or C-LF (default depends on rm)
        MATAsave(str)                      /// Save resulting mata object
        coefplot                           /// Coefficient  plot
        cached                             /// Use cached results
        colorspec(str asis)                /// special color handling
        ciopts(str)                        ///
        debug                              /// Print all the mata steps
        *                                  /// Options for coefplot
    ]

    if !inlist("`delta'", "rm", "sd", "") {
        disp as err "Only delta(rm) and delta(sd) available"
        exit 198
    }

    local rm = cond(inlist("`delta'", "rm", ""), "rm", "")
    local relativeMagnitudes = "`rm'" != ""
    if "`matasave'" == "" local results HonestEventStudy
    else local results: copy local matasave

    * Specify only ONE of
    *
    *     numpreperiods()
    *
    * OR
    *
    *     preperiodindices() and postperiodindices()
    *
    * If numpreperiods() is specified, then
    *
    *     preperiodindices(1 to numpreperiods)
    *     postperiodindices(numpreperiods+1 to length(e(b)))
    *
    * are assumed.

    local changegrid = 0
    if ( "`cached'" != "" ) {
        local results `s(HonestEventStudy)'
        cap mata mata desc `results'
        if ( _rc ) {
            disp as err "Cached results not found"
            exit 198
        }
        mata st_local("rmcached", `results'.options.relativeMagnitudes)
        if ( "`rm'" != "`rmcached'" ) {
            if ( "`rmcached'" != "" ) {
                disp as txt "{bf:note:} cached results uses relativeMagnitudes"
            }
            if ( "`rm'" != "" & "`delta'" != "" ) {
                disp as txt "{bf:warning:} cached results do not use relativeMagnitudes; option -rm- ignored"
            }
        }

        tempname gridcached gridlbcached gridubcached
        mata st_numscalar("`gridlbcached'", `results'.options.grid_lb)
        mata st_numscalar("`gridubcached'", `results'.options.grid_ub)
        mata st_numscalar("`gridcached'",   `results'.options.gridPoints)
        if "`grid_lb'"    != "" local changegrid = `changegrid' | ("`grid_lb'"    != "`=scalar(`gridlbcached')'")
        if "`grid_ub'"    != "" local changegrid = `changegrid' | ("`grid_ub'"    != "`=scalar(`gridubcached')'")
        if "`gridPoints'" != "" local changegrid = `changegrid' | ("`gridPoints'" != "`=scalar(`gridcached')'")
    }

    if "`grid_lb'"    == ""  local grid_lb .
    if "`grid_ub'"    == ""  local grid_ub .
    if "`gridPoints'" == ""  local gridPoints 1000

    if "`grid_lb'"    != "." confirm number `grid_lb'
    if "`grid_ub'"    != "." confirm number `grid_ub'
    if "`gridPoints'" != ""  confirm number `gridPoints'

    local dohonest = ("`b'`v'`l_vec'`alpha'`mvec'`method'`preperiodindices'`postperiodindices'" != "")
    local dohonest = `dohonest' | `changegrid' | (`numpreperiods' != 0)
    if ( `dohonest' & ("`cached'" != "") ) {
        disp as txt "{bf:warning:} cached results ignored if modifications are specified"
        local cached
    }

    if ( `dohonest' | ("`cached'" == "") ) {
        HonestSanityChecks, b(`b') vcov(`vcov') l_vec(`l_vec') mvec(`mvec') method(`method') `alpha' ///
            numpre(`numpreperiods') pre(`preperiodindices') post(`postperiodindices') `rm'
    }
    else {
        local alpha = 0.05
    }

    mata {
        if ( `dohonest' | ("`cached'" == "") ) {
            `results' = HonestDiD("`b'",                  ///
                                  "`vcov'",               ///
                                  `numpreperiods',        ///
                                  "`preperiodindices'",   ///
                                  "`postperiodindices'",  ///
                                  "`l_vec'",              ///
                                  "`mvec'",               ///
                                  `alpha',                ///
                                  "`method'",             ///
                                  "`debug'",              ///
                                  "`omit'",               ///
                                  `relativeMagnitudes',   ///
                                  `grid_lb',              ///
                                  `grid_ub',              ///
                                  `gridPoints')
            _honestPrintCI(`results')
        }
    }
    * `results'.timeVec = (2004, 2005, 2006, 2007, 2009, 2010, 2011, 2012)
    * `results'.referencePeriod = 2008

    tempname plotmatrix cimatrix dummycoef labels
    if ( "`coefplot'" != "" ) {
        cap which coefplot
        if ( _rc ) {
            disp as err "-coefplot- not found; required to plot CIs"
            exit _rc
        }
    }

    mata {
        if ( "`coefplot'" != "" ) {
            `plotmatrix' = `results'.CI
            `labels' = strofreal(`plotmatrix'[., 1]')
            `labels'[1] = "Original"

            st_matrix("`cimatrix'", `plotmatrix'[., 2::3]')
            st_matrix("`dummycoef'", ((`plotmatrix'[., 3] :+ `plotmatrix'[., 2])/2)')
            st_local("cimatlab", invtokens(`labels'))
        }
    }

    if ( "`coefplot'" != "" ) {
        matrix colnames `cimatrix'  = `cimatlab'
        matrix rownames `cimatrix'  = lb ub
        matrix colnames `dummycoef' = `cimatlab'

        if `"`colorspec'"' == "" local colorspec `""183 28 28" "13 71 161""'
        local optionsbak: copy local options
        local 0, `ciopts'
        syntax, [LColor(str) Color(str) *]
        local options: copy local optionsbak

        if ( `"`ciopts'"' == "" ) local ciopts recast(rcap)
        else local ciopts recast(rcap) `ciopts'

        tempname fullmatrix coefcall coefs
        matrix `fullmatrix' = `dummycoef' \ `cimatrix'
        mata `coefs' = J(1, `=colsof(`fullmatrix')', "")
        if ( "`lcolor'`color'" == "" ) {
            mata `coefcall' = `"(matrix(%s[1]), ci((2 3)) ciopt(lcolor("%s") `ciopts'))"'
        }
        else {
            mata `coefcall' = `"(matrix(%s[1]), ci((2 3)) ciopt(`ciopts'))"'
        }

        forvalues i = 1 / `=colsof(`fullmatrix')' {
            tempname x`i'
            local ix = min(`i', `:list sizeof colorspec')
            if ( "`lcolor'`color'" == "" ) {
                mata `coefs'[`i'] = sprintf(`coefcall', "`x`i''", tokens(st_local("colorspec"))[`ix'])
            }
            else {
                mata `coefs'[`i'] = sprintf(`coefcall', "`x`i''")
            }
            matrix `x`i'' = `fullmatrix'[1..rowsof(`fullmatrix'), `i']
        }
        mata st_local("matrices", invtokens(`coefs'))
        coefplot `matrices', vertical cionly yline(0) `options'

        * coefplot matrix(`dummycoef'), ci(`cimatrix') vertical cionly yline(0) `ciopts' `options'
        disp as err "({bf:warning:} horizontal distance in plot needn't be to scale)"
    }

    sreturn local HonestEventStudy = "`results'"
end

capture program drop HonestSanityChecks
program HonestSanityChecks
    syntax ,                        ///
    [                               ///
        b(str)                      /// name of coefficient matrix; default is e(b)
        vcov(str)                   /// name of vcov matrix; default is e(V)
        l_vec(str)                  /// Vector with parameters of interest (default is first period post event)
        mvec(str)                   /// Vector or list with with M-values
        alpha(real 0.05)            /// 1 - level of CI
        NUMPREperiods(int 0)        /// number of pre-treatment periods
        PREperiodindices(numlist)   /// pre-period indices
        POSTperiodindices(numlist)  /// post-period indices
        method(str)                 /// FLCI, Conditional, C-F or C-LF (default depends on rm)
        rm                          /// relative magnitudes
    ]

    local method = lower(`"`method'"')
    if ( `"`method'"' != "" ) {
             if ( `"`method'"' == "flci"        ) local method FLCI
        else if ( `"`method'"' == "conditional" ) local method Conditional
        else if ( `"`method'"' == "c-f"         ) local method C-F
        else if ( `"`method'"' == "c-lf"        ) local method C-LF
        else {
            disp as err "Method must equal one of: FLCI, Conditional, C-F or C-LF"
            exit 198
        }
    }

    if ( "`rm'" != "" ) {
        if "`method'" == "" local method C-LF
        if !inlist("`method'", "C-LF", "Conditional") {
            disp as err "Method '`method'' not allowed with RelativeMagnitudes"
            exit 198
        }
    }

    if ((`numpreperiods' != 0) & "`preperiodindices'`postperiodindices'" != "") {
        disp as err "Specify only one of numpre() or pre() and post()"
        exit 198
    }

    if ((`numpreperiods' == 0) & (("`preperiodindices'" == "") | ("`postperiodindices'" == "")) ) {
        disp as err "Specify either numpre() or both pre() and post()"
        exit 198
    }

    if ( `numpreperiods' < 0 ) {
        disp as err "numpre() must be positive"
        exit 198
    }

    tempname bb VV
    if ( "`b'" == "" ) {
        local b e(b)
        cap confirm matrix e(b)
        if ( _rc ) {
            disp as err "Last estimation coefficients not found; please specify vector."
            exit 198
        }
        matrix `bb' = e(b)
        local rowsb = rowsof(`bb')
        local colsb = colsof(`bb')
    }
    else {
        confirm matrix `b'
        local rowsb = rowsof(`b')
        local colsb = colsof(`b')
    }

    if ( "`vcov'" == "" ) {
        local vcov e(V)
        cap confirm matrix e(V)
        if ( _rc ) {
            disp as err "Last estimation vcov matrix not found; please specify matrix."
            exit 198
        }
        matrix `VV' = e(V)
        local rowsV = rowsof(`VV')
        local colsV = colsof(`VV')
    }
    else {
        confirm matrix `vcov'
        local rowsV = rowsof(`vcov')
        local colsV = colsof(`vcov')
    }

    if ( (`rowsb' != 1) & (`colsb' != 1) ) {
        disp as err "b() must be a vector; was `rowsb' x `colsb'"
        exit 198
    }

    if ( `rowsV' != `colsV' ) {
        disp as err "V() must be square"
        exit 198
    }

    if ( max(`rowsb', `colsb') != `rowsV' ) {
        disp as err "V() was `rowsV' x `colsV' not conformable with b() `rowsb' x `colsb'"
        exit 198
    }

    if ( (`numpreperiods' == 0) & ("`preperiodindices'" != "") & ("`postperiodindices'" != "") ) {
        local npre:  list sizeof preperiodindices
        local npost: list sizeof postperiodindices
        if ( max(`rowsb', `colsb') < (`npre' + `npost') ) {
            disp as err "Coefficient vector must be at least # pre + # post"
            exit 198
        }
    }

    if ( `numpreperiods' > 0 ) {
        if ( max(`rowsb', `colsb') <= `numpreperiods' ) {
            disp as err "Coefficient vector must be at least 1 + number of pre-treatment periods"
            exit 198
        }
    }

    if ( "`l_vec'" != "" ) confirm matrix `l_vec'
    if ( "`mvec'"  != "" ) {
        local 0, mvec(`mvec')
        cap syntax, mvec(numlist)
        if ( _rc == 0 ) c_local mvec: copy local mvec
    }

    c_local alpha:  copy local alpha
    c_local b:      copy local b
    c_local vcov:   copy local vcov
    c_local method: copy local method
end

if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
else local c_os_: di lower("`c(os)'")

cap program drop honestosqp_plugin
cap program honestosqp_plugin, plugin using("honestosqp_`c_os_'.plugin")

cap program drop honestecos_plugin
cap program honestecos_plugin, plugin using("honestecos_`c_os_'.plugin")
