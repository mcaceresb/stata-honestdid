*! version 0.1.2 17Jul2022 Mauricio Caceres Bravo, mauricio.caceres.bravo@gmail.com
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
        alpha(passthru)                    /// 1 - level of the CI
                                           ///
        REFERENCEperiodindex(int 0)        /// index for the reference period
        PREperiodindices(numlist)          /// pre-period indices
        POSTperiodindices(numlist)         /// post-period indices
        method(str)                        /// FLCI (xx default), Conditional, C-F or C-LF
        MATAsave(str)                      /// Save resulting mata object
        coefplot                           /// Coefficient  plot
        cached                             /// Use cached results
        ciopts(str)                        ///
        *                                  /// Options for coefplot
    ]

    if "`matasave'" == "" local results HonestEventStudy
    else local results: copy local matasave

    * Specify only ONE of
    *
    *     referenceperiodindex()
    *
    * OR
    *
    *     preperiodindices() and postperiodindices()
    *
    * If referenceperiodindex() is specified, then
    *
    *     preperiodindices(1 to referenceperiodindex)
    *     postperiodindices(referenceperiodindex+1 to length(e(b)))
    *
    * are assumed.

    if ( "`cached'" != "" ) {
        local results `s(HonestEventStudy)'
        cap mata mata desc `results'
        if ( _rc ) {
            disp as err "Cached results not found"
            exit 198
        }
    }

    local dohonest = ("`b'`v'`l_vec'`alpha'`mvec'`method'`preperiodindices'`postperiodindices'" != "")
    local dohonest = `dohonest' | (`referenceperiodindex' != 0)
    if ( `dohonest' & ("`cached'" != "") ) {
        disp as txt "{bf:warning:} cached results ignored if modifications are specified"
        local cached
    }

    if ( `dohonest' | ("`cached'" == "") ) {
        HonestSanityChecks, b(`b') vcov(`vcov') l_vec(`l_vec') mvec(`mvec') method(`method') `alpha' ///
            reference(`referenceperiodindex') pre(`preperiodindices') post(`postperiodindices')
    }
    else {
        local alpha = 0.05
    }

    mata {
        if ( `dohonest' | ("`cached'" == "") ) {
            `results' = HonestDiD("`b'",                  ///
                                  "`vcov'",               ///
                                  `referenceperiodindex', ///
                                  "`preperiodindices'",   ///
                                  "`postperiodindices'",  ///
                                  "`l_vec'",              ///
                                  "`mvec'",               ///
                                  `alpha',                ///
                                  "`method'")
            _honestPrintFLCI(`results')
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
            `plotmatrix' = `results'.FLCI
            `labels' = strofreal(`plotmatrix'[., 1]')
            `labels'[1] = "Original"

            st_matrix("`cimatrix'", `plotmatrix'[., 2::3]')
            st_matrix("`dummycoef'", ((`plotmatrix'[., 3] :+ `plotmatrix'[., 2])/2)')
            st_local("cimatlab", invtokens(`labels'))
        }
    }

    if ( "`coefplot'" != "" ) {
        if ( `"`ciopts'"' == "" ) local ciopts ciopts(recast(rcap))
        else local ciopts ciopts(recast(rcap) `ciopts')

        matrix colnames `cimatrix'  = `cimatlab'
        matrix rownames `cimatrix'  = lb ub
        matrix colnames `dummycoef' = `cimatlab'
        coefplot matrix(`dummycoef'), ci(`cimatrix') vertical cionly yline(0) `ciopts' `options'
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
        REFERENCEperiodindex(int 0) /// index for the reference period
        PREperiodindices(numlist)   /// pre-period indices
        POSTperiodindices(numlist)  /// post-period indices
        method(str)                 /// FLCI (xx default), Conditional, C-F or C-LF
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
        c_local method: copy local method
    }

    if ((`referenceperiodindex' != 0) & "`preperiodindices'`postperiodindices'" != "") {
        disp as err "Specify only one of reference() or pre() and post()"
        exit 198
    }

    if ((`referenceperiodindex' == 0) & (("`preperiodindices'" == "") | ("`postperiodindices'" == "")) ) {
        disp as err "Specify either reference() or both pre() and post()"
        exit 198
    }

    if ( `referenceperiodindex' < 0 ) {
        disp as err "reference() must be positive"
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

    if ( (`referenceperiodindex' == 0) & ("`preperiodindices'" != "") & ("`postperiodindices'" != "") ) {
        local npre:  list sizeof preperiodindices
        local npost: list sizeof postperiodindices
        if ( max(`rowsb', `colsb') < (`npre' + `npost') ) {
            disp as err "Coefficient vector must be at least # pre + # post"
            exit 198
        }
    }

    if ( `referenceperiodindex' > 0 ) {
        if ( max(`rowsb', `colsb') <= `referenceperiodindex' ) {
            disp as err "Coefficient vector must be at least 1 + reference period index"
            exit 198
        }
    }

    if ( "`l_vec'" != "" ) confirm matrix `l_vec'
    if ( "`mvec'"  != "" ) {
        local 0, mvec(`mvec')
        cap syntax, mvec(numlist)
        if ( _rc == 0 ) c_local mvec: copy local mvec
    }

    c_local alpha: copy local alpha
    c_local b: copy local b
    c_local vcov: copy local vcov
end

if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
else local c_os_: di lower("`c(os)'")

cap program drop honestosqp_plugin
cap program honestosqp_plugin, plugin using("honestosqp_`c_os_'.plugin")

cap program drop honestecos_plugin
cap program honestecos_plugin, plugin using("honestecos_`c_os_'.plugin")
