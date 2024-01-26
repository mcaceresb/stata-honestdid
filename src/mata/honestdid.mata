cap mata mata drop HonestDiD()
cap mata mata drop HonestDiDParse()
cap mata mata drop HonestSensitivityResults()
cap mata mata drop HonestSensitivityHelper()
cap mata mata drop HonestOriginalCS()
cap mata mata drop HonestOriginalCSHelper()
cap mata mata drop HonestEventStudy()
cap mata mata drop _honestPrintCI()

* Done
* ----
*
* constructOriginalCS()      <-> HonestOriginalCS()
* createSensitivityResults() <-> HonestSensitivityResults()
*
* Outsourced
* ----------
*
* createSensitivityPlot() <-> Done from Stata via coefplot
*
* TODO
* ----
*
* createEventStudyPlot()
* createSensitivityPlot_relativeMagnitudes()
* createSensitivityResults_relativeMagnitudes()

* b             = "`b'"
* V             = "`vcov'"
* numPrePeriods = `numpreperiods'
* pre           = "`preperiodindices'"
* post          = "`postperiodindices'"
* l_vec         = "`l_vec'"
* Mvec          = "`mvec'"
* alpha         = `alpha'
* method        = "`method'"
* debug         = "`debug'"
* omit          = "`omit'"
* rm            = `relativeMagnitudes'
* grid_lb       = `grid_lb'
* grid_ub       = `grid_ub'
* gridPoints    = `gridPoints'

mata
struct HonestEventStudy {
    struct _honestResults colvector Results
    struct _honestResults scalar OG
    struct _honestOptions scalar options

    real matrix CI
    real vector betahat
    real matrix sigma
    real vector timeVec
    real scalar numPrePeriods
    real scalar numPostPeriods
    real vector prePeriodIndices
    real vector postPeriodIndices
    real vector open
}

struct HonestEventStudy scalar HonestDiD(string scalar b,
                                         string scalar V,
                                         real scalar numPrePeriods,
                                         string scalar pre,
                                         string scalar post,
                                         string scalar l_vec,
                                         string scalar Mvec,
                                         real scalar alpha,
                                         string scalar method,
                                         string scalar debug,
                                         string scalar omit,
                                         real scalar rm,
                                         real scalar grid_lb,
                                         real scalar grid_ub,
                                         real scalar gridPoints)
{
    struct _honestOptions scalar options
    struct HonestEventStudy scalar results

    results = HonestDiDParse(b, V, numPrePeriods, pre, post, l_vec, Mvec, alpha,
                             method, debug, omit, rm, grid_lb, grid_ub, gridPoints)

    results.Results = HonestSensitivityResults(results, results.options)
    results.OG      = HonestOriginalCS(results, results.options)
    results.CI      = _honestSensitivityCIMatrix(results.Results, results.OG)
    results.open    = _honestSensitivityCIOpen(results.Results, results.OG)
    return(results)
}

struct HonestEventStudy scalar HonestDiDParse(string scalar b,
                                              string scalar V,
                                              real scalar numPrePeriods,
                                              string scalar pre,
                                              string scalar post,
                                              string scalar l_vec,
                                              string scalar Mvec,
                                              real scalar alpha,
                                              string scalar method,
                                              string scalar debug,
                                              string scalar omit,
                                              real scalar rm,
                                              real scalar grid_lb,
                                              real scalar grid_ub,
                                              real scalar gridPoints)
{
    struct _honestOptions scalar options
    struct HonestEventStudy scalar results
    real vector selomit, sel
    real scalar Mub

    stata(sprintf("_ms_omit_info %s", b))
    selomit = selectindex(!editvalue(st_matrix("r(omit)"), omit == "", 0))
    if ( omit != "" ) {
        if ( rows(st_matrix(b)) > cols(st_matrix(b)) ) {
            printf("-omit- requires b() to be a cow vector; option ignored\n")
            omit    = ""
            selomit = 1..length(st_matrix(b))
        }
    }

    if ( rows(st_matrix(V)) != cols(st_matrix(V)) ) {
        errprintf("vcov() is not a square matrix\n")
        _error(198)
    }

    if ( max((rows(st_matrix(b)), cols(st_matrix(b)))) != rows(st_matrix(V)) ) {
        errprintf("b() and vcov() not conformable\n")
        _error(198)
    }

    results = HonestEventStudy()
    if ( numPrePeriods > 0 ) {
        results.betahat = rowshape(st_matrix(b), 1)[selomit]
        results.sigma   = st_matrix(V)[selomit, selomit]
        results.numPrePeriods     = numPrePeriods
        results.numPostPeriods    = length(results.betahat)-numPrePeriods
        results.prePeriodIndices  = 1..numPrePeriods
        results.postPeriodIndices = (numPrePeriods+1)..length(results.betahat)
    }
    else {
        results.prePeriodIndices  = strtoreal(tokens(pre))
        results.postPeriodIndices = strtoreal(tokens(post))
        sel = results.prePeriodIndices, results.postPeriodIndices
        results.betahat = rowshape(st_matrix(b), 1)[selomit][sel]
        results.sigma   = st_matrix(V)[selomit, selomit][sel, sel]
        results.numPrePeriods  = length(results.prePeriodIndices)
        results.numPostPeriods = length(results.postPeriodIndices)
    }

    options.omit               = (omit != "")
    options.debug              = debug != ""
    options.rm                 = rm
    options.relativeMagnitudes = rm? "rm": ""
    options.method             = method
    options.alpha              = alpha
    options.grid_lb            = grid_lb
    options.grid_ub            = grid_ub
    options.gridPoints         = gridPoints
    if ( Mvec == "" ) {
        if ( rm ) {
            options.Mvec = _honestLinspace(0, 2, 10)[2..10]
        }
        else if ( results.numPrePeriods == 1) {
            options.Mvec = _honestLinspace(0, results.sigma[1, 1], 10)
        }
        else {
            Mub = _honestSDUpperBoundMpre(results.betahat, results.sigma, results.numPrePeriods, options.alpha)
            options.Mvec = _honestLinspace(0, Mub, 10)
        }
    }
    else {
        options.Mvec = strtoreal(tokens(Mvec))
        if ( missing(options.Mvec) ) {
            if ( length(st_matrix(Mvec)) == 0 ) {
                errprintf("HonestDiD(): ")
                errprintf("Unable to parse option mvec(); must be number list of vector name\n")
                _error(198)
            }
            else {
                options.Mvec = rowshape(st_matrix(Mvec), 1)
            }
        }
    }

    options.Mvec    = rowshape(sort(colshape(options.Mvec, 1), 1), 1)
    options.l_vec   = (l_vec == "")? _honestBasis(1, results.numPostPeriods): colshape(st_matrix(l_vec), 1)

    if ( min(options.Mvec) < 0 ) {
        errprintf("M must be greater than or equal to 0 (see mvec() option)\n")
        _error(198)
    }
    else if ( (min(options.Mvec) == 0) & rm ) {
        printf("Warning: M = 0 with Delta^RM imposes exact parallel trends in the\n")
        printf("post-treatment period, even if pre-treatment parallel trends is violated\n")
    }
    results.options = options
    return(results)
}

struct _honestResults colvector function HonestSensitivityResults(
    struct HonestEventStudy scalar EventStudy, struct _honestOptions scalar options)
{
    return(HonestSensitivityHelper(EventStudy.betahat,
                                   EventStudy.sigma,
                                   EventStudy.numPrePeriods,
                                   EventStudy.numPostPeriods,
                                   options))
}

struct _honestResults colvector function HonestSensitivityHelper(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    struct _honestOptions scalar options)
{
    struct _honestResults colvector Results
    struct _flciResults scalar temp_flci
    string scalar biasDirection, method, monotonicityDirection, Delta, hybrid_flag, bound
    real colvector l_vec
    real vector Mvec
    real scalar m, alpha, rm, debug, grid_lb, grid_ub, gridPoints, open
    real matrix temp_mat

    _honestWarnIfNotSymmPSD(sigma)

    alpha      = options.alpha
    Mvec       = options.Mvec
    l_vec      = options.l_vec
    method     = options.method
    rm         = options.rm
    debug      = options.debug
    grid_lb    = options.grid_lb
    grid_ub    = options.grid_ub
    gridPoints = options.gridPoints

    // TODO: xx hard-coded
    biasDirection         = ""
    monotonicityDirection = ""
    Results               = J(length(Mvec), 1, _honestResults())
    bound                 = "deviation from linear trend"
    bound                 = "deviation from parallel trends"

    Delta = ""
    if ( rm ) {
        if (method == "C-LF") {
            hybrid_flag = "LF"
        }
        else if (method == "Conditional") {
            hybrid_flag = "ARP"
        }
        else if (method == "") {
            hybrid_flag = "LF"
            method = "C-LF"
        }
        else {
            errprintf("HonestSensitivityHelper(): ")
            errprintf("Method with relativeMagnitudes must be either '' (blank), Conditional or C-LF.\n")
            _error(198)
        }
    }

    if ( biasDirection == "" & monotonicityDirection == "" ) {
        Delta = rm? "DeltaRM": "DeltaSD"
        if ( !rm & method == "" ) {
            method = "FLCI"
        }
        // computeConditionalCS_DeltaSD
    }
    else if ( (biasDirection != "") & (monotonicityDirection == "") ) {
        if ( !rm & method == "" ) {
            method = "C-F"
        }
        if (biasDirection == "positive") {
            Delta = rm? "DeltaRMPB": "DeltaSDPB"
        }
        else {
            Delta = rm? "DeltaRMNB": "DeltaSDNB"
        }
        // computeConditionalCS_DeltaSDB
        errprintf("HonestSensitivityHelper(): ")
        errprintf("Method %s with Delta = %s not implemented\n", method, Delta)
        _error(198)
    }
    else if ( (biasDirection == "") & (monotonicityDirection != "") ) {
        if ( !rm & method == "" ) {
            method = "C-F"
        }
        if ( monotonicityDirection == "increasing" ) {
            Delta = rm? "DeltaRMI": "DeltaSDI"
        }
        else {
            Delta = rm? "DeltaRMD": "DeltaSDD"
        }
        // computeConditionalCS_DeltaSDM
        errprintf("HonestSensitivityHelper(): ")
        errprintf("Method %s with Delta = %s not implemented\n", method, Delta)
        _error(198)
    }
    else {
        errprintf("HonestSensitivityHelper(): ")
        errprintf("Select either a shape restriction or sign restriction (not both).\n")
        _error(198)
    }

    if ( method == "FLCI" ) {
        if ( biasDirection != "" ) {
            errprintf("HonestSensitivityHelper(): ")
            errprintf("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!\n")
        }
        if ( monotonicityDirection != "" ) {
            errprintf("HonestSensitivityHelper(): ")
            errprintf("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!\n")
        }
    }
    else if ( method == "Conditional" ) {
        hybrid_flag = "ARP"
    }
    else if ( method == "C-F" ) {
        hybrid_flag = "FLCI"
    }
    else if ( method == "C-LF" ) {
        hybrid_flag = "LF"
    }
    else {
        errprintf("HonestSensitivityHelper(): ")
        errprintf("Method must equal one of: FLCI, Conditional, C-F or C-LF")
        _error(198)
    }

    open = (hybrid_flag == "FLCI")? numPostPeriods == 1: 1
    if ( method == "FLCI" ) {
        for (m = 1; m <= length(Mvec); m++) {
            temp_flci = _flciFindOptimal(betahat,
                                         sigma,
                                         numPrePeriods,
                                         numPostPeriods,
                                         Mvec[m],
                                         debug,
                                         l_vec,
                                         alpha)

            Results[m].lb     = temp_flci.FLCI[1]
            Results[m].ub     = temp_flci.FLCI[2]
            Results[m].method = method
            Results[m].Delta  = Delta
            Results[m].M      = Mvec[m]
            Results[m].lopen  = 0
            Results[m].uopen  = 0
        }
    }
    else if ( Delta == "DeltaSD" ) {
        for (m = 1; m <= length(Mvec); m++) {
            temp_mat = _honestSDConditionalCS(betahat,
                                              sigma,
                                              numPrePeriods,
                                              numPostPeriods,
                                              debug,
                                              l_vec,
                                              Mvec[m],
                                              alpha,
                                              hybrid_flag,
                                              grid_lb,
                                              grid_ub,
                                              gridPoints)

            Results[m].lb     = min(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].ub     = max(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].method = method
            Results[m].Delta  = Delta
            Results[m].M      = Mvec[m]
            Results[m].lopen  = open & (temp_mat[1, 2] == 1) + 4 * all(temp_mat[., 2] :== 0)
            Results[m].uopen  = open & (temp_mat[gridPoints, 2] == 1) + 4 * all(temp_mat[., 2] :== 0)
        }
    }
    else if ( Delta == "DeltaRM" ) {
        for (m = 1; m <= length(Mvec); m++) {
            temp_mat = _honestRMConditionalCS(betahat,
                                              sigma,
                                              numPrePeriods,
                                              numPostPeriods,
                                              debug,
                                              l_vec,
                                              Mvec[m],
                                              alpha,
                                              hybrid_flag,
                                              grid_lb,
                                              grid_ub,
                                              gridPoints)

            Results[m].lb     = min(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].ub     = max(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].method = method
            Results[m].Delta  = Delta
            Results[m].M      = Mvec[m]
            Results[m].lopen  = (open & (temp_mat[1, 2] == 1)) + 4 * all(temp_mat[., 2] :== 0)
            Results[m].uopen  = (open & (temp_mat[gridPoints, 2] == 1)) + 4 * all(temp_mat[., 2] :== 0)
        }
    }
    else {
        errprintf("HonestSensitivityHelper(): ")
        errprintf("Method %s with Delta = %s not implemented\n", method, Delta)
        _error(198)
    }

    options.method = method
    options.Delta  = Delta
    return(Results)
}

struct _honestResults scalar function HonestOriginalCS(
    struct HonestEventStudy scalar EventStudy, struct _honestOptions scalar options)
{
    return(HonestOriginalCSHelper(EventStudy.betahat,
                                  EventStudy.sigma,
                                  length(EventStudy.prePeriodIndices),
                                  length(EventStudy.postPeriodIndices),
                                  options))
}

struct _honestResults scalar function HonestOriginalCSHelper(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    struct _honestOptions scalar options)
{
    struct _honestResults scalar OG
    real vector sel, lb, ub, l_vec
    real scalar n, stdError, alpha

    alpha    = options.alpha
    l_vec    = options.l_vec
    n        = numPrePeriods + numPostPeriods
    sel      = (numPrePeriods+1)::n
    stdError = sqrt(l_vec' * sigma[sel, sel] * l_vec)

    lb = l_vec' * betahat[sel]' - invnormal(1-alpha/2)*stdError
    ub = l_vec' * betahat[sel]' + invnormal(1-alpha/2)*stdError

    OG.lb     = lb
    OG.ub     = ub
    OG.method = "Original"
    OG.Delta  = ""
    OG.M      = .
    OG.lopen  = 0
    OG.uopen  = 0

    return(OG)
}


void function _honestPrintCI(struct HonestEventStudy scalar EventStudy)
{
    real scalar i
    string vector asterisks
    asterisks = "", "(-)", "(+)", "(*)", J(1, 2, ("(!)", "", "", "")), "(!)"

    printf("\n")
    printf("|    M    |   lb   |   ub   |\n")
    printf("| ------- | ------ | ------ |\n")
    for (i = 1; i <= rows(EventStudy.CI); i++) {
        if ( i == 1 ) {
            printf("| %7.4f | %6.3f | %6.3f | (Original)\n",
                   EventStudy.CI[i, 1],
                   EventStudy.CI[i, 2],
                   EventStudy.CI[i, 3])
        }
        else {
            printf("| %7.4f | %6.3f | %6.3f | %s\n",
                   EventStudy.CI[i, 1],
                   EventStudy.CI[i, 2],
                   EventStudy.CI[i, 3],
                   asterisks[EventStudy.open[i]+1])
        }
    }
    printf("(method = %s, Delta = %s, alpha = %5.3f)\n",
           EventStudy.options.method,
           EventStudy.options.Delta,
           EventStudy.options.alpha)

    if ( any(EventStudy.open :== 1) ) {
        errprintf("(-) CI is open at lower endpoint; CI length may not be accurate.\n")
    }
    if ( any(EventStudy.open :== 2) ) {
        errprintf("(+) CI is open at upper endpoint; CI length may not be accurate.\n")
    }
    if ( any(EventStudy.open :== 3) ) {
        errprintf("(*) CI is open at both endpoints; CI length may not be accurate.\n")
    }
    if ( any(EventStudy.open :== 12) ) {
        errprintf("(!) Could not find upper and lower bounds using the default grid.\n")
    }
    if ( any(EventStudy.open ) ) {
        errprintf("Try expanding the grid using the grid_lb() and grid_ub() options\n")
    }
}
end
