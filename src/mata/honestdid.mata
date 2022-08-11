cap mata mata drop HonestDiD()
cap mata mata drop HonestSensitivityResults()
cap mata mata drop HonestSensitivityHelper()
cap mata mata drop HonestOriginalCS()
cap mata mata drop HonestOriginalCSHelper()
cap mata mata drop HonestEventStudy()
cap mata mata drop _honestPrintFLCI()

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

* b         = "`b'"
* V         = "`vcov'"
* reference = `referenceperiodindex'
* pre       = "`preperiodindices'"
* post      = "`postperiodindices'"
* l_vec     = "`l_vec'"
* Mvec      = "`mvec'"
* alpha     = 0.05

mata
struct HonestEventStudy {
    struct _honestResults colvector Results
    struct _honestResults scalar OG
    struct _honestOptions scalar options

    real matrix FLCI
    real vector betahat
    real matrix sigma
    real vector timeVec
    real scalar referencePeriod
    real vector prePeriodIndices
    real vector postPeriodIndices
}

struct HonestEventStudy scalar HonestDiD(string scalar b,
                                         string scalar V,
                                         real scalar reference,
                                         string scalar pre,
                                         string scalar post,
                                         string scalar l_vec,
                                         string scalar Mvec,
                                         real scalar alpha,
                                         string scalar method,
                                         string scalar debug,
                                         real scalar rm)
{
    struct _honestOptions scalar options
    struct HonestEventStudy scalar results
    real scalar Mub
    real vector sel

    results = HonestEventStudy()
    if ( reference > 0 ) {
        results.betahat = rowshape(st_matrix(b), 1)
        results.sigma   = st_matrix(V)
        results.prePeriodIndices  = 1..reference
        results.postPeriodIndices = (reference+1)..length(results.betahat)
    }
    else {
        results.prePeriodIndices  = strtoreal(tokens(pre))
        results.postPeriodIndices = strtoreal(tokens(post))
        sel = results.prePeriodIndices, results.postPeriodIndices
        results.betahat = rowshape(st_matrix(b), 1)[sel]
        results.sigma   = st_matrix(V)[sel, sel]
    }

    options.debug              = debug != ""
    options.rm                 = rm
    options.relativeMagnitudes = rm? "rm": ""
    options.method             = method
    options.alpha              = alpha
    if ( Mvec == "" ) {
        if ( rm ) {
            options.Mvec = _honestLinspace(0, 2, 10)
        }
        else if ( length(results.prePeriodIndices) == 1) {
            options.Mvec = _honestLinspace(0, results.sigma[1, 1], 10)
        }
        else {
            Mub = _honestSDUpperBoundMpre(results.betahat, results.sigma, length(results.prePeriodIndices), options.alpha)
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
    options.l_vec   = (l_vec == "")? _honestBasis(1, length(results.postPeriodIndices)): colshape(st_matrix(l_vec), 1)
    results.Results = HonestSensitivityResults(results, options)
    results.OG      = HonestOriginalCS(results, options)
    results.FLCI    = _honestSensitivityFLCIMatrix(results.Results, results.OG)
    results.options = options

    return(results)
}

struct _honestResults colvector function HonestSensitivityResults(
    struct HonestEventStudy scalar EventStudy, struct _honestOptions scalar options)
{
    return(HonestSensitivityHelper(EventStudy.betahat,
                                   EventStudy.sigma,
                                   length(EventStudy.prePeriodIndices),
                                   length(EventStudy.postPeriodIndices),
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
    real scalar m, alpha, rm, debug
    real matrix temp_mat

    alpha  = options.alpha
    Mvec   = options.Mvec
    l_vec  = options.l_vec
    method = options.method
    rm     = options.rm
    debug  = options.debug

    // TODO: xx hard-coded
    biasDirection         = ""
    monotonicityDirection = ""
    Results               = J(length(Mvec), 1, _honestResults())
    bound                 = "deviation from linear trend"
    bound                 = "deviation from parallel trends"

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
                                              hybrid_flag)

            Results[m].lb     = min(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].ub     = max(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].method = method
            Results[m].Delta  = Delta
            Results[m].M      = Mvec[m]
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
                                              hybrid_flag)

            Results[m].lb     = min(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].ub     = max(temp_mat[selectindex(temp_mat[., 2]), 1])
            Results[m].method = method
            Results[m].Delta  = Delta
            Results[m].M      = Mvec[m]
        }
    }
    else {
        errprintf("HonestSensitivityHelper(): ")
        errprintf("Method %s with Delta = %s not implemented\n", method, Delta)
        _error(198)
    }

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

    return(OG)
}


void function _honestPrintFLCI(struct HonestEventStudy scalar EventStudy)
{
    real scalar i
    printf("\n")
    printf("|    M    |   lb   |   ub   |\n")
    printf("| ------- | ------ | ------ |\n")
    for (i = 1; i <= rows(EventStudy.FLCI); i++) {
        if ( i == 1 ) {
            printf("| %7.4f | %6.3f | %6.3f | (Original)\n",
                   EventStudy.FLCI[i, 1], EventStudy.FLCI[i, 2], EventStudy.FLCI[i, 3])
        }
        else {
            printf("| %7.4f | %6.3f | %6.3f |\n",
                   EventStudy.FLCI[i, 1], EventStudy.FLCI[i, 2], EventStudy.FLCI[i, 3])
        }
    }
    printf("(method = %s, Delta = %s, alpha = %5.3f)\n",
           EventStudy.Results[1].method,
           EventStudy.Results[1].Delta,
           EventStudy.options.alpha)
}
end
