cap mata mata drop HonestDiD()
cap mata mata drop HonestSensitivityResults()
cap mata mata drop HonestSensitivityHelper()
cap mata mata drop HonestOriginalCS()
cap mata mata drop HonestOriginalCSHelper()
cap mata mata drop HonestEventStudy()

cap mata mata drop _honestSensitivityFLCIMatrix()
cap mata mata drop _honestOptions()
cap mata mata drop _honestResults()
cap mata mata drop _honestBasis()
cap mata mata drop _honestLinspace()
cap mata mata drop _honestInverseIndex()
cap mata mata drop _honestHelper()
cap mata mata drop _honestwToLFn()
cap mata mata drop _honestlToWFn()
cap mata mata drop _honestPrintFLCI()
cap mata mata drop _honestUpperBoundMpre()
cap mata mata drop _honestCreateASD()

* Done
* ----
*
* DeltaSD_upperBound_Mpre()  <-> _honestUpperBoundMpre()
* .create_A_SD()             <-> _honestCreateASD()
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
* alpha     = 0.01

mata
struct _honestOptions {
    real scalar alpha
    real vector l_vec
    real vector Mvec
}

struct _honestResults {
    real scalar lb
    real scalar ub
    string scalar method
    string scalar Delta
    real scalar M
}

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
                                         | string scalar l_vec,
                                         string scalar Mvec,
                                         real scalar alpha)
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

    options.alpha = alpha
    if ( Mvec == "" ) {
        if ( length(results.prePeriodIndices) == 1) {
            options.Mvec = _honestLinspace(0, results.sigma[1, 1], 10)
        }
        else {
            Mub = _honestUpperBoundMpre(results.betahat, results.sigma, length(results.prePeriodIndices), options.alpha)
            options.Mvec = _honestLinspace(0, Mub, 10)
        }
    }
    else {
        options.Mvec = strtoreal(tokens(Mvec))
        if ( missing(options.Mvec) ) {
            if ( length(st_matrix(Mvec)) == 0 ) {
                errprintf("Unable to parse option mvec(); must be number list of vector name\n")
                _error(198)
            }
            else {
                options.Mvec = rowshape(st_matrix(Mvec), 1)
            }
        }
    }

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
    struct _flciResults scalar temp
    string scalar biasDirection, method, monotonicityDirection, Delta, hybrid_flag
    real colvector l_vec
    real vector Mvec
    real scalar m, alpha

    alpha = options.alpha
    Mvec  = options.Mvec
    l_vec = options.l_vec

    // TODO: xx some hard-coded stuff here
    method                = "FLCI"
    biasDirection         = ""
    monotonicityDirection = ""
    Results               = J(length(Mvec), 1, _honestResults())

    if ( biasDirection == "" & monotonicityDirection == "" ) {
        Delta = "DeltaSD"
        if ( method == "" ) {
            method = "FLCI"
        }
        // computeConditionalCS_DeltaSD
    }
    else if ( (biasDirection != "") & (monotonicityDirection == "") ) {
        if ( method == "" ) {
            method = "C-F"
        }
        if (biasDirection == "positive") {
            Delta = "DeltaSDPB"
        }
        else {
            Delta = "DeltaSDNB"
        }
        // computeConditionalCS_DeltaSDB
    }
    else if ( (biasDirection == "") & (monotonicityDirection != "") ) {
        if ( method == "" ) {
            method = "C-F"
        }
        if ( monotonicityDirection == "increasing" ) {
            Delta = "DeltaSDI"
        }
        else {
            Delta = "DeltaSDD"
        }
        // computeConditionalCS_DeltaSDM
    }

    if ( method == "FLCI" ) {
        if ( biasDirection != "" ) {
            errprintf("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!\n")
        }
        if ( monotonicityDirection != "" ) {
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
        errprintf("Method must equal one of: FLCI, Conditional, C-F or C-LF")
        _error(198)
    }

    if ( method == "FLCI" ) {
        for (m = 1; m <= length(Mvec); m++) {
            temp = _flciFindOptimal(betahat,
                                    sigma,
                                    numPrePeriods,
                                    numPostPeriods,
                                    Mvec[m],
                                    l_vec,
                                    alpha)

            Results[m].lb     = temp.FLCI[1]
            Results[m].ub     = temp.FLCI[2]
            Results[m].method = method
            Results[m].Delta  = Delta
            Results[m].M      = Mvec[m]
        }
    }
    else {
        // computeConditionalCS_DeltaSD(betahat,
        //                              sigma,
        //                              numPrePeriods,
        //                              numPostPeriods,
        //                              Mvec[m],
        //                              l_vec,
        //                              alpha,
        //                              hybrid_flag = "ARP")
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

real matrix function _honestSensitivityFLCIMatrix(
    struct _honestResults colvector robustResults,
    struct _honestResults scalar originalResults,
    | real scalar rescaleFactor, real scalar maxM)
{
    if ( args() < 3 ) rescaleFactor = 1
    if ( args() < 4 ) maxM = .

    real colvector M, lb, ub, sel
    real scalar m

    M  = J(rows(robustResults), 1, .)
    lb = J(rows(robustResults), 1, .)
    ub = J(rows(robustResults), 1, .)
    for (m = 1; m <= length(M); m++) {
        M[m]  = robustResults[m].M
        lb[m] = robustResults[m].lb
        ub[m] = robustResults[m].ub
    }

    // Filter out observations above maxM (after rescaling)
    sel = selectindex(M :<= maxM)
    M   = M[sel]
    lb  = lb[sel]
    ub  = ub[sel]

    // Set M for OLS to be the min M in robust results minus the gap
    // between Ms in robust
    M  = (originalResults.M  \ M)  * rescaleFactor
    lb = (originalResults.lb \ lb) * rescaleFactor
    ub = (originalResults.ub \ ub) * rescaleFactor

    return(M, lb, ub)
}

void function _honestPrintFLCI(struct HonestEventStudy scalar EventStudy)
{
    real scalar i
    printf("\n")
    printf("|   M    |   lb   |   ub   |\n")
    printf("| ------ | ------ | ------ |\n")
    for (i = 1; i <= rows(EventStudy.FLCI); i++) {
        if ( i == 1 ) {
            printf("| %6.2f | %6.3f | %6.3f | (Original)\n",
                   EventStudy.FLCI[i, 1], EventStudy.FLCI[i, 2], EventStudy.FLCI[i, 3])
        }
        else {
            printf("| %6.2f | %6.3f | %6.3f |\n",
                   EventStudy.FLCI[i, 1], EventStudy.FLCI[i, 2], EventStudy.FLCI[i, 3])
        }
    }
    printf("(DeltaSD is hard-coded; alpha = %5.3f)\n", EventStudy.options.alpha)
}
end

mata
real colvector function _honestBasis(real scalar index, real scalar size)
{
    real colvector v
    v = J(size, 1, 0)
    v[index] = 1
    return(v)
}

real vector _honestLinspace(real scalar from, real scalar to, real scalar n) {
    real scalar step
    step = (to - from) / (n - 1)
    return(from :+ step :* (0..(n - 1)))
}

real matrix function _honestInverseIndex(real vector ix, real scalar n)
{
    real vector sel
    if (rows(ix) > cols(ix)) {
        sel = J(n, 1, 1)
        sel[ix] = J(length(ix), 1, 0)
    }
    else {
        sel = J(1, n, 1)
        sel[ix] = J(1, length(ix), 0)
    }
    return(selectindex(sel))
}

// Converts vector from w space to l space.
real vector function _honestwToLFn(real vector w) {
    real scalar numPrePeriods, col
    real matrix WtoLPreMat

    numPrePeriods = length(w)
    WtoLPreMat = I(numPrePeriods)
    if (numPrePeriods == 1) {
        WtoLPreMat = 1
    }
    else {
        for (col = 1; col < numPrePeriods; col++) {
            WtoLPreMat[col + 1, col] = -1
        }
    }
    return(WtoLPreMat * w)
}

// Converts vector from l space to w space
real vector function _honestlToWFn(real vector l_vec) {
    real scalar numPostPeriods, col
    real matrix lToWPostMat

    numPostPeriods = length(l_vec)
    lToWPostMat    = I(numPostPeriods)
    if (numPostPeriods == 1) {
        lToWPostMat = 1
    }
    else {
        for (col = 1; col < numPostPeriods; col++) {
            lToWPostMat[col + 1, col] = 1
        }
    }
    return(lToWPostMat * l_vec)
}

// TODO: xx _honestHelper has not been thoroughly tested
real scalar function _honestHelper(real scalar s, real colvector l_vec, real scalar numPostPeriods)
{
    return(abs((1..s) * l_vec[(numPostPeriods - s + 1)::numPostPeriods]))
}

// This function constructs an upper-bound for M at the 1-alpha
// level based on the observed pre-period coefficients.
real scalar function _honestUpperBoundMpre(real vector betahat,
                                           real matrix sigma,
                                           real scalar numPrePeriods,
                                           real scalar alpha) {

    real vector prePeriod_coef, prePeriodCoefDiffs, seDiffs, upperBoundVec
    real matrix prePeriod_sigma, A_SD, prePeriodSigmaDiffs

    prePeriod_coef  = betahat[1::numPrePeriods]'
    prePeriod_sigma = sigma[1::numPrePeriods, 1::numPrePeriods]

    A_SD                = _honestCreateASD(numPrePeriods, 0)
    prePeriodCoefDiffs  = A_SD * prePeriod_coef
    prePeriodSigmaDiffs = A_SD * prePeriod_sigma * A_SD'
    seDiffs             = sqrt(diagonal(prePeriodSigmaDiffs))

    upperBoundVec = prePeriodCoefDiffs + invnormal(1 - alpha) * seDiffs
    return(max(upperBoundVec))
}

real matrix function _honestCreateASD(real scalar numPrePeriods,
                                      real scalar numPostPeriods,
                                      | real scalar postPeriodMomentsOnly) {

    real vector postPeriodIndices, prePeriodOnlyRows
    real matrix Atilde
    real scalar r

    if ( args() < 3 ) postPeriodMomentsOnly = 0

    // This function creates a matrix for the linear constraints that
    // \delta \in Delta^SD(M).  It implements this using the general
    // characterization of A, NOT the sharp characterization of the
    // identified set.
    //
    // Inputs:
    //   numPrePeriods         = number of pre-periods. This is an element of resultsObjects.
    //   numPostPeriods        = number of post-periods. This is an element of resultsObjects.
    //   postPeriodMomentsOnly = whether to exlude moments relating only to pre-period (which don't affect ID set)
    //
    // First construct matrix Atilde -- (numPrePeriods+numPostPeriods-2) x (numPrePeriods+numPostPeriods+1)
    // Note Atilde is just the positive moments; is not related to Atilde, the rotate matrix, in the paper
    // Note: Atilde initially includes t = 0. We then drop it.

    Atilde = J(numPrePeriods + numPostPeriods - 1, numPrePeriods + numPostPeriods + 1, 0)
    for (r = 1; r <= (numPrePeriods + numPostPeriods - 1); r++) {
        Atilde[r, r..(r+2)] = (1, -2, 1)
    }
    Atilde = Atilde[., _honestInverseIndex(numPrePeriods+1, cols(Atilde))]

    // If postPeriodMomentsOnly == 1, exclude moments that only involve pre-periods
    if ( postPeriodMomentsOnly ) {
        postPeriodIndices = (numPrePeriods + 1)::cols(Atilde)
        prePeriodOnlyRows = selectindex(rowsum(Atilde[., postPeriodIndices] :!= 0) :== 0)
        Atilde = Atilde[_honestInverseIndex(prePeriodOnlyRows, cols(Atilde)), .]
    }

    // Construct A = [Atilde; -Atilde]
    return(Atilde \ -Atilde)
}
end
