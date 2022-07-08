cap mata mata drop HonestDiD()
cap mata mata drop HonestSensitivityResults()
cap mata mata drop HonestSensitivityHelper()
cap mata mata drop HonestOriginalCS()
cap mata mata drop HonestOriginalCSHelper()
cap mata mata drop HonestEventStudy()

cap mata mata drop _honestSensitivityFLCIMatrix()
cap mata mata drop _honestResults()
cap mata mata drop _honestBasis()
cap mata mata drop _honestLinspace()
cap mata mata drop _honestInverseIndex()
cap mata mata drop _honestHelper()
cap mata mata drop _honestwToLFn()
cap mata mata drop _honestPrintFLCI()

mata
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
                                         | string scalar l_vec)
{
    struct HonestEventStudy scalar results
    real vector l_vector, sel

    if ( args() < 6 ) l_vec = ""

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

    l_vector = (l_vec == "")? _honestBasis(1, length(results.postPeriodIndices)): colshape(st_matrix(l_vec), 1)
    results.Results = HonestSensitivityResults(results, l_vector)
    results.OG      = HonestOriginalCS(results, l_vector)
    results.FLCI    = _honestSensitivityFLCIMatrix(results.Results, results.OG)

    return(results)
}

struct _honestResults colvector function HonestSensitivityResults(
    struct HonestEventStudy scalar EventStudy, | real colvector l_vec)
{
    if ( args() < 2 ) l_vec = _honestBasis(1, length(EventStudy.postPeriodIndices))
    return(HonestSensitivityHelper(EventStudy.betahat,
                                   EventStudy.sigma,
                                   length(EventStudy.prePeriodIndices),
                                   length(EventStudy.postPeriodIndices),
                                   l_vec))
}

struct _honestResults colvector function HonestSensitivityHelper(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    | real colvector l_vec)
{
    struct _honestResults colvector Results
    struct _flciResults scalar temp

    real rowvector Mvec
    string scalar biasDirection, method, monotonicityDirection, Delta
    real scalar alpha, m

    if ( args() < 5 ) l_vec = _honestBasis(1, numPostPeriods)

    Mvec                  = (0..3) / 10
    biasDirection         = "negative"
    alpha                 = 0.05
    method                = "FLCI"
    monotonicityDirection = ""
    Delta                 = "DeltaSD"

    Results = J(length(Mvec), 1, _honestResults())
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

    return(Results)
}

struct _honestResults scalar function HonestOriginalCS(
    struct HonestEventStudy scalar EventStudy, | real colvector l_vec)
{
    if ( args() < 2 ) l_vec = _honestBasis(1, length(EventStudy.postPeriodIndices))
    return(HonestOriginalCSHelper(EventStudy.betahat,
                                  EventStudy.sigma,
                                  length(EventStudy.prePeriodIndices),
                                  length(EventStudy.postPeriodIndices),
                                  l_vec))
}

struct _honestResults scalar function HonestOriginalCSHelper(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    | real colvector l_vec,
    real scalar alpha)
{
    struct _honestResults scalar OG

    if ( args() < 5 ) l_vec = _honestBasis(1, numPostPeriods)
    if ( args() < 6 ) alpha = 0.05

    real vector sel, lb, ub
    real scalar n, stdError

    n = numPrePeriods + numPostPeriods
    sel = (numPrePeriods+1)::n
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

    real scalar m, Mgap, Mmin
    real colvector M, lb, ub, sel

    M  = J(rows(robustResults), 1, .)
    lb = J(rows(robustResults), 1, .)
    ub = J(rows(robustResults), 1, .)
    for (m = 1; m <= length(M); m++) {
        M[m]  = robustResults[m].M
        lb[m] = robustResults[m].lb
        ub[m] = robustResults[m].ub
    }

    // Set M for OLS to be the min M in robust results minus the gap
    // between Ms in robust
    Mgap = sort(M, 1)
    Mgap = min(Mgap[2::length(M)] :- Mgap[1::(length(M)-1)])
    Mmin = min(M)

    originalResults.M = Mmin - Mgap
    M  = (originalResults.M  \ M)  * rescaleFactor
    lb = (originalResults.lb \ lb) * rescaleFactor
    ub = (originalResults.ub \ ub) * rescaleFactor

    // Filter out observations above maxM (after rescaling)
    sel = selectindex(M :<= maxM)
    M   = M[sel]
    lb  = lb[sel]
    ub  = ub[sel]

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
    printf("(DeltaSD is hard-coded)\n")
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

// TODO: xx _honestHelper has not been thoroughly tested
real scalar function _honestHelper(real scalar s, real colvector l_vec, real scalar numPostPeriods)
{
    return(abs((1..s) * l_vec[(numPostPeriods - s + 1)::numPostPeriods]))
}
end
