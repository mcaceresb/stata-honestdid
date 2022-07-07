cap mata mata drop basisVector()
cap mata mata drop createSensitivityResults()
cap mata mata drop constructOriginalCS()
cap mata mata drop createSensitivityPlotMatrix()

cap mata mata drop HonestEventStudy()
cap mata mata drop HonestOptimalResults()
cap mata mata drop HonestSensitivityResults()

mata
struct HonestEventStudy {
    real vector betahat
    real matrix sigma
    real vector timeVec
    real scalar referencePeriod
    real vector prePeriodIndices
    real vector postPeriodIndices
}

struct HonestOptimalResults {
    real vector FCLI
    real vector optimalVec
    real scalar optimalHalfLength
    real scalar M
    string scalar status
}

struct HonestSensitivityResults {
    real scalar lb
    real scalar ub
    string scalar method
    string scalar Delta
    real scalar M
}
end

***********************************************************************
*                                                                     *
*                            Example Data                             *
*                                                                     *
***********************************************************************

mata
real colvector function basisVector(index, size)
{
    v = J(size, 1, 0)
    v[index] = 1
    return(v)
}

BCdata_EventStudy = HonestEventStudy()
BCdata_EventStudy.betahat = (
    0.006696352,
    0.029345034,
    -0.006472972,
    0.073014989,
    0.195961118,
    0.312063903,
    0.239541546,
    0.126042500)

BCdata_EventStudy.sigma = (
    ( 8.428358e-04,  4.768687e-04, 2.618051e-04, 0.0002354220, 0.0001676371, 0.0001128708, 1.992816e-05, -1.368265e-04) \
    ( 4.768687e-04,  6.425420e-04, 3.987425e-04, 0.0002435515, 0.0002201960, 0.0001804591, 3.843765e-05, -2.960422e-05) \
    ( 2.618051e-04,  3.987425e-04, 5.229950e-04, 0.0002117686, 0.0001840722, 0.0001458528, 7.005197e-05,  5.952995e-05) \
    ( 2.354220e-04,  2.435515e-04, 2.117686e-04, 0.0003089595, 0.0001197866, 0.0001334081, 1.016335e-04,  1.079052e-04) \
    ( 1.676371e-04,  2.201960e-04, 1.840722e-04, 0.0001197866, 0.0003599704, 0.0002478819, 1.749579e-04,  1.654257e-04) \
    ( 1.128708e-04,  1.804591e-04, 1.458528e-04, 0.0001334081, 0.0002478819, 0.0004263950, 2.171438e-04,  2.892748e-04) \
    ( 1.992816e-05,  3.843765e-05, 7.005197e-05, 0.0001016335, 0.0001749579, 0.0002171438, 4.886698e-04,  3.805322e-04) \
    (-1.368265e-04, -2.960422e-05, 5.952995e-05, 0.0001079052, 0.0001654257, 0.0002892748, 3.805322e-04,  7.617394e-04))

BCdata_EventStudy.timeVec = (2004, 2005, 2006, 2007, 2009, 2010, 2011, 2012)

BCdata_EventStudy.referencePeriod   = 2008
BCdata_EventStudy.prePeriodIndices  = (1, 2, 3, 4)
BCdata_EventStudy.postPeriodIndices = (5, 6, 7, 8)

BC_numPrePeriods  = length(BCdata_EventStudy.prePeriodIndices)
BC_numPostPeriods = length(BCdata_EventStudy.postPeriodIndices)
BC_l_vec          = basisVector(1, BC_numPostPeriods)

BC_DeltaSDNB_RobustResults = createSensitivityResults(BCdata_EventStudy.betahat,
                                                      BCdata_EventStudy.sigma,
                                                      BC_numPrePeriods,
                                                      BC_numPostPeriods,
                                                      BC_l_vec)

BC_OriginalResults = constructOriginalCS(BCdata_EventStudy.betahat,
                                         BCdata_EventStudy.sigma,
                                         BC_numPrePeriods,
                                         BC_numPostPeriods,
                                         BC_l_vec)

BC_plotMatrix = createSensitivityPlotMatrix(BC_DeltaSDNB_RobustResults, BC_OriginalResults)
st_matrix("_tmpci", BC_plotMatrix[., 2::3]')
st_local("_tmplab", invtokens(strofreal(BC_plotMatrix[., 1]')))
st_matrix("_tmpcoef", ((BC_plotMatrix[., 3] :+ BC_plotMatrix[., 2])/2)')
end

YAY!
CI slightly diff because of folded normal quantiles I think

matrix colnames _tmpci = `_tmplab'
matrix rownames _tmpci = lb95 ub95
matrix colnames _tmpcoef = `_tmplab'
coefplot matrix(_tmpcoef), ci(_tmpci) vertical cionly yline(0) ciopts(recast(rcap))
graph export xx.pdf, replace


sysuse auto, clear
matrix median = J(1,3,.)
matrix colnames median = mpg trunk turn 
matrix CI = J(2,3,.)
matrix colnames CI = mpg trunk turn 
matrix rownames CI = ll95 ul95
local i 0
foreach v of var mpg trunk turn {
    local ++ i
    quietly centile `v'
    matrix median[1, `i'] = r(c_1)
    matrix CI[1, `i'] = r(lb_1) \ r(ub_1)
}
matrix list median
matrix list CI
coefplot matrix(median), ci(CI)




mata:
struct HonestSensitivityResults colvector function createSensitivityResults(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    | real colvector l_vec)
{
    struct HonestSensitivityResults colvector Results
    struct _flciResults scalar temp

    real rowvector Mvec
    string scalar biasDirection, method, monotonicityDirection, Delta
    real scalar alpha

    if ( args() < 5 ) l_vec = basisVector(1, numPostPeriods)

    Mvec                  = (0..3) / 10
    biasDirection         = "negative"
    alpha                 = 0.05
    method                = "FLCI"
    monotonicityDirection = ""
    Delta                 = "DeltaSD"

    Results = J(length(Mvec), 1, HonestSensitivityResults())
    for (m = 1; m <= length(Mvec); m++) {
        temp = findOptimalFLCI(betahat,
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

struct HonestSensitivityResults scalar function constructOriginalCS(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    | real colvector l_vec,
    real scalar alpha)
{
    struct HonestSensitivityResults scalar OG

    if ( args() < 5 ) l_vec = basisVector(1, numPostPeriods)
    if ( args() < 6 ) alpha = 0.05

    real vector sel
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
x
    return(OG)
}

real matrix function createSensitivityPlotMatrix(
    struct HonestSensitivityResults colvector robustResults,
    struct HonestSensitivityResults scalar originalResults,
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
end

* ------------------------
* Ignore for now; switches
* ------------------------
* if ( (monotonicityDirection == "") & (biasDirection == "") ) {
*     Delta = "DeltaSD"
* }
* else if ( !(biasDirection == "") ) {
*     if (biasDirection == "positive") {
*         Delta = "DeltaSDPB"
*     } else {
*         Delta = "DeltaSDNB"
*     }
*     printf("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!\n")
* }
* else {
*     if (monotonicityDirection == "increasing") {
*         Delta = "DeltaSDI"
*     }
*     else {
*         Delta = "DeltaSDD"
*     }
*     printf("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!\n")
* }
* ------------------------
* Ignore for now; switches
* ------------------------
