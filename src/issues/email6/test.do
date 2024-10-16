/*

The issue is that ECOS sometimes tives a 10 code

    optimal subject to reduced feasibility

The R code allows this too, but the problem is that sometimes ECOS hits a
NaN and then they restore their "best" iteration, which can have an insane
tolerance. Example

    (within feastol=7.8e-316, reltol=6.6e-310, abstol=8.9e+252)

Those don't make sense, making me think there's some issue.  One idea is to
just discard results if optimal is code 10 but I'm worried this will cause
some other issues. At the moment I've modified the FLCI logic and I only
discard the reduced tolerance result if the CI is outside the pre-specified
grid, which I guess is already a dicey search so it's OK for that.

*/

* qui reghdfe Y F*event L*event if treatind==1 ,  a(i t) cluster(i)
import delimited ./beta_md.csv, clear
cap noi drop  ï
cap noi drop  v1
mata beta = st_data(., "*")
* mata st_matrix("beta", 10 * beta)
* mata st_matrix("beta", 5 * beta)
* mata st_matrix("beta", 7 * beta)
* mata st_matrix("beta", 3 * beta)
mata st_matrix("beta", beta)

import delimited ./vcov_md.csv, clear
cap noi drop  ï
cap noi drop  v1
mata vcov = st_data(., "*")
* mata st_matrix("vcov", 100 * vcov)
* mata st_matrix("vcov", 25 * vcov)
* mata st_matrix("vcov", 49 * vcov)
* mata st_matrix("vcov", 9 * vcov)
mata st_matrix("vcov", vcov)

matrix l_vec = 0.1429 \ 0.1429 \ 0.1429 \ 0.1429 \ 0.1429 \ 0.1429 \ 0.1429
honestdid, b(beta) vcov(vcov) l_vec(l_vec) pre(1/11) post(12/18) delta(sd) mvec(0(0.0008)0.004)

/*
* mata st_matrix("l_vec", J(1, 7, 1)/7)
local root .
local root ../../..
do `root'/src/ado/honestdid.ado
* qui {
*     do `root'/src/mata/osqp.mata
*     do `root'/src/mata/ecos.mata
*     do `root'/src/mata/utilities.mata
*     do `root'/src/mata/flci.mata
*     do `root'/src/mata/arp.mata
*     do `root'/src/mata/arpnn.mata
*     do `root'/src/mata/deltasd.mata
*     do `root'/src/mata/deltarm.mata
*     do `root'/src/mata/honestdid.mata
*     do `root'/src/mata/honestparallel.mata
* }
* qui do ../../..
local 0 , b(beta) vcov(vcov) l_vec(l_vec) pre(1/11) post(12/18) delta(sd) mvec(0(0.0008)0.004) parallel(0) gridPoints(100) 
* run honestdid pre
mata {
    b             = "`b'"
    V             = "`vcov'"
    numPrePeriods = `numpreperiods'
    pre           = "`preperiodindices'"
    post          = "`postperiodindices'"
    l_vec         = "`l_vec'"
    Mvec          = "`mvec'"
    alpha         = `alpha'
    method        = "`method'"
    debug         = "`debug'"
    omit          = "`omit'"
    rm            = `relativeMagnitudes'
    grid_lb       = `grid_lb'
    grid_ub       = `grid_ub'
    gridPoints    = `gridPoints'
    options = _honestOptions()
    results = HonestEventStudy()
    results = HonestDiDParse(b, V, numPrePeriods, pre, post, l_vec, Mvec, alpha,
                             method, debug, omit, rm, grid_lb, grid_ub, gridPoints)
    _honestSensitivityCIMatrix(HonestSensitivityHelper(results.betahat,
                                      results.sigma,
                                      results.numPrePeriods,
                                      results.numPostPeriods,
                                      results.options), results.OG)[3, .]

    results = HonestDiDParse(b, V, numPrePeriods, pre, post, l_vec, "0.0008", alpha,
                             method, debug, omit, rm, grid_lb, grid_ub, gridPoints)

    stata("log using test.log, replace")
    stuff = J(100, 3, .)
    for(i = 1; i <= rows(stuff); i++) {
        printf("------------------------------------------------------------------------")
        printf("%g", i)
        printf("------------------------------------------------------------------------")
        stuff[i, .] = _honestSensitivityCIMatrix(HonestSensitivityHelper(results.betahat,
                                          results.sigma,
                                          results.numPrePeriods,
                                          results.numPostPeriods,
                                          results.options), results.OG)[2, .]
    }
    stuff
    stata("log close _all")

    betahat               = results.betahat
    sigma                 = results.sigma
    numPrePeriods         = results.numPrePeriods
    numPostPeriods        = results.numPostPeriods
    alpha                 = results.options.alpha
    Mvec                  = results.options.Mvec
    l_vec                 = results.options.l_vec
    method                = results.options.method
    rm                    = results.options.rm
    debug                 = results.options.debug
    grid_lb               = results.options.grid_lb
    grid_ub               = results.options.grid_ub
    gridPoints            = results.options.gridPoints
    biasDirection         = ""
    monotonicityDirection = ""
    Results               = J(length(Mvec), 1, _honestResults())
    bound                 = "deviation from linear trend"
    bound                 = "deviation from parallel trends"
}
* run HonestSensitivityHelper pre
mata {
    m = 2
    M = Mvec[m]
    results   = _flciResults()
    numPoints = 10
    results   = _flciFindOptimalHelper(sigma, M, numPrePeriods, numPostPeriods, l_vec, alpha, numPoints)
    results.optimalVec
    results.optimalHalfLength
    (results.optimalVec' * betahat' - results.optimalHalfLength,
     results.optimalVec' * betahat' + results.optimalHalfLength)
}

mata mata drop _flciFindOptimalHelper()
mata
struct _flciResults scalar function _flciFindOptimalHelper(
    real matrix sigma,
    real scalar M,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    | real colvector l_vec,
    real scalar alpha,
    real scalar numPoints)
{
    struct _flciResults scalar results
    real rowvector hGrid
    real colvector sel, selmin, CI_halflength, optimal_l
    real scalar h0, hMin, n, nn, i, maxBias, diff, maxiter, xtol, iter, expected, step, first
    real matrix biasDF

    if ( args() < 5 ) l_vec     = _honestBasis(1, numPostPeriods)
    if ( args() < 6 ) alpha     = 0.05
    if ( args() < 7 ) numPoints = 10

    xtol    = sqrt(epsilon(1))
    maxiter = 100
    n       = numPrePeriods + numPrePeriods
    nn      = length((n/2 + 1)::n)
    h0      = _flciFindHMinBias(sigma, numPrePeriods, numPostPeriods, l_vec)
    hMin    = _flciFindLowestH(sigma, numPrePeriods, numPostPeriods, l_vec)
    diff    = 1
    iter    = 0
    printf("h0, hMin\n")
    h0, hMin

    hGrid  = _honestLinspace(hMin, h0, numPoints)
    printf("hGrid\n")
    hGrid
    biasDF = J(numPoints, 2 + n + 2 * nn, .)
    for (i = 1; i <= numPoints; i++) {
        biasDF[i, .] = _flciFindWorstCaseBiasGivenH(hGrid[i], sigma, numPrePeriods, numPostPeriods, l_vec)
    }
    sel     = selectindex((biasDF[., 1] :== 0) :| (biasDF[., 1] :== 10))
    biasDF  = biasDF[sel, .]
    hGrid   = hGrid'[sel, .]
    printf("hGrid\n")
    hGrid
    maxBias = biasDF[., 2] :* M
    printf("maxBias\n")
    maxBias
    CI_halflength = _flciFoldedNormalQuantiles(1-alpha, maxBias :/ hGrid) :* hGrid
    printf("CI_halflength\n")
    CI_halflength
    selmin = selectindex(min(CI_halflength) :== CI_halflength)[1]
    optimal_l = biasDF[selmin, (2 + n + nn + 1)..cols(biasDF)]'
    printf("optimal_l\n")
    optimal_l
    expected = floor(log((h0 - hMin) / xtol)/log((numPoints-1)/2))
    first = (selmin :== 1)
    step = hGrid[2] - hGrid[1]
    if ( (abs(M) < epsilon(1)) & first ) {
        diff = 0
        expected = 0
    }

    printf("step\n")
    step
    while ( ((diff > xtol) | (++iter <= expected)) & (iter < maxiter) ) {
        hGrid  = _honestLinspace(hGrid[selmin] - step, hGrid[selmin] + step, numPoints), hGrid[selmin]
        step   = 2 * step / (numPoints-1)
        printf("step hGrid\n")
        hGrid' \ . \ step
        biasDF = J(numPoints+1, 2 + n + 2 * nn, .)
        for (i = 1; i <= numPoints+1; i++) {
            printf("bias %g: %21.14f\n", i, hGrid[i])
            biasDF[i, .] = _flciFindWorstCaseBiasGivenH(hGrid[i], sigma, numPrePeriods, numPostPeriods, l_vec)
        }
biasDF'
_flciFindWorstCaseBiasGivenH(hGrid[1], sigma, numPrePeriods, numPostPeriods, l_vec)'
work = ECOS_workspace_abridged()
work = ECOS_solve("xx")
st_numscalar("__honestecos_obj")
work = ECOS_solve("__000002")
st_numscalar("__honestecos_obj")

mata buf1 = st_sdata(1, "xx")
mata buf2 = st_sdata(1, "__000002")

        sel     = selectindex(((biasDF[., 1] :== 0) :| (biasDF[., 1] :== 10)) :& (biasDF[., 2] :< .))
        biasDF  = biasDF[sel, .]
        hGrid   = hGrid'[sel, .]
        maxBias = biasDF[., 2] :* M

        printf("hGrid \ . \ min(hGrid) \ hMin \ maxBias\n")
        hGrid \ . \ min(hGrid) \ hMin \ maxBias
        CI_halflength = _flciFoldedNormalQuantiles(1-alpha, maxBias :/ hGrid) :* hGrid
        printf("CI_halflength\n")
        CI_halflength
        selmin = selectindex(min(CI_halflength) :== CI_halflength)[1]
        diff = max(reldif(optimal_l, biasDF[selmin, (2 + n + nn + 1)..cols(biasDF)]') \ (max(hGrid) - min(hGrid)))
        if ( first & (selmin :== 1) & (min(hGrid) < hMin) ) {
            diff = 0
            expected = 0
        }
        else {
            optimal_l = biasDF[selmin, (2 + n + nn + 1)..cols(biasDF)]'
        }
        first = (selmin :== 1)
    }

    results.optimalVec          = optimal_l \ l_vec
    results.optimalPrePeriodVec = optimal_l
    results.optimalHalfLength   = min(CI_halflength)
    results.M                   = M
    results.exitcode            = biasDF[selmin, 1]

    return(results)
}
end

mata mata drop _flciFindWorstCaseBiasGivenH()
mata
real rowvector function _flciFindWorstCaseBiasGivenH(real scalar hh,
                                                     real matrix sigma,
                                                     real scalar numPrePeriods,
                                                     real scalar numPostPeriods,
                                                     real colvector l_vec,
                                                     | real scalar M) {

    if ( args() < 6 ) M = 1

    struct _flciMatrices  scalar A_matrices
    real scalar threshold_sumweights, constant, n, b, u, s
    real vector c, A, h1, h21, h22, h, v, L
    real vector invPeriodIndices, optimal_x, optimal_l, optimal_w
    real matrix B, G1, G21, G22, G, W, C, _W0, _W1, _W2, _W3, W1, W2, W3
    struct ECOS_workspace_abridged scalar biasResult

    threshold_sumweights = (1..numPostPeriods) * l_vec
    constant = - threshold_sumweights
    for (s = 1; s <= numPostPeriods; s++) {
        constant = constant + _honestHelper(s, l_vec, numPostPeriods)
    }

    n  = numPrePeriods + numPrePeriods
    c  = (J(1, numPrePeriods, 1), J(1, numPrePeriods, 0))
    B  = _flciCreateConstraints(numPrePeriods)
    A  = B[1, .]
    b  = threshold_sumweights
    G1 = B[2::rows(B), .]
    h1 = J(n, 1, 0)

    A_matrices = _flciMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
    W  = A_matrices.A_quadratic_sd
    v  = A_matrices.A_linear_sd
    u  = A_matrices.A_constant_sd - hh^2

    invPeriodIndices = _honestInverseIndex(1..numPrePeriods, 2 * numPrePeriods)
    _W0 = ((W + W')/2)[invPeriodIndices, invPeriodIndices]
    symeigensystem(_W0, C=., L=.)
    _W1 = (C :* sqrt(L)) * C'
    _W2 = (C :/ sqrt(L)) * C'
    _W3 = (C :/     (L)) * C'
    _W1 = (_W1 + _W1')/2
    _W2 = (_W2 + _W2')/2
    _W3 = (_W3 + _W3')/2

    assert(max(reldif((C * diag(L) * C'), _W0))     < epsilon(1)^(3/4))
    assert(max(reldif(_W1 * _W1, _W0))              < epsilon(1)^(3/4))
    assert(max(reldif(_W1 * _W2, I(numPrePeriods))) < epsilon(1)^(3/4))
    assert(max(reldif(_W3 * _W0, I(numPrePeriods))) < epsilon(1)^(3/4))

    W1 = J(n, n, 0)
    W2 = J(n, n, 0)
    W3 = J(n, n, 0)
    W1[invPeriodIndices, invPeriodIndices] = _W1
    W2[invPeriodIndices, invPeriodIndices] = _W2
    W3[invPeriodIndices, invPeriodIndices] = _W3

    G21 = J(1, n, 0)
    G22 = -W1
    h21 = sqrt(0.25 * (v' * W3 * v) - u)
    h22 = 0.5 * W2 * v
    G   = G1 \ G21 \ G22
    h   = h1 \ h21 \ h22

    biasResult = ECOS(c, G, h, n, n + 1, 0, A, b, 0, 1)
    biasResult.info_obj_val = (biasResult.info_obj_val + constant) * M

    printf("ECOS output\n")
    biasResult.info_obj_val

    optimal_x = biasResult.solution_x'
    optimal_w = optimal_x[(length(optimal_x)/2 + 1)::length(optimal_x)]
    optimal_l = _honestwToLFn(optimal_w)

    return(biasResult.info_exitcode, biasResult.info_obj_val, optimal_x', optimal_w', optimal_l')
}
end

*/
