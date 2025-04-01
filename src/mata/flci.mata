cap mata mata drop _flciFindOptimal()
cap mata mata drop _flciResults()
cap mata mata drop _flciFindOptimalHelper()
cap mata mata drop _flciFoldedNormalQuantiles()
cap mata mata drop _flciMatrices()
cap mata mata drop _flciFindHMinBias()
cap mata mata drop _flciSigmaLFromW()
cap mata mata drop _flciMatricesForVarianceFromW()
cap mata mata drop _flciFindLowestH()
cap mata mata drop _flciCreateConstraints()
cap mata mata drop _flciFindWorstCaseBiasGivenH()

* Done
* ----
*
* .findOptimalFLCI()                 <-> _flciFindOptimal()
* .findOptimalFLCI_helper()          <-> _flciFindOptimalHelper()
* .findHForMinimumBias()             <-> _flciFindHMinBias()
* .computeSigmaLFromW()              <-> _flciSigmaLFromW()
* .findLowestH()                     <-> _flciFindLowestH()
* .createMatricesForVarianceFromW()  <-> _flciMatricesForVarianceFromW()
* .createConstraints_AbsoluteValue() <-> _flciCreateConstraints()
* .findWorstCaseBiasGivenH()         <-> _flciFindWorstCaseBiasGivenH()
* .qfoldednormal()                   <-> _flciFoldedNormalQuantiles()
* .wToLFn()                          <-> _honestwToLFn()
* .lToWFn()                          <-> _honestlToWFn()
*
* Not translated directly
* -----------------------
*
* .createObjectiveObject_MinimizeSD()    <-> Not needed given how OSQP is called
* .createConstraints_SumWeights()        <-> Not needed given how OSQP/ECOS are called
* .createConstraintsObject_SDLessThanH() <-> Not needed given how ECOS is called
* .createObjectiveObjectForBias()        <-> Not needed given how ECOS is called
*
* Extra internal Stata structures
* -------------------------------
*
* _flciResults()
* _flciMatrices()
*
* Notes
* -----
*
* .computeSigmaLFromW() passes optional arguments to
* .createMatricesForVarianceFromW() However, this is only called from
* .findHForMinimumBias() and the extra arguments are never used. Hence I
* drop them here.
*
* .createConstraintsObject_SDLessThanH() passes optional arguments
* to .createMatricesForVarianceFromW() However, this is only called
* from .findWorstCaseBiasGivenH() and the extra arguments are never
* used. Hence I drop them here.
*
* .createObjectiveObject_MinimizeSD() passes optional arguments to
* .createMatricesForVarianceFromW() However, this is only called from
* .findLowestH() and the extra arguments are never used. Hence I drop
* them here.

mata:
struct _flciResults {
    real vector FLCI
    real vector optimalVec
    real vector optimalPrePeriodVec
    real vector weights
    real scalar optimalHalfLength
    real scalar M
    real scalar exitcode
}

// findOptimalFLCI
struct _flciResults scalar function _flciFindOptimal(
    real vector betahat,
    real matrix sigma,
    real scalar numPrePeriods,
    real scalar numPostPeriods,
    real scalar M,
    real scalar debug,
    | real colvector l_vec,
    real scalar alpha,
    real scalar numPoints)
{
    struct _flciResults scalar results

    // NB: 5 is the smallest number of points for which this works.  While
    // in theory the smallest number of points (to minimize iterations) is
    // ~7.36, in practice it's not always smallest since the interval width
    // is not the sole stopping criteria. I do 10 and call it a day.

    if ( args() < 7 ) l_vec = _honestBasis(1, numPostPeriods)
    if ( args() < 8 ) alpha = 0.05
    if ( args() < 9 ) numPoints = 10

    if ( debug == 1 ) printf("\thonest debug: _flciFindOptimal()\n")

    results = _flciFindOptimalHelper(sigma, M, numPrePeriods, numPostPeriods, l_vec, alpha, numPoints)
    results.FLCI = (results.optimalVec' * betahat' - results.optimalHalfLength,
                    results.optimalVec' * betahat' + results.optimalHalfLength)
    results.weights = runningsum(results.optimalPrePeriodVec)

    return(results)
}
end

***********************************************************************
*                                                                     *
*                       .findOptimalFLCI_helper                       *
*                                                                     *
***********************************************************************

* In some early runs there was a slight numerical precision issue
* where the optimal value (the min) would occurr just beyond hMin or
* h0 (hMax). Hence I resolved to expand the grid beyond those values
* and continue the search. However, searching too far beyond the bounds
* results in errors.  Hence I break the loop if I expand the search in
* this way two consecutive times.

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

    hGrid  = _honestLinspace(hMin, h0, numPoints)
    biasDF = J(numPoints, 2 + n + 2 * nn, .)
    for (i = 1; i <= numPoints; i++) {
        biasDF[i, .] = _flciFindWorstCaseBiasGivenH(hGrid[i], sigma, numPrePeriods, numPostPeriods, l_vec)
    }
    sel     = selectindex((biasDF[., 1] :== 0) :| (biasDF[., 1] :== 10))
    biasDF  = biasDF[sel, .]
    hGrid   = hGrid'[sel, .]
    maxBias = biasDF[., 2] :* M
    CI_halflength = _flciFoldedNormalQuantiles(1-alpha, maxBias :/ hGrid) :* hGrid
    selmin = selectindex(min(CI_halflength) :== CI_halflength)[1]
    optimal_l = biasDF[selmin, (2 + n + nn + 1)..cols(biasDF)]'
    expected = floor(log((h0 - hMin) / xtol)/log((numPoints-1)/2))
    first = (selmin :== 1)
    step = hGrid[2] - hGrid[1]
    if ( (abs(M) < epsilon(1)) & first ) {
        diff = 0
        expected = 0
    }

    // NOTE: I implemented the switch here to hedge against the possibility of
    // the CI being outside the range in hGrid. I don't quite remember why it
    // happens. However, the problem here is that I update optimal_l even when
    // I break the loop because of this issue, which doesn't seem like a good
    // idea. So I'm moving the update to only happen if the estimation will
    // not get broken otherwise.
    //
    // Further note for these iterations, I've had problems allowing for
    // a success on code 10 (reduced feasibility) where ECOS can give bad
    // results (it can restore a run with tolerances that are too wide).

    while ( ((diff > xtol) | (++iter <= expected)) & (iter < maxiter) ) {
        hGrid  = _honestLinspace(hGrid[selmin] - step, hGrid[selmin] + step, numPoints), hGrid[selmin]
        step   = 2 * step / (numPoints-1)
        biasDF = J(numPoints+1, 2 + n + 2 * nn, .)
        for (i = 1; i <= numPoints+1; i++) {
            biasDF[i, .] = _flciFindWorstCaseBiasGivenH(hGrid[i], sigma, numPrePeriods, numPostPeriods, l_vec)
        }
        // sel     = selectindex(((biasDF[., 1] :== 0) :| (biasDF[., 1] :== 10)) :& (biasDF[., 2] :< .))
        sel     = selectindex((biasDF[., 1] :== 0) :& (biasDF[., 2] :< .))
        biasDF  = biasDF[sel, .]
        hGrid   = hGrid'[sel, .]
        maxBias = biasDF[., 2] :* M

        CI_halflength = _flciFoldedNormalQuantiles(1-alpha, maxBias :/ hGrid) :* hGrid
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

// Find exact quantiles of folded normal using bisective search
// ------------------------------------------------------------

real vector function _flciFoldedNormalQuantiles(real scalar p,
                                                real vector mu,
                                                | real scalar sd) {
    if ( args() < 3 ) sd = 1

    real scalar tol, diff, xdiff
    real vector left, right, exact, quant, sel, sell, selr, miss

    miss  = !(mu :< .)
    tol   = epsilon(1)^(3/4)
    left  = invnormal(p/2) * sd :+ mu
    right = invnormal((p + 9) / 10) * sd :+ abs(mu)
    sel   = selectindex(rowmin((reldif(left, mu), reldif(right, mu))) :< epsilon(1))
    if ( length(sel) ) {
        mu[sel] = J(length(sel), 1, .)
        miss    = !(mu :< .)
        left    = invnormal(p/2) * sd :+ mu
        right   = invnormal((p + 9) / 10) * sd :+ abs(mu)
    }

    assert(all(((normal(left  :- mu) - normal(-left  :- mu)) :<= p) :| miss))
    assert(all(((normal(right :- mu) - normal(-right :- mu)) :>= p) :| miss))

    quant = (left :+ right) / 2
    exact = normal(quant :- mu) - normal(-quant :- mu)
    diff  = max(abs(exact :- p))
    xdiff = 1
    while ( !missing(diff) & diff > tol & xdiff > tol ) {
        selr = selectindex(exact :> p)
        sell = selectindex(exact :< p)
        if ( length(selr) ) right[selr] = quant[selr]
        if ( length(sell) ) left[sell]  = quant[sell]
        xdiff = max(abs(quant :- (left :+ right) / 2))
        quant = (left :+ right) / 2
        exact = normal(quant :- mu) - normal(-quant :- mu)
        diff  = max(abs(exact :- p))
    }

    return(quant)
}
end

***********************************************************************
*                                                                     *
*                        .findHForMinimumBias                         *
*                                                                     *
***********************************************************************

mata
struct _flciMatrices {
    real matrix A_quadratic_sd
    real matrix A_linear_sd
    real matrix A_constant_sd
}

real scalar function _flciFindHMinBias(real matrix sigma,
                                       real scalar numPrePeriods,
                                       real scalar numPostPeriods,
                                       real colvector l_vec)
{
    real scalar hsquared
    real rowvector w
    w = (J(1, (numPrePeriods-1), 0), (1..numPostPeriods) * l_vec)
    hsquared = _flciSigmaLFromW(w, sigma, numPrePeriods, l_vec)
    return(sqrt(hsquared))
}

// Compute variance of affine estmiator with choice w
real scalar function _flciSigmaLFromW(real rowvector w,
                                      real matrix sigma,
                                      real scalar numPrePeriods,
                                      real colvector l_vec)
{
    struct _flciMatrices  scalar A_matrices
    real colvector UstackW
    real scalar varL

    A_matrices = _flciMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
    UstackW    = J(length(w), 1, 0) \ colshape(w, 1)
    varL = UstackW' * A_matrices.A_quadratic_sd * UstackW +
        A_matrices.A_linear_sd' * UstackW +
        A_matrices.A_constant_sd
    return(varL)
}

struct _flciMatrices scalar function _flciMatricesForVarianceFromW(real matrix sigma,
                                                                   real scalar numPrePeriods,
                                                                   real colvector l_vec,
                                                                   | real rowvector prePeriodIndices) {

    struct _flciMatrices  scalar A_matrices
    real rowvector postPeriodIndices
    real matrix SigmaPre, SigmaPrePost, WtoLPreMat
    real matrix UstackWtoLPreMat, A_quadratic_sd
    real colvector A_linear_sd
    real scalar SigmaPost, A_constant_sd, col

    if ( args() < 4 ) prePeriodIndices = 1..numPrePeriods
    postPeriodIndices = _honestInverseIndex(prePeriodIndices, cols(sigma))

    // Constructs matrix to compute the variance of the affine estimator for a choice of weights W.
    SigmaPre     = sigma[prePeriodIndices, prePeriodIndices]
    SigmaPrePost = sigma[prePeriodIndices, postPeriodIndices]
    SigmaPost    = l_vec' * sigma[postPeriodIndices, postPeriodIndices] * l_vec
    WtoLPreMat   = I(numPrePeriods)

    if (numPrePeriods == 1) {
        WtoLPreMat = 1
    }
    else {
        for (col = 1; col < numPrePeriods; col++) {
            WtoLPreMat[col + 1, col] = -1
        }
    }

    UstackWtoLPreMat = (J(numPrePeriods, numPrePeriods, 0), WtoLPreMat)
    A_quadratic_sd   = UstackWtoLPreMat' * SigmaPre * UstackWtoLPreMat
    A_linear_sd      = 2 * UstackWtoLPreMat' * SigmaPrePost * l_vec
    A_constant_sd    = SigmaPost

    A_matrices = _flciMatrices()
    A_matrices.A_quadratic_sd = A_quadratic_sd
    A_matrices.A_linear_sd    = A_linear_sd
    A_matrices.A_constant_sd  = A_constant_sd

    return(A_matrices)
}
end

***********************************************************************
*                                                                     *
*                            .findLowestH                             *
*                                                                     *
***********************************************************************

* Finds the minimum variance affine estimator.
mata
real scalar function _flciFindLowestH(real matrix sigma,
                                      real scalar numPrePeriods,
                                      real scalar numPostPeriods,
                                      real colvector l_vec) {

    struct _flciMatrices scalar A_matrices
    real matrix A, P
    struct OSQP_workspace_abridged scalar varResult

    real scalar threshold_sumweights, minimalVariance, minimalH
    real rowvector q, u, l

    threshold_sumweights = (1..numPostPeriods) * l_vec
    A_matrices = _flciMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)

    P = 2 * (A_matrices.A_quadratic_sd)
    q = A_matrices.A_linear_sd'
    A = _flciCreateConstraints(numPrePeriods)
    u = threshold_sumweights, J(1, 2 * numPrePeriods, 0)
    l = threshold_sumweights, J(1, 2 * numPrePeriods, .)
    varResult = OSQP(P, q, A, u, l)

    if ( (varResult.rc != 0) | (varResult.info_status != "solved") ) {
        errprintf("Error in optimization for hMin\n")
        _error(17390)
    }
    minimalVariance = varResult.info_obj_val + A_matrices.A_constant_sd
    minimalH = sqrt(minimalVariance)

    return(minimalH)
}
end

* _flciCreateConstraints
* ----------------------
*
* ### Inequality constraints
*
* This function creates linear constraints that help to minimize
* worst-case bias over w using the OSQP library, where we cast
* maximization of the absolute values as a linear problem.  To do this,
* we assume that the optimization is over a vector (U, W) where
*
* U_1 = |w_1|, U_2 = |w_1 + w_2|, etc.
*
*
* This function creates a matrix of linear constraints that imposes
* these equalities and can easily be passed to the solver,
*
* -U_1 + w_1 <= 0
* -U_1 - w_1 <= 0, and so on.
*
* ### Equality constraints
*
* In additionl we have equality constraints of the  form
*
*     thresh = sum_t t * i_vec[t]
*     sum_t W = thresh
*
* ### Stack
*
* The solver only accepts one constraint, so we stack these
* two constraints together. They can be expressed as
*
*    (-infty, thresh) <= [A, S] * (U, W) <= (0, thresh)
*
* with S a vector of numPrePeriods 0s and numPrePeriods 1s

mata
real matrix function _flciCreateConstraints(real scalar numPrePeriods) {

    // Dense version
    // -------------

    real matrix lowerTriMat, A_absolutevalue, A
    real vector A_sumweights

    lowerTriMat = J(numPrePeriods, numPrePeriods, 1)
    _lowertriangle(lowerTriMat)
    A_absolutevalue = ((-I(numPrePeriods),  lowerTriMat ) \
                       (-I(numPrePeriods), -lowerTriMat))
    A_sumweights = (J(1, numPrePeriods, 0), J(1, numPrePeriods, 1))
    A = A_sumweights \ A_absolutevalue

    // Sparse version
    // --------------
    //
    // struct OSQP_csc_matrix scalar A
    // real vector helper
    // real scalar i
    //
    // A.data    = J(1, 2 * numPrePeriods, -1)
    // A.indptr  = 0, (1..numPrePeriods) * 2
    // A.indices = rowshape(((1::numPrePeriods), (1::numPrePeriods) :+ numPrePeriods), 1)
    // for (i = 1; i <= numPrePeriods; i++) {
    //     helper    = 0, i..(numPrePeriods), (i+numPrePeriods)..(numPrePeriods+numPrePeriods)
    //     A.indices = A.indices, helper
    //     A.indptr  = A.indptr, A.indptr[i + numPrePeriods] + length(helper)
    //     A.data    = A.data, 1, J(1, numPrePeriods-i+1, 1), J(1, numPrePeriods-i+1, -1)
    // }

    return(A)
}
end

***********************************************************************
*                                                                     *
*                      .findWorstCaseBiasGivenH                       *
*                                                                     *
***********************************************************************

* A second-order cone program (SOCP) is a convex optimization problem of
* the form
*
*    min  c' x
*    s.t. Ax = b
*         Gx <(K) h
*
* where c is a vector conformable with x, A, G are matrices conformable
* with x, and h is a vector of the same dimension as Gx. In this case,
* <(K) indicates that h - Gx belongs to the cone K, and any number of
* first and second order cones may be specified.
*
* Concretely, say G' = [G1' G2'] and h' = [h1' h2'] with h1-G1x and
* h2-G2x belonging to the first and second-order cones, resp. The
* first-order cone is the positive ortahnt R^n_+, so the constraint is
*
*    G1x <= h1
*
* Now a vector (u_1, u_2) in R x R^{n - 1} belongs to the second-order
* cone if
*
*    ||u_2||_2 <= u_1
*
* with ||.||_2 the L_2 or Eucledean norm. So the constraint is
*
*    L_2(h22 - G22x) <= h21 - G21x
*
* For G2' = [G21' G22'] and h2' = [h21 h22'], noting G21 is a row
* vector, G22 is a matrix, h21 is a scalar, and h22 is a column vector.
* Now we're in business: For our worst-case bias, we are specifically
* interested in the problem with form
*
*    min  c' x
*    s.t. Ax = b
*         Dx <= 0
*         x' W x + v' x + u <= 0
*
* With some modifications we can fit this into the SOCP! The objective
* and equality constraints are unmodified. G1 = D and h1 = 0. So far so
* good. For the quadratic constraint, we can show
*
*    x' W x + v' x + u <= 0
*
* iff
*
*    L_2(W^{1/2} x + 0.5 W^{-1/2} v) <= (0.25 v' W^{-1} v - u)^0.5
*
* Set
*
*     G21 = 0
*     G22 = -W^{1/2}
*     h21 = (0.25 v' W^{-1} v - u)^0.5
*     h22 = 0.5 W^{-1/2} v
*
* and we're done.

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

    // This function minimizes worst-case bias over Delta^{SD}(M)
    // subject to the constraint that the SD of the estimator is <= hh.
    // Note: this function assumes M = 1 unless specified otherwise.
    //
    // NB: It's a bit unfortunate that here we use h and the same letter
    // is used in the SOCP notation. To aboid confusion, I name our h "hh"

    threshold_sumweights = (1..numPostPeriods) * l_vec
    constant = - threshold_sumweights
    for (s = 1; s <= numPostPeriods; s++) {
        constant = constant + _honestHelper(s, l_vec, numPostPeriods)
    }

    // First the equality and linear inequality constraints
    // n  = numPrePeriods + numPostPeriods
    n  = numPrePeriods + numPrePeriods
    c  = (J(1, numPrePeriods, 1), J(1, numPrePeriods, 0))
    B  = _flciCreateConstraints(numPrePeriods)
    A  = B[1, .]
    b  = threshold_sumweights
    G1 = B[2::rows(B), .]
    h1 = J(n, 1, 0)

    // Now the second-order constraints
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

    // manually check the eigendecomposition succeeded
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

    // Run ECOS
    //     min
    //         c * x'
    //     s.t.
    //         A * x' - threshold_sumweights
    //         G1 * x'
    //         x * W * x' + x * v + u
    biasResult = ECOS(c, G, h, n, n + 1, 0, A, b)

    // Multiply objective by M (note that solution otherwise doesn't
    // depend on M, so no need to run many times with many different Ms)
    biasResult.info_obj_val = (biasResult.info_obj_val + constant) * M

    // NB: ECOS returns missing if solver fails; I think that should be
    // equivalent to filtering if < Ifn.
    //
    // work.info_exitcode == 0 is Optimal solution found, best
    // work.info_exitcode == 10 is Optimal solution found subject to reduced tolerances, OK

    // Compute the implied w and l
    optimal_x = biasResult.solution_x'
    optimal_w = optimal_x[(length(optimal_x)/2 + 1)::length(optimal_x)]
    optimal_l = _honestwToLFn(optimal_w)

    return(biasResult.info_exitcode, biasResult.info_obj_val, optimal_x', optimal_w', optimal_l')
}
end
