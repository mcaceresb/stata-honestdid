mata:
void function findOptimalFLCI()
{
    FLCI_Results = _findOptimalFLCI_helper()
    FLCI_Results.FLCI = (FLCI_Results.optimalVec' * betahat - FLCI_Results.optimalHalfLength,
                         FLCI_Results.optimalVec' * betahat + FLCI_Results.optimalHalfLength)

    return(FLCI_Results)
}

void function _findOptimalFLCI_helper()
{
    h0 = _findHForMinimumBias(sigma,
                              numPrePeriods,
                              numPostPeriods,
                              l_vec)

    hMin = _findLowestH(sigma,
                        numPrePeriods,
                        numPostPeriods,
                        l_vec)

    numPoints = 100
    hGrid = _honest_linspace(hMin, h0, numPoints)
    biasDF = J(1, numPoints, .)
    for (i = 1; i <= numPoints; i++) {
        biasDF[i] = _findWorstCaseBiasGivenH(hGrid[i], sigma, numPrePeriods, numPostPeriods, l_vec)
    }

    // Ok, return missing if solver fails. This takes care of the < Inf
    // and status == "optimal" etc cases, I think.
    //
    // Then just select the min
    maxBias = biasDF[selectindex(biasDF :< .)] :* Mvec[m]
    CI_halflength = _qfoldednormal(1-alpha, maxBias :/ hGrid) :* hGrid
    optimalHalfLength = min(CI_halflength)
    xx optimal_l = xx[selectindex(optimalHalfLength :== CI_halflength), .]

    results.optimalVec = (optimal_l, l_vec)
    results.optimalPrePeriodVec = (optimal_l, l_vec)
    results.optimalHalfLength = optimalHalfLength
    results.M = Mvec[m]

    return(results)
}

function _qfoldednormal(p, mu) {
    sd = 1
    numSims = 1e6
    seed = 0

    // Computes the pth quantile of the folded normal distribution with
    // mean mu and sd = sd Vectorized over mu

    // set.seed(seed)
    normDraws  = rnormal(1, numSims, 0, sd)
    pQuantiles = purrr::map_dbl(.x = mu, .f = function(mu){quantile(abs(normDraws + mu), probs = p)})
    return(pQuantiles)
}

end

* ---------------------------------------------------------------------
* Sparse matrix

// K  = 4
// Ab = 2 * A_quadratic_sd[K:, K:]
// Ab[np.tril_indices_from(Ab, -1)] = 0
// P = sparse.block_diag([sparse.csc_matrix((K, K)), sparse.triu(Ab)], format='csc')
// # P  = sparse.block_diag([sparse.csc_matrix((K, K)), Ab], format='csc')
// q  = A_lienar_sd
// A  = sparse.csc_matrix(np.vstack([A_sumweights, A_absolutevalue]))
// u  = np.array([1] + [0] * 8)
// l  = np.array([1] + [-np.inf] * 8)

* ---------------------------------------------------------------------
* _findHForMinimumBias

mata
struct _honestMatrices {
    real matrix A_quadratic_sd
    real matrix A_linear_sd
    real matrix A_constant_sd
}

void function _findHForMinimumBias(sigma, numPrePeriods, numPostPeriods, l_vec)
{
    hsquared = _computeSigmaLFromW((J(1, (numPrePeriods-1), 0), (1..numPostPeriods) * l_vec),
                                   sigma,
                                   numPrePeriods,
                                   numPostPeriods,
                                   l_vec)
    h = sqrt(hsquared)
    return(h)
}

// Compute variance of affine estmiator with choice w
void function _computeSigmaLFromW(w, sigma, numPrePeriods, numPostPeriods, l_vec) {
      A_matrices = _createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
      UstackW    = J(length(w), 1, 0) \ colshape(w, 1)
      varL = UstackW' * A_matrices.A_quadratic_sd * UstackW +
          A_matrices.A_linear_sd' * UstackW +
          A_matrices.A_constant_sd
      return(varL)
}

real scalar _honestMatrices function _createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec) {
    // if ( args() < 5 ) prePeriodIndices = 1..numPrePeriods
    prePeriodIndices  = 1..numPrePeriods
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

    A_matrices = _honestMatrices()
    A_matrices.A_quadratic_sd = A_quadratic_sd
    A_matrices.A_linear_sd    = A_linear_sd
    A_matrices.A_constant_sd  = A_constant_sd

    return(A_matrices)
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

end

* ---------------------------------------------------------------------
* _findHForMinimumBias

// function _createObjective_MinimizeSD(sigma, numPrePeriods, numPostPeriods, UstackW, l_vec) {
//     // Create objective function to minimize the standard deviation.
//     objective = CVXR::Minimize(CVXR::quad_form(UstackW, A_matrices$A_quadratic_sd)
//                                + t(A_matrices$A_linear_sd) %*% UstackW
//                                + A_matrices$A_constant_sd)
// }

mata
void function _findLowestH(sigma, numPrePeriods, numPostPeriods, l_vec) {
    struct OSQP_csc_matrix scalar A, P
    struct OSQP_workspace_abridged scalar varResult
    string scalar fname
    // Finds the minimum variance affine estimator.

    threshold_sumweights = (1..numPostPeriods) * l_vec
    A_matrices = _createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
    P = OSQP_csc_convert(2 * uppertriangle(A_matrices.A_quadratic_sd))
    q = A_matrices.A_linear_sd'
    A = OSQP_csc_convert(_createConstraints(numPrePeriods))
    u = threshold_sumweights, J(1, 2 * numPrePeriods, 0)
    l = threshold_sumweights, J(1, 2 * numPrePeriods, .)

    fname = st_tempfilename()
    OSQP_setup(fname, P, q, A, u, b)
    varResult = OSQP_solve(fname)

    if( (varResult.rc != 0) | (varResult.status != "solved") ) {
        errprintf("Error in optimization for h0\n")
    }
    minimalVariance = varResult.info_obj_val + A_matrices.A_constant_sd
    minimalH = sqrt(minimalVariance)

    return(minimalH)
}

// _createConstraints
// ------------------
//
// ### Inequality constraints
//
// This function creates linear constraints that help to minimize
// worst-case bias over w using the OSQP library, where we cast
// maximization of the absolute values as a linear problem.  To do this,
// we assume that the optimization is over a vector (U, W) where
//
// U_1 = |w_1|, U_2 = |w_1 + w_2|, etc.
//
//
// This function creates a matrix of linear constraints that imposes
// these equalities and can easily be passed to the solver,
//
// -U_1 + w_1 <= 0
// -U_1 - w_1 <= 0, and so on.
//
// ### Equality constraints
//
// In additionl we have equality constraints of the  form
//
//     thresh = sum_t t * i_vec[t]
//     sum_t W = thresh
//
// ### Stack
//
// The solver only accepts one constraint, so we stack these
// two constraints together. They can be expressed as
//
//    (-infty, thresh) <= [A, S] * (U, W) <= (0, thresh)
//
// with S a vector of numPrePeriods 0s and numPrePeriods 1s

struct OSQP_csc_matrix scalar function _createConstraints(real scalar numPrePeriods) {

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

real vector _honest_linspace(real scalar from, real scalar to, real scalar n) {
    real scalar step
    step = (to - from) / (n - 1)
    return(from :+ step :* (0..(n - 1)))
}
end

* ---------------------------------------------------------------------
* _findWorstCaseBiasGivenH

real scalar function _honest_helper(s, l_vec, numPostPeriods)
{
    return(abs((1..s) * l_vec[(numPostPeriods - s + 1)::numPostPeriods]))
}

_findWorstCaseBiasGivenH <- function(h, sigma, numPrePeriods, numPostPeriods, l_vec) {
    M = 1
    returnDF = 1

    // This function minimizes worst-case bias over Delta^{SD}(M)
    // subject to the constraint that the SD of the estimator is <= h.
    // Note: this function assumes M = 1 unless specified otherwise.

alglib?
q = (J(1, numPrePeriods, 1), J(1, numPrePeriods, 0))
min q * x
// x' * A_matrices.A_quadratic_sd * x + A_matrices.A_linear_sd' * x <= h^2 - A_matrices$A_constant_sd
// Ax <= c(threshold_sumweights, 0)
cholesky(makesymmetric(A_matrices.A_quadratic_sd))

    A_matrices = _createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
    threshold_sumweights = (1..numPostPeriods) * l_vec
    A = _createConstraints(numPrePeriods)
    constant = - threshold_sumweights
    for (s = 1; s <= numPostPeriods; s++) {
        constant = constant + _honest_helper(s, l_vec, numPostPeriods)
    }

    objectiveBias <- _createObjectiveObjectForBias(numPrePeriods  = numPrePeriods,
                                                   UstackW        = UstackW,
                                                   numPostPeriods = numPostPeriods,
                                                   l_vec          = l_vec)
    quad_constraint <- .createConstraintsObject_SDLessThanH(sigma = sigma, numPrePeriods = numPrePeriods, l_vec = l_vec, UstackW = UstackW, h = h)
    biasProblem = CVXR::Problem(objectiveBias, constraints = list(abs_constraint, sum_constraint, quad_constraint))
    biasResult <- CVXR::psolve(biasProblem)


    // Multiply objective by M (note that solution otherwise doesn't
    // depend on M, so no need to run many times with many different Ms)
    biasResult.info_obj_val = biasResult.info_obj_val * M

    // NB: xx Return missing if solver fails; I think that should be
    // equivalent to filtering if < Ifn

    return(biasResult.info_obj_val)

    // What is this for?
    // xx apparently optimal.l is used @_@
    //
    // // Compute the implied w and l
    // optimal_x <- biasResult.solution_x
    // optimal_w <- optimal_x[(length(optimal_x)/2 + 1)::length(optimal_x)]
    // optimal_l <- _wToLFn(optimal_w)
    //
    // if ( returnDF ) {
    //     temp = list(status = biasResult$status,
    //                 value   = biasResult$value,
    //                 optimal.x = I(list(unlist(optimal.x))),
    //                 optimal.w = I(list(unlist(optimal.w))),
    //                 optimal.l = I(list(unlist(optimal.l)))
    //     )
    //     return(as.data.frame(temp, stringsAsFactors = FALSE))
    // }
    // else {
    //     return(list(status = biasResult$status,
    //                 value  = biasResult$value,
    //                 optimal.x = optimal.x,
    //                 optimal.w = optimal.w,
    //                 optimal.l = optimal.l))
    // }
}

real vector function _wToLFn(w) {
    // Converts vector from w space to l space.
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

.createObjectiveObjectForBias <- function(numPrePeriods, numPostPeriods, l_vec, UstackW){
    objective.UstackW = Minimize()
    return(objective.UstackW)
}

.createConstraintsObject_SDLessThanH <- function(sigma, numPrePeriods, l_vec, UstackW, h, ...){
}
