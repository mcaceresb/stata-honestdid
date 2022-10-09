cap mata mata drop _honestARPComputeCI()
cap mata mata drop _honestARPConstructGamma()
cap mata mata drop _honestReducedRowEchelon()
cap mata mata drop _honestARPLeastFavorableCV()
cap mata mata drop _honestARPDeltaLPWrapper()
cap mata mata drop _honestARPConditionalTest()
cap mata mata drop _honestARPDualLP()
cap mata mata drop _honestARPDeltaTest()
cap mata mata drop _honestARPDualVLO()
cap mata mata drop _honestARPCheckSolHelper()
cap mata mata drop _honestARPMaxProgram()
cap mata mata drop _honestARPFLCIVloVup()
cap mata mata drop _honestARPGeneralizedNorm()

* Done
* ----
*
* .ARP_computeCI()              <-> _honestARPComputeCI()
* .construct_Gamma()            <-> _honestARPConstructGamma()
* .leading_one()                <-> [deprecated]
* .compute_eta()                <-> [deprecated]
* .compute_least_favorable_cv() <-> _honestARPLeastFavorableCV()
* .test_delta_lp_fn_wrapper()   <-> _honestARPDeltaLPWrapper()
* .lp_conditional_test_fn()     <-> _honestARPConditionalTest()
* .lp_dual_fn()                 <-> _honestARPDualLP()
* .test_delta_lp_fn()           <-> _honestARPDeltaTest()
* .vlo_vup_dual_fn()            <-> _honestARPDualVLO()
* .check_if_solution_helper()   <-> _honestARPCheckSolHelper()
* .max_program()                <-> _honestARPMaxProgram()
* .FLCI_computeVloVup()         <-> _honestARPFLCIVloVup()
* .norminvp_generalized()       <-> _honestARPGeneralizedNorm()

* TODO
* ----
* xx ALL THE COMMENTS

mata
real matrix function _honestARPComputeCI(real rowvector betahat,
                                         real matrix sigma,
                                         real scalar numPrePeriods,
                                         real scalar numPostPeriods,
                                         real matrix A,
                                         real colvector d,
                                         real scalar debug,
                                         real colvector l_vec,
                                         real scalar alpha,
                                         string scalar hybrid_flag,
                                         struct _honestHybridList scalar hybrid_list,
                                         real scalar returnLength,
                                         real scalar grid_lb,
                                         real scalar grid_ub,
                                         real scalar gridPoints,
                                         real vector rowsForARP)
{
    real scalar i, theta, left, right, reject
    real matrix Gamma, AGammaInv, AGammaInv_one, AGammaInv_minusOne
    real matrix sigmaY, resultsGrid
    real vector thetaGrid, testResultsGrid, thetaDiff, gridLength, Y

    // This function computes the ARP confidence interval for delta \in Delta = {A delta <= d}.
    //
    // Inputs:
    //   betahat             = vector of estimated event study coefficients
    //   sigma               = covariance matrix of estimated event study coefficients
    //   l_vec               = vector that defines parameter of interest, l'beta.
    //   numPrePeriods       = number of pre-periods
    //   numPostPeriods      = number of post-periods
    //   A                   = matrix A that defines Delta.
    //   d                   = vector d that defines Delta
    //   alpha               = size of CI
    //   hybrid_flag         = flag for hybrid
    //   hybrid_list         = list of objects needed for hybrid
    //   returnLength        = True only returns the length of the CI, False returns the full grid of values
    //   grid.lb/ub          = upper and lower bound of the grid.
    //   gridPoints          = number of grid points to test over.
    //
    // Outputs:
    //   testResultsGrid

    if ( debug == 1 ) {
        printf("\thonest debug: _honestARPComputeCI()\n")
        printf("\thonest debug: \thybrid_flag = %s\n", hybrid_flag)
    }

    // Construct grid of theta values to test over.
    thetaGrid = _honestLinspace(grid_lb, grid_ub, gridPoints)

    // Construct matrix Gamma and A*(Gamma)^{-1}
    Gamma = _honestARPConstructGamma(l_vec)
    AGammaInv = A[., (numPrePeriods+1)::(numPrePeriods+numPostPeriods)] * qrinv(Gamma)
    AGammaInv_one = AGammaInv[., 1]
    // This will be the matrix X_T in ARP functions
    AGammaInv_minusOne = AGammaInv[., _honestInverseIndex(1, cols(AGammaInv))]

    // Compute Y = A betahat - d and its variance A*sigma*A'
    Y = A * betahat' - d
    sigmaY = A * sigma * A'

    // HYBRID: Check if hybrid is least favorable, then compute
    // least-favorable CV and add it to hybrid list
    if ( hybrid_flag == "LF" ) {
        hybrid_list.lf_cv = _honestARPLeastFavorableCV(
            sigmaY,
            hybrid_list.hybrid_kappa,
            AGammaInv_minusOne,
            rowsForARP)
    }

    // Loop over theta grid and store acceptances
    testResultsGrid = J(gridPoints, 1, .)
    for (i = 1; i <= gridPoints; i++) {
        theta = thetaGrid[i]

        // If FLCI Hybrid, compute dbar
        if (hybrid_flag == "FLCI") {
            left  = (hybrid_list.vbar' * d) + (1 - hybrid_list.vbar' * AGammaInv_one) * theta
            right = (hybrid_list.vbar' * d) - (1 - hybrid_list.vbar' * AGammaInv_one) * theta
            hybrid_list.dbar = (hybrid_list.flci_halflength - left) \ (hybrid_list.flci_halflength + right)
        }

        // For the FLCI, we compute the modified FLCI vlo, vup inside
        // this function now.
        reject = _honestARPDeltaLPWrapper(
            theta,
            Y-AGammaInv_one*theta,
            AGammaInv_minusOne,
            sigmaY,
            alpha,
            hybrid_flag,
            hybrid_list,
            rowsForARP)

        testResultsGrid[i] = 1 - reject
    }
    resultsGrid = thetaGrid', testResultsGrid

    // Compute length, else return grid
    if ( returnLength ) {
        if ( debug == 1 ) printf("\thonest debug: \treturnLength = 1\n")
        thetaDiff = thetaGrid[2::gridPoints] :- thetaGrid[1::(gridPoints-1)]
        gridLength = ((0, thetaDiff) + (thetaDiff, 0)) :/ 2
        return(gridLength * resultsGrid[., 2])
    }
    else {
        if ( debug == 1 ) printf("\thonest debug: \treturnLength = 0\n")
        return(resultsGrid)
    }
}

real matrix function _honestARPConstructGamma(real vector l) {
    // This function constructs the invertible matrix Gamma that
    // has the 1 x \bar{T} vector l' as its first row. To do so, it uses
    // the following linear algebra trick:
    //
    // Construct matrix B = [l e_1 ... e_{\bar{T}}]. Convert it to reduced
    // row echelon form and find the columns containing the pivots of
    // each row. The colums in the original matrix associated with the pivots
    // form a basis of R^{\bar{T}} and we use this to construct \Gamma.
    //
    // This function is used in the change of basis to transform Delta into
    // a MI of the form that can be used by the ARP method.

    real matrix R, B, Gamma
    real vector leading_ones
    R = B = colshape(l, 1), I(length(l))

    // Drop collinear columns
    leading_ones = _honestReducedRowEchelon(R)
    Gamma = B[., leading_ones]'
    if (rank(Gamma) < max((cols(Gamma), rows(Gamma)))) {
        errprintf("Something went wrong in QR algorithm.\n")
        _error(198)
    }
    return(Gamma)
}

real rowvector function _honestReducedRowEchelon(real matrix M)
{
    real scalar i, r, pivot, nrow, ncol
    real rowvector row, pivots

    pivots = pivot = 1
    nrow = rows(M)
    ncol = cols(M)

    if ( min((nrow, ncol)) == 0 ) {
        return
    }
    else {
        for(r = 1; r <= nrow; r++) {
            if ( pivot >= ncol ) {
                pivots = (r == 1)? pivot: (pivots, pivot)
                return(pivots)
            }
            i = r
            while ( M[i, pivot] == 0 ) {
                if ( i == nrow ) {
                    i = r
                    pivot = pivot + 1
                    if ( ncol == pivot ) {
                        pivots = (r == 1)? pivot: (pivots, pivot)
                        return(pivots)
                    }
                }
                else {
                    i = i + 1
                }
            }
            row     = M[i, .]
            M[i, .] = M[r, .]
            M[r, .] = row
            M[r, .] = M[r, .] / M[r, pivot]
            for (i = 1; i <= nrow; i++) {
                if (i != r) {
                    M[i, .] = M[i, .] :- M[r, .] :* M[i, pivot]
                }
            }
            pivots = (r == 1)? pivot: (pivots, pivot)
            pivot = pivot + 1
        }
    }
    return(pivots)
}

real scalar function _honestARPLeastFavorableCV(real matrix sigma,
                                                real scalar hybrid_kappa,
                                                | real matrix X_T,
                                                real vector rowsForARP,
                                                real scalar sims) {

    real scalar i, dimDelta
    real matrix xi_draws, Cons, C0
    real colvector eta_vec, sdVec
    real rowvector f, f0, A0
    real matrix _X_T, _sigma
    string scalar rseedcache

    // Computes the least favorable critical value following the algorithm in
    // Section 6.2 of Andrews, Roth, Pakes (2019).
    //
    // Inputs:
    //   X_T   = matrix
    //   sigma = covariance matrix of y_T
    //   hybrid_kappa = desired size
    //   sims  = number of simulations, default = 10000

    // NB: I use Halton sequences here with the intent of getting good
    // coverage for the distribution. However, for high dimensions
    // Halton sequences are highly correlated. Hence I use a modified
    // algorithm where I jumble the Halton sequences, so the coverage
    // is preserved and the correlation disappears.

    rseedcache = rseed()
    rseed(0)
    if ( args() < 5 ) sims = 1000
    if ( args() < 3 ) {
        // no nuisance parameter case; the issue is that sigma
        // may not be full-rank, so the quantile is simulated
        xi_draws = (invnormal(_honestHaltonJumble(sims, cols(sigma))) * _honestMatrixPow(sigma, 0.5)) :/ sqrt(diagonal(sigma))'
        eta_vec  = rowmax(xi_draws)
        return(_honestQuantiles(eta_vec, 1 - hybrid_kappa))
    }
    else {
        if ( args() >= 4 ) {
            _X_T   = X_T[rowsForARP, .]
            _sigma = sigma[rowsForARP, rowsForARP]
        }
        else {
            _X_T   = X_T
            _sigma = sigma
        }

        // Nuisance parameter case
        // Compute elements for LP that will be re-used
        sdVec    = sqrt(diagonal(_sigma)) // standard deviation of y_T elements
        dimDelta = cols(_X_T)             // dimension of delta
        f        = 1, J(1, dimDelta, 0)   // Define objective function
        Cons     = (sdVec, _X_T)          // Define linear constraint
        xi_draws = invnormal(_honestHaltonJumble(sims, cols(_sigma))) * _honestMatrixPow(_sigma, 0.5)

        // This is doing
        //
        //     min_x
        //         f * x
        //     s.t.
        //         sigma^{1/2} * Z <= Cons * x
        //
        // with Z a multivariate standard normal. The issue is, again, that
        // sigma needn't be full-rank, so we do a simulation.
        f0 = f, 0
        A0 = (J(1, length(f), 0), 1)
        C0 = -Cons, J(rows(Cons), 1, 0)

        eta_vec = J(sims, 1, .)
        for (i = 1; i <= sims; i++) {
            eta_vec[i] = ECOS_obj(f0, C0, -xi_draws[i, .]', rows(C0), 0, 0, A0, 1)
        }

        // We compute the 1-kappa quantile of eta_vec and return this value
        return(_honestQuantiles(eta_vec, 1-hybrid_kappa))
    }
    rseed(rseedcache)
}

// NB: theta passed but not used
real scalar function _honestARPDeltaLPWrapper(real scalar theta,
                                              real colvector y_T,
                                              real matrix X_T,
                                              real matrix sigma,
                                              real scalar alpha,
                                              string scalar hybrid_flag,
                                              struct _honestHybridList scalar hybrid_list,
                                              | real vector rowsForARP) {

    real scalar results
    if ( args() < 8 ) rowsForARP = 1::length(y_T)

    // This function providers a wrapper for the function test_delta_lp_fn.
    // it only returns an indicator for whether the test rejects.  This is
    // used by other functions to construct the ARP CIs.

    results = _honestARPConditionalTest(theta,
                                        y_T,
                                        X_T,
                                        sigma,
                                        alpha,
                                        hybrid_flag,
                                        hybrid_list,
                                        rowsForARP)
    return(results)
}

// NB: In R this also gives back eta_star, delta_star, lambda but here I
// only return reject so the wrapper is deprecated.
real scalar function _honestARPConditionalTest(real scalar theta,
                                               real colvector y_T,
                                               real matrix X_T,
                                               real matrix sigma,
                                               real scalar alpha,
                                               string scalar hybrid_flag,
                                               struct _honestHybridList scalar hybrid_list,
                                               | real vector rowsForARP,
                                               real scalar tol_lambda) {

    real scalar M, k, eta, mod_size, degenerate_flag, fullRank_flag, sigma_gamma
    real scalar vlo, vup, sigma_B_dual, maxstat, zlo_dual, zup_dual, cval, zlo, zup
    real scalar size_B, sigma2_B, sigma_B

    real colvector y_T_ARP, gamma_tilde, gamma_full, S
    real rowvector B_index, Bc_index, lpDualSoln, vFLCI
    real colvector sdVec, sdVec_B, sdVec_Bc, e1, v_B, rho, maximand_or_minimand

    real matrix X_T_ARP, sigma_ARP, VbarMat, X_TB
    real matrix X_TBc, S_B, S_Bc, Gamma_B

    struct ECOS_workspace_abridged scalar linSoln

    if ( args() < 8 ) rowsForARP = 1::length(y_T)
    if ( args() < 9 ) tol_lambda = 1e-6

    // Performs ARP test of moment inequality E[y_T - X_T \delta] <= 0.
    // It returns an indicator for whether the test rejects, as
    // as the maximum statistic and delta.
    //
    // Inputs:
    //   theta        = value of theta to test
    //   y_T          = vector
    //   X_T          = matrix
    //   sigma        = covariance matrix of y_T
    //   alpha        = desired size of test
    //   hybrid_flag  = Flag for whether implement hybrid test.
    //                   "LF"   = least-favorable hybrid,
    //                   "FLCI" = FLCI hybrid
    //                   "ARP"  = ARP
    //   hybrid_list   = list object with following elements:
    //                     hybrid_kappa    = size of first-stage test
    //                     lf_cv           = least-favorable cv (only if LF hybrid)
    //                     flci_vlo        = value of vlo for FLCI (only if FLCI hybrid)
    //                     flci_vup        = value of vup for FLCI (only if FLCI hybrid)
    //                     flci_pe         = FLCI point estimate (only if FLCI hybrid)
    //                     flci_halflength = half-length of FLCI (only if FLCI hybrid)
    // rowsForARP = an index of which rows of Y_T and X_T should be used in the ARP problem.
    // This allows for ARP to be based on a subset of moments used for the FLCI hybrid (e.g. excluding postperiods)
    //
    // Outputs: list with elements
    //   reject = indicator for wheter test rejects
    //   eta = eta value
    //   delta = delta value
    //   lambda = value of lagrange multipliers

    y_T_ARP   = y_T[rowsForARP]
    X_T_ARP   = X_T[rowsForARP, .]
    sigma_ARP = sigma[rowsForARP, rowsForARP]

// Debug: Use this in place of *_ARP to switch to primal
//     y_T = 1::3
//     X_T = (1 \ -2 \ 1), (-1 \ 1 \ 0)
//     sigma = I(3)
//
// Debug: Use this in place of *_ARP with +/- y_T to test out the VLO/VUP function
//     rseed(1234)
//     y_T = -(1::5)
//     X_T = (1 \ -2 \ 1 \ 1 \ 0), (-1 \ 1 \ 0 \ 1 \ 1), (0 \ 0 \ 0 \ 0 \ 1)
//     sigma = makesymmetric(runiform(5, 5))

    // Dimensions of objects
    M = rows(sigma_ARP)
    k = cols(X_T_ARP)

    // Compute eta, and the argmin delta
    linSoln = _honestARPDeltaTest(y_T_ARP, X_T_ARP, sigma_ARP)
    eta     = linSoln.info_obj_val

    // Check to see if valid solution to LP obtained
    if ( linSoln.success == 0 ) {
        errprintf("_honestARPConditionalTest(): ")
        errprintf("LP for eta did not converge properly. Not rejecting\n")
        return(0)
    }

    // HYBRID: Implement first-stage test for hybrid.
    if (hybrid_flag == "LF") {
        // Least-favorable hybrid: check if estimated eta is larger than
        // least favorable critical value. If so, reject immediately.
        mod_size = (alpha - hybrid_list.hybrid_kappa) / (1 - hybrid_list.hybrid_kappa)
        if (eta > hybrid_list.lf_cv) {
            return(1)
        }
    }
    else if (hybrid_flag == "FLCI") {
        // FLCI hybrid: Test whether parameter of interest falls within
        // FLCI. If not, reject
        mod_size = (alpha - hybrid_list.hybrid_kappa)/(1 - hybrid_list.hybrid_kappa)
        VbarMat  = (hybrid_list.vbar, -hybrid_list.vbar)'
        if ( max(VbarMat * y_T :- hybrid_list.dbar) > 0 ) {
            return(1)
        }
    }
    else if (hybrid_flag == "ARP") {
        // No hybrid selected
        mod_size = alpha
    }
    else {
        // if hybrid flag does not equal LF, FLCI or ARP, return error
        errprintf("_honestARPConditionalTest(): ")
        errprintf("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'\n")
        _error(198)
    }

    // We now check for conditions under which the primal and dual
    // solutions are equal to one another. If these conditions don't hold,
    // we will switch to the dual.
    degenerate_flag = (sum(linSoln.solution_z :> tol_lambda) != (k+1))

    // Store which moments are binding: Do so using lambda rather than the moments.
    B_index  = (linSoln.solution_z :> tol_lambda)
    Bc_index = (1 :- B_index)
    X_TB     = X_T_ARP[selectindex(B_index), .]  // select binding moments
    // Check whether binding moments have full rank
    fullRank_flag = (rank(X_TB) == min((rows(X_TB), cols(X_TB))))

    // If degenerate or binding moments don't have full rank, switch to dual
    if ( !fullRank_flag | degenerate_flag ) {
        // -------------
        // Dual approach
        // -------------
        // warning('Using dual approach')

        // We calculate vlo and vup using the bisection approach that
        // conditions on having a gamma_tilde - a vertex of the dual. This
        // is returned as lambda, from the functioon test_delta_lp_fn.

        gamma_tilde = edittozerotol(linSoln.solution_z', tol_lambda) // numerical precision foo
        lpDualSoln  = _honestARPDualLP(y_T_ARP, X_T_ARP, eta, gamma_tilde, sigma_ARP)

        vlo = lpDualSoln[1]
        vup = lpDualSoln[2]
        sigma_B_dual = sqrt(gamma_tilde' * sigma_ARP * gamma_tilde)

        // If sigma_B_dual is 0 to numerical precision, reject iff eta > 0
        if ( sigma_B_dual < 1e-10 ) {
            return(eta > 0)
        }
        maxstat = eta/sigma_B_dual

        // HYBRID: Modify vlo, vup for the hybrid test
        if (hybrid_flag == "LF") {
            // Modify only vup using least-favorable CV for the
            // least-favorable hybrid
            zlo_dual = vlo/sigma_B_dual
            zup_dual = min((vup, hybrid_list.lf_cv))/sigma_B_dual
        }
        else if (hybrid_flag == "FLCI") {
            // Compute vlo_FLCI, vup_FLCI.
            gamma_full = J(length(y_T), 1, 0)
            gamma_full[rowsForARP] = gamma_tilde
            sigma_gamma = (sigma * gamma_full) * invsym(gamma_full' * sigma * gamma_full)

            S = y_T - sigma_gamma * (gamma_full' * y_T)
            vFLCI = _honestARPFLCIVloVup(hybrid_list.vbar, hybrid_list.dbar, S, sigma_gamma)

            // Modify vlo, vup using the FLCI vlo, vup values
            zlo_dual = max((vlo, vFLCI[1]))/sigma_B_dual
            zup_dual = min((vup, vFLCI[2]))/sigma_B_dual
        }
        else if (hybrid_flag == "ARP") {
            // If ARP, don't modify
            zlo_dual = vlo/sigma_B_dual
            zup_dual = vup/sigma_B_dual
        }
        else {
            // if hybrid flag does not equal LF, FLCI or ARP, return error
            errprintf("_honestARPConditionalTest(): ")
            errprintf("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'\n")
            _error(198)
        }

        if (!((missing(zlo_dual) | zlo_dual <= maxstat) & (maxstat <= zup_dual | missing(zup_dual)))) {
            return(0)
        }
        else {
            // Per ARP (2021), CV = max(0, c_{1-alpha}), where
            // c_{1-alpha} is the 1-alpha quantile of truncated normal.
            cval = max((0, _honestARPGeneralizedNorm(1 - mod_size, zlo_dual, zup_dual)))
            return(maxstat > cval)
        }
    }
    else {

        // ---------------
        // primal approach
        // ---------------
        size_B   = sum(B_index)
        sdVec    = sqrt(diagonal(sigma_ARP))
        sdVec_B  = sdVec[selectindex(B_index)]
        sdVec_Bc = sdVec[selectindex(Bc_index)]

        X_TBc = X_T_ARP[selectindex(Bc_index), .]
        S_B   = I(M)
        S_B   = S_B[selectindex(B_index), .]
        S_Bc  = I(M)
        S_Bc  = S_Bc[selectindex(Bc_index), .]

        Gamma_B  = (sdVec_Bc, X_TBc) * qrinv((sdVec_B, X_TB)) * S_B - S_Bc
        e1       = (1 \ J(size_B-1, 1, 0))
        v_B      = (e1' * qrinv((sdVec_B, X_TB)) * S_B)'
        sigma2_B = v_B' * sigma_ARP * v_B
        sigma_B  = sqrt(sigma2_B)
        rho      = (Gamma_B * sigma_ARP * v_B) :/ sigma2_B

        // Compute truncation values
        maximand_or_minimand = ((-Gamma_B * y_T_ARP) :/ rho) :+ (v_B' * y_T_ARP)

        vlo = any(rho :> 0)? max(maximand_or_minimand[selectindex(rho :> 0)]): .
        vup = any(rho :< 0)? min(maximand_or_minimand[selectindex(rho :< 0)]): .

        // HYBRID: Modify vlo, vup for the hybrid
        if (hybrid_flag == "LF") {
            // if LF, modify vup.
            zlo = vlo/sigma_B
            zup = min((vup, hybrid_list.lf_cv))/sigma_B
        }
        else if (hybrid_flag == "FLCI") {
            // Compute vlo_FLCI, vup_FLCI.
            gamma_full = J(length(y_T), 1, 0)
            gamma_full[rowsForARP] = v_B

            sigma_gamma = (sigma * gamma_full) * invsym((gamma_full' * sigma * gamma_full))
            S = y_T :- sigma_gamma * (gamma_full' * y_T)
            vFLCI = _honestARPFLCIVloVup(hybrid_list.vbar, hybrid_list.dbar, S, sigma_gamma)

            // if FLCI, modify both vlo, vup
            zlo = max((vlo, vFLCI[1]))/sigma_B
            zup = min((vup, vFLCI[2]))/sigma_B
        }
        else if (hybrid_flag == "ARP") {
            // If not hybrid, don't modify
            zlo = vlo/sigma_B
            zup = vup/sigma_B
        }
        else {
            // if hybrid flag does not equal LF, FLCI or ARP, return error
            errprintf("_honestARPConditionalTest(): ")
            errprintf("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'\n")
            _error(198)
        }

        // Compute test statistic
        maxstat = eta/sigma_B
        if ( !((missing(zlo) | zlo <= maxstat) & (maxstat <= zup | missing(zup))) ) {
            // errprintf("max stat (%3f) is not between z_lo (%3f) and z_up (%3f) in the primal approach",
            //           maxstat, zlo, zup)
            return(0)
        }
        else {
            // Per ARP (2021), CV = max(0, c_{1-alpha}), where c_{1-alpha}
            // is the 1-alpha quantile of truncated normal.
            cval = max((0, _honestARPGeneralizedNorm(1 - mod_size, zlo, zup)))
            return(maxstat > cval)
        }
    }
}

struct ECOS_workspace_abridged scalar function _honestARPDeltaTest(
    real vector y_T,
    real matrix X_T,
    real matrix sigma) {

    struct ECOS_workspace_abridged scalar res
    real matrix C
    real rowvector f
    real colvector b, sdVec
    real scalar dimDelta

    // Returns the value of eta that solves
    //   min_{eta, delta} eta
    //     s.t. y_T - x_T delta <= eta * diag(sigma)
    //
    // Inputs:
    //   y_T = vector
    //   X_T = matrix
    //   sigma = covariance matrix of y_T
    // Outputs: list with elements
    //   etaStar   = minimum of linear program
    //   deltaStar = minimizer of linear program
    //   lambda    = lagrange multipliers on primal i.e. solution to dual
    //   error_flag = whether the linear program was solved.
    //
    // Note:
    //   To solve the lin. program, we need to transform the
    //   linear program such that it is in the form
    //     max f' X
    //     s.t. C X \leq b

    dimDelta = cols(X_T)             // dimension of delta
    sdVec    = sqrt(diagonal(sigma)) // standard deviation of y_T elements

    // Define objective function
    f = 1, J(1, dimDelta+1, 0)

    // Define linear constraint
    C = -(sdVec, X_T, J(rows(X_T), 1, 0))
    b = -y_T

    // Define linear program using lpSolveAPI
    res = ECOS(f, C, b, rows(C), 0, 0, (J(1, length(f)-1, 0), 1), 1)
    res.solution_x = res.solution_x[1::(length(f)-1)]

    return(res)
}

real rowvector function _honestARPDualLP(real colvector y_T,
                                         real matrix X_T,
                                         real scalar eta,
                                         real colvector gamma_tilde,
                                         real matrix sigma) {

    real colvector sdVec
    real matrix W_T, s_T
    real scalar inv

    // Wrapper function to calculate vlo, vlup using the bisection approach
    //
    // Inputs:
    //   y_T = vector
    //   X_T = matrix
    //   eta = solution to LP from test_delta_lp_fn
    //   gamma_tilde = vertex of the dual, the output lambda from test_delta_lp_fn
    //   sigma = covariance matrix of y_T.

    sdVec = sqrt(diagonal(sigma))
    W_T   = (sdVec, X_T)
    inv   = edittozerotol(gamma_tilde' * sigma * gamma_tilde, epsilon(1)) // numerical precision
    s_T   = (I(length(y_T)) - (1/inv) * (sigma * (gamma_tilde * gamma_tilde'))) * y_T
    return(_honestARPDualVLO(eta, s_T, gamma_tilde, sigma, W_T))
}

real rowvector function _honestARPDualVLO(real scalar eta,
                                          real colvector s_T,
                                          real colvector gamma_tilde,
                                          real matrix sigma,
                                          real matrix W_T) {

    struct ECOS_workspace_abridged scalar linprog
    real scalar tol_c, tol_eq, sigma_B, low_initial, high_initial, maxiters, switchiters
    real scalar iters, vlo, vup, checksol, dif, low, high, mid
    real colvector b

    // This function computes vlo and vup for the dual linear program for the
    // conditional values using the bisection method described in Appendix
    // D of ARP (2021) in Algorithm 1.
    //
    // Inputs:
    //   eta         = solution to LP from test_delta_lp_fn
    //   s_T         =
    //   gamma_tilde = dual solution to LP from test_delta_lp_fn, given in output lambda.
    //   sigma       = covariance matrix of y_T
    //   W_T         =
    //
    // Output: list with elements
    //   vlo
    //   vup

    // Options for bisection algorithm
    tol_c        = 1e-6
    tol_eq       = 1e-6
    sigma_B      = sqrt(gamma_tilde' * sigma * gamma_tilde)
    low_initial  = min((-100, eta - 20 * sigma_B))
    high_initial = max(( 100, eta + 20 * sigma_B))
    maxiters     = 10000
    switchiters  = 10
    checksol     = _honestARPCheckSolHelper(eta, tol_eq, s_T, gamma_tilde, sigma, W_T)
    if ( missing(checksol) | (checksol == 0) ) {
        // errprintf("User-supplied eta is not a solution. Not rejecting automatically\n")
        return((eta, .))
    }

    // Compute vup
    linprog = ECOS_workspace_abridged()
    if ( _honestARPCheckSolHelper(high_initial, tol_eq, s_T, gamma_tilde, sigma, W_T, linprog) ) {
        vup = .
    }
    else {
        // Shortcut instead of bisection search
        dif   = 0
        iters = 0
        b     = (1/(gamma_tilde' * sigma * gamma_tilde)) * (sigma * gamma_tilde)
        mid   = (linprog.solution_x * s_T) / (1 - linprog.solution_x * b)
        while ( (++iters < maxiters) & !_honestARPCheckSolHelper(mid, tol_eq, s_T, gamma_tilde, sigma, W_T, linprog) ) {
            if ( iters >= switchiters ) {
                dif = tol_c + 1
                break
            }
            mid = (linprog.solution_x * s_T) / (1 - linprog.solution_x * b)
        }

        // Switch back to bijection if shortcut taking too long
        low   = eta
        high  = mid
        while ( (dif > tol_c) & (++iters < maxiters) ) {
            mid = (high + low)/2
            if ( _honestARPCheckSolHelper(mid, tol_eq, s_T, gamma_tilde, sigma, W_T, linprog) ) {
                low = mid
            }
            else {
                high = mid
            }
            dif = high - low
        }
        vup = mid
    }

    // Compute vlo
    if ( _honestARPCheckSolHelper(low_initial, tol_eq, s_T, gamma_tilde, sigma, W_T, linprog) ) {
        // NB: -Inf, Inf are not actual concepts in Stata (missing is
        // technically equivalent to Inf but it's not always the best
        // idea to use it in that way. In any case, if either bound is
        // missing the condition evaluates to True downstream.
        vlo = .
    }
    else {
        dif   = 0
        iters = 0
        b     = (1/(gamma_tilde' * sigma * gamma_tilde)) * (sigma * gamma_tilde)
        mid   = (linprog.solution_x * s_T) / (1 - linprog.solution_x * b)
        while ( (++iters < maxiters) & !_honestARPCheckSolHelper(mid, tol_eq, s_T, gamma_tilde, sigma, W_T, linprog) ) {
            if ( iters >= switchiters ) {
                dif = tol_c + 1
                break
            }
            mid = (linprog.solution_x * s_T) / (1 - linprog.solution_x * b)
        }

        // Switch back to bijection if shortcut taking too long
        low   = mid
        high  = eta
        while ( (dif > tol_c) & (++iters < maxiters) ) {
            mid = (low + high)/2
            if ( _honestARPCheckSolHelper(mid, tol_eq, s_T, gamma_tilde, sigma, W_T, linprog) ) {
                high = mid
            }
            else {
                low = mid
            }
            dif = high - low
        }
        vlo = mid
    }

    return((vlo, vup))
}

real scalar function _honestARPCheckSolHelper(real scalar c,
                                              real scalar tol,
                                              real colvector s_T,
                                              real colvector gamma_tilde,
                                              real matrix sigma,
                                              real matrix W_T,
                                              | struct ECOS_workspace_abridged linprog) {

    if ( args() < 7 ) linprog = ECOS_workspace_abridged()
    real scalar min_value
    min_value = _honestARPMaxProgram(s_T, gamma_tilde, sigma, W_T, c, linprog)
    return((missing(min_value)? .: abs(c - min_value) <= tol))
}

real scalar function _honestARPMaxProgram(real colvector s_T,
                                          real colvector gamma_tilde,
                                          real matrix sigma,
                                          real matrix W_T,
                                          real scalar c,
                                          | struct ECOS_workspace_abridged linprog) {

    if ( args() < 6 ) linprog = ECOS_workspace_abridged()

    real scalar n
    real colvector f, beq
    real matrix Aeq

    // Define objective and constraints
    f   = s_T :+ (1/(gamma_tilde' * sigma * gamma_tilde)) * (sigma * gamma_tilde) * c
    n   = length(f)
    Aeq = W_T'
    beq = 1 \ J(rows(Aeq)-1, 1, 0)

    // Set up linear program. Provide a lower bound of 0, which apparently
    // is the default lower bound in the ROI package.  Solve linear program
    // and return negative of objective because we want the max.
    linprog = ECOS(-f, -I(n), J(n, 1, 0), n, 0, 0, Aeq, beq)

    return(-linprog.info_obj_val)
}

real rowvector function _honestARPFLCIVloVup(real colvector vbar,
                                             real colvector dbar,
                                             real colvector S,
                                             real colvector c) {
    real matrix VbarMat
    real colvector max_or_min
    real scalar vlo, vup
    // This function computes the values of Vlo, Vup modified for the FLCI hybrid.
    // Outputs:
    //   flci_vlo = value of vlo associated with FLCI
    //   flci_vup = value of vup associated with FLCI
    // Compute vbarMat
    VbarMat = (vbar, -vbar)'
    // Comptute max_min
    max_or_min = (dbar :- (VbarMat * S)) :/ (VbarMat * c)
    // Compute Vlo, Vup for the FLCI
    // NB: If empty this returns missing for either bound
    vlo = max(max_or_min[selectindex((VbarMat * c) :< 0)])
    vup = min(max_or_min[selectindex((VbarMat * c) :> 0)])
    // Return Vlo, Vup
    return((vlo, vup))
}

// TODO: xx Precise way to compute inverse truncated normal (current
// method falls back on numerical integration approx).
real scalar function _honestARPGeneralizedNorm(real scalar p,
                                               real scalar l,
                                               real scalar u,
                                               | real scalar mu,
                                               real scalar sd) {

    real scalar left, right, quant
    if ( args() < 4 ) mu = 0
    if ( args() < 5 ) sd = 1
    right = missing(u)? 1: normal((u-mu)/sd)
    left  = missing(l)? 0: normal((l-mu)/sd)
    quant = mu + sd * invnormal(p * (right - left) + left)
    return(missing(quant)? _honestTruncNormInv(p, l, u, mu, sd): quant)
}
end
