cap mata mata drop _honestSDConditionalCS()
cap mata mata drop _honestSDUpperBoundMpre()
cap mata mata drop _honestSDComputeIDSet()
cap mata mata drop _honestSDCreateA()
cap mata mata drop _honestSDCreated()

* Done
* ----
*
* DeltaSD_upperBound_Mpre()      <-> _honestSDUpperBoundMpre()
* .create_A_SD()                 <-> _honestSDCreateA()
* .create_d_SD()                 <-> _honestSDCreated()
* .compute_IDset_DeltaSD()       <-> _honestSDComputeIDSet()
* computeConditionalCS_DeltaSD() <-> _honestSDConditionalCS()
* .ARP_computeCI()               <-> _honestARPComputeCI()
* .APR_computeCI_NoNuis()        <-> _honestARPComputeCINoNuis()

* TODO
* ----
* xx ALL THE COMMENTS

mata
real matrix function _honestSDConditionalCS(real rowvector betahat,
                                            real matrix sigma,
                                            real scalar numPrePeriods,
                                            real scalar numPostPeriods,
                                            real scalar debug,
                                            | real colvector l_vec,
                                            real scalar M,
                                            real scalar alpha,
                                            string scalar hybrid_flag,
                                            real scalar grid_lb,
                                            real scalar grid_ub,
                                            real scalar gridPoints)
{

    struct _honestHybridList scalar hybrid_list
    real rowvector IDset
    real colvector sel
    struct _flciResults scalar flci
    struct OSQP_workspace_abridged scalar result

    // NB: If inputs are assigned the input object is changed; 
    // hence grid_lb and grid_ub need to be local variables.
    real scalar hybrid_kappa, returnLength, postPeriodMomentsOnly, sdTheta, gridlb, gridub
    real vector d_SD, rowsForARP, postPeriodIndices, q
    real matrix A_SD, P, CI

    // This function computes the ARP CI that includes nuisance parameters
    // for Delta^{SD}(M). This functions uses ARP_computeCI for all
    // of its computations.
    //
    // Inputs:
    //   betahat             = vector of estimated event study coefficients
    //   sigma               = covariance matrix of estimated event study coefficients
    //   numPrePeriods       = number of pre-periods
    //   numPostPeriods      = number of post-periods
    //   l_vec               = vector that defines parameter of interest
    //   M                   = tuning parameter for Delta^SD(M), default M = 0.
    //   alpha               = desired size of CI, default alpha = 0.05
    //   hybrid_flag         = flag for hybrid, default = "FLCI"
    //   hybrid_kappa        = desired size of first-stage hybrid test, default = NULL
    //   returnLength        = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
    //   gridPoints            = number of gridpoints to test over, default = 1000
    //   postPeriodMomentsOnly = exclude moments for delta^SD that only include pre-period coefs
    //
    //  Outputs:
    //   data_frame containing upper and lower bounds of CI.

    if ( args() < 6  ) l_vec = _honestBasis(1, numPostPeriods)
    if ( args() < 7  ) M = 0
    if ( args() < 8  ) alpha = 0.05
    if ( args() < 9  ) hybrid_flag = "FLCI"
    if ( args() < 10 ) grid_lb     = .
    if ( args() < 11 ) grid_ub     = .
    if ( args() < 12 ) gridPoints  = 1e3

    // TODO xx hard-coded
    hybrid_kappa          = alpha/10
    returnLength          = 0
    postPeriodMomentsOnly = 1

    // Construct A_SD, d_SD
    A_SD = _honestSDCreateA(numPrePeriods, numPostPeriods, 0)
    d_SD = _honestSDCreated(numPrePeriods, numPostPeriods, M, 0)

    if ( postPeriodMomentsOnly & (numPostPeriods > 1) ) {
        postPeriodIndices = (numPrePeriods +1)..cols(A_SD)
        rowsForARP = selectindex(rowsum(A_SD[. , postPeriodIndices] :!= 0) :> 0)
    }
    else {
        rowsForARP = 1::rows(A_SD)
    }

    // Create hybrid_list object
    hybrid_list.hybrid_kappa = hybrid_kappa

    // If there is only one post-period, we use the no-nuisance parameter functions
    if ( debug == 1 ) printf("\thonest debug: _honestSDConditionalCS()\n")
    if (numPostPeriods == 1) {
        if ( debug == 1 ) printf("\thonest debug: \tnumPostPeriods = 1\n")

        if ( hybrid_flag == "FLCI" ) {
            if ( debug == 1 ) printf("\thonest debug: \thybrid_flag = FLCI\n")
            // Compute FLCI
            flci = _flciFindOptimalHelper(sigma, M, numPrePeriods, numPostPeriods, l_vec, hybrid_kappa)

            // Add objects to hybrid_list: flci l vector
            hybrid_list.flci_l = flci.optimalVec

            // Add objects to hybrid_list: flci half-length
            hybrid_list.flci_halflength = flci.optimalHalfLength

            // compute FLCI ub and FLCI lb
            gridlb = missing(grid_lb)? (flci.optimalVec' * betahat') - flci.optimalHalfLength: grid_lb
            gridub = missing(grid_ub)? (flci.optimalVec' * betahat') + flci.optimalHalfLength: grid_ub
        }
        else if ( hybrid_flag == "LF" ) {
            if ( debug == 1 ) printf("\thonest debug: \thybrid_flag = LF\n")
            // Compute LF CV
            hybrid_list.lf_cv = _honestARPLeastFavorableCV(A_SD * sigma * A_SD', hybrid_kappa)

            // construct theta grid
            if ( missing(grid_lb) | missing(grid_ub) ) {
                // Compute identified set under parallel trends
                sel = (numPrePeriods+1)::(numPrePeriods+numPostPeriods)
                IDset = _honestSDComputeIDSet(M, J(1, numPrePeriods + numPostPeriods, 0),
                                              l_vec, numPrePeriods, numPostPeriods)
                sdTheta = sqrt(l_vec' * sigma[sel, sel] * l_vec)
            }
            gridlb = missing(grid_lb)? IDset[1] - 20 * sdTheta: grid_lb
            gridub = missing(grid_ub)? IDset[1] + 20 * sdTheta: grid_ub
        }
        else if ( hybrid_flag == "ARP" ) {
            if ( debug == 1 ) printf("\thonest debug: \thybrid_flag = ARP\n")
            // construct theta grid
            if ( missing(grid_lb) | missing(grid_ub) ) {
                // Compute identified set under parallel trends
                sel = (numPrePeriods+1)::(numPrePeriods+numPostPeriods)
                IDset = _honestSDComputeIDSet(M, J(1, numPrePeriods + numPostPeriods, 0),
                                              l_vec, numPrePeriods, numPostPeriods)
                sdTheta = sqrt(l_vec' * sigma[sel, sel] * l_vec)
            }
            gridlb = missing(grid_lb)? IDset[1] - 20 * sdTheta: grid_lb
            gridub = missing(grid_ub)? IDset[2] + 20 * sdTheta: grid_ub
        }
        else {
            errprintf("_honestSDConditionalCS(): ")
            errprintf("hybrid_flag must equal 'APR' or 'FLCI' or 'LF'\n")
            _error(198)
        }

        CI = _honestARPComputeCINoNuis(betahat,
                                       sigma,
                                       numPrePeriods,
                                       numPostPeriods,
                                       A_SD,
                                       d_SD,
                                       debug,
                                       l_vec,
                                       alpha,
                                       hybrid_flag,
                                       hybrid_list,
                                       returnLength,
                                       gridlb,
                                       gridub,
                                       gridPoints)

    }
    else {
        // CASE: NumPostPeriods > 1
        // HYBRID: If hybrid, compute FLCI
        if (hybrid_flag == "FLCI") {
            flci = _flciFindOptimalHelper(sigma, M, numPrePeriods, numPostPeriods, l_vec, hybrid_kappa)

            // Add objects to hybrid_list: flci l vector
            hybrid_list.flci_l = flci.optimalVec

            // Add vbar to flci l vector
            P = 2 * A_SD * A_SD'
            q = 2 * A_SD * flci.optimalVec

            // NB: This isn't just invsym(P) * q because P is not full-rank.
            result = OSQP(P, -q, I(rows(P)), J(1, rows(P), .), J(1, rows(P), .))
            hybrid_list.vbar = result.solution_x'

            // Add objects to hybrid_list: flci half-length
            hybrid_list.flci_halflength = flci.optimalHalfLength

            // compute FLCI ub and FLCI lb
            gridlb = missing(grid_lb)? (flci.optimalVec' * betahat') - flci.optimalHalfLength: grid_lb
            gridub = missing(grid_ub)? (flci.optimalVec' * betahat') + flci.optimalHalfLength: grid_ub
        }
        else {
            if ( missing(grid_lb) | missing(grid_ub) ) {
                // Compute identified set under parallel trends
                sel = (numPrePeriods+1)::(numPrePeriods+numPostPeriods)
                IDset = _honestSDComputeIDSet(M, J(1, numPrePeriods + numPostPeriods, 0),
                                              l_vec, numPrePeriods, numPostPeriods)
                sdTheta = sqrt(l_vec' * sigma[sel, sel] * l_vec)
            }
            gridlb = missing(grid_lb)? IDset[1] - 20 * sdTheta: grid_lb
            gridub = missing(grid_ub)? IDset[2] + 20 * sdTheta: grid_ub
        }

        // Compute ARP CI for l'beta using Delta^SD
        CI = _honestARPComputeCI(betahat,
                                 sigma,
                                 numPrePeriods,
                                 numPostPeriods,
                                 A_SD,
                                 d_SD,
                                 debug,
                                 l_vec,
                                 alpha,
                                 hybrid_flag,
                                 hybrid_list,
                                 returnLength,
                                 gridlb,
                                 gridub,
                                 gridPoints,
                                 rowsForARP)

    }
    // Returns CI
    return(CI)
}

// This function constructs an upper-bound for M at the 1-alpha
// level based on the observed pre-period coefficients.
real scalar function _honestSDUpperBoundMpre(real vector betahat,
                                             real matrix sigma,
                                             real scalar numPrePeriods,
                                             real scalar alpha) {

    real vector prePeriod_coef, prePeriodCoefDiffs, seDiffs, upperBoundVec
    real matrix prePeriod_sigma, A_SD, prePeriodSigmaDiffs

    prePeriod_coef  = betahat[1::numPrePeriods]'
    prePeriod_sigma = sigma[1::numPrePeriods, 1::numPrePeriods]

    A_SD                = _honestSDCreateA(numPrePeriods, 0)
    prePeriodCoefDiffs  = A_SD * prePeriod_coef
    prePeriodSigmaDiffs = A_SD * prePeriod_sigma * A_SD'
    seDiffs             = sqrt(diagonal(prePeriodSigmaDiffs))

    upperBoundVec = prePeriodCoefDiffs + invnormal(1 - alpha) * seDiffs
    return(max(upperBoundVec))
}

real matrix function _honestSDCreateA(real scalar numPrePeriods,
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

real colvector function _honestSDCreated(real scalar numPrePeriods,
                                         real scalar numPostPeriods,
                                         real scalar M,
                                         | real scalar postPeriodMomentsOnly) {

    real matrix A_SD

    if ( args() < 4 ) postPeriodMomentsOnly = 0

    // This function creates a vector for the linear constraints that \delta \in Delta^SD(M).
    // It implements this using the general characterization of d, NOT the sharp
    // characterization of the identified set.
    //
    // Inputs:
    //   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
    //   numPostPeriods = number of post-periods. This is an element of resultsObjects.
    //   M              = smoothness parameter of Delta^SD(M).
    //   postPeriodMomentsOnly = whether to exlude moments relating only to pre-period (which don't affect ID set)

    A_SD = _honestSDCreateA(numPrePeriods, numPostPeriods, postPeriodMomentsOnly)
    return(J(rows(A_SD), 1, M))
}


real rowvector function _honestSDComputeIDSet(real scalar M,
                                              real rowvector trueBeta,
                                              real colvector l_vec,
                                              real scalar numPrePeriods,
                                              real scalar numPostPeriods) {

    struct ECOS_workspace_abridged scalar result_min, result_max

    real matrix A_SD_eq, A_SD_in
    real colvector fDelta
    real rowvector d_SD_eq, d_SD_in
    real scalar id_ub, id_lb

    // This function computes the upper and lower bound of the identified set
    // given the event study coefficients, lvec and M.
    //
    // Note: lvec is assumed to be non-negative.
    //
    // Inputs:
    //   M              = smoothness param of Delta^SD
    //   trueBeta       = vector of population event study coefficients
    //   l_vec          = vector l defining parameter of interest
    //   numPrePeriods  = number of pre-periods
    //   numPostPeriods = number of post-periods
    //
    // Outputs:
    //   dataframe with columns
    //     id.ub = upper bound of ID set
    //     id.lb = lower bound of ID set
    //     M     = M value passed

    // Create objective function: Wish to min/max l'delta_post
    fDelta = J(numPrePeriods, 1, 0) \ l_vec

    // Create A, d that define Delta^SDPB(M)
    A_SD_in = _honestSDCreateA(numPrePeriods, numPostPeriods)
    d_SD_in = _honestSDCreated(numPrePeriods, numPostPeriods, M)

    // Create equality constraint that delta_pre = beta_pre
    A_SD_eq  = I(numPrePeriods), J(numPrePeriods, numPostPeriods, 0)
    d_SD_eq  = trueBeta[1..numPrePeriods]

    result_max = ECOS(-fDelta, A_SD_in, d_SD_in, rows(A_SD_in), 0, 0, A_SD_eq, d_SD_eq)
    result_min = ECOS( fDelta, A_SD_in, d_SD_in, rows(A_SD_in), 0, 0, A_SD_eq, d_SD_eq)

    if ( !(result_max.success & result_min.success) ) {
        errprintf("_honestSDComputeIDSet(): Solver did not find an optimum\n")
    }

    if ( result_max.success ) {
        id_lb = l_vec' * trueBeta[(numPrePeriods+1)..(numPrePeriods+numPostPeriods)]' + result_max.info_obj_val
    }
    else {
        id_lb = l_vec' * trueBeta[(numPrePeriods+1)..(numPrePeriods+numPostPeriods)]'
    }

    if ( result_min.success ) {
        id_ub = l_vec' * trueBeta[(numPrePeriods+1)..(numPrePeriods+numPostPeriods)]' - result_min.info_obj_val
    }
    else {
        id_ub = l_vec' * trueBeta[(numPrePeriods+1)..(numPrePeriods+numPostPeriods)]'
    }

    // Return identified set
    return((id_lb, id_ub))
}
end
