cap mata mata drop _honestRMConditionalCS()
cap mata mata drop _honestRMConditionalCSFixedS()
cap mata mata drop _honestRMCreateA()
cap mata mata drop _honestRMCreated()

* Done
* ----
*
* computeConditionalCS_DeltaRM() <-> _honestRMConditionalCS()
* .computeConditionalCS_DeltaRM_fixedS() <-> _honestRMConditionalCSFixedS()
* .create_A_RM() <-> _honestRMCreateA()
* .create_d_RM() <-> _honestRMCreated()

* TODO
* ----
* xx ALL THE COMMENTS

mata
real matrix function _honestRMConditionalCS(real rowvector betahat,
                                            real matrix sigma,
                                            real scalar numPrePeriods,
                                            real scalar numPostPeriods,
                                            real scalar debug,
                                            | real colvector l_vec,
                                            real scalar Mbar,
                                            real scalar alpha,
                                            string scalar hybrid_flag,
                                            real scalar grid_lb,
                                            real scalar grid_ub,
                                            real scalar gridPoints)
{

    real matrix CIs_RM_plus_allS, CIs_RM_minus_allS, CI_s_minus, CI_s_plus, resultsGrid
    real scalar hybrid_kappa, returnLength, postPeriodMomentsOnly, sdTheta, s_i
    real colvector s_indices, sel, CIs_RM_plus_maxS, CIs_RM_minus_maxS
    real rowvector thetaGrid, thetaDiff, gridLength, gridlb, gridub

    // This function computes the ARP CI that includes nuisance parameters
    // for Delta^{RM}(Mbar). This functions uses ARP_computeCI for all
    // of its computations.
    //
    // Inputs:
    //   betahat               = vector of estimated event study coefficients
    //   sigma                 = covariance matrix of estimated event study coefficients
    //   numPrePeriods         = number of pre-periods
    //   numPostPeriods        = number of post-periods
    //   l_vec                 = vector that defines parameter of interest
    //   Mbar                  = tuning parameter for Delta^RM(Mbar), default Mbar = 0.
    //   alpha                 = desired size of CI, default alpha = 0.05
    //   hybrid_flag           = flag for hybrid, default = "LF". Must be either "LF" or "ARP"
    //   hybrid_kappa          = desired size of first-stage hybrid test, default = NULL
    //   returnLength          = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
    //   gridPoints            = number of gridpoints to test over, default = 1000
    //   postPeriodMomentsOnly = exclude moments for delta^MB that only include pre-period coefs
    //
    //  Outputs:
    //   data_frame containing upper and lower bounds of CI.

    if ( args() < 6  ) l_vec       = _honestBasis(1, numPostPeriods)
    if ( args() < 7  ) Mbar        = 0
    if ( args() < 8  ) alpha       = 0.05
    if ( args() < 9  ) hybrid_flag = "LF"

    // TODO xx hard-coded
    hybrid_kappa          = alpha/10
    returnLength          = 0
    postPeriodMomentsOnly = 1

    // Create minimal s index for looping.
    s_indices = (-(numPrePeriods - 1))::0

    // If grid.ub, grid.lb is not specified, we set these bounds to be
    // equal to the id set under parallel trends
    //     {0} +- 20*sdTheta (i.e. [-20*sdTheta, 20*sdTheta].
    sel = (numPrePeriods+1)::(numPrePeriods+numPostPeriods)
    sdTheta = sqrt(l_vec' * sigma[sel, sel] * l_vec)

    if ( args() < 10 ) grid_lb     = .
    if ( args() < 11 ) grid_ub     = .
    if ( args() < 12 ) gridPoints  = 1e3

    gridlb = missing(grid_lb)? ((betahat * l_vec) - 20*sdTheta): grid_lb
    gridub = missing(grid_ub)? ((betahat * l_vec) + 20*sdTheta): grid_ub

    // Loop over s values for (+), (-), left join the resulting CIs based on the grid
    CIs_RM_plus_allS  = J(gridPoints, length(s_indices), 0)
    CIs_RM_minus_allS = J(gridPoints, length(s_indices), 0)

    if ( debug == 1 ) {
        printf("\thonest debug: _honestRMConditionalCS()\n")
        printf("\thonest debug: \tnumPostPeriods = %g\n", numPostPeriods)
        printf("\thonest debug: \thybrid_flag    = %s\n", hybrid_flag)
    }

    for (s_i = 1; s_i <= length(s_indices); s_i++) {
        // Compute CI for s, (+) and bind it to all CI's for (+)
        CI_s_plus = _honestRMConditionalCSFixedS(s_indices[s_i],
                                                 1,
                                                 Mbar,
                                                 betahat,
                                                 sigma,
                                                 numPrePeriods,
                                                 numPostPeriods,
                                                 l_vec,
                                                 alpha,
                                                 hybrid_flag,
                                                 hybrid_kappa,
                                                 postPeriodMomentsOnly,
                                                 gridlb,
                                                 gridub,
                                                 gridPoints)

        // Compute CI for s, (-) and bind it to all CI's for (-)
        CI_s_minus = _honestRMConditionalCSFixedS(s_indices[s_i],
                                                  0,
                                                  Mbar,
                                                  betahat,
                                                  sigma,
                                                  numPrePeriods,
                                                  numPostPeriods,
                                                  l_vec,
                                                  alpha,
                                                  hybrid_flag,
                                                  hybrid_kappa,
                                                  postPeriodMomentsOnly,
                                                  gridlb,
                                                  gridub,
                                                  gridPoints)

        CIs_RM_plus_allS[., s_i]  = CI_s_plus[., 2]
        CIs_RM_minus_allS[., s_i] = CI_s_minus[., 2]
    }

    CIs_RM_plus_maxS  = rowmax(CIs_RM_plus_allS)
    CIs_RM_minus_maxS = rowmax(CIs_RM_minus_allS)

    // Take the max between (+), (-) and Construct grid containing theta
    // points and whether any CI accepted
    thetaGrid   = _honestLinspace(gridlb, gridub, gridPoints)
    resultsGrid = thetaGrid', rowmax((CIs_RM_plus_maxS, CIs_RM_minus_maxS))
    if ( returnLength ) {
        thetaDiff  = thetaGrid[2::gridPoints] :- thetaGrid[1::(gridPoints-1)]
        gridLength = ((0, thetaDiff) + (thetaDiff, 0)) :/ 2
        return(gridLength * resultsGrid[., 2])
    }
    else {
        return(resultsGrid)
    }
}

real matrix function _honestRMConditionalCSFixedS(real scalar s,
                                                  real scalar max_positive,
                                                  real scalar Mbar,
                                                  real rowvector betahat,
                                                  real matrix sigma,
                                                  real scalar numPrePeriods,
                                                  real scalar numPostPeriods,
                                                  real colvector l_vec,
                                                  real scalar alpha,
                                                  string scalar hybrid_flag,
                                                  real scalar hybrid_kappa,
                                                  real scalar postPeriodMomentsOnly,
                                                  real scalar grid_lb,
                                                  real scalar grid_ub,
                                                  real scalar gridPoints)
{

    struct _honestHybridList scalar hybrid_list
    real matrix A_RM_s, CI
    real vector d_RM, rowsForARP, postPeriodIndices

    // This function computes the ARP CI that includes nuisance parameters
    // for Delta^{RM}(Mbar) for a fixed s and (+),(-). This functions uses ARP_computeCI for all
    // of its computations. It is used as a helper function in computeConditionalCS_DeltaRM below.

    // Check that hybrid_flag equals LF or ARP
    if ( hybrid_flag != "LF" & hybrid_flag != "ARP" ) {
        errprintf("_honestRMConditionalCSFixedS(): ")
        errprintf("hybrid_flag must equal 'ARP' or 'LF'.\n")
        _error(198)
    }

    // Create hybrid_list object
    hybrid_list.hybrid_kappa = hybrid_kappa

    // Create matrix A_RM_s, and vector d_RM
    A_RM_s = _honestRMCreateA(numPrePeriods, numPostPeriods, s, Mbar, max_positive)
    d_RM   = _honestRMCreated(numPrePeriods, numPostPeriods)

    // If only use post period moments, construct indices for the post period moments only.
    if ( postPeriodMomentsOnly & (numPostPeriods > 1) ) {
        postPeriodIndices = (numPrePeriods +1)..cols(A_RM_s)
        rowsForARP = selectindex(rowsum(A_RM_s[., postPeriodIndices] :!= 0) :> 0)
    }
    else {
        rowsForARP = 1::rows(A_RM_s)
    }

    // if there is only one post-period, we use the no-nuisance parameter functions
    if (numPostPeriods == 1) {
        if (hybrid_flag == "LF") {
            // Compute LF CV and store it in hybrid_list
            hybrid_list.lf_cv = _honestARPLeastFavorableCV(A_RM_s * sigma * A_RM_s', hybrid_kappa)
        }

        // Compute confidence set
        CI = _honestARPComputeCINoNuis(betahat,
                                       sigma,
                                       numPrePeriods,
                                       numPostPeriods,
                                       A_RM_s,
                                       d_RM,
                                       0,
                                       l_vec,
                                       alpha,
                                       hybrid_flag,
                                       hybrid_list,
                                       0,
                                       grid_lb,
                                       grid_ub,
                                       gridPoints)
    }
    else {
        // CASE: NumPostPeriods > 1, use the nuisance parameter functions.
        // Compute ARP CI for l'beta using Delta^RM
        CI = _honestARPComputeCI(betahat,
                                 sigma,
                                 numPrePeriods,
                                 numPostPeriods,
                                 A_RM_s,
                                 d_RM,
                                 0,
                                 l_vec,
                                 alpha,
                                 hybrid_flag,
                                 hybrid_list,
                                 0,
                                 grid_lb,
                                 grid_ub,
                                 gridPoints,
                                 rowsForARP)
    }

    return(CI)
}

real matrix function _honestRMCreateA(real scalar numPrePeriods,
                                      real scalar numPostPeriods,
                                      real scalar s,
                                      | real scalar Mbar,
                                      real scalar max_positive,
                                      real scalar dropZero)
{

    real matrix A, Atilde, A_UB
    real rowvector v_max_dif
    real scalar r

    if ( args() < 4 ) Mbar = 1
    if ( args() < 5 ) max_positive = 1
    if ( args() < 6 ) dropZero = 1

    // This function creates a matrix for the linear constraints that
    // \delta \in Delta^RM_{s,(.)}(Mbar), where (.) is + if max_positve = T and (-) if max_positive = F.
    //
    // Inputs:
    //   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
    //   numPostPeriods = number of post-periods. This is an element of resultsObjects.

    // First construct matrix Atilde that takes first diffrences --
    // (numPrePeriods+numPostPeriods) x (numPrePeriods+numPostPeriods+1)
    //
    // Note Atilde is just the positive moments; is not related to Atilde, the rotate matrix, in the paper
    // Note: Atilde initially includes t = 0. We then drop it.

    Atilde = J(numPrePeriods+numPostPeriods, numPrePeriods+numPostPeriods+1, 0)
    for (r = 1; r <= (numPrePeriods+numPostPeriods); r++) {
        Atilde[r, r..(r+1)] = (-1, 1)
    }

    // Create a vector to extract the max first dif, which corresponds
    // with the first dif for period s, or minus this if max_positive == F
    v_max_dif = J(1, numPrePeriods + numPostPeriods + 1, 0)
    v_max_dif[(numPrePeriods+s)..(numPrePeriods+1+s)] = (-1, 1)
    if ( max_positive == 0 ) {
        v_max_dif = -v_max_dif
    }

    // The bounds for the first dif starting with period t are
    // 1*v_max_dif if t<=0 and M*v_max_dif if t>0
    A_UB = J(numPrePeriods, 1, v_max_dif) \ J(numPostPeriods, 1, Mbar :* v_max_dif)

    // Construct A that imposes |Atilde * delta | <= A_UB * delta
    A = (Atilde - A_UB \ -Atilde - A_UB)

    // Remove all-zero rows of the matrix Atilde, corresponding with the constraint
    // (delta_s - delta_s-1) - (delta_s - delta_s-1) <= (delta_s delta_s-1) - (delta_s - delta_s-1)
    A = A[selectindex(rowsum(A:^2) :> 1e-10), ]

    // Remove the period corresponding with t=0
    if ( dropZero ) {
        A = A[., _honestInverseIndex(numPrePeriods+1, numPrePeriods+numPostPeriods+1)]
    }

    return(A)
}

real colvector function _honestRMCreated(real scalar numPrePeriods,
                                         real scalar numPostPeriods,
                                         | real scalar dropZero)
{

    real matrix A_RM

    if ( args() < 3 ) dropZero = 1

    // This function creates a vector for the linear constraints that
    // delta is in Delta^RM_{s,(.)}(Mbar), where (.) is + if max_positve = T and - if max_positive = F.
    // It implements this using the general characterization of d, NOT the sharp
    // characterization of the identified set.
    //
    // Inputs:
    //   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
    //   numPostPeriods = number of post-periods. This is an element of resultsObjects.

    // d doesn't depend on Mbar or s; we just use this to get the dims right
    A_RM = _honestRMCreateA(numPrePeriods, numPostPeriods, 0, 0, dropZero)
    return(J(rows(A_RM), 1, 0))
}
end
