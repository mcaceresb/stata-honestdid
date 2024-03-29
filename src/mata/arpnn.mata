cap mata mata drop _honestARPComputeCINoNuis()
cap mata mata drop _honestARPTestOverGrid()
cap mata mata drop _honestTestInIdSet()
cap mata mata drop _honestTestInIdSet_FLCI_Hybrid()
cap mata mata drop _honestTestInIdSet_LF_Hybrid()

* Done
* ----
*
* .APR_computeCI_NoNuis()            <-> _honestARPComputeCINoNuis()
* .testOverThetaGrid()               <-> _honestARPTestOverGrid()
* .testInIdentifiedSet()             <-> _honestTestInIdSet()
* .testInIdentifiedSet_FLCI_Hybrid() <-> _honestTestInIdSet_FLCI_Hybrid()
* .testInIdentifiedSet_LF_Hybrid()   <-> _honestTestInIdSet_LF_Hybrid()

* TODO
* ----
* xx ALL THE COMMENTS

mata
real matrix function _honestARPComputeCINoNuis(real rowvector betahat,
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
                                               real scalar gridPoints)
{
    real matrix resultsGrid
    real rowvector thetaGrid, thetaDiff, gridLength

    // This function computes the ARP confidence interval for Delta^SD(M) for a given event study.
    // It takes the following inputs:
    //   betahat = vector of estimated event study coefficients
    //   sigma   = covariance matrix of estimated event study coefficients
    //   A       = matrix defining the set Delta
    //   d       = vector defining the set Delta
    //   numPrePeriods = number of pre-periods
    //   numPostPeriods = number of post-periods
    //   M = tuning parameter of Delta^SD(M), default M = 0
    //   postPeriod = post-period of interest
    //   alpha = size of CI, default 0.05.

    // Construct grid of tau values to test and test which values of tau lie in ID set.
    thetaGrid = _honestLinspace(grid_lb, grid_ub, gridPoints)
    if ( debug == 1 ) printf("\thonest debug: _honestARPComputeCINoNuis()\n")
    if (hybrid_flag == "ARP") {
        if ( debug == 1 ) printf("\thonest debug: \thybrid_flag = ARP\n")
        resultsGrid = _honestARPTestOverGrid(betahat,
                                             sigma,
                                             thetaGrid,
                                             A,
                                             d,
                                             alpha,
                                             numPrePeriods,
                                             &_honestTestInIdSet(),
                                             hybrid_list)
    }
    else if (hybrid_flag == "FLCI") {
        if ( debug == 1 ) printf("\thonest debug: \thybrid_flag = FLCI\n")
        resultsGrid = _honestARPTestOverGrid(betahat,
                                             sigma,
                                             thetaGrid,
                                             A,
                                             d,
                                             alpha,
                                             numPrePeriods,
                                             &_honestTestInIdSet_FLCI_Hybrid(),
                                             hybrid_list)
    }
    else if (hybrid_flag == "LF") {
        if ( debug == 1 ) printf("\thonest debug: \thybrid_flag = LF\n")
        resultsGrid = _honestARPTestOverGrid(betahat,
                                             sigma,
                                             thetaGrid,
                                             A,
                                             d,
                                             alpha,
                                             numPrePeriods,
                                             &_honestTestInIdSet_LF_Hybrid(),
                                             hybrid_list)
    }
    else {
        errprintf("_honestARPComputeCINoNuis(): ")
        errprintf("hybrid_flag must equal 'LF', 'FLCI' or 'ARP'\n")
        _error(198)
    }

    // Compute length, else return grid
    if ( returnLength ) {
        if ( debug == 1 ) printf("\thonest debug: \treturnLength = 1\n")
        thetaDiff = thetaGrid[2..gridPoints] :- thetaGrid[1..(gridPoints-1)]
        gridLength = ((0, thetaDiff) + (thetaDiff, 0)) :/ 2
        return(gridLength * resultsGrid[., 2])
    }
    else {
        if ( debug == 1 ) printf("\thonest debug: \treturnLength = 0\n")
        return(resultsGrid)
    }
}

real matrix function _honestARPTestOverGrid(real rowvector betahat,
                                            real matrix sigma,
                                            real rowvector thetaGrid,
                                            real matrix A,
                                            real colvector d,
                                            real scalar alpha,
                                            real scalar numPrePeriods,
                                            pointer(function) testFn,
                                            struct _honestHybridList scalar hybrid_list)
{
    real scalar i, theta, reject
    real colvector y, testResultsGrid

    // Tests whether values in a grid lie in the identified set.
    testResultsGrid = J(length(thetaGrid), 1, .)
    for (i = 1; i <= length(thetaGrid); i++) {
        theta  = thetaGrid[i]
        y      = betahat' - _honestBasis(numPrePeriods + 1, length(betahat)) * theta
        reject = (*testFn)(y, sigma, A, d, alpha, hybrid_list)
        testResultsGrid[i] = 1 - reject
    }
    return((thetaGrid', testResultsGrid))
}

real scalar function _honestTestInIdSet(real colvector y,
                                        real matrix sigma,
                                        real matrix A,
                                        real colvector d,
                                        real scalar alpha,
                                        struct _honestHybridList scalar hybrid_list,
                                        | real matrix Abar_additional,
                                        real colvector dbar_additional)
{
    real scalar reject, maxLocation, maxMoment, criticalVal, sigmabar
    real matrix sigmaTilde, Atilde, T_B, Abar
    real colvector dtilde, normalizedMoments, iota, gamma, dbar, c, z
    real vector i, VLoVUpVec

    if ( args() < 7 ) Abar_additional = .
    if ( args() < 8 ) dbar_additional = .

    // Runs ARP test of the moments E[AY] - d <= 0, where Y ~ N(mu, sigma)
    // and mu <= 0 under the null.  The ARP test conditions on the location
    // of the binding moment, which we write as Abar Y <= dbar.  This can be
    // used, e.g. for hybrid values with the FLCI

    sigmaTilde = sqrt(diagonal(A * sigma * A'))
    Atilde = diag(1 :/ sigmaTilde) * A
    dtilde = diag(1 :/ sigmaTilde) * d

    normalizedMoments = Atilde * y - dtilde
    maxindex(normalizedMoments, 1, i=., .)
    maxLocation = i[1]
    maxMoment   = normalizedMoments[maxLocation]

    T_B   = _honestSelectionMat(maxLocation, rows(Atilde), 1)
    iota  = J(rows(Atilde), 1, 1)
    gamma = (T_B * Atilde)'
    Abar  = Atilde - iota * T_B * Atilde
    dbar  = (I(rows(dtilde)) - iota * T_B) * dtilde

    // If statement, modifies Abar for the FLCI hybrid
    if ( !missing(Abar_additional) ) {
        Abar = Abar \ Abar_additional
    }

    if ( !missing(dbar_additional) ) {
        dbar = dbar \ dbar_additional
    }

    sigmabar  = sqrt(gamma' * sigma * gamma)
    c         = sigma * gamma :/ (gamma' * sigma * gamma)
    z         = (I(rows(y)) - c * gamma') * y
    VLoVUpVec = _honestVLoVUpFN(gamma, sigma, Abar, dbar, z)

    // Per ARP (2021), CV = max(0, c_{1-alpha}), where c_{1-alpha} is
    // the 1-alpha quantile of truncated normal.
    criticalVal = max((0, _honestARPGeneralizedNorm(1-alpha, VLoVUpVec[1], VLoVUpVec[2], T_B * dtilde, sigmabar)))
    reject = ((maxMoment + T_B * dtilde) > criticalVal)

    return(reject)
}

real scalar function _honestTestInIdSet_FLCI_Hybrid(real colvector y,
                                                    real matrix sigma,
                                                    real matrix A,
                                                    real colvector d,
                                                    real scalar alpha,
                                                    struct _honestHybridList scalar hybrid_list)
{
    real scalar reject, alphatilde
    real matrix A_firststage
    real colvector d_firststage

    // This function does a hybrid test where we first test if |l'y| > halflength
    // where l and halflength are for the FLCI of size beta
    // If so, we reject in the first stage
    // If not, we add the event that |l'y| <= halflength to the conditioning event,
    // and we adjust second stage size accordingly

    // Note that if y = (betahat - tau), then can derive that $tau \in l'\betahat \pm \chi$ iff |l'y| \leq \chi.
    // Note that this assume l places weight of 1 on \tau

    // Create A_firststage and d_firststage to capture the absolutevalue constraints
    // We reject if A_firststage %*%y - d_firststage has any positive elements
    // otherwise we add these to the constraints

    A_firststage = (hybrid_list.flci_l, -hybrid_list.flci_l)'
    d_firststage = (hybrid_list.flci_halflength \ hybrid_list.flci_halflength)

    // Run the first-stage test
    if ( max(A_firststage * y - d_firststage) > 0 ) {
        reject = 1
    } else {
        // Per ARP (2021), CV = max(0, c_{1-alpha-tilde}), where alpha-tilde = (alpha - kappa)/(1-kappa)
        // quantile of truncated normal that accounts for failing to reject in the first stage.
        alphatilde = (alpha - hybrid_list.hybrid_kappa) / (1 - hybrid_list.hybrid_kappa)
        reject = _honestTestInIdSet(y, sigma, A, d, alphatilde, hybrid_list, A_firststage, d_firststage)
    }

    return(reject)
}

real scalar function _honestTestInIdSet_LF_Hybrid(real colvector y,
                                                  real matrix sigma,
                                                  real matrix A,
                                                  real colvector d,
                                                  real scalar alpha,
                                                  struct _honestHybridList scalar hybrid_list)
{
    real scalar reject, maxLocation, maxMoment, alphatilde, criticalVal, sigmabar
    real matrix sigmaTilde, Atilde, T_B, Abar
    real colvector dtilde, normalizedMoments, iota, gamma, dbar, c, z
    real vector i, VLoVUpVec

    sigmaTilde = sqrt(diagonal(A * sigma * A'))
    Atilde = diag(1 :/ sigmaTilde) * A
    dtilde = diag(1 :/ sigmaTilde) * d

    normalizedMoments = Atilde * y - dtilde
    maxindex(normalizedMoments, 1, i=., .)
    maxLocation = i[1]
    maxMoment   = normalizedMoments[maxLocation]

    if ( maxMoment > hybrid_list.lf_cv ) {
        reject = 1
    }
    else {
        T_B   = _honestSelectionMat(maxLocation, rows(Atilde), 1)
        iota  = J(rows(Atilde), 1, 1)
        gamma = (T_B * Atilde)'
        Abar  = Atilde - iota * T_B * Atilde
        dbar  = (I(rows(dtilde)) - iota * T_B) * dtilde

        sigmabar  = sqrt(gamma' * sigma * gamma)
        c         = sigma * gamma :/ (gamma' * sigma * gamma)
        z         = (I(rows(y)) - c * gamma') * y
        VLoVUpVec = _honestVLoVUpFN(gamma, sigma, Abar, dbar, z)

        // Per ARP (2021), CV = max(0, c_{1-alpha-tilde}), where alpha-tilde = (alpha - kappa)/(1-kappa)
        // quantile of truncated normal that accounts for failing to reject in the first stage.

        alphatilde  = (alpha - hybrid_list.hybrid_kappa) / (1 - hybrid_list.hybrid_kappa)
        criticalVal = max((0, _honestARPGeneralizedNorm(1-alphatilde, VLoVUpVec[1], VLoVUpVec[2], T_B * dtilde, sigmabar)))
        reject      = ((maxMoment + T_B * dtilde) > criticalVal)
    }
    return(reject)
}
end
