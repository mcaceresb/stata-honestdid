mata mata drop basisVector()
mata mata drop createSensitivityResults()
mata mata drop constructOriginalCS()
mata mata drop createSensitivityPlot()

* glpk-5.0-1
*
* This command should work like a post-estimation command that looks for
* e(b) and e(V) by default

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

// ------------
// Example data
// ------------

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
end

mata:
real colvector function basisVector(index, size)
{
    v = J(size, 1, 0)
    v[index] = 1
    return(v)
}

void function createSensitivityResults()
{
    betahat        = BCdata_EventStudy.betahat
    sigma          = BCdata_EventStudy.sigma
    numPrePeriods  = BC_numPrePeriods
    numPostPeriods = BC_numPostPeriods
    l_vec          = BC_l_vec
    Mvec           = (0..3) / 10
    biasDirection  = "negative"
    alpha          = 0.05

    method                = "FLCI"
    monotonicityDirection = ""
    parallel              = 0

    // ------------------------
    // Ignore for now; switches
    // ------------------------
    if ( (monotonicityDirection == "") & (biasDirection == "") ) {
        Delta = "DeltaSD"
    }
    else if ( !(biasDirection == "") ) {
        if (biasDirection == "positive") {
            Delta = "DeltaSDPB"
        } else {
            Delta = "DeltaSDNB"
        }
        printf("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!\n")
    }
    else {
        if (monotonicityDirection == "increasing") {
            Delta = "DeltaSDI"
        }
        else {
            Delta = "DeltaSDD"
        }
        printf("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!\n")
    }
    // ------------------------
    // Ignore for now; switches
    // ------------------------

    Results = J(length(Mvec), 1, HonestSensitivityResults())
    for (m = 1; m <= length(Mvec); m++) {
        temp = findOptimalFLCI(betahat,
                               sigma,
                               numPrePeriods,
                               numPostPeriods,
                               l_vec,
                               Mvec[m],
                               alpha)

        Results[m].lb     = temp.FLCI[1]
        Results[m].ub     = temp.FLCI[2]
        Results[m].method = method
        Results[m].Delta  = Delta
        Results[m].M      = Mvec[m]
    }
}

void function constructOriginalCS()
{
}

void function createSensitivityPlot()
{
}
end
