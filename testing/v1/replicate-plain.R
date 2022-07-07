# ---------------------------------------------------------------------
# basisVector

BC_DeltaSDNB_RobustResults <-
    createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                             sigma          = BCdata_EventStudy$sigma,
                             numPrePeriods  = BC_numPrePeriods-1,
                             numPostPeriods = BC_numPostPeriods+1,
                             l_vec          = basisVector(index = 1, size = BC_numPostPeriods+1),
                             Mvec           = seq(from=0, to=0.3, by=0.1),
                             biasDirection  = "negative")

# ---------------------------------------------------------------------
# createSensitivityResults

betahat        <- BCdata_EventStudy$betahat
sigma          <- BCdata_EventStudy$sigma
numPrePeriods  <- BC_numPrePeriods
numPostPeriods <- BC_numPostPeriods
l_vec          <- BC_l_vec
Mvec           <- seq(from=0, to=0.3, by=0.1)
biasDirection  <- "negative"
alpha          <- 0.05

method = "FLCI"
monotonicityDirection = NULL
parallel = FALSE

# xx DeltaSD_upperBound_Mpre
# xx findOptimalFLCI

# Recycling array of length 1 in vector-array arithmetic is deprecated.
# Use c() or as.vector() instead.

createSensitivityResults <- function(betahat, sigma,
                                     numPrePeriods, numPostPeriods,
                                     method = NULL,
                                     Mvec = NULL,
                                     l_vec = .basisVector(index = 1, size = numPostPeriods),
                                     monotonicityDirection = NULL,
                                     biasDirection = NULL,
                                     alpha = 0.05,
                                     parallel = FALSE) {

    # If Mvec is null, construct default Mvec
    if (is.null(Mvec)) {
        if (numPrePeriods == 1) {
            Mvec <- seq(from = 0, to = sqrt(sigma[1, 1]), length.out = 10)
        } else {
            Mub <- DeltaSD_upperBound_Mpre(betahat       = betahat,
                                           sigma         = sigma,
                                           numPrePeriods = numPrePeriods,
                                           alpha         = 0.05)
            Mvec <- seq(from = 0, to = Mub, length.out = 10)
        }
    }

    # If Delta = Delta^{SD}, construct confidence intervals
    method <- "FLCI"
    if (is.null(monotonicityDirection) & is.null(biasDirection)) {
        Delta <- "DeltaSD"
    } else if (!is.null(biasDirection)) {
        if (biasDirection == "positive") {
            Delta <- "DeltaSDPB"
        } else {
            Delta <- "DeltaSDNB"
        }
        warning("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!")
    } else {
        if (monotonicityDirection == "increasing") {
            Delta <- "DeltaSDI"
        } else {
            Delta <- "DeltaSDD"
        }
        warning("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!")
    }

    # There are other methods, ways to specifcy Delta. Let's start here.
    # Basically the basic instance of the function is a wrapper for
    # findOptimalFLCI.

    Results <- foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
        temp <- findOptimalFLCI(betahat        = betahat,
                                sigma          = sigma,
                                numPrePeriods  = numPrePeriods,
                                numPostPeriods = numPostPeriods,
                                l_vec          = l_vec,
                                M              = Mvec[m],
                                alpha          = alpha)

        tibble(lb     = temp$FLCI[1],
               ub     = temp$FLCI[2],
               method = "FLCI",
               Delta  = Delta,
               M      = Mvec[m])
    }

    return(Results)
}

# ---------------------------------------------------------------------
# constructOriginalCS



# ---------------------------------------------------------------------
# createSensitivityPlot
