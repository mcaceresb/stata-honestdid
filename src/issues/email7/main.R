library(HonestDiD)
data(BCdata_EventStudy)
data(LWdata_EventStudy)

BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec          <- basisVector(index = 1, size = BC_numPostPeriods)
BC_l_vec          <- cbind(c(1, 0, 0, 0))

scale <- 100
og <- 
constructOriginalCS(betahat        = scale *   BCdata_EventStudy$betahat,
                    sigma          = scale^2 * BCdata_EventStudy$sigma,
                    numPrePeriods  = BC_numPrePeriods,
                    numPostPeriods = BC_numPostPeriods,
                    l_vec          = BC_l_vec)

step <- as.numeric(og[2] - og[1])
x <- 
createSensitivityResults_relativeMagnitudes(betahat        = scale *   BCdata_EventStudy$betahat,
                                            sigma          = scale^2 * BCdata_EventStudy$sigma,
                                            numPrePeriods  = BC_numPrePeriods,
                                            numPostPeriods = BC_numPostPeriods,
                                            l_vec          = BC_l_vec,
                                            gridPoints     = 10000,
                                            grid.ub        = as.numeric(og[2]+step),
                                            grid.lb        = as.numeric(og[1]-step),
                                            bound          = "deviation from parallel trends",
                                            Mbarvec        = 0)

2 * 3 * step / 1000
abs(cbind(og)[1:2]-cbind(x)[1:2]) / c(max(abs(cbind(og)[1]), 1), max(abs(cbind(og)[2]), 1))
library(TruncatedNormal)
help(norminvp)

  createSensitivityResults(betahat        = scale *   scale *   BCdata_EventStudy$betahat,
                           sigma          = scale^2 * scale^2 * BCdata_EventStudy$sigma,
                           numPrePeriods  = BC_numPrePeriods,
                           numPostPeriods = BC_numPostPeriods,
                           l_vec          = BC_l_vec,
                           method         = "FLCI",
                           Mvec           = 0)
