# install.packages("remotes")
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# remotes::install_github("asheshrambachan/HonestDiD")
library(HonestDiD)
data(BCdata_EventStudy)

# createSensitivityResults and createSensitivityPlot
# Just have it for delta SD

BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec          <- basisVector(index = 1, size = BC_numPostPeriods)

BC_DeltaSDNB_RobustResults <-
    createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                             sigma          = BCdata_EventStudy$sigma,
                             numPrePeriods  = BC_numPrePeriods,
                             numPostPeriods = BC_numPostPeriods,
                             l_vec          = BC_l_vec,
                             method         = "FLCI",
                             Mvec           = seq(from=0, to=0.3, by=0.1),
                             biasDirection  = "negative")

BC_OriginalResults <- 
    constructOriginalCS(betahat        = BCdata_EventStudy$betahat,
                        sigma          = BCdata_EventStudy$sigma,
                        numPrePeriods  = BC_numPrePeriods,
                        numPostPeriods = BC_numPostPeriods,
                        l_vec          = BC_l_vec)

BC_DeltaSDNB_SensitivityPlot <-
    createSensitivityPlot(robustResults   = BC_DeltaSDNB_RobustResults,
                          originalResults = BC_OriginalResults)
BC_DeltaSDNB_SensitivityPlot

