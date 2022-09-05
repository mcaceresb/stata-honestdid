library(lfe)
library(HonestDiD)
data(BCdata_EventStudy)
data(LWdata_EventStudy)
# Rscript --no-save --no-restore --verbose test/unit-tests-consistency.R > test/unit-tests-consistency.R.log 2>&1

BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec          <- basisVector(index = 1, size = BC_numPostPeriods)
BC_l_vec          <- cbind(c(1, 0, 0, 0))
BC_l_alt          <- cbind(c(0, 0, 1, 0))
LW_numPrePeriods  <- length(LWdata_EventStudy$prePeriodIndices)
LW_numPostPeriods <- length(LWdata_EventStudy$postPeriodIndices)
LW_l_vec          <- basisVector(15 - (-2), LW_numPostPeriods)
LW_l_alt          <- basisVector(15 - (-1), LW_numPostPeriods)

cat("\n")
cat("Base relative magnitudes comparison\n")
cat("-----------------------------------\n")
cat("\n")

timer = Sys.time()
results <- createSensitivityResults_relativeMagnitudes(betahat = BCdata_EventStudy$betahat,
                                            sigma              = BCdata_EventStudy$sigma,
                                            numPrePeriods      = BC_numPrePeriods,
                                            numPostPeriods     = BC_numPostPeriods,
                                            bound              = "deviation from parallel trends")
print(results)
print(Sys.time() - timer)

timer = Sys.time()
results <- createSensitivityResults_relativeMagnitudes(betahat        = BCdata_EventStudy$betahat,
                                                       sigma          = BCdata_EventStudy$sigma,
                                                       numPrePeriods  = BC_numPrePeriods,
                                                       numPostPeriods = BC_numPostPeriods,
                                                       l_vec          = BC_l_vec,
                                                       gridPoints     = 100,
                                                       grid.ub        = 1,
                                                       grid.lb        = -1,
                                                       bound          = "deviation from parallel trends",
                                                       Mbarvec        = seq(from=0, to=2, by=0.5))
print(results)
print(Sys.time() - timer)

timer = Sys.time()
results <- createSensitivityResults_relativeMagnitudes(betahat        = BCdata_EventStudy$betahat,
                                                       sigma          = BCdata_EventStudy$sigma,
                                                       numPrePeriods  = BC_numPrePeriods,
                                                       numPostPeriods = BC_numPostPeriods,
                                                       l_vec          = BC_l_vec,
                                                       gridPoints     = 100,
                                                       grid.ub        = 1,
                                                       grid.lb        = -1,
                                                       method         = "Conditional",
                                                       bound          = "deviation from parallel trends",
                                                       Mbarvec        = seq(from=0, to=2, by=0.5))
print(results)
print(Sys.time() - timer)

timer = Sys.time()
results <- createSensitivityResults_relativeMagnitudes(betahat        = BCdata_EventStudy$betahat,
                                                       sigma          = BCdata_EventStudy$sigma,
                                                       numPrePeriods  = BC_numPrePeriods,
                                                       numPostPeriods = BC_numPostPeriods,
                                                       l_vec          = BC_l_alt,
                                                       gridPoints     = 100,
                                                       grid.ub        = 1,
                                                       grid.lb        = -1,
                                                       bound          = "deviation from parallel trends",
                                                       Mbarvec        = seq(from=0, to=2, by=0.5))
print(results)
print(Sys.time() - timer)

timer = Sys.time()
results <- createSensitivityResults_relativeMagnitudes(betahat        = BCdata_EventStudy$betahat,
                                                       sigma          = BCdata_EventStudy$sigma,
                                                       numPrePeriods  = BC_numPrePeriods,
                                                       numPostPeriods = BC_numPostPeriods,
                                                       l_vec          = BC_l_alt,
                                                       gridPoints     = 100,
                                                       grid.ub        = 1,
                                                       grid.lb        = -1,
                                                       method         = "Conditional",
                                                       bound          = "deviation from parallel trends",
                                                       Mbarvec        = seq(from=0, to=2, by=0.5))
print(results)
print(Sys.time() - timer)

cat("\n")
cat("Base non-rm comparison\n")
cat("----------------------\n")
cat("\n")

timer = Sys.time()
results <- createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                    sigma          = BCdata_EventStudy$sigma,
                                    numPrePeriods  = BC_numPrePeriods,
                                    numPostPeriods = BC_numPostPeriods)
print(results)
print(Sys.time() - timer)

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    results <- createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                        sigma          = BCdata_EventStudy$sigma,
                                        numPrePeriods  = BC_numPrePeriods,
                                        numPostPeriods = BC_numPostPeriods,
                                        l_vec          = BC_l_vec,
                                        method         = meth,
                                        Mvec           = seq(from=0, to=0.3, by=0.1))
    print(results)
    print(Sys.time() - timer)
}

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    results <- createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                        sigma          = BCdata_EventStudy$sigma,
                                        numPrePeriods  = BC_numPrePeriods,
                                        numPostPeriods = BC_numPostPeriods,
                                        l_vec          = BC_l_alt,
                                        method         = meth,
                                        Mvec           = seq(from=0, to=0.3, by=0.1))
    print(results)
    print(Sys.time() - timer)
}

cat("\n")
cat("Large non-rm comparison\n")
cat("-----------------------\n")
cat("\n")

timer = Sys.time()
results <- createSensitivityResults(betahat        = LWdata_EventStudy$betahat,
                                    sigma          = LWdata_EventStudy$sigma,
                                    numPrePeriods  = LW_numPrePeriods,
                                    numPostPeriods = LW_numPostPeriods)
print(results)
print(Sys.time() - timer)

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    results <- createSensitivityResults(betahat        = LWdata_EventStudy$betahat,
                                        sigma          = LWdata_EventStudy$sigma,
                                        numPrePeriods  = LW_numPrePeriods,
                                        numPostPeriods = LW_numPostPeriods,
                                        l_vec          = LW_l_vec,
                                        method         = meth,
                                        Mvec           = seq(from=0, to=0.04, by=0.005))
    print(results)
    print(Sys.time() - timer)
}

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    results <- createSensitivityResults(betahat        = LWdata_EventStudy$betahat,
                                        sigma          = LWdata_EventStudy$sigma,
                                        numPrePeriods  = LW_numPrePeriods,
                                        numPostPeriods = LW_numPostPeriods,
                                        l_vec          = LW_l_alt,
                                        method         = meth,
                                        Mvec           = seq(from=0, to=0.04, by=0.005))
    print(results)
    print(Sys.time() - timer)
}

# timer = Sys.time()
# createSensitivityResults_relativeMagnitudes(betahat        = LWdata_EventStudy$betahat,
#                                             sigma          = LWdata_EventStudy$sigma,
#                                             numPrePeriods  = LW_numPrePeriods,
#                                             numPostPeriods = LW_numPostPeriods,
#                                             l_vec          = LW_l_vec,
#                                             Mbarvec        = seq(from=0, to=0.04, by=0.005))
# print(Sys.time() - timer)

cat("\n")
cat("One post-period\n")
cat("---------------\n")
cat("\n")

opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods + BC_numPostPeriods - 1,
             numPostPeriods = 1,
             Mvec           = seq(from=0, to=0.3, by=0.1))
for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    print(do.call(createSensitivityResults, c(opts, list(method=meth))))
    print(Sys.time() - timer)
}

opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods + BC_numPostPeriods - 1,
             numPostPeriods = 1,
             Mbarvec        = seq(from=0, to=2, by=0.5))
for (meth in c("Conditional", "C-LF")) {
    timer = Sys.time()
    print(do.call(createSensitivityResults_relativeMagnitudes, c(opts, list(method=meth))))
    print(Sys.time() - timer)
}
