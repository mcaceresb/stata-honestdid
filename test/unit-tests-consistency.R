# remove.packages("HonestDiD")
# install.packages(".", repos=NULL, type="source")
# quit(save="no")
# R
# source('~/projects/ra/roth/honestDiD/references/HonestDiD/R/arp-nuisance.R')
# Rscript --no-save --no-restore --verbose test/unit-tests-consistency.R > test/unit-tests-consistency.R.log 2>&1

library(HonestDiD)
library(haven)
library(lfe)
data(BCdata_EventStudy)

# Write results to CSV for automated Stata comparison.
# Output goes to test/r-results/; re-run this script to refresh.
OUT <- "test/r-results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_results <- function(results, filename, outdir=OUT) {
  df <- data.frame(method = results$method,
                   M      = if ("M" %in% colnames(results)) results$M else results$Mbar,
                   lb     = results$lb,
                   ub     = results$ub)
  write.csv(df, file.path(outdir, filename), row.names = FALSE)
  cat(sprintf("  -> %s (%d rows)\n", filename, nrow(df)))
}

BC_numPrePeriods     <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods    <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec             <- basisVector(index = 1, size = BC_numPostPeriods)
BC_l_vec             <- cbind(c(1, 0, 0, 0))
BC_l_alt             <- cbind(rep(1, BC_numPostPeriods)/BC_numPostPeriods)

LWdata_RawData       <- haven::read_dta("test/LWdata_RawData.dta")
LWdata_RawData$nobs  <- round(LWdata_RawData$nobs)
LWdata_Results       <- lfe::felm(emp ~
    rtESV13  + rtESV14  + rtESV15  +
    rtESV16  + rtESV17  + rtESV18  +
    rtESV19  + rtESV110 + rtESV111 + # End Pre-periods
    rtESV113 + rtESV114 + rtESV115 +
    rtESV116 + rtESV117 + rtESV118 +
    rtESV119 + rtESV120 + rtESV121 +
    rtESV122 + rtESV123 + rtESV124 +
    rtESV125 + rtESV126 + rtESV127 +
    rtESV128 + rtESV129 + rtESV130 +
    rtESV131 + rtESV132 + rtESV133 +
    rtESV134 + rtESV135 + # End post-periods
    yearsfcor + yearsflr + aveitc + fscontrol +
    asian + black + hispanic + other |
    factor(PUS_SURVEY_YEAR)*factor(BIRTHYEAR) +
    factor(PUS_SURVEY_YEAR) + factor(BIRTHSTATE) |
    0 | BIRTHSTATE,
    data = LWdata_RawData,
    weights = LWdata_RawData$nobs)

# Extract coefficients for the event-study variables (first 32: 9 pre + 23 post)
LW_betahat    <- 100 * coef(LWdata_Results)[1:32]
LW_sigma      <- 100^2 * LWdata_Results$clustervcv[1:32, 1:32]
LW_numPrePeriods  <- 9
LW_numPostPeriods <- 23
LW_l_vec          <- basisVector(15 - (-2), LW_numPostPeriods)
LW_l_alt          <- cbind(rep(1, LW_numPostPeriods) / LW_numPostPeriods)

# Export betahat and sigma so Stata can use identical inputs
write.csv(data.frame(betahat=LW_betahat), file.path(OUT, "lw_betahat.csv"), row.names = FALSE)
write.csv(as.data.frame(LW_sigma),        file.path(OUT, "lw_sigma.csv"),   row.names = FALSE)

cat("\n")
cat("Base relative magnitudes comparison\n")
cat("-----------------------------------\n")
cat("\n")

timer = Sys.time()
results <- createSensitivityResults_relativeMagnitudes(betahat        = BCdata_EventStudy$betahat,
                                                       sigma          = BCdata_EventStudy$sigma,
                                                       numPrePeriods  = BC_numPrePeriods,
                                                       numPostPeriods = BC_numPostPeriods,
                                                       gridPoints     = 1000,
                                                       grid.ub        = -0.4,
                                                       grid.lb        =  0.8,
                                                       bound          = "deviation from parallel trends")
print(results)
write_results(results, "bc_rm_default.csv")
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
write_results(results, "bc_rm_lvec_clf.csv")
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
write_results(results, "bc_rm_lvec_conditional.csv")
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
write_results(results, "bc_rm_lalt_clf.csv")
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
write_results(results, "bc_rm_lalt_conditional.csv")
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
write_results(results, "bc_sd_default.csv")
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
    tag <- tolower(gsub("-", "", meth))
    write_results(results, sprintf("bc_sd_lvec_%s.csv", tag))
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
    tag <- tolower(gsub("-", "", meth))
    write_results(results, sprintf("bc_sd_lalt_%s.csv", tag))
    print(Sys.time() - timer)
}

cat("\n")
cat("Large non-rm comparison\n")
cat("-----------------------\n")
cat("\n")

timer = Sys.time()
results <- createSensitivityResults(betahat        = LW_betahat,
                                    sigma          = LW_sigma,
                                    numPrePeriods  = LW_numPrePeriods,
                                    numPostPeriods = LW_numPostPeriods)
print(results)
write_results(results, "lw_sd_default.csv")
print(Sys.time() - timer)

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    results <- createSensitivityResults(betahat        = LW_betahat,
                                        sigma          = LW_sigma,
                                        numPrePeriods  = LW_numPrePeriods,
                                        numPostPeriods = LW_numPostPeriods,
                                        l_vec          = LW_l_vec,
                                        method         = meth,
                                        Mvec           = seq(from=0, to=0.04, by=0.005))
    print(results)
    tag <- tolower(gsub("-", "", meth))
    write_results(results, sprintf("lw_sd_lvec_%s.csv", tag))
    print(Sys.time() - timer)
}

# TODO: As of April, 2026 doing 0.005 leads to a discrepency with Stata at
# TODO: 0.005 for FLCI.  Determine if numerical precision or deep issue.

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    timer = Sys.time()
    results <- createSensitivityResults(betahat        = LW_betahat,
                                        sigma          = LW_sigma,
                                        numPrePeriods  = LW_numPrePeriods,
                                        numPostPeriods = LW_numPostPeriods,
                                        l_vec          = LW_l_alt,
                                        method         = meth,
                                        Mvec           = seq(from=0, to=0.04, by=0.01))
    print(results)
    tag <- tolower(gsub("-", "", meth))
    write_results(results, sprintf("lw_sd_lalt_%s.csv", tag))
    print(Sys.time() - timer)
}

# timer = Sys.time()
# createSensitivityResults_relativeMagnitudes(betahat        = LW_betahat,
#                                             sigma          = LW_sigma,
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
    results <- do.call(createSensitivityResults, c(opts, list(method=meth)))
    print(results)
    tag <- tolower(gsub("-", "", meth))
    write_results(results, sprintf("onepost_sd_%s.csv", tag))
    print(Sys.time() - timer)
}

opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods + BC_numPostPeriods - 1,
             numPostPeriods = 1,
             Mbarvec        = seq(from=0, to=2, by=0.5))
for (meth in c("Conditional", "C-LF")) {
    timer = Sys.time()
    results <- do.call(createSensitivityResults_relativeMagnitudes, c(opts, list(method=meth)))
    print(results)
    tag <- tolower(gsub("-", "", meth))
    write_results(results, sprintf("onepost_rm_%s.csv", tag))
    print(Sys.time() - timer)
}
