library(HonestDiD)

beta  <- read.csv("beta_md.csv")
beta  <- as.vector(as.matrix(beta[2:length(beta)]))
vcov  <- read.csv("vcov_md.csv")
vcov  <- unname(as.matrix(vcov[,2:ncol(vcov)]))
res   <- createSensitivityResults(betahat        = beta,
                                  sigma          = vcov,
                                  numPrePeriods  = 11,
                                  numPostPeriods = 7,
                                  l_vec          = rep(0.1429, 7),
                                  Mvec           = seq(0, 0.004, 0.0008),
                                  method         = "FLCI")
res
