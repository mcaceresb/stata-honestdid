load("debug.RData")
ls()

post = num_ann+num_post
delta_sd_results <-
HonestDiD::createSensitivityResults(betahat = 1e2 * betahat,
sigma = 1e2 * sigma,
numPrePeriods = num_pre,
numPostPeriods = post, # ann+treat
Mvec = seq(from = 0, to = 0.04, by =0.01),
l_vec = c(rep(0,num_ann),rep(1/num_post,num_post)),
parallel = FALSE)
