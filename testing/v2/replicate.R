# install.packages("remotes")
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# remotes::install_github("asheshrambachan/HonestDiD")

#######################################################################
#                                                                     #
#                               Startup                               #
#                                                                     #
#######################################################################

library(HonestDiD)
data(BCdata_EventStudy)

BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec          <- basisVector(index = 1, size = BC_numPostPeriods)

# createSensitivityResults
# ------------------------

betahat        <- BCdata_EventStudy$betahat
sigma          <- BCdata_EventStudy$sigma
numPrePeriods  <- BC_numPrePeriods
numPostPeriods <- BC_numPostPeriods
l_vec          <- BC_l_vec
Mvec           <- seq(from=0, to=0.3, by=0.1)
biasDirection  <- "negative"

alpha                 <- 0.05
method                <- "FLCI"
monotonicityDirection <- NULL
parallel              <- FALSE

# xx DeltaSD_upperBound_Mpre
# xx findOptimalFLCI

# If Delta = Delta^{SD}, construct confidence intervals
method <- "FLCI"
Delta  <- "DeltaSD"

# There are other methods, ways to specifcy Delta. Let's start here.
# Basically the basic instance of the function is a wrapper for
# findOptimalFLCI

# xx Results <- foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
# xx     temp <- findOptimalFLCI(betahat        = betahat,
# xx                             sigma          = sigma,
# xx                             numPrePeriods  = numPrePeriods,
# xx                             numPostPeriods = numPostPeriods,
# xx                             l_vec          = l_vec,
# xx                             M              = Mvec[m],
# xx                             alpha          = alpha)
# xx 
# xx     tibble(lb     = temp$FLCI[1],
# xx            ub     = temp$FLCI[2],
# xx            method = "FLCI",
# xx            Delta  = Delta,
# xx            M      = Mvec[m])
# xx }

# constructOriginalCS
# -------------------

# createSensitivityPlot
# ---------------------

# originalResults <- BC_OriginalResults

#######################################################################
#                                                                     #
#                           findOptimalFLCI                           #
#                                                                     #
#######################################################################

# FLCI_Results <- .findOptimalFLCI_helper
#
# FLCI <- c(t(FLCI_Results$optimalVec) %*% betahat - FLCI_Results$optimalHalfLength,
#           t(FLCI_Results$optimalVec) %*% betahat + FLCI_Results$optimalHalfLength)
#
# Results <- list(
#     FLCI              = FLCI,
#     optimalVec        = FLCI_Results$optimalVec,
#     optimalHalfLength = FLCI_Results$optimalHalfLength,
#     M                 = FLCI_Results$M,
#     status            = FLCI_Results$status
# )


#######################################################################
#                                                                     #
#                       .findOptimalFLCI_helper                       #
#                                                                     #
#######################################################################

numPoints <- 100

# .findHForMinimumBiasn
# ---------------------

# xx h0 <- .findHForMinimumBias(sigma          = sigma,
# xx                            numPrePeriods  = numPrePeriods,
# xx                            numPostPeriods = numPostPeriods,
# xx                            l_vec          = l_vec)

# .findLowestH
# ------------

# xx hMin <- .findLowestH(sigma          = sigma,
# xx                      numPrePeriods  = numPrePeriods,
# xx                      numPostPeriods = numPostPeriods,
# xx                      l_vec          = l_vec)

# .findWorstCaseBiasGivenH
# ------------------------

# xx hGrid  <- seq(from = hMin, to = c(h0), length.out = numPoints)
# xx biasDF <- purrr::map_dfr(.x = hGrid,
# xx                          .f = function(h) {
# xx                              .findWorstCaseBiasGivenH(h,
# xx                                                       sigma          = sigma,
# xx                                                       numPrePeriods  = numPrePeriods,
# xx                                                       numPostPeriods = numPostPeriods,
# xx                                                       l_vec          = l_vec,
# xx                                                       returnDF       = T) %>% mutate(h = h)})

# Misc manipulation
# -----------------

# xx biasDF <- biasDF %>% rename(bias = value)
# xx biasDF <-
# xx     left_join(biasDF %>% mutate(id = 1),
# xx               data.frame(m = M, id = 1),
# xx               by = "id") %>% dplyr::select(-id)
# xx biasDF <- biasDF %>%
# xx     rename(maxBias = bias) %>% filter(maxBias < Inf)

# .qfoldednormal
# --------------

# xx biasDF <- biasDF %>%
# xx     mutate(maxBias = maxBias * m) %>%
# xx     mutate(CI.halflength = .qfoldednormal(p = 1-alpha, mu = maxBias/h) * h)
# xx optimalCIDF <- biasDF %>%
# xx     group_by(m) %>%
# xx     filter(status == "optimal" | status == "optimal_inaccurate") %>%
# xx     filter(CI.halflength == min(CI.halflength))
# xx results <- list(optimalVec = c(unlist(optimalCIDF$optimal.l), l_vec),
# xx                optimalPrePeriodVec = unlist(optimalCIDF$optimal.l),
# xx                optimalHalfLength = optimalCIDF$CI.halflength,
# xx                M = optimalCIDF$m,
# xx                status = optimalCIDF$status)
# xx .qfoldednormal <- function(p, mu = 0, sd = 1, numSims = 10^6, seed = 0) {
# xx     set.seed(seed)
# xx     normDraws  <- rnorm(n = numSims, sd = sd)
# xx     pQuantiles <- purrr::map_dbl(.x = mu, .f = function(mu){quantile(abs(normDraws + mu), probs = p)})
# xx     return(pQuantiles)
# xx }

#######################################################################
#                                                                     #
#                        .findHForMinimumBias                         #
#                                                                     #
#######################################################################

w <- c(rep(0,numPrePeriods-1), t(1:numPostPeriods) %*% l_vec)

# xx hsquared <- .computeSigmaLFromW(w,
# xx                                 sigma          = sigma,
# xx                                 numPrePeriods  = numPrePeriods,
# xx                                 numPostPeriods = numPostPeriods,
# xx                                 l_vec          = l_vec)

# .createMatricesForVarianceFromW
# -------------------------------

prePeriodIndices <- 1:numPrePeriods

# Constructs matrix to compute the variance of the affine estimator for
# a choice of weights W.

SigmaPre     <- sigma[prePeriodIndices, prePeriodIndices]
SigmaPrePost <- sigma[prePeriodIndices, -prePeriodIndices]
SigmaPost    <- t(l_vec) %*% sigma[-prePeriodIndices, -prePeriodIndices] %*% l_vec
WtoLPreMat   <- diag(numPrePeriods)
if (numPrePeriods == 1) {
    WtoLPreMat = 1
} else {
    for(col in 1:(numPrePeriods-1) ) {
        WtoLPreMat[col+1, col] <- -1
    }
}

UstackWtoLPreMat <- cbind(matrix(0, nrow = numPrePeriods, ncol = numPrePeriods), WtoLPreMat)
A_quadratic_sd   <- t(UstackWtoLPreMat) %*% SigmaPre %*% UstackWtoLPreMat
A_linear_sd      <- 2 * t(UstackWtoLPreMat) %*% SigmaPrePost %*% l_vec
A_constant_sd    <- SigmaPost
A_matrices       <- list(A_quadratic_sd = A_quadratic_sd,
                         A_linear_sd    = A_linear_sd,
                         A_constant_sd  = A_constant_sd)

# .computeSigmaLFromW
# -------------------

# xx A_matrices <- .createMatricesForVarianceFromW(sigma = sigma, numPrePeriods = numPrePeriods, l_vec = l_vec)
UstackW <- matrix(c( rep(0, length(w)), w ), ncol = 1)
varL <- t(UstackW) %*% A_matrices$A_quadratic_sd %*% UstackW +
        t(A_matrices$A_linear_sd) %*% UstackW + A_matrices$A_constant_sd
hsquared <- varL
h0 <- h <- sqrt(hsquared)

#######################################################################
#                                                                     #
#                            .findLowestH                             #
#                                                                     #
#######################################################################

UstackW <- Variable(numPrePeriods + numPrePeriods)

# .createConstraints_AbsoluteValue
# --------------------------------

K <- numPrePeriods
lowerTriMat <- diag(K)
lowerTriMat[lower.tri(lowerTriMat)] <- 1
A_absolutevalue <-
    rbind(cbind( -diag(K),  lowerTriMat ),
          cbind( -diag(K), -lowerTriMat))
threshold_absolutevalue <- rep(0, NROW(A_absolutevalue))
direction_absolutevalue <- rep("<=", NROW(A_absolutevalue))
constraint <- A_absolutevalue %*% UstackW <= threshold_absolutevalue
abs_constraint <- constraint

# .createConstraints_SumWeights
# -----------------------------

# Creates constraints on the sum of the weights.
numPostPeriods = length(l_vec)
A_sumweights <- c(rep(0,numPrePeriods), rep(1,numPrePeriods))
threshold_sumweights <- t((1:numPostPeriods)) %*% l_vec
constraint <- t(A_sumweights) %*% UstackW == threshold_sumweights
sum_constraint <- constraint

# .createMatricesForVarianceFromW
# -------------------------------

# Constructs matrix to compute the variance of the affine estimator for
# a choice of weights W.
prePeriodIndices <- 1:numPrePeriods
SigmaPre <- sigma[prePeriodIndices, prePeriodIndices]
SigmaPrePost <- sigma[prePeriodIndices, -prePeriodIndices]
SigmaPost <- t(l_vec) %*% sigma[-prePeriodIndices, -prePeriodIndices] %*% l_vec
WtoLPreMat <- diag(numPrePeriods)
if (numPrePeriods == 1) {
    WtoLPreMat = 1
} else {
    for(col in 1:(numPrePeriods-1) ){
        WtoLPreMat[col+1, col] <- -1
    }
}
UstackWtoLPreMat <- cbind(matrix(0, nrow = numPrePeriods, ncol = numPrePeriods), WtoLPreMat)
A_quadratic_sd <-  t(UstackWtoLPreMat) %*% SigmaPre %*% UstackWtoLPreMat
A_linear_sd <- 2 * t(UstackWtoLPreMat) %*% SigmaPrePost %*% l_vec
A_constant_sd <- SigmaPost

A_matrices <- list(A_quadratic_sd = A_quadratic_sd,
                   A_linear_sd = A_linear_sd,
                   A_constant_sd = A_constant_sd)

# .createObjectiveObject_MinimizeSD
# ---------------------------------

# Create objective function to minimize the standard deviation.
objective <- CVXR::Minimize(CVXR::quad_form(UstackW, A_matrices$A_quadratic_sd)
                            + t(A_matrices$A_linear_sd) %*% UstackW
                            + A_matrices$A_constant_sd)
objectiveVariance <- objective


# .findLowestH
# ------------

varProblem <- CVXR::Problem(objectiveVariance, constraints = list(abs_constraint, sum_constraint))
varResult  <- CVXR::psolve(varProblem)
if(varResult$status != "optimal" & varResult$status != "optimal_inaccurate"){
    warning("Error in optimization for h0")
}
minimalVariance <- varResult$value
minimalH <- sqrt(minimalVariance)
hMin <- minimalH

#######################################################################
#                                                                     #
#                      .findWorstCaseBiasGivenH                       #
#                                                                     #
#######################################################################

UstackW <- CVXR::Variable(numPrePeriods + numPrePeriods)

# .wToLFn
# -------

.wToLFn <- function(w){
    # Converts vector from w space to l space.
    numPrePeriods <- length(w)
    WtoLPreMat <- diag(numPrePeriods)
    if (numPrePeriods == 1) {
        WtoLPreMat = 1
    } else {
        for(col in 1:(numPrePeriods-1) ){
            WtoLPreMat[col+1, col] <- -1
        }
    }
    l <- c(WtoLPreMat %*% w)
    return(l)
}

# .createObjectiveObjectForBias
# -----------------------------

# Constructs the objective function for the worst-case bias.
constant = sum(sapply(1:numPostPeriods, FUN = function(s) { abs(t(1:s) %*% l_vec[(numPostPeriods - s + 1):numPostPeriods]) })) - t((1:numPostPeriods)) %*% l_vec
objective.UstackW = CVXR::Minimize( constant +  t(c(rep(1,numPrePeriods), rep(0,numPrePeriods))) %*% UstackW )
objectiveBias <- objective.UstackW

# .createConstraintsObject_SDLessThanH
# ------------------------------------

# Creates constraint that the variance of the affine estimator must be
# less than some level h.
threshold_sd_constraint <- h^2
constraint <- CVXR::quad_form(UstackW, q_matrices$A_quadratic_sd) + t(A_matrices$A_linear_sd) %*% UstackW + A_matrices$A_constant_sd <= threshold_sd_constraint
quad_constraint <- constraint
xx A_matrices <- .createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)

# .findWorstCaseBiasGivenH
# ------------------------

M <- 1
returnDF <- TRUE

xx abs_constraint <- .createConstraints_AbsoluteValue(sigma = sigma, numPrePeriods = numPrePeriods, UstackW = UstackW)
xx sum_constraint <- .createConstraints_SumWeights(numPrePeriods = numPrePeriods, l_vec = l_vec, UstackW = UstackW)

biasProblem <- CVXR::Problem(objectiveBias, constraints = list(abs_constraint, sum_constraint, quad_constraint))
biasResult <- CVXR::psolve(biasProblem)

# Multiply objective by M (note that solution otherwise doesn't depend on M,
# so no need to run many times with many different Ms)
biasResult$value <- biasResult$value * M

# Compute the implied w and l
optimal.x <- biasResult$getValue(UstackW)
optimal.w <- optimal.x[ (length(optimal.x)/2 + 1):length(optimal.x) ]
optimal.l <- .wToLFn(optimal.w)

temp <- list(status = biasResult$status,
             value = biasResult$value,
             optimal.x = I(list(unlist(optimal.x))),
             optimal.w = I(list(unlist(optimal.w))),
             optimal.l = I(list(unlist(optimal.l))))
temp <- as.data.frame(temp, stringsAsFactors = FALSE)
