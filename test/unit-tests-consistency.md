R vs Stata Comparison
=====================

See [`unit-tests-consistency.do`](unit-tests-consistency.do) and [`unit-tests-consistency.R`](unit-tests-consistency.R) for details.

Benzarti and Carloni (2019)
---------------------------

### Compare Relative Magnitudes

Stata call:

```stata
tempname beta sigma
mata st_matrix(st_local("beta"),  _honestExampleBCBeta())
mata st_matrix(st_local("sigma"), _honestExampleBCSigma())
honestdid, numpre(4) b(`beta') vcov(`sigma')
```

R call:

```R
library(HonestDiD)
data(BCdata_EventStudy)
BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods,
             numPostPeriods = BC_numPostPeriods,
             bound          = "deviation from parallel trends")
do.call(createSensitivityResults_relativeMagnitudes, opts)
```

Results:

| Stata M | Stata lb | Stata ub | R lb     | R ub   | R M    |
| ------- | -------- | -------- | -------- | ------ | ------ |
|  0.0000 |    0.159 |    0.233 |  0.159   | 0.233  | 0      |
|  0.2222 |    0.143 |    0.248 |  0.143   | 0.248  | 0.222  |
|  0.4444 |    0.124 |    0.266 |  0.123   | 0.267  | 0.444  |
|  0.6667 |    0.101 |    0.287 |  0.101   | 0.288  | 0.667  |
|  0.8889 |    0.079 |    0.308 |  0.0786  | 0.309  | 0.889  |
|  1.1111 |    0.056 |    0.329 |  0.0551  | 0.330  | 1.11   |
|  1.3333 |    0.032 |    0.351 |  0.0308  | 0.352  | 1.33   |
|  1.5556 |    0.008 |    0.375 |  0.00722 | 0.376  | 1.56   |
|  1.7778 |   -0.016 |    0.379 | -0.0171  | 0.379  | 1.78   |
|  2.0000 |   -0.041 |    0.379 | -0.0414  | 0.379  | 2      |

### Compare Relative Magnitudes (change grid, search, method)

Stata call:

```stata
mata st_matrix("l_vec", _honestBasis(1, 4))
local opts mvec(0(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1) l_vec(l_vec)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts'
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' method(Conditional)
```

R call:

```R
BC_l_vec <- basisVector(index = 1, size = BC_numPostPeriods)
opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods,
             numPostPeriods = BC_numPostPeriods,
             l_vec          = BC_l_vec,
             gridPoints     = 100,
             grid.ub        = 1,
             grid.lb        = -1,
             bound          = "deviation from parallel trends",
             Mbarvec        = seq(from=0, to=2, by=0.5))
do.call(createSensitivityResults_relativeMagnitudes, opts)
do.call(createSensitivityResults_relativeMagnitudes, c(opts, list(method="Conditional")))
```

Results (same for both methods):

| Stata M | Stata lb | Stata ub | R lb     | R ub   | R M  |
| ------- | -------- | -------- | -------- | ------ | ---- |
|  0.0000 |    0.172 |    0.232 |  0.172   | 0.232  | 0    |
|  0.5000 |    0.131 |    0.253 |  0.131   | 0.253  | 0.5  |
|  1.0000 |    0.071 |    0.313 |  0.0707  | 0.313  | 1    |
|  1.5000 |    0.030 |    0.354 |  0.0303  | 0.354  | 1.5  |
|  2.0000 |   -0.030 |    0.414 | -0.0303  | 0.414  | 2    |

Stata call:

```stata
mata st_matrix("l_alt", _honestBasis(3, 4))
local opts mvec(0(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1) l_vec(l_alt)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts'
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' method(Conditional)
```

R call:

```R
opts$l_vec <- basisVector(index = 3, size = BC_numPostPeriods)
do.call(createSensitivityResults_relativeMagnitudes, opts)
do.call(createSensitivityResults_relativeMagnitudes, c(opts, list(method="Conditional")))
```

Results (same for both methods):

| Stata M | Stata lb | Stata ub | R lb     | R ub   | R M  |
| ------- | -------- | -------- | -------- | ------ | ---- |
|  0.0000 |    0.212 |    0.273 |  0.212   | 0.273  | 0    |
|  0.5000 |    0.071 |    0.414 |  0.0707  | 0.414  | 0.5  |
|  1.0000 |   -0.091 |    0.576 | -0.0909  | 0.576  | 1    |
|  1.5000 |   -0.253 |    0.737 | -0.253   | 0.737  | 1.5  |
|  2.0000 |   -0.434 |    0.919 | -0.434   | 0.919  | 2    |


### Compare non-Relative Magnitudes

Stata call:

```stata
honestdid, numpre(4) b(`beta') vcov(`sigma') norelmag
```

R call:

```R
opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods,
             numPostPeriods = BC_numPostPeriods)
do.call(createSensitivityResults, opts)
```

Results:

| Stata M | Stata lb | Stata ub | R lb     | R ub   | R M    |
| ------- | -------- | -------- | -------- | ------ | ------ |
|  0.0000 |    0.131 |    0.216 |  0.131   | 0.216  | 0      |
|  0.0225 |    0.183 |    0.327 |  0.184   | 0.328  | 0.0225 |
|  0.0449 |    0.174 |    0.363 |  0.175   | 0.363  | 0.0449 |
|  0.0674 |    0.152 |    0.386 |  0.152   | 0.386  | 0.0674 |
|  0.0899 |    0.130 |    0.408 |  0.130   | 0.408  | 0.0899 |
|  0.1123 |    0.107 |    0.431 |  0.107   | 0.431  | 0.112  |
|  0.1348 |    0.085 |    0.453 |  0.0847  | 0.453  | 0.135  |
|  0.1572 |    0.062 |    0.476 |  0.0622  | 0.476  | 0.157  |
|  0.1797 |    0.040 |    0.498 |  0.0397  | 0.498  | 0.180  |
|  0.2022 |    0.017 |    0.521 |  0.0173  | 0.521  | 0.202  |

### Compare non-Relative Magnitudes (change grid, search, method)

Stata call:

```stata
local opts mvec(0(0.1)0.3) l_vec(l_vec)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' norelmag method(`meth')
}

local opts mvec(0(0.1)0.3) l_vec(l_alt)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' norelmag method(`meth')
}
```

R call:

```R
opts <- list(betahat        = BCdata_EventStudy$betahat,
             sigma          = BCdata_EventStudy$sigma,
             numPrePeriods  = BC_numPrePeriods,
             numPostPeriods = BC_numPostPeriods,
             l_vec          = BC_l_vec,
             Mvec           = seq(from=0, to=0.3, by=0.1))
for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    print(do.call(createSensitivityResults, c(opts, list(method = meth))))
}

opts$l_vec <- basisVector(index = 3, size = BC_numPostPeriods)
for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    print(do.call(createSensitivityResults, c(opts, list(method = meth))))
}
```

Results for `l_vec` (Conditional, C-F, C-LF nearly identical):

| Method | Stata M | Stata lb | Stata ub | R lb     | R ub   | R M |
| ------ | ------- | -------- | -------- | -------- | ------ | --- |
|   FLCI |  0.0000 |    0.131 |    0.216 |   0.131  | 0.216  | 0   |
|   FLCI |  0.1000 |    0.119 |    0.419 |   0.119  | 0.419  | 0.1 |
|   FLCI |  0.2000 |    0.019 |    0.519 |   0.0194 | 0.519  | 0.2 |
|   FLCI |  0.3000 |   -0.081 |    0.619 |  -0.0806 | 0.619  | 0.3 |
|   C-LF |  0.0000 |    0.210 |    0.328 |   0.209  | 0.329  | 0   |
|   C-LF |  0.1000 |    0.120 |    0.419 |   0.119  | 0.419  | 0.1 |
|   C-LF |  0.2000 |    0.019 |    0.519 |   0.0191 | 0.519  | 0.2 |
|   C-LF |  0.3000 |   -0.081 |    0.618 |  -0.0809 | 0.620  | 0.3 |

Results for `l_alt`:

| Method | Stata M | Stata lb | Stata ub | R lb     | R ub   | R M |
| ------ | ------- | -------- | -------- | -------- | ------ | --- |
|   FLCI |  0.0000 |    0.144 |    0.260 |   0.144  | 0.260  | 0   |
|   FLCI |  0.1000 |   -0.244 |    1.161 |  -0.244  | 1.16   | 0.1 |
|   FLCI |  0.2000 |   -0.844 |    1.761 |  -0.844  | 1.76   | 0.2 |
|   FLCI |  0.3000 |   -1.444 |    2.361 |  -1.44   | 2.36   | 0.3 |
|   C-LF |  0.0000 |    0.339 |    0.442 |   0.339  | 0.442  | 0   |
|   C-LF |  0.1000 |   -0.243 |    1.042 |  -0.243  | 1.04   | 0.1 |
|   C-LF |  0.2000 |   -0.843 |    1.642 |  -0.843  | 1.64   | 0.2 |
|   C-LF |  0.3000 |   -1.443 |    2.242 |  -1.44   | 2.24   | 0.3 |

Lovenheim and Willen (2019)
---------------------------

### Compare non-Relative Magnitudes

Stata call:

```stata
use test/LWdata_RawData.dta, clear
qui mata stata(_honestExampleLWCall())
matrix b = 100 * e(b)
matrix V = 100^2 * e(V) * 1.038349
honestdid, pre(1/9) post(10/32) b(b) vcov(V) norelmag
```

R call:

```R
data(LWdata_EventStudy)
LW_numPrePeriods  <- length(LWdata_EventStudy$prePeriodIndices)
LW_numPostPeriods <- length(LWdata_EventStudy$postPeriodIndices)
opts <- list(betahat        = LWdata_EventStudy$betahat,
             sigma          = LWdata_EventStudy$sigma,
             numPrePeriods  = LW_numPrePeriods,
             numPostPeriods = LW_numPostPeriods)
do.call(createSensitivityResults, opts)
```

Results (FLCI):

| Stata M | Stata lb | Stata ub | R lb   | R ub   | R M    |
| ------- | -------- | -------- | ------ | ------ | ------ |
|  0.0000 |    0.794 |    1.762 |  0.795 | 1.76   | 0      |
|  0.4432 |   -2.253 |    0.193 | -2.25  | 0.194  | 0.443  |
|  0.8863 |   -2.678 |    0.759 | -2.68  | 0.758  | 0.886  |
|  1.3295 |   -3.107 |    1.227 | -3.11  | 1.23   | 1.33   |
|  1.7726 |   -3.551 |    1.670 | -3.55  | 1.67   | 1.77   |
|  2.2158 |   -3.994 |    2.113 | -3.99  | 2.11   | 2.22   |
|  2.6590 |   -4.437 |    2.557 | -4.44  | 2.56   | 2.66   |
|  3.1021 |   -4.880 |    3.000 | -4.88  | 3.00   | 3.10   |
|  3.5453 |   -5.323 |    3.443 | -5.32  | 3.44   | 3.55   |
|  3.9884 |   -5.766 |    3.886 | -5.77  | 3.89   | 3.99   |

### Compare non-Relative Magnitudes (change grid, search, method)

Stata call:

```stata
mata st_matrix("l_vec", _honestBasis(15 - (-2), 23))
mata st_matrix("l_alt", _honestBasis(15 - (-1), 23))

local opts norelmag mvec(0(0.005)0.04) l_vec(l_vec)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(`meth')
}

local opts norelmag mvec(0(0.005)0.04) l_vec(l_alt)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(`meth')
}
```

R call:

```R
LW_l_vec <- basisVector(15 - (-2), LW_numPostPeriods)
LW_l_alt <- basisVector(15 - (-1), LW_numPostPeriods)

opts <- list(betahat        = LWdata_EventStudy$betahat,
             sigma          = LWdata_EventStudy$sigma,
             numPrePeriods  = LW_numPrePeriods,
             numPostPeriods = LW_numPostPeriods,
             l_vec          = LW_l_vec,
             Mvec           = seq(from=0, to=0.04, by=0.005))

for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    print(do.call(createSensitivityResults, c(opts, list(method = meth))))
}

opts$l_vec <- LW_l_alt
for (meth in c("FLCI", "Conditional", "C-F", "C-LF")) {
    print(do.call(createSensitivityResults, c(opts, list(method = meth))))
}
```

Results for `l_vec` (Conditional and C-LF nearly identical):

| Method | Stata M | Stata lb | Stata ub | R lb    | R ub    | R M   |
| ------ | ------- | -------- | -------- | ------- | ------- | ----- |
|   FLCI |  0.0000 |   3.878  |    8.342 |   3.91  | 8.37    | 0     |
|   FLCI |  0.0050 |   2.124  |    8.178 |   1.63  | 7.69    | 0.005 |
|   FLCI |  0.0100 |  -0.083  |    8.194 |  -0.108 | 8.17    | 0.01  |
|   FLCI |  0.0150 |  -2.382  |    8.036 |  -2.46  | 7.95    | 0.015 |
|   FLCI |  0.0200 |  -4.942  |    7.514 |  -5.01  | 7.44    | 0.02  |
|   FLCI |  0.0250 |  -8.120  |    6.236 |  -8.12  | 6.24    | 0.025 |
|   FLCI |  0.0300 |  -10.555 |    5.558 |  -10.5  | 5.62    | 0.03  |
|   FLCI |  0.0350 |  -12.772 |    5.047 |  -12.7  | 5.10    | 0.035 |
|   FLCI |  0.0400 |  -14.104 |    5.390 |  -14.1  | 5.39    | 0.04  |
|   C-F  |  0.0000 |   2.920  |    3.208 |   2.95  |  3.23   | 0     |
|   C-F  |  0.0050 |   1.418  |    1.811 |   2.03  |  2.36   | 0.005 |
|   C-F  |  0.0100 |  -0.374  |    0.202 |  -0.581 |  0.0375 | 0.01  |
|   C-F  |  0.0150 |  -2.173  |   -1.234 |  -2.32  | -1.33   | 0.015 |
|   C-F  |  0.0200 |  -4.003  |   -1.300 |  -3.89  | -1.30   | 0.02  |
|   C-F  |  0.0250 |  -5.894  |   -0.636 |  -6.08  | -0.644  | 0.025 |
|   C-F  |  0.0300 |  -7.892  |    0.031 |  -8.01  |  0.0387 | 0.03  |
|   C-F  |  0.0350 |  -10.082 |    0.732 |  -10.0  |  0.729  | 0.035 |
|   C-F  |  0.0400 |  -12.632 |    1.419 |  -12.6  |  1.44   | 0.04  |
|   C-LF |  0.0000 |  -12.759 |   -3.895 |  -12.8  | -3.97   | 0     |
|   C-LF |  0.0050 |  -13.524 |   -3.290 |  -13.5  | -3.37   | 0.005 |
|   C-LF |  0.0100 |  -14.289 |   -2.675 |  -14.3  | -2.73   | 0.01  |
|   C-LF |  0.0150 |  -15.054 |   -2.034 |  -15.1  | -2.09   | 0.015 |
|   C-LF |  0.0200 |  -15.819 |   -1.378 |  -15.8  | -1.44   | 0.02  |
|   C-LF |  0.0250 |  -16.584 |   -0.714 |  -16.6  | -0.780  | 0.025 |
|   C-LF |  0.0300 |  -17.349 |   -0.052 |  -17.3  | -0.0868 | 0.03  |
|   C-LF |  0.0350 |  -18.114 |    0.671 |  -18.1  |  0.598  | 0.035 |
|   C-LF |  0.0400 |  -18.879 |    1.342 |  -18.9  |  1.30   | 0.04  |

Results for `l_alt` (Conditional and C-LF nearly identical):

| Method | Stata M | Stata lb | Stata ub | R lb    | R ub    | R M   |
| ------ | ------- | -------- | -------- | ------- | ------- | ----- |
|   FLCI |  0.0000 |   4.686  |    9.223 |   4.70  | 9.24    | 0     |
|   FLCI |  0.0050 |   3.014  |    8.949 |   2.56  | 8.50    | 0.005 |
|   FLCI |  0.0100 |   0.869  |    8.837 |   0.936 | 8.90    | 0.01  |
|   FLCI |  0.0150 |  -1.358  |    8.565 |  -1.28  | 8.64    | 0.015 |
|   FLCI |  0.0200 |  -3.848  |    7.923 |  -3.94  | 7.83    | 0.02  |
|   FLCI |  0.0250 |  -6.967  |    6.513 |  -7.02  | 6.46    | 0.025 |
|   FLCI |  0.0300 |  -9.085  |    5.979 |  -9.16  | 5.90    | 0.03  |
|   FLCI |  0.0350 |  -11.254 |    5.345 |  -11.2  | 5.41    | 0.035 |
|   FLCI |  0.0400 |  -12.416 |    5.686 |  -12.4  | 5.73    | 0.04  |
|   C-F  |  0.0000 |   3.705  |    4.017 |   3.72  | 4.03    | 0     |
|   C-F  |  0.0050 |   2.283  |    2.689 |   2.90  | 3.25    | 0.005 |
|   C-F  |  0.0100 |   0.572  |    1.153 |   0.427 | 1.04    | 0.01  |
|   C-F  |  0.0150 |  -1.146  |   -0.229 |  -1.20  | 0.255   | 0.015 |
|   C-F  |  0.0200 |  -2.895  |   -0.735 |  -3.07  | 0.734   | 0.02  |
|   C-F  |  0.0250 |  -4.706  |   -0.158 |  -4.72  | 0.151   | 0.025 |
|   C-F  |  0.0300 |  -6.626  |    0.439 |  -6.52  | 0.442   | 0.03  |
|   C-F  |  0.0350 |  -8.741  |    1.052 |  -8.65  | 1.06    | 0.035 |
|   C-F  |  0.0400 |  -11.146 |    1.676 |  -11.2  | 1.66    | 0.04  |
|   C-LF |  0.0000 |  -14.493 |   -3.032 |  -14.5  | -3.15   | 0     |
|   C-LF |  0.0050 |  -15.173 |   -2.506 |  -15.2  | -2.60   | 0.005 |
|   C-LF |  0.0100 |  -15.853 |   -1.952 |  -15.9  | -2.05   | 0.01  |
|   C-LF |  0.0150 |  -16.533 |   -1.407 |  -16.5  | -1.47   | 0.015 |
|   C-LF |  0.0200 |  -17.213 |   -0.810 |  -17.2  | -0.913  | 0.02  |
|   C-LF |  0.0250 |  -17.893 |   -0.233 |  -17.9  | -0.304  | 0.025 |
|   C-LF |  0.0300 |  -18.573 |    0.353 |  -18.6  |  0.279  | 0.03  |
|   C-LF |  0.0350 |  -19.253 |    0.983 |  -19.3  |  0.906  | 0.035 |
|   C-LF |  0.0400 |  -19.933 |    1.616 |  -19.9  |  1.54   | 0.04  |
