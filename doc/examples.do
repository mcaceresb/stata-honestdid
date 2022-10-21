* Example 1: Benzarti and Carloni (2019)
* --------------------------------------

* Section 6.2 of {browse :Rambachan and Roth (2021)} explains the
* underlying event study specification for the coefficients and
* variance-covariance matrix used below. We load the stored estimates
* since the underlying microdata is not available; further note the
* referenced {cmd:beta} and {cmd:sigma} objects only contain the entries
* corresponding to the event study coefficients. Now we construct robust
* confidence intervals for DeltaRM(Mbar). In the R Vignette, the test
* uses option bound set to "deviation from linear trend"; here only the
* analogue to "deviation from parallel trends" has been implemented.

tempname beta sigma
mata {
    st_matrix(st_local("beta"),  _honestExampleBCBeta())
    st_matrix(st_local("sigma"), _honestExampleBCSigma())
}
local opts mvec(0.5(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts'

* The results are printed to the Stata console and saved in a mata
* object. The user can specify the name of the mata object via the 
* {opt mata()} option. The default name is HonestEventStudy and the
* object's name is saved in {cmd:s(HonestEventStudy)}.  In addition 
* to the CI, all the inputs and options are saved:

mata `s(HonestEventStudy)'.CI
mata `s(HonestEventStudy)'.betahat
mata `s(HonestEventStudy)'.sigma
mata `s(HonestEventStudy)'.numPrePeriods
mata `s(HonestEventStudy)'.numPostPeriods
mata `s(HonestEventStudy)'.prePeriodIndices
mata `s(HonestEventStudy)'.postPeriodIndices
mata `s(HonestEventStudy)'.open

mata `s(HonestEventStudy)'.options.alpha
mata `s(HonestEventStudy)'.options.l_vec
mata `s(HonestEventStudy)'.options.Mvec
mata `s(HonestEventStudy)'.options.rm
mata `s(HonestEventStudy)'.options.method
mata `s(HonestEventStudy)'.options.Delta
mata `s(HonestEventStudy)'.options.grid_lb
mata `s(HonestEventStudy)'.options.grid_ub
mata `s(HonestEventStudy)'.options.gridPoints

mata _honestPrintCI(`s(HonestEventStudy)')

* For ease of use, the package also provides a way to plot the CIs
* using the {cmd:coefplot} package. This can be done when the CIs
* are computed or using the results cached in mata.

honestdid, coefplot cached
honestdid, coefplot cached xtitle(Mbar) ytitle(95% Robust CI)
graph export coefplot.pdf, replace

* Example 2: Lovenheim and Willen (2019)
* --------------------------------------

* Section 6.3 of {browse :Rambachan and Roth (2021)} explains the
* underlying event study specification and rationale for the regression
* below, which is based on Equation (20) (also see p. 11 of the 
* {browse :HonestDiD Vignette}). The data is the same one provided in the
* HonestDiD R package and can be downloaded by installing {cmd:honestdid}
* with the {cmd:all} option. Since the data is available, we run the
* regression using {cmd:reghdfe} and then use the resulting estimates.

use test/LWdata_RawData.dta, clear           
mata stata(_honestExampleLWCall())      
honestdid, pre(1/9) post(10/32) coefplot

* (As noted above, this requires the user package {cmd:reghdfe}.) Note
* we we did not need to specify the coefficient vector or the
* variance-covariance vector, and instead the function took the stored
* results {cmd e(b)} and {cmd e(V)} from {cmd:reghdfe} to do the
* computations. Further, the coefficient plot was created using the
* {cmd:coefplot} package.

* Now to mirror the Vignette:
matrix b = 100 * e(b)
matrix V = 100^2 * e(V)
mata st_matrix("l_vec", _honestBasis(15 - (-2), 23))
local opts delta(sd) mvec(0(0.005)0.04) l_vec(l_vec)
local plot coefplot xtitle(M) ytitle(95% Robust CI)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' `plot'
graph export coefplot.pdf, replace
