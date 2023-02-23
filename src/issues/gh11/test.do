* https://github.com/mcaceresb/stata-honestdid/issues/11
mata
b = (.00705319, .01688509, .00335541, .01051399, .0143172, 0, -.32077168, -.31881965)
V = (6.804e-07, 0,          0,         0,         0,         0, 0,         0 ) \
    (5.542e-07, 5.614e-07,  0,         0,         0,         0, 0,         0) \
    (4.199e-07, 3.536e-07,  4.390e-07, 0,         0,         0, 0,         0) \
    (4.180e-07, 3.498e-07,  4.224e-07, 4.402e-07, 0,         0, 0,         0) \
    (3.587e-07, 2.996e-07,  3.264e-07, 3.286e-07, 3.236e-07, 0, 0,         0) \
    (        0, 0,          0,         0,         0,         0, 0,         0) \
    (1.521e-07, 9.880e-08,  1.520e-07, 1.470e-07, 1.098e-07, 0, 3.529e-07, 0) \
    (4.001e-08, -6.690e-08, 5.505e-08, 5.128e-08, 1.767e-08, 0, 2.317e-07, 4.327e-07)
V = makesymmetric(V)
st_matrix("b", b)
st_matrix("V", V)
end
honestdid, pre(1/5) post(7/8) b(b) vcov(V) mvec(0.5(0.5)2)
honestdid, pre(1/5) post(7/8) b(b) vcov(V)
honestdid, pre(1/5) post(7/8) b(b) vcov(V) delta(sd)
honestdid, pre(1/5) post(7/8) b(b) vcov(V) grid_lb(-0.5) grid_ub(0.5)

* do src/install
* exit, clear
* stata14-mp
* cd src/issues/gh11
* cd ../../..

* mata st_matrix("e(V)")[7, 7]
* mata sigma[7, 7]
* mata sqrt(V[7, 7])
* mata 20 * .0005941
* mata numPrePeriods = 5
* mata numPostPeriods = 2
* mata betahat = b
* mata sdTheta = 20 * .0005941
* mata maxpre = max(abs((betahat, 0)[2..(numPrePeriods+1)] :- (betahat, 0)[1..numPrePeriods]))
*
* mata sel = (numPrePeriods+1)::(numPrePeriods+numPostPeriods)
* mata sdTheta = sqrt(l_vec' * sigma[sel, sel] * l_vec)
* mata gridlb = missing(grid_lb)? -20*sdTheta: grid_lb
*
* findfile honestdid.ado
* qui do `r(fn)'
*
* local 0 , pre(1/5) post(7/8) b(b) vcov(V)
* syntax, [ b(str) vcov(str) l_vec(str) mvec(str) grid_lb(str) grid_ub(str) gridPoints(str) alpha(passthru) omit delta(str) NUMPREperiods(int 0) PREperiodindices(numlist) POSTperiodindices(numlist) method(str) MATAsave(str) parallel(str) coefplot cached colorspec(str asis) ciopts(str) debug * ]
*
* local rm = cond(inlist("`delta'", "rm", ""), "rm", "")
* local relativeMagnitudes = "`rm'" != ""
* if "`matasave'" == "" local results HonestEventStudy
* else local results: copy local matasave
*
* if "`grid_lb'"    == ""  local grid_lb .
* if "`grid_ub'"    == ""  local grid_ub .
* if "`gridPoints'" == ""  local gridPoints 1000
*
* if "`grid_lb'"    != "." confirm number `grid_lb'
* if "`grid_ub'"    != "." confirm number `grid_ub'
* if "`gridPoints'" != ""  confirm number `gridPoints'
*
* local dohonest = ("`b'`v'`l_vec'`alpha'`mvec'`method'`preperiodindices'`postperiodindices'" != "")
* local dohonest = `dohonest' | `changegrid' | (`numpreperiods' != 0)
* HonestSanityChecks, b(`b') vcov(`vcov') l_vec(`l_vec') mvec(`mvec') method(`method') `alpha' numpre(`numpreperiods') pre(`preperiodindices') post(`postperiodindices') `rm'
*
* tempfile honestfile
* local parallel 0
*
* mata
* b             = "`b'"
* V             = "`vcov'"
* numPrePeriods = `numpreperiods'
* pre           = "`preperiodindices'"
* post          = "`postperiodindices'"
* l_vec         = "`l_vec'"
* Mvec          = "`mvec'"
* alpha         = `alpha'
* method        = "`method'"
* debug         = "`debug'"
* omit          = "`omit'"
* rm            = `relativeMagnitudes'
* grid_lb       = `grid_lb'
* grid_ub       = `grid_ub'
* gridPoints    = `gridPoints'
*
* results        = HonestDiDParse(b, V, numPrePeriods, pre, post, l_vec, Mvec, alpha, method, debug, omit, rm, grid_lb, grid_ub, gridPoints)
* betahat        = results.betahat
* sigma          = results.sigma
* numPrePeriods  = results.numPrePeriods
* numPostPeriods = results.numPostPeriods
* alpha          = results.options.alpha
* Mvec           = results.options.Mvec
* l_vec          = results.options.l_vec
* method         = results.options.method
* rm             = results.options.rm
* debug          = results.options.debug
* grid_lb        = results.options.grid_lb
* grid_ub        = results.options.grid_ub
* gridPoints     = results.options.gridPoints
* biasDirection         = ""
* monotonicityDirection = ""
* Results               = J(length(Mvec), 1, _honestResults())
* bound                 = "deviation from linear trend"
* bound                 = "deviation from parallel trends"
* hybrid_flag           = "LF"
* Delta                 = "DeltaRM"
* open                  = (hybrid_flag == "FLCI")? numPostPeriods == 1: 1
* temp_mat = _honestRMConditionalCS(betahat,
*                                   sigma,
*                                   numPrePeriods,
*                                   numPostPeriods,
*                                   debug,
*                                   l_vec,
*                                   1,
*                                   alpha,
*                                   hybrid_flag,
*                                   -0.5,
*                                   0.5,
*                                   gridPoints)
* min(temp_mat[selectindex(temp_mat[., 2]), 1])
* max(temp_mat[selectindex(temp_mat[., 2]), 1])
*     hybrid_kappa          = alpha/10
*     returnLength          = 0
*     postPeriodMomentsOnly = 1
*     s_indices = (-(numPrePeriods - 1))::0
*     sel = (numPrePeriods+1)::(numPrePeriods+numPostPeriods)
*     sdTheta = sqrt(l_vec' * sigma[sel, sel] * l_vec)
*     gridlb = missing(grid_lb)? -20*sdTheta: grid_lb
*     gridub = missing(grid_ub)?  20*sdTheta: grid_ub
*
*     // Loop over s values for (+), (-), left join the resulting CIs based on the grid
*     CIs_RM_plus_allS  = J(gridPoints, length(s_indices), 0)
*     CIs_RM_minus_allS = J(gridPoints, length(s_indices), 0)
*
* results.Results = HonestSensitivityHelper(results.betahat, results.sigma, results.numPrePeriods, results.numPostPeriods, results.options)
* results.OG      = HonestOriginalCS(results, results.options)
* results.CI      = _honestSensitivityCIMatrix(results.Results, results.OG)
* results.open    = _honestSensitivityCIOpen(results.Results, results.OG)
* _honestPrintCI(results)
* OSQP_cleanup()
* ECOS_cleanup()
* end
