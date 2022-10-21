* stata14-mp -b do test/unit-tests-consistency.do
clear all
log close _all
log using test/unit-tests-consistency.do.log, replace
set rmsg on
tempname beta sigma
mata {
    st_matrix(st_local("beta"),  _honestExampleBCBeta())
    st_matrix(st_local("sigma"), _honestExampleBCSigma())
}

* Base relative magnitudes comparison
honestdid, numpre(4) b(`beta') vcov(`sigma')

mata st_matrix("l_vec", _honestBasis(1, 4))
mata st_matrix("l_alt", J(4, 1, 1)/4)
local opts mvec(0(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1) l_vec(l_vec)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts'
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' method(Conditional)

local opts mvec(0(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1) l_vec(l_alt)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts'
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' method(Conditional)

* Base non-rm comparison
honestdid, numpre(4) b(`beta') vcov(`sigma') delta(sd)

local opts mvec(0(0.1)0.3) l_vec(l_vec)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(FLCI)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(Conditional)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(C-F)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(C-LF)

local opts mvec(0(0.1)0.3) l_vec(l_alt)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(FLCI)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(Conditional)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(C-F)
honestdid, numpre(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(C-LF)

* Large non-rm comparison
* -----------------------

use test/LWdata_RawData.dta, clear
mata stata(_honestExampleLWCall())
matrix b = 100 * e(b)
matrix V = 100^2 * e(V) * 1.038349
mata st_matrix("l_vec", _honestBasis(15 - (-2), 23))
mata st_matrix("l_alt", J(23, 1, 1)/23)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) delta(sd)

local opts delta(sd) mvec(0(0.005)0.04) l_vec(l_vec)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(FLCI)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(Conditional)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(C-F)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(C-LF)

local opts delta(sd) mvec(0(0.005)0.04) l_vec(l_alt)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(FLCI)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(Conditional)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(C-F)
honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(C-LF)

* One post-period
* ---------------

local opts mvec(0(0.1)0.3) delta(sd) numpre(7)
honestdid, b(`beta') vcov(`sigma') `opts' method(FLCI)
honestdid, b(`beta') vcov(`sigma') `opts' method(Conditional)
honestdid, b(`beta') vcov(`sigma') `opts' method(C-F)
honestdid, b(`beta') vcov(`sigma') `opts' method(C-LF)

local opts mvec(0(0.5)2) delta(rm) numpre(7)
honestdid, b(`beta') vcov(`sigma') `opts' method(FLCI)
honestdid, b(`beta') vcov(`sigma') `opts' method(Conditional)
honestdid, b(`beta') vcov(`sigma') `opts' method(C-F)
honestdid, b(`beta') vcov(`sigma') `opts' method(C-LF)
