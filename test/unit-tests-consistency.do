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
honestdid, reference(4) b(`beta') vcov(`sigma')

mata st_matrix("l_vec", _honestBasis(1, 4))
mata st_matrix("l_alt", _honestBasis(3, 4))
local opts mvec(-2(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1) l_vec(l_vec)
honestdid, reference(4) b(`beta') vcov(`sigma') `opts'
honestdid, reference(4) b(`beta') vcov(`sigma') `opts' method(Conditional)

local opts mvec(-2(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1) l_vec(l_alt)
honestdid, reference(4) b(`beta') vcov(`sigma') `opts'
honestdid, reference(4) b(`beta') vcov(`sigma') `opts' method(Conditional)

* Base non-rm comparison
honestdid, reference(4) b(`beta') vcov(`sigma') delta(sd)

local opts mvec(-0.3(0.1)0.3) l_vec(l_vec)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, reference(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(`meth')
}

local opts mvec(-0.3(0.1)0.3) l_vec(l_alt)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, reference(4) b(`beta') vcov(`sigma') `opts' delta(sd) method(`meth')
}

* Large non-rm comparison
* -----------------------

use test/LWdata_RawData.dta, clear
mata stata(_honestExampleLWCall())
matrix b = 100 * e(b)
matrix V = 100^2 * e(V) * 1.038349
mata st_matrix("l_vec", _honestBasis(15 - (-2), 23))
mata st_matrix("l_alt", _honestBasis(15 - (-1), 23))
honestdid, pre(1/9) post(10/32) b(b) vcov(V) delta(sd)

local opts delta(sd) mvec(0(0.005)0.04) l_vec(l_vec)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(`meth')
}

local opts delta(sd) mvec(0(0.005)0.04) l_vec(l_alt)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' method(`meth')
}

* One post-period
* ---------------

local opts mvec(-0.3(0.1)0.3) delta(sd) reference(7)
foreach meth in FLCI Conditional C-F C-LF {
    honestdid, b(`beta') vcov(`sigma') `opts' method(`meth')
}

local opts mvec(-2(0.5)2) delta(rm) reference(7)
foreach meth in Conditional C-LF {
    honestdid, b(`beta') vcov(`sigma') `opts' method(`meth')
}
