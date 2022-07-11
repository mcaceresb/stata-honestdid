* mata: mata clear
* mata: mata set matastrict on
* mata: mata set mataoptimize on
* cd src/build
* do ../mata/osqp.mata
* do ../mata/ecos.mata
* do ../mata/flci.mata
* do ../mata/honestdid.mata
* do ../ado/honestdid.ado

* Run these if Stata doesn't find the functions:
* ----------------------------------------------
* do osqp.mata
* do ecos.mata
* do flci.mata
* do honestdid.mata
* do honestdid.ado

***********************************************************************
*                                                                     *
*                           Matrix Example                            *
*                                                                     *
***********************************************************************

mata
betahat = (
    0.006696352,
    0.029345034,
    -0.006472972,
    0.073014989,
    0.195961118,
    0.312063903,
    0.239541546,
    0.126042500)

sigma = (
    ( 8.428358e-04,  4.768687e-04, 2.618051e-04, 0.0002354220, 0.0001676371, 0.0001128708, 1.992816e-05, -1.368265e-04) \
    ( 4.768687e-04,  6.425420e-04, 3.987425e-04, 0.0002435515, 0.0002201960, 0.0001804591, 3.843765e-05, -2.960422e-05) \
    ( 2.618051e-04,  3.987425e-04, 5.229950e-04, 0.0002117686, 0.0001840722, 0.0001458528, 7.005197e-05,  5.952995e-05) \
    ( 2.354220e-04,  2.435515e-04, 2.117686e-04, 0.0003089595, 0.0001197866, 0.0001334081, 1.016335e-04,  1.079052e-04) \
    ( 1.676371e-04,  2.201960e-04, 1.840722e-04, 0.0001197866, 0.0003599704, 0.0002478819, 1.749579e-04,  1.654257e-04) \
    ( 1.128708e-04,  1.804591e-04, 1.458528e-04, 0.0001334081, 0.0002478819, 0.0004263950, 2.171438e-04,  2.892748e-04) \
    ( 1.992816e-05,  3.843765e-05, 7.005197e-05, 0.0001016335, 0.0001749579, 0.0002171438, 4.886698e-04,  3.805322e-04) \
    (-1.368265e-04, -2.960422e-05, 5.952995e-05, 0.0001079052, 0.0001654257, 0.0002892748, 3.805322e-04,  7.617394e-04))

st_matrix("betahat", betahat)
st_matrix("sigma", sigma)
end

matrix mvec = 0.3, 0.4, 0.5
honestdid, b(betahat) vcov(sigma) reference(4) alpha(0.01)
honestdid, b(betahat) vcov(sigma) reference(4) alpha(0.01) coefplot xlabel(,angle(60))
honestdid, b(betahat) vcov(sigma) reference(4) mvec(0 0.2 0.4)
honestdid, b(betahat) vcov(sigma) reference(4) mvec(mvec)
honestdid, b(betahat) vcov(sigma) reference(4) mvec(0(0.1)0.3) coefplot
* graph export ../../test/coefplot.pdf, replace

***********************************************************************
*                                                                     *
*                            Data Exmaple                             *
*                                                                     *
***********************************************************************

* use ../../test/LWdata_RawData.dta, clear
use LWdata_RawData.dta, clear
replace nobs = round(nobs, 1)
reghdfe emp rtESV13 rtESV14 rtESV15 rtESV16 rtESV17 rtESV18 rtESV19 rtESV110 rtESV111 rtESV113 rtESV114 rtESV115 rtESV116 rtESV117 rtESV118 rtESV119 rtESV120 rtESV121 rtESV122 rtESV123 rtESV124 rtESV125 rtESV126 rtESV127 rtESV128 rtESV129 rtESV130 rtESV131 rtESV132 rtESV133 rtESV134 rtESV135 yearsfcor yearsflr aveitc fscontrol asian black hispanic other [fw = nobs], absorb(PUS_SURVEY_YEAR BIRTHSTATE PUS_SURVEY_YEAR#BIRTHYEAR) cluster(BIRTHSTATE) noconstant
honestdid, pre(1/9) post(10/32)
matrix b = 100 * e(b)
matrix V = 100^2 * e(V)
honestdid, b(b) vcov(V) pre(1/9) post(10/32) coefplot mvec(0(0.1)0.3)
* graph export ../../test/coefplot2.pdf, replace
