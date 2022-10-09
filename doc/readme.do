* Install here coefplot, ftools, reghdfe, plot scheme
* local github https://raw.githubusercontent.com
* ssc install coefplot,      replace
* ssc install ftools,        replace
* ssc install reghdfe,       replace
* net install scheme-modern, replace from(`github'/mdroste/stata-scheme-modern/master)
* set scheme modern

* Load data
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
l in 1/5

* Keep years before 2016. Drop the 2016 cohort
keep if (year < 2016) & (missing(yexp2) | (yexp2 != 2015))

* Create a treatment dummy
gen byte D = (yexp2 == 2014)
gen `:type year' Dyear = cond(D, year, 2013)

* Run the TWFE spec
reghdfe dins b2013.Dyear, absorb(stfips year) cluster(stfips) noconstant

local plotopts ytitle("Estimate and 95% Conf. Int.") title("Effect on dins")
coefplot, vertical yline(0) ciopts(recast(rcap)) xlabel(,angle(45)) `plotopts'
graph export doc/readme_coefplot.png, replace width(1600)

honestdid, pre(1/5) post(7/8) mvec(0.5(0.5)2)
honestdid, pre(1 2 3 4 5) post(7 8) mvec(0.5(0.5)2)

reghdfe dins b2013.year##D, absorb(stfips year) cluster(stfips) noconstant
matrix list e(b)
honestdid, pre(1/5) post(6/7) mvec(0.5(0.5)2) omit

reghdfe dins b2013.Dyear, absorb(stfips year) cluster(stfips) noconstant
honestdid, numpre(5) mvec(0.5(0.5)2) omit

mata index = 1..5, 7..8
mata st_matrix("b", st_matrix("e(b)")[index])
mata st_matrix("V", st_matrix("e(V)")[index, index])
matrix list b
matrix list V
honestdid, b(b) vcov(V) numpre(5) mvec(0.5(0.5)2)
honestdid, coefplot cached

local plotopts xtitle(Mbar) ytitle(95% Robust CI)
honestdid, cached coefplot `plotopts'
graph export doc/readme_deltarm_ex1.png, replace width(1600)

local plotopts xtitle(M) ytitle(95% Robust CI)
honestdid, pre(1/5) post(6/7) mvec(0(0.01)0.05) delta(sd) omit coefplot `plotopts'
graph export doc/readme_deltasd_ex1.png, replace width(1600)

matrix l_vec = 0.5 \ 0.5
local plotopts xtitle(Mbar) ytitle(95% Robust CI)
honestdid, pre(1/5) post(6/7) mvec(0(0.5)2) l_vec(l_vec) omit coefplot `plotopts'
graph export doc/readme_deltarm_ex2.png, replace width(1600)

* CS
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
qui sum year, meanonly
replace yexp2 = cond(mi(yexp2), r(max) + 1, yexp2)
qui csdid dins, time(year) ivar(stfips) gvar(yexp2) long2 notyet
csdid_estat event, window(-4 5) estore(csdid)
estimates restore csdid

local plotopts xtitle(Mbar) ytitle(95% Robust CI)
honestdid, pre(3/6) post(7/12) mvec(0.5(0.5)2) coefplot `plotopts'
graph export doc/readme_deltarm_csdid.png, replace width(1600)
