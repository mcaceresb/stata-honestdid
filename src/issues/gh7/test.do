local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
qui sum year, meanonly
replace yexp2 = cond(mi(yexp2), r(max) + 1, yexp2)
qui csdid dins, time(year) ivar(stfips) gvar(yexp2) long2 notyet
csdid_estat event, window(-4 5) estore(csdid)
*csdid_plot, name(es, replace)
estimates restore csdid

// Sensitivity analysis - relative magnitude
honestdid, pre(3/6) post(7/12) mvec(0.5(0.5)2) coefplot // default relative time (Tp0)

// Other post-periods (e.g. Tp1, Tp2), changing post(./12)
honestdid, pre(3/6) post(8/12) mvec(0.5(0.5)2) coefplot
honestdid, pre(3/6) post(9/12) mvec(0.5(0.5)2) coefplot

// Sensitivity analysis using smoothness restriction
honestdid, pre(3/6) post(7/12) mvec(0(0.01)0.05) delta(sd) omit coefplot // default period

// Other periods, e.g. Tp1:
matrix l_vec= (0\1\0\0\0\0)
honestdid, pre(3/6) post(7/12) mvec(0(0.01)0.05) l_vec(l_vec) delta(sd) omit coefplot
local fails = 0
forvalues i = 1 / 20 {
    cap honestdid, pre(3/6) post(7/12) mvec(0(0.01)0.05) l_vec(l_vec) delta(sd) omit coefplot
    if _rc local ++fails
}
disp `fails'
mata st_matrix("e(b)")'
mata st_matrixcolstripe("e(b)")
