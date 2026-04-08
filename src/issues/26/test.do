use ~/bulk/lib/HetEventStudies/output/df.dta, clear
gen t2 = t + 16
gen g2 = g + 16
jwdid y, ivar(i) time(t2) gvar(g2) cluster(g2) never
estat event, post
cap frame drop coefs
frame create coefs
frame coefs {
    mata: stata(sprintf("set obs %g", length(st_matrix("e(b)"))))
    gen b = .
    gen t = _n - 16
    mata: st_store(., "b", st_matrix("e(b)")')
    twoway (scatter b t), xlabel(-15(1)9)
    graph export jwdid.pdf, replace
}
shell convert jwdid.pdf jwdid.png

use ~/bulk/lib/HetEventStudies/output/df.dta, clear
did_imputation y i t g, allhorizons pretrends(15)
event_plot, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") title("Borusyak et al. (2021) imputation estimator") xlabel(-15(1)9))
graph export did.pdf, replace
