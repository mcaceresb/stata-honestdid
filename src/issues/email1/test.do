do src/issues/email1/data.do
did_multiplegt Y id t D, robust_dynamic cluster(id) breps(20)
did_multiplegt Y id t D, robust_dynamic dynamic(10) placebo(10) breps(20) cluster(id)
mata st_matrix("vcov", variance(st_matrix("bootstrap")))
mata st_matrix("b", colsum(st_matrix("bootstrap"))/rows(st_matrix("bootstrap")))
honestdid, b(b) vcov(vcov) pre(12/21) post(2/11) delta(sd)
* event_plot e(estimates)#e(variances), default_look ///
*     graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") ///
*     title("did_multiplegt") xlabel(-10(1)10)) stub_lag(Effect_#) stub_lead(Placebo_#) together
