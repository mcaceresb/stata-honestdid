* Run the regression
* ------------------
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
keep if (year < 2016) & (missing(yexp2) | (yexp2 != 2015))
gen byte D = (yexp2 == 2014)
gen `:type year' Dyear = cond(D, year, 2013)
reghdfe dins b2013.Dyear, absorb(stfips year) cluster(stfips) noconstant

* Run HonestDiD and put results into a matrix
* -------------------------------------------
honestdid, pre(1/5) post(7/8) mvec(0.5(0.5)2) nocoefplot
mata st_matrix("honesttest", `s(HonestEventStudy)'.CI)
matrix colnames honesttest = M lb ub
matlist honesttest, names(columns)

* Export just the matrix to excel
* -------------------------------
putexcel set honesttest.xlsx, modify sheet(honesttest)
putexcel A1 = matrix(honesttest), colnames
putexcel save
putexcel clear

* Export with additional info
* ---------------------------
mata st_numscalar("honestalpha",  `s(HonestEventStudy)'.options.alpha)
mata st_strscalar("honestmethod", `s(HonestEventStudy)'.options.method)
mata st_strscalar("honestDelta",  `s(HonestEventStudy)'.options.Delta)

putexcel set honesttest.xlsx, modify sheet(honesttest)
putexcel A1 = matrix(honesttest), colnames
putexcel A`=`:rowsof honesttest'+2' = "alpha"
putexcel A`=`:rowsof honesttest'+3' = "method"
putexcel A`=`:rowsof honesttest'+4' = "Delta"
putexcel B`=`:rowsof honesttest'+2' = honestalpha
putexcel B`=`:rowsof honesttest'+3' = honestmethod
putexcel B`=`:rowsof honesttest'+4' = honestDelta
putexcel save
putexcel clear
