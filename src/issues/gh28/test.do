import excel using con.xlsx, clear firstrow

reg delay index 交互项 ROA lev PE BM asset turnover i.year i.sfdm
mat list e(b)
matrix l_vec = (1/3) \ (1/3) \ (1/3)	
local plotopts ylabel(0(0.02)0.22)
honestdid, l_vec(l_vec) pre(1/13) post(14/16) mvec(0) omit delta(sd)
