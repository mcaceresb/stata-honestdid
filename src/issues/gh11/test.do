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
matrix l_vec = 1 \ 1
honestdid, pre(1/5) post(7/8) b(b) vcov(V) mvec(0.5) l_vec(l_vec) parallel(0)
matrix l_vec = -0.5 \ -2
honestdid, pre(1/5) post(7/8) b(b) vcov(V) mvec(0.5) l_vec(l_vec) parallel(0)
honestdid, pre(1/5) post(7/8) b(b) vcov(V)
honestdid, pre(1/5) post(7/8) b(b) vcov(V) delta(sd)
honestdid, pre(1/5) post(7/8) b(b) vcov(V) grid_lb(-0.5) grid_ub(0.5)

* do src/install
* exit, clear
* stata14-mp
* cd src/issues/gh11
* cd ../../..

* local c_os_: di lower("`c(os)'")
* cap program drop honestosqp_plugin
* cap program drop honestecos_plugin
* program honestosqp_plugin, plugin using(`c(pwd)'/src/build/honestosqp_`c_os_'.plugin)
* program honestecos_plugin, plugin using(`c(pwd)'/src/build/honestecos_`c_os_'.plugin)
*
* mata
* sdVec = .006107373 \ .005716205 \ .006923872 \ .006813589
* _X_T  = -1 \ 2 \ 1 \ -2
* Cons  = (sdVec, _X_T)
* C0 = -Cons, J(rows(Cons), 1, 0)
* x0 = (-.000889745 \ .0006105562 \ .000035116 \ -.0014651852)
* f0 = (1 \ 0 \ 0)
* A0 = (J(1, length(f0)-1, 0), 1)
* ECOS_obj(f0, C0/100, -x0/100, rows(C0), 0, 0, A0, 1, 1, 1)
*
* st_numscalar("__honestecos_maxit", 10^round(-log10(min(1e-2 \ min(abs(sdVec))))))
* ECOS_obj(f0, C0, -x0, rows(C0), 0, 0, A0, 1, 1, 1)
* st_numscalar("__honestecos_maxit", J(0, 0, .))
* end
