capture program drop basic_checks
program basic_checks
    syntax, [*]
    tempname beta sigma
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }
    matrix mvec  = 0.3, 0.2, 0, 0.1
    honestdid, b(`beta') vcov(`sigma') `options'
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0) alpha(0.00123) coefplot xlabel(,angle(60))
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0.3 0.2 0 0.1)
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0.01(0.025)0.05)
    honestdid, b(`beta') vcov(`sigma') `options' mvec(mvec)
    honestdid, b(`beta') vcov(`sigma') `options' mvec(10(10)50)
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0(0.1)0.3) alpha(0.01)
    honestdid, coefplot cached
    honestdid, coefplot cached xlabel(,angle(60))
end

capture program drop reg_checks
program reg_checks
    syntax, [*]
    use test/LWdata_RawData.dta, clear
    replace nobs = round(nobs, 1)
    qui reghdfe emp rtESV13 rtESV14 rtESV15 rtESV16 rtESV17 rtESV18 rtESV19 rtESV110 rtESV111 rtESV113 rtESV114 rtESV115 rtESV116 rtESV117 rtESV118 rtESV119 rtESV120 rtESV121 rtESV122 rtESV123 rtESV124 rtESV125 rtESV126 rtESV127 rtESV128 rtESV129 rtESV130 rtESV131 rtESV132 rtESV133 rtESV134 rtESV135 yearsfcor yearsflr aveitc fscontrol asian black hispanic other [fw = nobs], absorb(PUS_SURVEY_YEAR BIRTHSTATE PUS_SURVEY_YEAR#BIRTHYEAR) cluster(BIRTHSTATE) noconstant

    honestdid, pre(1/9) post(10/32) `options'
    honestdid, pre(1/9) post(10/32) `options' mvec(0(0.1)0.3) alpha(0.01)
    honestdid, coefplot cached
    honestdid, coefplot cached xlabel(,angle(60))
    * honestdid, pre(1/9) post(10/32) `options' mvec(0(10)50)

    matrix b = 100 * e(b)
    matrix V = 100^2 * e(V)
    honestdid, b(b) vcov(V) pre(1/9) post(10/32) `options' mvec(0.01(0.1)0.3) alpha(0.01)
    honestdid, coefplot cached
    honestdid, coefplot cached xlabel(,angle(60))
    * honestdid, b(b) vcov(V) pre(1/9) post(10/32) `options' mvec(0(10)50)
end

capture program drop basic_failures
program basic_failures
    tempname beta sigma
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }
    cap honestdid, delta(sd) b(`beta') vcov(`sigma')
    assert _rc == 198
    cap honestdid, delta(sd) b(`beta') vcov(`sigma') numpre(0)
    assert _rc == 198
    cap honestdid, delta(sd) b(`beta') vcov(`sigma') numpre(8)
    assert _rc == 198
    cap honestdid, delta(sd) b(`beta') vcov(`sigma') numpre(4) l_vec(xxx)
    assert _rc == 111
    cap honestdid, delta(sd) b(xxx) vcov(`sigma') numpre(4)
    assert _rc == 111
    cap honestdid, delta(sd) b(`beta') vcov(xxx) numpre(4)
    assert _rc == 111
    cap honestdid, delta(sd) numpre(4)
    assert _rc == 198
    cap honestdid, b(`beta') vcov(`sigma') numpre(4) delta(rm) method(FLCI)
    assert _rc == 198
    cap honestdid, b(`beta') vcov(`sigma') numpre(4) delta(rm) method(C-F)
    assert _rc == 198
    cap honestdid, b(`beta') vcov(`sigma') numpre(4) mvec(-1)
    assert _rc == 198
end

* qui do src/install.do
* qui do src/ado/honestdid.ado
* local 0 basic_checks, numpre(4) method(FLCI)
* syntax namelist(min=1 max=1), [NOIsily tab(int 4) *]
* local 0 , `options'
* tempname beta sigma
* mata {
*     st_matrix(st_local("beta"),  _honestExampleBCBeta())
*     st_matrix(st_local("sigma"), _honestExampleBCSigma())
* }
* timer clear
* timer on 1
* matrix mvec = 0, 0.3, 0.4, 0.5
* honestdid, b(`beta') vcov(`sigma') mvec(0.1) numpre(4) method(Conditional) delta(rm) debug
* * matrix lvec = 1 \ 0 \ 0 \ 0
* * honestdid, b(`beta') vcov(`sigma') mvec(0) numpre(4) l_vec(lvec) method(FLCI) delta(sd) debug
* timer off 1
* timer list
* timer clear
