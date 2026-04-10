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

* Each matrix is a separate linear combination of post-period coefficients.
* The call below is equivalent to calling honestdid three times with
* l_vec(lvec0), l_vec(lvec1), l_vec(lvec2) separately, but returns
* all results in a single s(HonestEventStudy) mata object with a
* 4-column CI matrix [lvec_idx, M, lb, ub] and matching open matrix.
*
* do src/install.do
* exit, clear
* stata14-mp

capture program drop basic_multilvec
program basic_multilvec
    tempname beta sigma r0 r1 r2 rmulti
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }
    matrix l_vec0 = 1 \ 0 \ 0 \ 0
    matrix l_vec1 = 0 \ 1 \ 0 \ 0
    matrix l_vec2 = 0 \ 0 \ 1 \ 0
    local opts b(`beta') vcov(`sigma') numpre(4) delta(sd) method(FLCI) mvec(0(0.1)0.3)

    * Individual calls (sequential, reference)
    honestdid, `opts' l_vec(l_vec0) parallel(0) matasave(`r0')
    honestdid, `opts' l_vec(l_vec1) parallel(0) matasave(`r1')
    honestdid, `opts' l_vec(l_vec2) parallel(0) matasave(`r2')

    * Multi-l_vec sequential call
    honestdid, `opts' l_vecs(l_vec0 l_vec1 l_vec2) parallel(0) matasave(`rmulti')
    * mata `rmulti'.CI
    * mata `rmulti'.options.l_vec_names
    * mata `rmulti'.options.l_vec
    * mata rows(`rmulti'.OG)

    * Assert multi-l_vec matches individual sequential calls (skip original CS rows where M=.)
    mata {
        m = `rmulti'.CI[selectindex(`rmulti'.CI[.,1]:==1), 2..4]
        s = `r0'.CI
        assert(mreldif(select(m, m[.,1]:!=.), select(s, s[.,1]:!=.)) < 1e-10)
        m = `rmulti'.CI[selectindex(`rmulti'.CI[.,1]:==2), 2..4]
        s = `r1'.CI
        assert(mreldif(select(m, m[.,1]:!=.), select(s, s[.,1]:!=.)) < 1e-10)
        m = `rmulti'.CI[selectindex(`rmulti'.CI[.,1]:==3), 2..4]
        s = `r2'.CI
        assert(mreldif(select(m, m[.,1]:!=.), select(s, s[.,1]:!=.)) < 1e-10)
    }

    * Test with parallel if available
    cap parallel numprocessors
    if ( !_rc ) {
        tempname r0p r1p r2p rmultip

        * Individual calls (parallel)
        honestdid, `opts' l_vec(l_vec0) matasave(`r0p')
        honestdid, `opts' l_vec(l_vec1) matasave(`r1p')
        honestdid, `opts' l_vec(l_vec2) matasave(`r2p')

        * Multi-l_vec parallel call
        honestdid, `opts' l_vecs(l_vec0 l_vec1 l_vec2) matasave(`rmultip')

        * Parallel individual should match sequential individual
        mata {
            m = `r0p'.CI
            s = `r0'.CI
            assert(mreldif(select(m, m[.,1]:!=.), select(s, s[.,1]:!=.)) < 1e-10)
            m = `r1p'.CI
            s = `r1'.CI
            assert(mreldif(select(m, m[.,1]:!=.), select(s, s[.,1]:!=.)) < 1e-10)
            m = `r2p'.CI
            s = `r2'.CI
            assert(mreldif(select(m, m[.,1]:!=.), select(s, s[.,1]:!=.)) < 1e-10)
        }

        * Parallel multi-l_vec should match sequential multi-l_vec
        mata {
            m = `rmultip'.CI
            s = `rmulti'.CI
            assert(mreldif(select(m, m[.,2]:!=.), select(s, s[.,2]:!=.)) < 1e-10)
        }
    }
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
