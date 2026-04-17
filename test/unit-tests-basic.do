* Structural assertion: check CI dimensions, lb <= ub, struct field consistency.
* Call after each non-cached honestdid invocation.
capture mata mata drop _test_assert_ci()
mata
void function _test_assert_ci(struct HonestEventStudy scalar res) {
    real scalar c, K, L, n
    real matrix ci

    K = cols(res.options.l_vec)
    L = length(res.options.Mvec)
    n = res.numPrePeriods + res.numPostPeriods

    // CI column count: 3 for single l_vec, 4 for multi
    c = cols(res.CI)
    if ( K == 1 ) assert(c == 3)
    else          assert(c == 4)

    // CI row count: K * (L + 1) where +1 is original CS row
    assert(rows(res.CI) == K * (L + 1))

    // lb <= ub for non-missing rows
    ci = select(res.CI, rowmissing(res.CI) :< c)
    assert(rows(ci) > 0)
    if ( c == 3 ) assert(all(ci[.,2] :<= ci[.,3]))
    else          assert(all(ci[.,3] :<= ci[.,4]))

    // l_vec rows must equal numPostPeriods
    assert(rows(res.options.l_vec) == res.numPostPeriods)

    // betahat and sigma dimensions
    assert(length(res.betahat) == n)
    assert(rows(res.sigma) == n)
    assert(cols(res.sigma) == n)

    // alpha is valid
    assert((res.options.alpha > 0) & (res.options.alpha < 1))

    // OG has K entries (one original CS per l_vec)
    assert(rows(res.OG) == K)

    // open: K columns for multi (idx + value), 1 column for single
    if ( K == 1 ) assert(cols(res.open) == 1)
    else          assert(cols(res.open) == 2)
}
end

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
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0) alpha(0.00123) coefplot xlabel(,angle(60))
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0.3 0.2 0 0.1)
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0.01(0.025)0.05)
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, b(`beta') vcov(`sigma') `options' mvec(mvec)
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, b(`beta') vcov(`sigma') `options' mvec(10(10)50)
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, b(`beta') vcov(`sigma') `options' mvec(0(0.1)0.3) alpha(0.01)
    mata _test_assert_ci(`s(HonestEventStudy)')
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
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, pre(1/9) post(10/32) `options' mvec(0(0.1)0.3) alpha(0.01)
    mata _test_assert_ci(`s(HonestEventStudy)')
    honestdid, coefplot cached
    honestdid, coefplot cached xlabel(,angle(60))
    * honestdid, pre(1/9) post(10/32) `options' mvec(0(10)50)

    matrix b = 100 * e(b)
    matrix V = 100^2 * e(V)
    honestdid, b(b) vcov(V) pre(1/9) post(10/32) `options' mvec(0.01(0.1)0.3) alpha(0.01)
    mata _test_assert_ci(`s(HonestEventStudy)')
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

* Insert a zero column into beta/sigma to simulate an omitted coefficient,
* set the Stata column stripe so _ms_omit_info detects it, run with omit,
* and compare against the clean (no extra column) result.
capture program drop basic_omit
program basic_omit
    tempname beta sigma b_clean V_clean b_omit V_omit r_clean r_omit
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }

    * Clean version: 8 coefficients, numpre=4
    matrix `b_clean' = `beta'
    matrix `V_clean' = `sigma'
    local opts numpre(4) delta(sd) method(FLCI) mvec(0(0.5)1)
    honestdid, b(`b_clean') vcov(`V_clean') `opts' matasave(`r_clean')

    * Insert a zero at position 5 (first post-period), making 9 coefficients.
    * The 5th column is the omitted one.
    tempname b8
    tempname V8 b9 V9
    mata {
        `b8' = st_matrix(st_local("b_clean"))
        `V8' = st_matrix(st_local("V_clean"))
        `b9' = `b8'[1..4], 0, `b8'[5..8]
        `V9' = (`V8'[.,1..4],  J(8,1,0),  `V8'[.,5..8])
        `V9' = (`V9'[1..4,.] \ J(1,9,0) \ `V9'[5..8,.])
        st_matrix(st_local("b_omit"), `b9')
        st_matrix(st_local("V_omit"), `V9')
    }

    * Label the 5th column as omitted via Stata's o. convention
    matrix colnames `b_omit' = x1 x2 x3 x4 o.x5 x6 x7 x8 x9
    matrix colnames `V_omit' = x1 x2 x3 x4 o.x5 x6 x7 x8 x9
    matrix rownames `V_omit' = x1 x2 x3 x4 o.x5 x6 x7 x8 x9

    honestdid, b(`b_omit') vcov(`V_omit') `opts' omit matasave(`r_omit')

    * Assert omit result matches clean result
    mata assert(mreldif(`r_clean'.CI, `r_omit'.CI) < 1e-10)
end

* matasave field correctness (option-specific checks beyond _test_assert_ci)
capture program drop basic_matasave
program basic_matasave
    tempname beta sigma r1 r2
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }

    * Single l_vec, 4 M values
    local opts numpre(4) delta(sd) method(FLCI) mvec(0 0.1 0.2)
    honestdid, b(`beta') vcov(`sigma') `opts' matasave(`r1')
    mata {
        _test_assert_ci(`r1')
        assert(`r1'.options.alpha  == 0.05)
        assert(`r1'.options.method == "FLCI")
        assert(`r1'.options.Mvec   == (0, 0.1, 0.2))
    }

    * Multi l_vec, 3 M values: check option-specific fields
    matrix l_vecA = 1 \ 0 \ 0 \ 0
    matrix l_vecB = 0 \ 1 \ 0 \ 0
    honestdid, b(`beta') vcov(`sigma') `opts' l_vecs(l_vecA l_vecB) matasave(`r2')
    mata {
        _test_assert_ci(`r2')
        assert(`r2'.options.l_vec_names[1] == "l_vecA")
        assert(`r2'.options.l_vec_names[2] == "l_vecB")
        assert(`r2'.options.Mvec == (0, 0.1, 0.2))
    }
end

* edge cases
capture program drop basic_edgecases
program basic_edgecases
    tempname beta sigma r1p r1 rdup1 rdup2
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }

    * 1 post-period (numpre=7 with 8-element beta)
    local opts delta(sd) method(FLCI) mvec(0 0.1)
    honestdid, b(`beta') vcov(`sigma') numpre(7) `opts' matasave(`r1p')
    mata {
        _test_assert_ci(`r1p')
        assert(`r1p'.numPostPeriods == 1)
        assert(rows(`r1p'.options.l_vec) == 1)
    }

    * 1 post-period with l_vecs (single scalar l_vec)
    matrix l_s = 1
    honestdid, b(`beta') vcov(`sigma') numpre(7) `opts' l_vec(l_s) matasave(`r1')
    mata {
        _test_assert_ci(`r1')
        assert(mreldif(`r1p'.CI, `r1'.CI) < 1e-10)
    }

    * duplicate l_vecs should produce identical rows
    matrix l_vecX = 1 \ 0 \ 0 \ 0
    honestdid, b(`beta') vcov(`sigma') numpre(4) `opts' matasave(`rdup1') l_vecs(l_vecX l_vecX)
    honestdid, b(`beta') vcov(`sigma') numpre(4) `opts' matasave(`rdup2') l_vec(l_vecX)
    tempname d1 d2
    mata {
        _test_assert_ci(`rdup1')
        `d1' = `rdup1'.CI[selectindex(`rdup1'.CI[.,1]:==1), 2..4]
        `d2' = `rdup1'.CI[selectindex(`rdup1'.CI[.,1]:==2), 2..4]
        assert(mreldif(select(`d1', `d1'[.,1]:!=.), select(`d2', `d2'[.,1]:!=.)) < 1e-10)
        assert(mreldif(select(`d1', `d1'[.,1]:!=.), select(`rdup2'.CI, `rdup2'.CI[.,1]:!=.)) < 1e-10)
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
