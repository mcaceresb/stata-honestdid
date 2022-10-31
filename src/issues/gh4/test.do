set seed 1729
capture program drop test
program test, eclass
    tempname b V
    mata {
        _b = (
            .22417568,
            -.37208402,
            -.6068593,
            -.47421248,
            -.50601109,
            -.54184265,
            -.68990523,
            -.87048126,
            -1.4504416,
            -.11910901,
            -.28366669,
            .49299721,
            .53201816,
            .76635512,
            .63850458,
            .79138165,
            1.3078362,
            .92813858,
            2.6659446,
            2.3487094,
            1.0566333,
            .59196122,
            .64191653,
            1.1827365,
            -4.3573762
        )
        _V = variance(rnormal(length(_b), length(_b), 0, 1)) * 0.4
        st_eclear()
        st_matrix("`b'", _b)
        st_matrix("`V'", _V)

        alpha    = 0.05
        l_vec    = _honestBasis(5, 13)
        n        = length(_b)
        sel      = 13::n
        stdError = sqrt(l_vec' * _V[sel, sel] * l_vec)

        lb = l_vec' * _b[sel]' - invnormal(1-alpha/2)*stdError
        ub = l_vec' * _b[sel]' + invnormal(1-alpha/2)*stdError

    }
    matrix rownames `V' = `:colnames `V''
    ereturn post `b' `V'
end
test

matrix list e(b)
mata rows(st_matrix("e(V)")), cols(st_matrix("e(V)"))
matrix l_vec=(0\0\0\0\1\0\0\0\0\0\0\0\0)
honestdid, l_vec(l_vec) pre(1/12) post(13/25) mvec(0(0.025)0.25)
mata `s(HonestEventStudy)'.CI
mata `s(HonestEventStudy)'.betahat
mata `s(HonestEventStudy)'.numPrePeriods
mata `s(HonestEventStudy)'.numPostPeriods
mata `s(HonestEventStudy)'.options.alpha
mata `s(HonestEventStudy)'.options.l_vec
mata `s(HonestEventStudy)'.options.Mvec
mata `s(HonestEventStudy)'.options.grid_lb
mata `s(HonestEventStudy)'.options.grid_ub
honestdid, cached coefplot `plotopts'
graph export test.pdf, replace

* honestdid, l_vec(l_vec) pre(1/12) post(13/25) mvec(0 0.25)
* local plotopts xtitle(Mbar) ytitle(95% Robust CI)
* honestdid, cached coefplot `plotopts'
* graph export test.pdf, replace
