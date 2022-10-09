cap mata mata drop _honestOptions()
cap mata mata drop _honestResults()
cap mata mata drop _honestHybridList()

cap mata mata drop _honestSensitivityCIMatrix()
cap mata mata drop _honestSensitivityCIOpen()
cap mata mata drop _honestBasis()
cap mata mata drop _honestLinspace()
cap mata mata drop _honestInverseIndex()
cap mata mata drop _honestHelper()
cap mata mata drop _honestwToLFn()
cap mata mata drop _honestlToWFn()
cap mata mata drop _honestQuantiles()
cap mata mata drop _honestMatrixPow()
cap mata mata drop _honestPrimes()
cap mata mata drop _honestPrimesn()
cap mata mata drop _honestHalton()
cap mata mata drop _honestHaltonVector()
cap mata mata drop _honestHaltonJumble()
cap mata mata drop _honestHaltonJumbleVector()
cap mata mata drop _honestSummary()
cap mata mata drop _honestVLoVUpFN()
cap mata mata drop _honestLeeCFN()
cap mata mata drop _honestSelectionMat()
cap mata mata drop _honestTruncNormInv()
cap mata mata drop _honestTruncNormCDF()
cap mata mata drop _honestExampleBCBeta()
cap mata mata drop _honestExampleBCSigma()
cap mata mata drop _honestExampleLWCall()

mata
struct _honestOptions {
    real scalar alpha
    real vector l_vec
    real vector Mvec
    real scalar rm
    real scalar omit
    real scalar debug
    real scalar grid_lb
    real scalar grid_ub
    real scalar gridPoints
    string scalar Delta
    string scalar method
    string scalar relativeMagnitudes
}

struct _honestResults {
    real scalar lb
    real scalar ub
    string scalar method
    string scalar Delta
    real scalar M
    real scalar lopen
    real scalar uopen
}

struct _honestHybridList {
    real scalar hybrid_kappa
    real scalar flci_halflength
    real scalar lf_cv
    real colvector flci_l
    real colvector vbar
    real colvector dbar
}


real matrix function _honestSensitivityCIMatrix(
    struct _honestResults colvector robustResults,
    struct _honestResults scalar originalResults,
    | real scalar rescaleFactor, real scalar maxM)
{
    if ( args() < 3 ) rescaleFactor = 1
    if ( args() < 4 ) maxM = .

    real colvector M, lb, ub, sel
    real scalar m

    M  = J(rows(robustResults), 1, .)
    lb = J(rows(robustResults), 1, .)
    ub = J(rows(robustResults), 1, .)
    for (m = 1; m <= length(M); m++) {
        M[m]  = robustResults[m].M
        lb[m] = robustResults[m].lb
        ub[m] = robustResults[m].ub
    }

    // Filter out observations above maxM (after rescaling)
    sel = selectindex(M :<= maxM)
    M   = M[sel]
    lb  = lb[sel]
    ub  = ub[sel]

    // Set M for OLS to be the min M in robust results minus the gap
    // between Ms in robust
    M  = (originalResults.M  \ M)  * rescaleFactor
    lb = (originalResults.lb \ lb) * rescaleFactor
    ub = (originalResults.ub \ ub) * rescaleFactor

    return(M, lb, ub)
}

real vector function _honestSensitivityCIOpen(
    struct _honestResults colvector robustResults,
    struct _honestResults scalar originalResults)
{
    real scalar m
    real vector open
    open = (originalResults.lopen + 2 * originalResults.uopen) \ J(rows(robustResults), 1, .)
    for (m = 1; m <= rows(robustResults); m++) {
        open[m+1] = robustResults[m].lopen + 2 * robustResults[m].uopen
    }
    return(open)
}

real colvector function _honestBasis(real scalar index, real scalar size)
{
    real colvector v
    v = J(size, 1, 0)
    v[index] = 1
    return(v)
}

real vector _honestLinspace(real scalar from, real scalar to, real scalar n) {
    real scalar step
    step = (to - from) / (n - 1)
    return(from :+ step :* (0..(n - 1)))
}

real matrix function _honestInverseIndex(real vector ix, real scalar n)
{
    real vector sel
    if (rows(ix) > cols(ix)) {
        sel = J(n, 1, 1)
        sel[ix] = J(length(ix), 1, 0)
    }
    else {
        sel = J(1, n, 1)
        sel[ix] = J(1, length(ix), 0)
    }
    return(selectindex(sel))
}

// Converts vector from w space to l space.
real vector function _honestwToLFn(real vector w) {
    real scalar numPrePeriods, col
    real matrix WtoLPreMat

    numPrePeriods = length(w)
    WtoLPreMat = I(numPrePeriods)
    if (numPrePeriods == 1) {
        WtoLPreMat = 1
    }
    else {
        for (col = 1; col < numPrePeriods; col++) {
            WtoLPreMat[col + 1, col] = -1
        }
    }
    return(WtoLPreMat * w)
}

// Converts vector from l space to w space
real vector function _honestlToWFn(real vector l_vec) {
    real scalar numPostPeriods, col
    real matrix lToWPostMat

    numPostPeriods = length(l_vec)
    lToWPostMat    = I(numPostPeriods)
    if (numPostPeriods == 1) {
        lToWPostMat = 1
    }
    else {
        for (col = 1; col < numPostPeriods; col++) {
            lToWPostMat[col + 1, col] = 1
        }
    }
    return(lToWPostMat * l_vec)
}

real scalar function _honestHelper(real scalar s, real colvector l_vec, real scalar numPostPeriods)
{
    return(abs((1..s) * l_vec[(numPostPeriods - s + 1)::numPostPeriods]))
}

real vector function _honestQuantiles(real vector x, real vector p) {
    real scalar n
    real colvector m, i, g, sel
    real matrix orderx
    n   = length(x)
    m   = 1 :- colshape(p, 1)
    i   = floor(n * colshape(p, 1) + m)
    sel = i \ i :+ 1
    g   = n * colshape(p, 1) + m - i
    if ( any(i :< 1) ) {
        sel[selectindex(sel :< 1)] = 1
    }
    if ( any(i :>= n) ) {
        sel[selectindex(sel :> n)] = n
    }
    orderx = rowshape(sort(colshape(x, 1), 1)[sel], 2)'
    return(((1 :- g) :* orderx[., 1] + g :* orderx[., 2]))
}

real matrix function _honestMatrixPow(real matrix X, real scalar pow) {
    real matrix C
    real vector L
    if ( rows(X) != cols(X) ) {
        errprintf("_honestMatrixPow() requires a square matrix\n")
        _error(198)
    }
    if ( (rank(X) != cols(X)) & (pow < 0) ) {
        errprintf("_honestMatrixPow() requires a full-rank matrix for inverses\n")
        _error(198)
    }
    eigensystem(X, C=., L=.)
    return(Re((C :* (L:^pow)) * C'))
}

real matrix function _honestHaltonJumble(real scalar n, real scalar d)
{
    real matrix r
    real colvector p
    real scalar j
    r = J(n, d, 0)
    p = _honestPrimes(d)
    for (j = 1; j <= d; j++) {
        r[., j] = _honestHaltonJumbleVector(n, p[j])
    }
    return(r)
}

real colvector function _honestHaltonJumbleVector(real scalar n, real scalar p)
{
    real colvector r
    real scalar i, k, f
    r = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        k = i
        f = 1
        while ( k > 0 ) {
            f    = f / p
            r[i] = r[i] + (f * mod(k, p))
            k    = floor(k / p)
        }
    }
    return(jumble(r))
}

real matrix function _honestHalton(real scalar n, real scalar d)
{
    real matrix r
    real colvector p
    real scalar j
    r = J(n, d, 0)
    p = _honestPrimes(d)
    for (j = 1; j <= d; j++) {
        r[., j] = _honestHaltonVector(n, p[j])
    }
    return(r)
}

real colvector function _honestHaltonVector(real scalar n, real scalar p)
{
    real colvector r
    real scalar i, k, f
    r = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        k = i
        f = 1
        while ( k > 0 ) {
            f    = f / p
            r[i] = r[i] + (f * mod(k, p))
            k    = floor(k / p)
        }
    }
    return(r)
}

real matrix function _honestPrimes(real scalar k)
{
    real scalar n, nk
    n  = 2
    nk = n / log(n)
    while (nk < (k + 1)) {
        ++n
        nk = n / log(n)
    }
    return(_honestPrimesn(n)[1::k])
}

real matrix function _honestPrimesn(real scalar n)
{
    real scalar start, m
    real colvector s, sel
    s = 2 :* (1::(trunc(n/2)-1)) :+ 1
    for (m = 3; m <= trunc(sqrt(n)); m++, m++) {
        if ( s[trunc((m - 3) / 2) + 1] ) {
            start  = trunc((m * m - 3) / 2) + 1
            sel    = m :* (0::trunc((length(s) - start) / m)) :+ start
            if ( start <= length(sel) ) {
                s[sel] = J(length(sel), 1, 0)
            }
        }
    }
    return((2 \ s[selectindex(s :> 0)]))
}

real matrix function _honestSummary(real matrix x)
{
    real scalar i
    real matrix qres
    real vector qpts
    qpts = (0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
    qres = J(length(qpts), cols(x), .)
    for (i = 1; i <= cols(x); i++) {
        qres[., i] = _honestQuantiles(x[., i], qpts)
    }
    return(qres)
}

real matrix function _honestSelectionMat(real vector selection, real scalar size, | real scalar srows)
{
    real matrix m
    if ( args() < 3 ) srows = 0
    if ( srows ) {
        m = J(length(selection), size, 0)
        m[., selection] = I(length(selection))
    }
    else {
        m = J(size, length(selection), 0)
        m[selection, .] = I(length(selection))
    }
    return(m)
}

real vector function _honestVLoVUpFN(real colvector eta, real matrix Sigma, real matrix A, real colvector b, real colvector z){
    real colvector c, objective, ACNegativeIndex, ACPositiveIndex
    real scalar VLo, VUp
    c = _honestLeeCFN(eta, Sigma)
    objective = (b - A * z) :/ (A * c)
    ACNegativeIndex = selectindex((A * c) :< 0)
    ACPositiveIndex = selectindex((A * c) :> 0)
    VLo = max(objective[ACNegativeIndex])
    VUp = min(objective[ACPositiveIndex])
    return((VLo, VUp))
}

real colvector function _honestLeeCFN(real colvector eta, real matrix Sigma)
{
    return((Sigma * eta) :/ (eta' * Sigma * eta))
}

real scalar function _honestTruncNormCDF(real scalar z,
                                         real scalar l,
                                         real scalar u,
                                         | real scalar mu,
                                         real scalar sd,
                                         real scalar n) {

    real vector grid, width, area, odd, even
    real scalar n1, n2, ub, lb, zn
    if ( args() < 4 ) mu = 0
    if ( args() < 5 ) sd = 1
    if ( args() < 6 ) n  = 1e3

    zn = (z - mu)/sd
    lb = missing(l)? -37: (l-mu)/sd
    ub = missing(u)?  37: (u-mu)/sd
    n1 = floor(n * ((zn - lb) / (ub - lb)))
    n2 = n - n1 + 1

    grid  = (_honestLinspace(lb, zn, n1), _honestLinspace(zn, ub, n2)[|2 \ n2|])
    width = grid[|2 \ n|] :- grid[|1 \ n-1|]
    area  = width :* (normalden(grid[|1 \ n-1|]) :+ normalden(grid[|2 \ n|])) :/ 2
    return(quadsum(area[|1 \ n1|]) :/ quadsum(area))
}

real scalar function _honestTruncNormInv(real scalar p,
                                         real scalar l,
                                         real scalar u,
                                         | real scalar mu,
                                         real scalar sd,
                                         real scalar n) {
    real vector grid, width, area
    real scalar ub, lb, index, pscale, a
    if ( args() < 4 ) mu = 0
    if ( args() < 5 ) sd = 1
    if ( args() < 6 ) n  = 1e3

    lb = missing(l)? -37: max((min(((l-mu)/sd,  37)), -37))
    ub = missing(u)?  37: min((max(((u-mu)/sd, -37)),  37))

    if (ub <= lb) {
        return(.)
    }
    else if ( lb >= 1 ) {
        grid = exp(_honestLinspace(log(lb), log(ub), n))
    }
    else if ( ub <= -1 ) {
        grid = -exp(_honestLinspace(log(-lb), log(-ub), n))
    }
    else {
        return(mu + sd * invnormal(p * (normal(ub) - normal(lb)) + normal(lb)))
    }

    width  = grid[|2 \ n|] :- grid[|1 \ n-1|]
    area   = quadrunningsum(width :* (normalden(grid[|1 \ n-1|]) :+ normalden(grid[|2 \ n|])) :/ 2)
    index  = 0
    pscale = p * area[n-1]

    while ( (++index < n) & (pscale > area[index]) ) {}
    a = (pscale - area[index-1]) / (area[index] - area[index-1])
    return(index < n? mu + sd * (a * grid[index] + (1 - a) * grid[index-1]): .)
}

real rowvector function _honestExampleBCBeta() {
    return((0.006696352,
            0.029345034,
           -0.006472972,
            0.073014989,
            0.195961118,
            0.312063903,
            0.239541546,
            0.126042500))
}

real matrix function _honestExampleBCSigma() {
    return(( 8.428358e-04,  4.768687e-04, 2.618051e-04, 0.0002354220, 0.0001676371, 0.0001128708, 1.992816e-05, -1.368265e-04) \
           ( 4.768687e-04,  6.425420e-04, 3.987425e-04, 0.0002435515, 0.0002201960, 0.0001804591, 3.843765e-05, -2.960422e-05) \
           ( 2.618051e-04,  3.987425e-04, 5.229950e-04, 0.0002117686, 0.0001840722, 0.0001458528, 7.005197e-05,  5.952995e-05) \
           ( 2.354220e-04,  2.435515e-04, 2.117686e-04, 0.0003089595, 0.0001197866, 0.0001334081, 1.016335e-04,  1.079052e-04) \
           ( 1.676371e-04,  2.201960e-04, 1.840722e-04, 0.0001197866, 0.0003599704, 0.0002478819, 1.749579e-04,  1.654257e-04) \
           ( 1.128708e-04,  1.804591e-04, 1.458528e-04, 0.0001334081, 0.0002478819, 0.0004263950, 2.171438e-04,  2.892748e-04) \
           ( 1.992816e-05,  3.843765e-05, 7.005197e-05, 0.0001016335, 0.0001749579, 0.0002171438, 4.886698e-04,  3.805322e-04) \
           (-1.368265e-04, -2.960422e-05, 5.952995e-05, 0.0001079052, 0.0001654257, 0.0002892748, 3.805322e-04,  7.617394e-04))
}

string scalar function _honestExampleLWCall() {
    return("reghdfe emp rtESV13 rtESV14 rtESV15 rtESV16 rtESV17 rtESV18 rtESV19 rtESV110 rtESV111 rtESV113 rtESV114 rtESV115 rtESV116 rtESV117 rtESV118 rtESV119 rtESV120 rtESV121 rtESV122 rtESV123 rtESV124 rtESV125 rtESV126 rtESV127 rtESV128 rtESV129 rtESV130 rtESV131 rtESV132 rtESV133 rtESV134 rtESV135 yearsfcor yearsflr aveitc fscontrol asian black hispanic other [fw = nobs], absorb(PUS_SURVEY_YEAR BIRTHSTATE PUS_SURVEY_YEAR#BIRTHYEAR) cluster(BIRTHSTATE) noconstant")
}
end
