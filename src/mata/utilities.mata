cap mata mata drop _honestOptions()
cap mata mata drop _honestResults()
cap mata mata drop _honestHybridList()

cap mata mata drop _honestSensitivityFLCIMatrix()
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
cap mata mata drop _honestSummary()
cap mata mata drop _honestVLoVUpFN()
cap mata mata drop _honestLeeCFN()
cap mata mata drop _honestSelectionMat()

mata
struct _honestOptions {
    real scalar alpha
    real vector l_vec
    real vector Mvec
    real scalar rm
    string scalar method
    string scalar relativeMagnitudes
}

struct _honestResults {
    real scalar lb
    real scalar ub
    string scalar method
    string scalar Delta
    real scalar M
}

struct _honestHybridList {
    real scalar hybrid_kappa
    real scalar flci_halflength
    real scalar lf_cv
    real colvector flci_l
    real colvector vbar
    real colvector dbar
}


real matrix function _honestSensitivityFLCIMatrix(
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
end
