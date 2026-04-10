cap mata mata drop HonestDiDPLL()
cap mata mata drop _honestPLLAppendReplace()
cap mata mata drop _honestPLLRun()
cap mata mata drop _honestPLLFinish()
cap mata mata drop _honestPLLFinishMultiK()
cap mata mata drop _honestPLLGetKvecVec()
cap mata mata drop _honestPLLSave()
cap mata mata drop _honestPLLLoad()

mata
// Single l_vec: select M-values subset and compute sensitivity.
// Multi-l_vec: pass kvec_vec = vector of l_vec indices for this worker's rows.
//   All Mvec values are computed for each l_vec; l_vec_names encodes kvec order for assembly.
void HonestDiDPLL(struct HonestEventStudy scalar EventStudy, real vector mindex, | real vector kvec_vec)
{
    real scalar ki, k
    real vector orig_Mvec
    real matrix orig_l_vec
    struct _honestOptions scalar opts_k
    struct _honestResults colvector R_k

    if ( args() < 3 ) {
        EventStudy.options.Mvec = EventStudy.options.Mvec[mindex]
        EventStudy.Results      = HonestSensitivityResults(EventStudy, EventStudy.options)
        return
    }

    // Multi-l_vec mode: mindex ignored; compute full Mvec for each kvec in kvec_vec.
    orig_Mvec  = EventStudy.options.Mvec
    orig_l_vec = EventStudy.options.l_vec

    EventStudy.Results             = J(0, 1, _honestResults())
    EventStudy.options.Mvec        = J(1, 0, .)
    EventStudy.options.l_vec       = J(rows(orig_l_vec), 0, .)
    EventStudy.options.l_vec_names = J(1, length(kvec_vec), "")

    for (ki = 1; ki <= length(kvec_vec); ki++) {
        k            = kvec_vec[ki]
        opts_k       = EventStudy.options
        opts_k.l_vec = orig_l_vec[., k]
        opts_k.Mvec  = orig_Mvec
        R_k          = HonestSensitivityResults(EventStudy, opts_k)
        EventStudy.Results                 = EventStudy.Results \ R_k
        EventStudy.options.Mvec            = EventStudy.options.Mvec, orig_Mvec
        EventStudy.options.l_vec           = EventStudy.options.l_vec, orig_l_vec[., k]
        EventStudy.options.l_vec_names[ki] = strofreal(k)
    }
}

// Build kvecVec contribution from one loaded parfile: J(L,1,k) for each k in l_vec_names.
real colvector function _honestPLLGetKvecVec(struct HonestEventStudy scalar pf, real scalar L)
{
    real scalar ki
    real colvector kvecVec
    kvecVec = J(0, 1, .)
    for (ki = 1; ki <= length(pf.options.l_vec_names); ki++) {
        kvecVec = kvecVec \ J(L, 1, strtoreal(pf.options.l_vec_names[ki]))
    }
    return(kvecVec)
}

void function _honestPLLAppendReplace(struct HonestEventStudy scalar EventStudy, struct HonestEventStudy scalar AppendReplace)
{
    real vector Mvec
    Mvec = EventStudy.options.Mvec, AppendReplace.options.Mvec

    EventStudy.CI                = EventStudy.CI \ AppendReplace.CI
    EventStudy.betahat           = AppendReplace.betahat
    EventStudy.sigma             = AppendReplace.sigma
    EventStudy.prePeriodIndices  = AppendReplace.prePeriodIndices
    EventStudy.postPeriodIndices = AppendReplace.postPeriodIndices
    EventStudy.open              = EventStudy.open \ AppendReplace.open
    EventStudy.options           = AppendReplace.options
    EventStudy.OG                = AppendReplace.OG
    EventStudy.Results           = EventStudy.Results \ AppendReplace.Results
    EventStudy.options.Mvec      = Mvec
}

real scalar function _honestPLLRun(string scalar fname, struct HonestEventStudy scalar EventStudy, real scalar ncores)
{
    return(_stata(sprintf(`"HonestParallel `"%s"' %g %g"', fname, length(EventStudy.options.Mvec), ncores)))
}

void function _honestPLLFinish(struct HonestEventStudy scalar EventStudy)
{
    EventStudy.OG   = HonestOriginalCS(EventStudy, EventStudy.options)
    EventStudy.CI   = _honestSensitivityCIMatrix(EventStudy.Results, EventStudy.OG)
    EventStudy.open = _honestSensitivityCIOpen(EventStudy.Results, EventStudy.OG)
}

void function _honestPLLSave(string scalar fname, struct HonestEventStudy scalar EventStudy)
{
    colvector C
    real scalar m, fh, ki
    fh = fopen(fname, "rw")
    C = bufio()

    fbufput(C, fh, "%4bu", rows(EventStudy.CI))
    fbufput(C, fh, "%4bu", cols(EventStudy.CI))
    fbufput(C, fh, "%8z",  EventStudy.CI)
    fbufput(C, fh, "%4bu", length(EventStudy.betahat))
    fbufput(C, fh, "%8z",  EventStudy.betahat)
    fbufput(C, fh, "%4bu", rows(EventStudy.sigma))
    fbufput(C, fh, "%4bu", cols(EventStudy.sigma))
    fbufput(C, fh, "%8z",  EventStudy.sigma)
    fbufput(C, fh, "%8z",  EventStudy.numPrePeriods)
    fbufput(C, fh, "%8z",  EventStudy.numPostPeriods)
    fbufput(C, fh, "%4bu", length(EventStudy.prePeriodIndices))
    fbufput(C, fh, "%8z",  EventStudy.prePeriodIndices)
    fbufput(C, fh, "%4bu", length(EventStudy.postPeriodIndices))
    fbufput(C, fh, "%8z",  EventStudy.postPeriodIndices)
    fbufput(C, fh, "%4bu", length(EventStudy.open))
    fbufput(C, fh, "%8z",  EventStudy.open)

    fbufput(C, fh, "%8z",  EventStudy.options.alpha)
    fbufput(C, fh, "%4bu", length(EventStudy.options.l_vec))
    fbufput(C, fh, "%8z",  EventStudy.options.l_vec)
    fbufput(C, fh, "%4bu", length(EventStudy.options.l_vec_names))
    for (ki = 1; ki <= length(EventStudy.options.l_vec_names); ki++) {
        fbufput(C, fh, "%4bu", strlen(EventStudy.options.l_vec_names[ki]))
        fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.options.l_vec_names[ki])), EventStudy.options.l_vec_names[ki])
    }
    fbufput(C, fh, "%4bu", length(EventStudy.options.Mvec))
    fbufput(C, fh, "%8z",  EventStudy.options.Mvec)
    fbufput(C, fh, "%8z",  EventStudy.options.rm)
    fbufput(C, fh, "%8z",  EventStudy.options.omit)
    fbufput(C, fh, "%8z",  EventStudy.options.debug)
    fbufput(C, fh, "%8z",  EventStudy.options.grid_lb)
    fbufput(C, fh, "%8z",  EventStudy.options.grid_ub)
    fbufput(C, fh, "%8z",  EventStudy.options.gridPoints)
    fbufput(C, fh, "%4bu", strlen(EventStudy.options.Delta))
    fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.options.Delta)), EventStudy.options.Delta)
    fbufput(C, fh, "%4bu", strlen(EventStudy.options.method))
    fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.options.method)), EventStudy.options.method)
    fbufput(C, fh, "%4bu", strlen(EventStudy.options.relativeMagnitudes))
    fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.options.relativeMagnitudes)), EventStudy.options.relativeMagnitudes)

    fbufput(C, fh, "%4bu", length(EventStudy.Results))
    for(m = 1; m <= length(EventStudy.Results); m++) {
        fbufput(C, fh, "%8z", EventStudy.Results[m].lb)
        fbufput(C, fh, "%8z", EventStudy.Results[m].ub)
        fbufput(C, fh, "%8z", EventStudy.Results[m].M)
        fbufput(C, fh, "%8z", EventStudy.Results[m].lopen)
        fbufput(C, fh, "%8z", EventStudy.Results[m].uopen)
        fbufput(C, fh, "%4bu", strlen(EventStudy.Results[m].method))
        fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.Results[m].method)), EventStudy.Results[m].method)
        fbufput(C, fh, "%4bu", strlen(EventStudy.Results[m].Delta))
        fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.Results[m].Delta)),  EventStudy.Results[m].Delta)
    }

    fbufput(C, fh, "%4bu", rows(EventStudy.OG))
    if ( rows(EventStudy.OG) > 0 ) {
        for(m = 1; m <= rows(EventStudy.OG); m++) {
            fbufput(C, fh, "%8z", EventStudy.OG[m].lb)
            fbufput(C, fh, "%8z", EventStudy.OG[m].ub)
            fbufput(C, fh, "%8z", EventStudy.OG[m].M)
            fbufput(C, fh, "%8z", EventStudy.OG[m].lopen)
            fbufput(C, fh, "%8z", EventStudy.OG[m].uopen)
            fbufput(C, fh, "%4bu", strlen(EventStudy.OG[m].method))
            fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.OG[m].method)), EventStudy.OG[m].method)
            fbufput(C, fh, "%4bu", strlen(EventStudy.OG[m].Delta))
            fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.OG[m].Delta)),  EventStudy.OG[m].Delta)
        }
    }

    fclose(fh)
}

struct HonestEventStudy scalar function _honestPLLLoad(string scalar fname)
{
    struct HonestEventStudy scalar results

    colvector C
    real scalar m, fh, rows, cols, k_lvnames, k_og, ki

    results         = HonestEventStudy()
    results.options = _honestOptions()
    fh              = fopen(fname, "r")
    C               = bufio()

    rows                       = fbufget(C, fh, "%4bu")
    cols                       = fbufget(C, fh, "%4bu")
    results.CI                 = fbufget(C, fh, "%8z", rows, cols)
    cols                       = fbufget(C, fh, "%4bu")
    results.betahat            = fbufget(C, fh, "%8z", 1, cols)
    rows                       = fbufget(C, fh, "%4bu")
    cols                       = fbufget(C, fh, "%4bu")
    results.sigma              = fbufget(C, fh, "%8z", rows, cols)
    results.numPrePeriods      = fbufget(C, fh, "%8z", 1)
    results.numPostPeriods     = fbufget(C, fh, "%8z", 1)
    cols                       = fbufget(C, fh, "%4bu")
    results.prePeriodIndices   = fbufget(C, fh, "%8z", 1, cols)
    cols                       = fbufget(C, fh, "%4bu")
    results.postPeriodIndices  = fbufget(C, fh, "%8z", 1, cols)
    rows                       = fbufget(C, fh, "%4bu")
    results.open               = fbufget(C, fh, "%8z", rows, 1)

    results.options.alpha      = fbufget(C, fh, "%8z")
    rows                       = fbufget(C, fh, "%4bu")
    results.options.l_vec      = fbufget(C, fh, "%8z", rows, 1)
    k_lvnames                  = fbufget(C, fh, "%4bu")
    if ( k_lvnames > 0 ) {
        results.options.l_vec       = colshape(results.options.l_vec, k_lvnames)
        results.open                = colshape(results.open, 2)
        results.options.l_vec_names = J(1, k_lvnames, "")
        for (ki = 1; ki <= k_lvnames; ki++) {
            results.options.l_vec_names[ki] = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
        }
    }
    else {
        results.options.l_vec_names = J(1, 0, "")
    }
    cols                       = fbufget(C, fh, "%4bu")
    results.options.Mvec       = fbufget(C, fh, "%8z", 1, cols)
    results.options.rm         = fbufget(C, fh, "%8z")
    results.options.omit       = fbufget(C, fh, "%8z")
    results.options.debug      = fbufget(C, fh, "%8z")
    results.options.grid_lb    = fbufget(C, fh, "%8z")
    results.options.grid_ub    = fbufget(C, fh, "%8z")
    results.options.gridPoints = fbufget(C, fh, "%8z")

    results.options.Delta              = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
    results.options.method             = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
    results.options.relativeMagnitudes = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))

    rows = fbufget(C, fh, "%4bu")
    results.Results = J(rows, 1, _honestResults())
    for(m = 1; m <= rows; m++) {
        results.Results[m].lb     = fbufget(C, fh, "%8z")
        results.Results[m].ub     = fbufget(C, fh, "%8z")
        results.Results[m].M      = fbufget(C, fh, "%8z")
        results.Results[m].lopen  = fbufget(C, fh, "%8z")
        results.Results[m].uopen  = fbufget(C, fh, "%8z")
        results.Results[m].method = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
        results.Results[m].Delta  = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
    }

    k_og       = fbufget(C, fh, "%4bu")
    results.OG = J(k_og, 1, _honestResults())
    for(m = 1; m <= k_og; m++) {
        results.OG[m].lb     = fbufget(C, fh, "%8z")
        results.OG[m].ub     = fbufget(C, fh, "%8z")
        results.OG[m].M      = fbufget(C, fh, "%8z")
        results.OG[m].lopen  = fbufget(C, fh, "%8z")
        results.OG[m].uopen  = fbufget(C, fh, "%8z")
        results.OG[m].method = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
        results.OG[m].Delta  = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
    }

    fclose(fh)

    return(results)
}

// Build 4-column CI for multi-l_vec after parallel assembly.
// kvecVec: J(L,1,k) blocks in Results order (one block per l_vec batch).
// orig_l_vec/orig_l_vec_names: the original K-column l_vec matrix and user-supplied names.
// L: number of M values (mveclen).
void function _honestPLLFinishMultiK(
    struct HonestEventStudy scalar EventStudy,
    real colvector kvecVec,
    real matrix orig_l_vec,
    string rowvector orig_l_vec_names,
    real scalar L)
{
    real scalar k, i, result_idx, nresult_batches
    struct _honestOptions scalar opts_k
    struct _honestResults colvector R_k
    struct _honestResults scalar OG_k
    real matrix CI_k
    real colvector open_k

    nresult_batches = rows(kvecVec) / L

    EventStudy.CI   = J(0, 4, .)
    EventStudy.open = J(0, 2, .)
    EventStudy.OG   = J(0, 1, _honestResults())

    // Restore original l_vec so HonestOriginalCS uses correct data
    EventStudy.options.l_vec       = orig_l_vec
    EventStudy.options.l_vec_names = orig_l_vec_names

    result_idx = 1
    for (i = 1; i <= nresult_batches; i++) {
        k            = kvecVec[result_idx]
        R_k          = EventStudy.Results[result_idx .. result_idx+L-1]
        opts_k       = EventStudy.options
        opts_k.l_vec = orig_l_vec[., k]
        OG_k         = HonestOriginalCS(EventStudy, opts_k)
        CI_k         = _honestSensitivityCIMatrix(R_k, OG_k)
        open_k       = colshape(_honestSensitivityCIOpen(R_k, OG_k), 1)
        EventStudy.CI   = EventStudy.CI   \ (J(rows(CI_k),   1, k), CI_k)
        EventStudy.open = EventStudy.open \ (J(rows(open_k), 1, k), open_k)
        EventStudy.OG   = EventStudy.OG   \ OG_k
        result_idx      = result_idx + L
    }
}
end
