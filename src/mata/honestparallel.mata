cap mata mata drop HonestDiDPLL()
cap mata mata drop _honestPLLAppendReplace()
cap mata mata drop _honestPLLRun()
cap mata mata drop _honestPLLFinish()
cap mata mata drop _honestPLLSave()
cap mata mata drop _honestPLLLoad()

mata
void HonestDiDPLL(struct HonestEventStudy scalar EventStudy, real vector mindex)
{
    EventStudy.options.Mvec = EventStudy.options.Mvec[mindex]
    EventStudy.Results      = HonestSensitivityResults(EventStudy, EventStudy.options)
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
    real scalar m, fh
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
    fbufput(C, fh, "%4bu", length(EventStudy.prePeriodIndices))
    fbufput(C, fh, "%8z",  EventStudy.prePeriodIndices)
    fbufput(C, fh, "%4bu", length(EventStudy.postPeriodIndices))
    fbufput(C, fh, "%8z",  EventStudy.postPeriodIndices)
    fbufput(C, fh, "%4bu", length(EventStudy.open))
    fbufput(C, fh, "%8z",  EventStudy.open)

    fbufput(C, fh, "%8z",  EventStudy.options.alpha)
    fbufput(C, fh, "%4bu", length(EventStudy.options.l_vec))
    fbufput(C, fh, "%8z",  EventStudy.options.l_vec)
    fbufput(C, fh, "%4bu", length(EventStudy.options.Mvec))
    fbufput(C, fh, "%8z",  EventStudy.options.Mvec)
    fbufput(C, fh, "%8z",  EventStudy.options.rm)
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

    fbufput(C, fh, "%8z", EventStudy.OG.lb)
    fbufput(C, fh, "%8z", EventStudy.OG.ub)
    fbufput(C, fh, "%8z", EventStudy.OG.M)
    fbufput(C, fh, "%8z", EventStudy.OG.lopen)
    fbufput(C, fh, "%8z", EventStudy.OG.uopen)
    fbufput(C, fh, "%4bu", strlen(EventStudy.OG.method))
    fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.OG.method)), EventStudy.OG.method)
    fbufput(C, fh, "%4bu", strlen(EventStudy.OG.Delta))
    fbufput(C, fh, sprintf("%%%gs", strlen(EventStudy.OG.Delta)),  EventStudy.OG.Delta)

    fclose(fh)
}

struct HonestEventStudy scalar function _honestPLLLoad(string scalar fname)
{
    struct HonestEventStudy scalar results

    colvector C
    real scalar m, fh, rows, cols

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
    cols                       = fbufget(C, fh, "%4bu")
    results.prePeriodIndices   = fbufget(C, fh, "%8z", 1, cols)
    cols                       = fbufget(C, fh, "%4bu")
    results.postPeriodIndices  = fbufget(C, fh, "%8z", 1, cols)
    rows                       = fbufget(C, fh, "%4bu")
    results.open               = fbufget(C, fh, "%8z", rows, 1)

    results.options.alpha      = fbufget(C, fh, "%8z")
    rows                       = fbufget(C, fh, "%4bu")
    results.options.l_vec      = fbufget(C, fh, "%8z", rows, 1)
    cols                       = fbufget(C, fh, "%4bu")
    results.options.Mvec       = fbufget(C, fh, "%8z", 1, cols)
    results.options.rm         = fbufget(C, fh, "%8z")
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

    results.OG        = _honestResults()
    results.OG.lb     = fbufget(C, fh, "%8z")
    results.OG.ub     = fbufget(C, fh, "%8z")
    results.OG.M      = fbufget(C, fh, "%8z")
    results.OG.lopen  = fbufget(C, fh, "%8z")
    results.OG.uopen  = fbufget(C, fh, "%8z")
    results.OG.method = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))
    results.OG.Delta  = fbufget(C, fh, sprintf("%%%gs", fbufget(C, fh, "%4bu")))

    fclose(fh)

    return(results)
}
end
