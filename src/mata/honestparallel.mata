mata
struct HonestEventStudy scalar HonestDiDPLL(struct HonestEventStudy scalar EventStudy, real vector mindex)
{
    EventStudy.options.Mvec = EventStudy.options.Mvec[mindex]
    EventStudy.Results      = HonestSensitivityResults(results, EventStudy.options)
}


void function _honestPLLAppendReplace(struct HonestEventStudy scalar EventStudy, struct HonestEventStudy scalar AppendReplace)
{
    // xx
}

void function _honestPLLFinish(struct HonestEventStudy scalar EventStudy)
{
    EventStudy.OG   = HonestOriginalCS(results, EventStudy.options)
    EventStudy.CI   = _honestSensitivityCIMatrix(EventStudy.Results, EventStudy.OG)
    EventStudy.open = _honestSensitivityCIOpen(EventStudy.Results, EventStudy.OG)
}

void function _honestPLLSave(string scalar fname, struct HonestEventStudy scalar EventStudy)
{
    // xx save results
}

struct HonestEventStudy scalar EventStudy function _honestPLLLoad(string scalar fname)
{
    // xx load results
}
end
