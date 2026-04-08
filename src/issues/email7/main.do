* forvalues i = 1(1)`numpost'{
*     mat lcomb[1, `i'] = 1
*     honestdid, b(B) vcov(V) pre(1/19) post(20/40) l_vec(lcomb) mvec(1) delta(rm)
*     mata: st_matrix("StataMatrix", `s(HonestEventStudy)'.CI)
*     mat res[`numpre' + `i', 5] = StataMatrix[2, 3]
*     mat res[`numpre' + `i', 6] = StataMatrix[2, 2]
*     mat lcomb[1, `i'] = 0
* }


* tempname beta sigma
* mata {
*     st_matrix(st_local("beta"),  _honestExampleBCBeta())
*     st_matrix(st_local("sigma"), _honestExampleBCSigma())
* }
* local opts mvec(0) gridPoints(1000) grid_lb(0.1) grid_ub(0.3)
* honestdid, numpre(4) b(`beta') vcov(`sigma') `opts'



tempname beta sigma
mata {
    st_matrix(st_local("beta"),  100 * _honestExampleBCBeta())
    st_matrix(st_local("sigma"), 100^2 * _honestExampleBCSigma())
}
local lower = 100 * 0.13353
local upper = 100 * 0.25839
local grid  = 100000
local step  = (`upper' - `lower')
local opts gridPoints(`grid') grid_lb(`=`lower'-`step'/10') grid_ub(`=`upper'+`step'/10')
honestdid, numpre(4) b(`beta') vcov(`sigma') mvec(0) delta(rm) `opts' debug parallel(0) alpha(0.05) method(Conditional)
mata `s(HonestEventStudy)'.CI
mata (`s(HonestEventStudy)'.CI[1, 2] - `s(HonestEventStudy)'.CI[2, 2])
mata (`s(HonestEventStudy)'.CI[1, 3] - `s(HonestEventStudy)'.CI[2, 3])

mata b  = 100 * _honestExampleBCBeta()[5]
mata se = 100 * sqrt(_honestExampleBCSigma()[5, 5])
mata a  = 0.05
mata b + se * invnormal(1-a/2)
mata b - se * invnormal(1-a/2)
mata b + se * invnormal(1-((a - a/10)/(1 - a/10))/2)
mata b - se * invnormal(1-((a - a/10)/(1 - a/10))/2)
mata b + se * invnormal(1-((a - a/10)/(1 - a/10)))
mata b - se * invnormal(1-((a - a/10)/(1 - a/10)))
mata ((a - a/10)/(1 - a/10))
mata _honestARPGeneralizedNorm(((a - a/10)/(1 - a/10))/2, ., ., b, se)
mata _honestARPGeneralizedNorm(1-((a - a/10)/(1 - a/10)), ., ., b, se)

disp `step'
disp 2 * `step'/1000
mata `s(HonestEventStudy)'.CI
mata `s(HonestEventStudy)'.Results[1].lb
mata `s(HonestEventStudy)'.Results[1].ub
mata `s(HonestEventStudy)'.options.gridPoints






tempname beta sigma
mata {
    st_matrix(st_local("beta"),  _honestExampleBCBeta())
    st_matrix(st_local("sigma"), _honestExampleBCSigma())
}
mat B = `beta'
mat V = `sigma'


parallel setclusters 4
capture program drop honestes
program honestes, rclass
    syntax anything(everything), es_i(int) es_numpost(int) [*]
    tempname StataMatrix
    mat lcomb = J(1, `es_numpost', 0)
    mat lcomb[1, `es_i'] = 1
    honestdid `anything', `options'
    mata: st_matrix("`StataMatrix'", `s(HonestEventStudy)'.CI)
    mat lcomb[1, `es_i'] = 0
    return matrix StataMatrix = `StataMatrix'
end
parallel honestes

mata {
    string scalar function honestblob(real vector v)
    {
        colvector C
        string scalar buf
        C = bufio()
        buf = (4 + 4 * length(v)) * char(0);
        bufput(C, buf, 0, "%4bu", length(v));
        bufput(C, buf, 4, "%4bu", v);
        return(buf)
    }
}
end

parallel setclusters 4, f
local numpre = 4
local numpost = 4
mat lcomb = J(1, `numpost', 0)
frame create honestes
frame honestes: set obs `numpost'
frame honestes: gen strL StataMatrix = ""
capture program drop honestes
program honestes
    tempname beta sigma
    mata {
        st_matrix(st_local("beta"),  _honestExampleBCBeta())
        st_matrix(st_local("sigma"), _honestExampleBCSigma())
    }
    mat B = `beta'
    mat V = `sigma'
    local opts b(B) vcov(V) pre(1/4) post(5/8) mvec(1) delta(rm) parallel(0)
    mat lcomb = J(1, 4, 0)
    forvalues i = 1 / `=_N'{
        mat lcomb[1, `i'] = 1
        honestdid, l_vec(lcomb) `opts'
        mata {
            out = `s(HonestEventStudy)'.CI[rows(`s(HonestEventStudy)'.CI), 2..3]
            C = bufio()
            buf = (4 + 4 * length(out)) * char(0);
            bufput(C, buf, 0, "%4bu", length(out));
            bufput(C, buf, 4, "%4bu", out);
            (void) st_sstore(`i', "StataMatrix", buf)
        }
        mat lcomb[1, `i'] = 0
    }
end
frame honestes: parallel, prog(honestes): honestes
frame honestes: format %20s StataMatrix
frame honestes: mata C = bufio()
frame honestes: mata bufget(C, st_sdata(1, "StataMatrix"), 0, "%4bu")
frame honestes: mata bufget(C, st_sdata(1, "StataMatrix"), 4, "%4bu", 1, 2)

parallel printlog 1
set rmsg on












mat res   = J(`numpre' + `numpost', 6, .)
mat res[`numpre' + `i', 5] = StataMatrix[2, 3]
mat res[`numpre' + `i', 6] = StataMatrix[2, 2]
