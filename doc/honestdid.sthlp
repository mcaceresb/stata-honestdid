{smcl}
{* *! version 0.4.2 29Aug2022}{...}
{viewerdialog honestdid "dialog honestdid"}{...}
{vieweralsosee "[R] honestdid" "mansection R honestdid"}{...}
{viewerjumpto "Syntax" "honestdid##syntax"}{...}
{viewerjumpto "Description" "honestdid##description"}{...}
{viewerjumpto "Options" "honestdid##options"}{...}
{viewerjumpto "Examples" "honestdid##examples"}{...}
{title:Title}

{p2colset 5 18 18 2}{...}
{p2col :{cmd:honestdid} {hline 2}}Stata implementation of the HonestDiD R package{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
Stata version of the HonestDiD R package, which implements robust inference for difference-in-differences and event study design methods developed in {browse "https://asheshrambachan.github.io/assets/files/hpt-draft.pdf":Rambachan and Roth (2021)}. See the {help honestdid##examples:examples} section below for an illustration of how users can use this package to conduct sensitivity analysis on the parallel trends assumption in difference-in-differences and event study designs.

{p 8 15 2}
{cmd:honestdid}
{cmd:,}
[{it:{help honestdid##table_options:options}} {it:{help coefplot:coefplot_options}}]

{pstd}
Typically at least one of {opt reference()} or {opt pre()} and {opt post()} are required options.

{synoptset 27 tabbed}{...}
{marker table_options}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{opth reference:periodindex(int)}} reference period; assumed omitted from vector (required or specify pre()/post()){p_end}
{synopt :{opth pre:periodindex(numlist)}} pre-period indices (required or specify reference()){p_end}
{synopt :{opth post:periodindex(numlist)}} post-period indices (required or specify reference()){p_end}
{synopt :{opt no:relmag}} (do not) use relative magnitudes (rel. mag. is default){p_end}
{synopt :{opt b(str)}} name of coefficient matrix; default is e(b){p_end}
{synopt :{opt vcov(str)}} name of vcov matrix; default is e(V){p_end}
{synopt :{opt l_vec(str)}} name of vector with parameters of interest (default is first period post event){p_end}
{synopt :{opt mvec(str)}} either name of vector with M values or number list of M values{p_end}
{synopt :{opth grid_lb(real)}} lower bound for grid search (ignored with FLCI); default is selected internally based on estimates {p_end}
{synopt :{opth grid_ub(real)}} upper bound for grid search (ignored with FLCI); default is selected internally based on estimates {p_end}
{synopt :{opt gridPoints(str)}} number of grid points for search; default 1000 {p_end}
{synopt :{opth alpha(real)}} 1 - confidence level; default 0.05{p_end}
{synopt :{opt method(str)}} C-LF (default unless norelmag is specified), FLCI (overrides relmag option; default with norelmag), Conditional, C-F  {p_end}
{synopt :{opt mata:save(str)}} save resulting mata object (default: HonestEventStudy){p_end}
{synopt :{opt coefplot}} coefficient plot{p_end}
{synopt :{opt cached}} use cached results for coefficient plot{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:honestdid} See the {browse "https://github.com/asheshrambachan/HonestDiD/blob/master/doc/HonestDiD_Example.pdf":HonestDiD Vignette}. While not all the functionality is mirrored, the base functionality has been implemented (see examples below).

{marker example}{...}
{title:Example 1: Benzarti and Carloni (2019)}

{pstd}
Construct robust confidence intervals for DeltaRM(Mbar). In the R Vignette,
the test uses option bound set to "deviation from linear trend"; here only
the analogue to "deviation from parallel trends" has been implemented.

{phang2}{cmd:. tempname beta sigma                                            }{p_end}
{phang2}{cmd:. mata {c -(}                                                    }{p_end}
{phang2}{cmd:      st_matrix(st_local("beta"),  _honestExampleBCBeta())       }{p_end}
{phang2}{cmd:      st_matrix(st_local("sigma"), _honestExampleBCSigma())      }{p_end}
{phang2}{cmd:  {c )-}                                                         }{p_end}
{phang2}{cmd:. local opts mvec(0(0.5)2) gridPoints(100) grid_lb(-1) grid_ub(1)}{p_end}
{phang2}{cmd:. honestdid, reference(4) b(`beta') vcov(`sigma') `opts'         }{p_end}

{pstd}
The results are printed to the Stata console and saved in a mata object. The user
can specify the name of the mata object via the {opt mata()} option. The default
name is HonestEventStudy and the object's name is saved in {cmd:s(HonestEventStudy)}.
In addition to the CI, all the inputs and options are saved:

{phang2}{cmd:. mata `s(HonestEventStudy)'.CI                }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.betahat           }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.sigma             }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.referencePeriod   }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.prePeriodIndices  }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.postPeriodIndices }{p_end}

{phang2}{cmd:. mata `s(HonestEventStudy)'.options.alpha     }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.l_vec     }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.Mvec      }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.rm        }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.method    }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.Delta     }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.grid_lb   }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.grid_ub   }{p_end}
{phang2}{cmd:. mata `s(HonestEventStudy)'.options.gridPoints}{p_end}

{pstd}
For ease of use, the package also provides a way to plot the CIs
using the {cmd:coefplot} package. This can be done when the CIs
are computed or using the results cached in mata.

{phang2}{cmd:. honestdid, coefplot cached                                   }{p_end}
{phang2}{cmd:. honestdid, coefplot cached xtitle(Mbar) ytitle(95% Robust CI)}{p_end}
{phang2}{cmd:. graph export coefplot.pdf, replace                           }{p_end}

{title:Example 2: Lovenheim and Willen (2019)}

{phang2}{cmd:. use test/LWdata_RawData.dta, clear      }{p_end}
{phang2}{cmd:. mata stata(_honestExampleLWCall())      }{p_end}
{phang2}{cmd:. honestdid, pre(1/9) post(10/32) coefplot}{p_end}

{pstd}
This required the user package {cmd:reghdfe}. First, you can see
that we did not need to specify the coefficient vector or the
variance-covariance vector, and instead the function took the stored
results {cmd:e(b)} and {cmd:e(V)} from {cmd:reghdfe} to do the
computations. Further, the coefficient plot was created using the
{cmd:coefplot} package.

{pstd}
Now to mirror the Vignette:

{phang2}{cmd:. matrix b = 100 * e(b)                                     }{p_end}
{phang2}{cmd:. matrix V = 100^2 * e(V)                                   }{p_end}
{phang2}{cmd:. mata st_matrix("l_vec", _honestBasis(15 - (-2), 23))      }{p_end}
{phang2}{cmd:. local opts norelmag mvec(0(0.005)0.04) l_vec(l_vec)       }{p_end}
{phang2}{cmd:. local plot coefplot xtitle(Mbar) ytitle(95% Robust CI)    }{p_end}
{phang2}{cmd:. honestdid, pre(1/9) post(10/32) b(b) vcov(V) `opts' `plot'}{p_end}
{phang2}{cmd:. graph export coefplot.pdf, replace                        }{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:honestdid} stores the following in {cmd:s()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Macros}{p_end}
{synopt:{cmd:s(HonestEventStudy)}}{cmd:honestdid} results mata object{p_end}

{marker mata}{...}
{pstd}
The following data are available in {cmd:s(HonestEventStudy)} (default name: HonestEventStudy):

        real matrix CI
            Matrix with M in first column, original CI and robust CIs in columns 2, 3

        real vector betahat
            coefficient vector (pre and post period coefficients only)

        real matrix sigma
            vcov matrix (pre and post period vcov matrix only)

        real vector prePeriodIndices
            pre period indices of coef vector

        real vector postPeriodIndices
            post period indices of coef vector

        struct _honestOptions scalar options
            structure with options used in internal computations

The following data are available in {cmd:s(HonestEventStudy).options} (default name: HonestEventStudy.options):

        real scalar alpha
            1 - confidence level
            
        real vector l_vec
            vector with parameters of interest (default is first period post event)

        real vector Mvec
            vector with M values

        real scalar rm
            whether relative magnitudes were requested

        real scalar grid_lb
            lower bound for grid search (ignored with FLCI)

        real scalar grid_ub
            upper bound for grid search (ignored with FLCI)

        real scalar gridPoints
            number of grid points

        string scalar Delta
            type of Delta (DeltaSD or DeltaRM)

        string scalar method
            method used (C-LF, FLCI, Conditional, C-F)

{marker references}{...}
{title:References}

{pstd}
See the {browse "https://github.com/asheshrambachan/HonestDiD/blob/master/doc/HonestDiD_Example.pdf":HonestDiD Vignette} and the paper by {browse "https://asheshrambachan.github.io/assets/files/hpt-draft.pdf":Rambachan and Roth (2021)}.

