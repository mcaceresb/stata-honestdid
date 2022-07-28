{smcl}
{* *! version 0.2.0 28Jul2022}{...}
{viewerdialog honestdid "dialog honestdid"}{...}
{vieweralsosee "[R] honestdid" "mansection R honestdid"}{...}
{viewerjumpto "Syntax" "honestdid##syntax"}{...}
{viewerjumpto "Description" "honestdid##description"}{...}
{viewerjumpto "Options" "honestdid##options"}{...}
{viewerjumpto "Examples" "honestdid##examples"}{...}
{title:Title}

{p2colset 5 14 14 2}{...}
{p2col :{cmd:honestdid} {hline 2}}HonestDiD{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
HonestDiD alpha release

{p 8 15 2}
{cmd:honestdid}
{cmd:,}
[{it:{help honestdid##table_options:options}} {it:{help coefplot:coefplot_options}}]

{synoptset 18 tabbed}{...}
{marker table_options}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{opt b(str)}} name of coefficient matrix; default is e(b)                   {p_end}
{synopt :{opt vcov(str)}} name of vcov matrix; default is e(V)                       {p_end}
{synopt :{opt l_vec(str)}} Vector with parameters of interest (default is first period post event){p_end}
{synopt :{opt mvec(str)}} Vector name with M values or list of M values{p_end}
{synopt :{opt alpha(real)}} 1 - confience level; default 0.05{p_end}
{synopt :{opt reference:periodindex(str)}} index for the reference period            {p_end}
{synopt :{opt pre:periodindex(str)}} pre-period indices                              {p_end}
{synopt :{opt post:periodindex(str)}} post-period indices                            {p_end}
{synopt :{opt method(str)}} FLCI (xx default), Conditional, C-F or C-LF{p_end}
{synopt :{opt mata:save(str)}} Save resulting mata object (default: HonestEventStudy){p_end}
{synopt :{opt coefplot}} Coefficient plot                                            {p_end}
{synopt :{opt cached}} Use cached results for coefficient plot                       {p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:honestdid} xx

{marker options}{...}
{title:Options}

{dlgtab:Command Options}

{phang}{opth b(str)} xx

{marker example}{...}
{title:Example 1: xx}

{phang2}{cmd:. xx}{p_end}

{title:Example 2: xx}

{phang2}{cmd:. xx}{p_end}

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

        real matrix FLCI
            Matrix with M in first clumn, original CI and FLCIs in columns 2, 3

        real vector betahat
            coefficient vector

        real matrix sigma
            vcov matrix

        real vector prePeriodIndices
            pre period indices of coef vector

        real vector postPeriodIndices
            post period indices of coef vector

        struct _honestResults colvector Results
            xx

        struct _honestResults scalar OG
            xx

{marker references}{...}
{title:References}

{pstd}
HonestDiD xx
