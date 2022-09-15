HonestDiD
=========

The HonestDiD package implements the tools for robust inference and
sensitivity analysis for differences-in-differences and event study
designs developed in [Rambachan and Roth (2022)](https://asheshrambachan.github.io/assets/files/hpt-draft.pdf).

`version 0.4.5 30Aug2022` | [Background](#background) | [Installation](#package-installation) | [Examples](#example-usage-medicaid-expansions) | [Acknowledgements](#acknowledgements)

## Background

The robust inference approach in Rambachan and Roth formalizes the
intuition that pre-trends are informative about violations of parallel
trends. They provide a few different ways of formalizing what this
means.

**Bounds on relative magnitudes.** One way of formalizing this idea is
to say that the violations of parallel trends in the post-treatment
period cannot be much bigger than those in the pre-treatment period.
This can be formalized by imposing that the post-treatment violation of
parallel trends is no more than some constant $\bar{M}$ larger than the
maximum violation of parallel trends in the pre-treatment period. The
value of $\bar{M} = 1$, for instance, imposes that the post-treatment
violation of parallel trends is no longer than the worst pre-treatment
violation of parallel trends (between consecutive periods). Likewise,
setting $\bar{M} = 2$ implies that the post-treatment violation of
parallel trends is no more than twice that in the pre-treatment period.

**Smoothness restrictions.** A second way of formalizing this is to say
that the post-treatment violations of parallel trends cannot deviate too
much from a linear extrapolation of the pre-trend. In particular, we can
impose that the slope of the pre-trend can change by no more than *M*
across consecutive periods, as shown in the figure below for an example
with three periods.

![diagram-smoothness-restriction](src/assets/deltaSD.png)

Thus, imposing a smoothness restriction with $M = 0$ implies that the
counterfactual difference in trends is exactly linear, whereas larger
values of $M$ allow for more non-linearity.

**Other restrictions**. The Rambachan and Roth framework allows for a
variety of other restrictions on the differences in trends as well. For
example, the smoothness restrictions and relative magnitudes ideas can
be combined to impose that the non-linearity in the post-treatment
period is no more than $\bar{M}$ times larger than that in the
pre-treatment periods. The researcher can also impose monotonicity or
sign restrictions on the differences in trends as well.

**Robust confidence intervals**. Given restrictions of the type
described above, Rambachan and Roth provide methods for creating robust
confidence intervals that are guaranteed to include the true parameter
at least 95% of the time when the imposed restrictions on satisfied.
These confidence intervals account for the fact that there is estimation
error both in the treatment effects estimates and our estimates of the
pre-trends.

**Sensitivity analysis**. The approach described above naturally lends
itself to sensitivity analysis. That is, the researcher can report
confidence intervals under different assumptions about how bad the
post-treatment violation of parallel trends can be (e.g., different
values of $\bar{M}$ or $M$.) They can also report the "breakdown value"
of $\bar{M}$ (or $M$) for a particular conclusion---e.g. the largest
value of $\bar{M}$ for which the effect is still significant.

## Package installation

The package may be installed by using `net install`:

```stata
local github "https://raw.githubusercontent.com"
net install multe, from(`github'/mcaceresb/stata-honestdid/main) replace
```

## Example usage -- Medicaid expansions

As an illustration of the package, we will examine the effects of
Medicaid expansions on insurance coverage using publicly-available data
derived from the ACS. We first load the data and packages relevant for
the analysis.

```stata

* Install here coefplot, ftools, reghdfe
ssc install coefplot, replace
ssc install ftools,   replace
ssc install reghdfe,  replace

* Load data
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear
l in 1/5
```

```
     +--------------------------------------------+
     |  stfips   year       dins   yexp2        W |
     |--------------------------------------------|
  1. | alabama   2008   .6814122       .   613156 |
  2. | alabama   2009   .6580621       .   613156 |
  3. | alabama   2010   .6313651       .   613156 |
  4. | alabama   2011   .6563886       .   613156 |
  5. | alabama   2012   .6708115       .   613156 |
     +--------------------------------------------+
```

The data is a state-level panel with information on health insurance
coverage and Medicaid expansion. The variable `dins` shows the share of
low-income childless adults with health insurance in the state. The
variable `yexp2` gives the year that a state expanded Medicaid coverage
under the Affordable Care Act, and is missing if the state never
expanded.

### Estimate the baseline DiD

For simplicity, we will first focus on assessing sensitivity to
violations of parallel trends in a non-staggered DiD (see below
regarding methods for staggered timing). We therefore restrict the
sample to the years 2015 and earlier, and drop the small number of
states who are first treated in 2015. We are now left with a panel
dataset where some units are first treated in 2014 and the remaining
units are not treated during the sample period. We can then estimate the
effects of Medicaid expansion using a canonical two-way fixed effects
event-study specification,

$$
Y_{it} = \alpha_i + \lambda_t + \sum_{s \ne 2013} 1[s = t] \times D_i \times \beta_s + u_{it}
$$

where $D$ is 1 if a unit is first treated in 2014 and 0 otherwise.


```stata
local mixtape https://raw.githubusercontent.com/Mixtape-Sessions
use `mixtape'/Advanced-DID/main/Exercises/Data/ehec_data.dta, clear

* Keep years before 2016. Drop the 2016 cohort
keep if (year < 2016) & (missing(yexp2) | (yexp2 != 2015))

* Create a treatment dummy
gen byte D = yexp2 == 2014

* Run the TWFE spec
reghdfe dins b2013.year##D, absorb(stfips year) cluster(stfips) noconstant

local plotopts ytitle("Estimate and 95% Conf. Int.") title("Effect on dins")
coefplot, vertical yline(0) ciopts(recast(rcap)) xlabel(,angle(45)) `plotopts'
```

## Sensitivity analysis using relative magnitudes restrictions

xx

## Additional options and resources

See the [vignette](https://github.com/asheshrambachan/HonestDiD/blob/master/doc/HonestDiD_Example.pdf) for the R package. You can also view a video presentation about this paper [here](https://www.youtube.com/watch?v=6-NkiA2jN7U).

## Acknowledgements

This software package is based upon work supported by the National
Science Foundation Graduate Research Fellowship under Grant DGE1745303
(Rambachan) and Grant DGE1144152 (Roth). We thank Mauricio Cáceres Bravo
for his help in developing the package.
