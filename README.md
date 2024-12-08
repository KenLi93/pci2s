
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pci2s: Regression-based proximal causal inference using two-stage-least-squares

<!-- badges: start -->

[![R-CMD-check](https://github.com/r-lib/usethis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-lib/usethis/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

We provide estimation and inference for regression-based proximal causal
inference with a variety data types of outcome variables and negative
control outcomes. Details can be found at Liu et
al. ([2024](#ref-p2sls-glm)) and Li et al. ([2024](#ref-p2sls-ah)) An
introduction of proximal causal inference can be found at Shi et
al. ([2020](#ref-review)) and Tchetgen Tchetgen ([2024](#ref-intro)).

## Installation

Install from GitHub:

``` r
# install.packages("pak")
pak::pak("KenLi93/pci2s")
```

## References

<div id="refs" class="references">

<div id="ref-p2sls-glm">

Liu, J., Park, C., Li, K. and Tchetgen Tchetgen, E.J., 2024.
“Regression-based proximal causal inference.” *American Journal of
Epidemiology* p.kwae370.

</div>

<div id="ref-p2sls-ah">

Li, K., Linderman, G.C., Shi, X. and Tchetgen, E.J.T., 2024.
“Regression-based proximal causal inference for right-censored
time-to-event data.” *arXiv preprint* arXiv:2409.08924.

</div>

<div id="ref-review">

Shi, X., Miao, W. and Tchetgen, E.T., 2020. “A selective review of
negative control methods in epidemiology.” *Current epidemiology
reports* 7, pp.190-202.

</div>

<div id="ref-intro">

Tchetgen Tchetgen, E.J., Ying, A., Cui, Y., Shi, X. and Miao, W., 2024.
“An introduction to proximal causal inference.” *Statistical Science*
39(3), pp.375-390.

</div>

</div>
