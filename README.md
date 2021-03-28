
 [![R-CMD-check](https://github.com/kabarigou/StanMoMo/workflows/R-CMD-check/badge.svg)](https://github.com/kabarigou/StanMoMo/actions)

StanMoMo
================

The `StanMoMo` package performs Bayesian **Mo**rtality **Mo**deling with
**Stan**, which is a C++ package for performing full Bayesian inference
(see <https://mc-stan.org/>). The current package supports a variety of
popular mortality models: the Lee-Carter (LC) model, the
Renshaw-Haberman model (LC with cohort effect), the Age-Period-Cohort
(APC) model, the Cairns-Blake-Dowd (CBD) model and the M6 model (CBD
with cohort effect). By a simple function call, the user obtains the
MCMC simulations for each parameter, the log likelihoods and death rates
predictions. Moreover, the package includes tools for model selection
and Bayesian model averaging by future-out cross-validation.

## Installation

To install `StanMoMo` from GitHub you will need `devtools`:

``` r
install.packages("devtools")
devtools::install_github('kabarigou/StanMoMo')
```

The installation may take a few minutes (models need to be compiled).
However, once installed, you can perform Bayesian mortality forecasting
in a matter of seconds.

After installing the package, you have to load the package via:

``` r
library(StanMoMo)
```

## Important note

The main purpose of this package is to provide users high-level
functions for estimating and forecasting mortality models in a Bayesian
setting without requiring any knowledge of the Stan modeling language.
This package depends on the [rstan](https://mc-stan.org/rstan/) package,
which translates the Stan model to C++ code, which is compiled into a
dynamic shared object (DSO). During the package installation, the
mortality models are pre-compiled and therefore you need a C++ compiler
on your machine (for instance, Rtools for Windows or Xcode on Mac, see
[here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for
more details). Once the package will be released on CRAN (date not
defined yet), all mortality models will be already pre-compiled in order
to run immediately when called, avoiding compilation time.

If you have any comments or suggestions about the package, feel free to
email <karim.barigou@univ-lyon1.fr>
