[![R-CMD-check](https://github.com/kabarigou/StanMoMo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kabarigou/StanMoMo/actions/workflows/R-CMD-check.yaml)


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
and Bayesian model averaging by leave-future-out validation.

## Installation

If you have R 4.0.0 or later and use Windows or Mac, you can directly
install the binary package via

``` r
install.packages("StanMoMo",repos=c("https://cloud.r-project.org",
"https://kabarigou.github.io/drat"),type = "binary",dependencies = TRUE)
```

Otherwise, you can also install the source package from Github via
`devtools`:

``` r
install.packages("devtools")
devtools::install_github('kabarigou/StanMoMo')
```

The installation of the source package may take a few minutes (models
need to be compiled). For this reason, we recommend to install the
binary package instead. Once the package is installed, you can perform
Bayesian mortality forecasting in a matter of seconds.

After installing the package, you have to load the package via:

``` r
library(StanMoMo)
```

## Important note if you install from source (install\_github)

The main purpose of this package is to provide users high-level
functions for estimating and forecasting mortality models in a Bayesian
setting without requiring any knowledge of the Stan modeling language.
This package depends on the [rstan](https://mc-stan.org/rstan/) package,
which translates the Stan model to C++ code, which is compiled into a
dynamic shared object (DSO). If you install the binary package (what we
recommend), the mortality models are already pre-compiled. If you
install the source package, the models are compiled during the
installation and therefore you need a C++ compiler on your machine (for
instance, Rtools for Windows or Xcode on Mac, see
[here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for
more details).

If you have any comments or suggestions about the package, feel free to
email <karim.barigou@univ-lyon1.fr>
