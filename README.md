StanMoMo
================

The `StanMoMo` package performs Bayesian Mortality Modeling with Stan for a variety of popular mortality models. The current package supports the Lee-Carter (LC) model, the Renshaw-Haberman model (LC with cohort effect), the Age-Period-Cohort (APC) model, the Cairns-Blake-Dowd (CBD) model and the M6 model (CBD with cohort effect). By a simple function call, the user obtains the MCMC simulations for each parameter, the log likelihoods and death rates predictions. Moreover, the package includes tools for model selection and Bayesian model averaging by future-out cross-validation.

Installation
------------

To install `StanMoMo` from GitHub you will need `devtools`:

``` r
install.packages("devtools")
devtools::install_github('kabarigou/StanMoMo')
```

After installing the package, you have to load the package via:

``` r
library(StanMoMo)
```
