## Resubmission

Comments from Kurt Hornik:
```
  Dear maintainer,

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_StanMoMo.html>.

Please correct before 2022-07-04 to safely retain your package on CRAN.

Best,
-k
```

The error was caused by the function load_HMD_data which depends on a website that was recently updated. We have now removed this problematic function and the updated package has no more errors.

## Test environments
* macOS-latest, R release
* windows-latest, R release
* ubuntu-latest, R release

## R CMD check results

0 errors | 0 warnings | 2 notes

❯ checking installed package size ... NOTE
    installed size is  6.9Mb
    sub-directories of 1Mb or more:
      libs   5.8Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
  

* This is an update of the StanMoMo package.



