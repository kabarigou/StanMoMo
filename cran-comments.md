## Resubmission

Comments from reviewer Gregor Seyer:
```
\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.)


Please fix and resubmit.

Best,
Gregor Seyer

```
In this resubmission, we have unwrapped the examples and simplified them so that they are executable in < 5 sec. Only two of them were wrapped with \donttest because they took more than 5 sec to be executed.


## Resubmission

Comments from reviewer Uwe Higges:
```
 Check: Overall checktime, Result: NOTE
   Overall checktime 11 min > 10 min

mainly from

* checking re-building of vignette outputs ... [488s] OK

Please reduce the vignette build timings by using
  - small toy data only
  - few iterations
  - or by providing precomputed results for the most lengthy parts. 
Otherwise we cannot afford checking the vignette regularly on CRAN.

Please fix and resubmit.

Best,
Uwe Ligges


```
In this resubmission, we have - used smaller data and - used fewer iterations. 
Now, the vignette checks below 10 minutes: 
```
checking re-building of vignette outputs ... [35s]
```
Moreover, we have via check_win_devel: 
```
check time in seconds: 194
```

## Resubmission

Comments from reviewer Julia Haider:
```
  Please do not start the description with "This package", package name, title or similar.
Besides that: Please always write package names, software names and API (application programming interface) names in single quotes in title and description. e.g: --> 'StanMoMo', 'Stan'

Please always explain all acronyms in the description text.

If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> authors (year) <arXiv:...> authors (year, ISBN:...) or if those are not available: <https:...> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or 'F' as vector names.), e.g.:
   man/fit_mo_mo.Rd

Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. 
\value{No return value, called for side effects} or similar) Missing Rd-tags:
      forecasting_plot.Rd: \value

Please fix and resubmit.

Best,
Julia Haider

```
In this resubmission, we have - revised the description following your comments. 
- Wrote references in the description of DESCRIPTION file in the form: 
"authors (year) <doi:...> " - Wrote FALSE instead of F - add \value to missing Rd-tags. 
  
## Resubmission

Comments from reviewer Uwe Ligges:
```
   Found the following (possibly) invalid URLs:
     URL: http://www.mortality.org/ (moved to https://www.mortality.org/)
       From: man/FRMaleData.Rd
       Status: 301
       Message: Moved Permanently

Please change http --> https, add trailing slashes, or follow moved content as appropriate.

Please fix and resubmit.

Best,
Uwe Ligges
```
In this resubmission, we have changed http --> https as requested.

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
  

* This is a new release.



