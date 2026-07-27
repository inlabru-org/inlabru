## Submission notes

* Feature release 2.15.0

## R CMD check results and comments

* R CMD check results for 2.15.1 are clean on all platforms, with no new
  warnings or errors compared to 2.14.1.

* The non-CRAN Suggested package INLA has been extensively tested with inlabru
  locally and in github actions for both Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:
```
Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/testing
``` 

## revdepcheck results

We checked all 9 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
