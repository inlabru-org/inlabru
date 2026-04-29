## Submission notes

* Bugfix release 2.14.1

## R CMD check results and comments

* CRAN tests reporting
  "Error in `sort.int(x, na.last = na.last, decreasing = decreasing, ...)`: 'x' must be atomic"
  were due to overly narrow error detection, and have been fixed.
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

We checked all 8 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
