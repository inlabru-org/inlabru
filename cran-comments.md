## Submission notes

* Bugfix release 2.14.1

## R CMD check results and comments

* R CMD check results for 2.14.1 are clean on all platforms, with no new
  warnings or errors compared to 2.14.0.
* CRAN tests for 2.14.0 reported for R-devel:
  "Error in `sort.int(x, na.last = na.last, decreasing = decreasing, ...)`: 'x' must be atomic"
  which was due to overly narrow error detection, and have been fixed.
* CRAN messages in the tests for 2.14.0 for r-devel detecting possible > 2
  thread requests are due to an issue in the INLA package. The inlabru tests
  have been adjusted to avoid triggering the extra threads on CRAN:

  > test-aggregate.R: OMP: Warning #96: Cannot form a team with 3 threads, using 2 instead.
  > test-aggregate.R: OMP: Hint Consider unsetting KMP_DEVICE_THREAD_LIMIT (KMP_ALL_THREADS), KMP_TEAMS_THREAD_LIMIT, and OMP_THREAD_LIMIT (if any are set).

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
