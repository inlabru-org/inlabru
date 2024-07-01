## Submission notes

* Resubmission of 2.11.0 as 2.11.1, with additional bugfixes
* Fixes issue with package cross references in documentation; these had gone
  undetected by all local and github testing, but were caught by CRAN checks.
  I was unable to get any R version check (including the very latest devel version)
  to generate the actual NOTE, but the problem was clear and quickly fixed.

## R CMD check results and comments

* No change in RMD CMD check results
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

We checked 6 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
