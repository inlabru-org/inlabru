## Submission notes

* Minor release 2.10.0, including new features and bug fixes
* CRAN checks for old version, 2.9.0:
  NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
* Checks for new version, 2.10.0 (with latest INLA, 23.10.28):
  - Spurious error message about potentially invalid doi, see below

## R CMD check results and comments

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

We checked 7 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
