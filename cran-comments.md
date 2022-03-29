## Test environments
* ubuntu 20.04 (local, with INLA), R 4.0.5, R 4.1.3, R devel
* ubuntu 20.04 (on github, with INLA), R 4.0.5, R 4.1.3, R devel
* macOS-latest (on github, with INLA), R 4.1.3
* windows-latest (on github, with INLA), R 4.1.3
* win-builder, R devel
* For the github platforms, separate tests were also
  done without installing packages in Suggests.

## Submission notes
* Bugfix release 2.5.1
* NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
* CRAN checks for old version, 2.5.0:
  ERROR: In v2.5.0, `...names()` had been introduced without noting it was only
  introduced in R4.1. Pre-submission testing had been turned off for the
  old R release due to technical issues in the past that are now resolved.
  The 2.5.1 code has been reverted to the `names(list(...))` pattern,
  to support use on R 4.0.5.
* Checks for new version, 2.5.1 (with latest INLA, 22.03.28):
  - Spurious error message about potentially invalid doi, see below

## R CMD check results

Comments:

* The non-CRAN Suggested package INLA has been extensively tested with inlabru
  locally and in github actions for both Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:

  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/testing
    
* Spurious doi problems noted by win-builder.
  Manually accessing https://doi.org/10.1111/2041-210X.13168 correctly leads to
  https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13168
  Manually accessing https://doi.org/10.1214/17-AOAS1078 correctly leads to
  https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/
    issue-4/Point-process-models-for-spatio-temporal-distance-sampling-data-from/10.1214/17-AOAS1078.full

  The checking messages:
  
  Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1111/2041-210X.13168
    From: README.md
    Status: 503
    Message: Service Unavailable

  Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.13168
    From: DESCRIPTION
          inst/CITATION
    Status: Service Unavailable
    Message: 503
  DOI: 10.1214/17-AOAS1078
    From: inst/CITATION
    Status: Internal Server Error
    Message: 500


## Downstream dependencies
inlabru does not have any reverse dependencies

## Failure reports for previous inlabru version 2.5.1

checking R code for possible problems ... [13s/13s] NOTE
bru_mapper.default: no visible global function definition for
  ‘...names’
gg.inla.mesh: no visible global function definition for ‘...names’
Undefined global functions or variables:
  ...names
