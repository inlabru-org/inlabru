## Test environments
* ubuntu 20.04 (local, with INLA), R 3.6.2, R 4.0.4, R devel
* ubuntu 20.04 (on github, with INLA), R 4.0.4, R devel
* macOS-latest (on github, with INLA), R 4.0.4
* windows-latest (on github, with INLA), R 4.0.4
* win-builder, R 4.0.4, R devel
* R-hub (Ubuntu 20.04.1) R 4.0.4

## Submission notes
* Bugfix release 2.3.1
* NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
* CRAN checks for old version, 2.3.0:
  - ERROR: Examples for fm_CRS failed on pre-PROJ6 system.
           Examples that relied on PROJ6 or later now check if PROJ6 is available
  - NOTE: Rd xrefs; remaining unavailable xrefs have been removed
  - WARN: Error in readRDS; package-unrelated problem with the CRAN system
          for Check: Rd cross-references
          Spurious issue for r-devel-linux-x86_64-debian-gcc
* Checks for new version, 2.3.1:
  - Spurious note about orphaned 'ggmap' on windows systems, see below
  - Spurious error message about potentially invalid doi, see below

With the new version, 2.3.1, there were no ERRORs or WARNINGs on the test systems,
except for a spurious "Suggests orphaned package: ggmap" on windows systems,
and a spurious doi error; see below for details.

## R CMD check results

Comments:

* The non-CRAN Suggested package INLA has been extensively tested with inlabru
  locally and in github actions for both Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:

  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/testing
    
* Spurious doi problem. The 10.1214/17-AOAS1078 doi is listed on
    https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-4/Point-process-models-for-spatio-temporal-distance-sampling-data-from/10.1214/17-AOAS1078.full
  which works, but https://doi.org/10.1214/17-AOAS1078 gives an invalid redirect
  response with a malformed version of the target URL.
  Found the following (possibly) invalid DOIs:
    DOI: 10.1214/17-AOAS1078
      From: inst/CITATION
      Status: libcurl error code 6:
      	Could not resolve host: projecteuclid.org%5Cjournals%5Cannals-of-applied-statistics%5Cvolume-11%5Cissue-4%5CPoint-process-models-for-spatio-temporal-distance-sampling-data-from%5C10.1214
      Message: Error

* When checking on Windows on R-hub and github, it says that
     Suggests orphaned package: 'ggmap'
  To the best of my knowledge, this is a problem with how the system itself
  tracks orphanded packages.
  ggmap _was_ orphaned at some point, but is no longer orphaned.
  rhub and others appear to be confused about this, based on the discussion here:
  https://community.rstudio.com/t/orphaned-package-on-windows-build/84165/4

## Downstream dependencies
inlabru does not have any reverse dependencies
