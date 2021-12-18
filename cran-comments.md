## Test environments
* ubuntu 20.04 (local, with INLA), R 4.1.2, R devel
* ubuntu 20.04 (on github, with INLA), R 4.1.2, R devel
* macOS-latest (on github, with INLA), R 4.1.2
* windows-latest (on github, with INLA), R 4.1.2
* win-builder, R devel
* R-hub (Fedora Linux, clang, gfortran) R devel

## Submission notes
* Bugfix and minor feature release 2.4.0
* NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
* CRAN checks for old version, 2.3.1:
  - M1 test errors: "... did not throw the expected warning"; this was
    due to changes in rgdal, and the test for this warning has been removed for
    broader compatibility with rgal versions.
  - M1 test error in
    "  ── Error (test_ipoints.R:141:3): Polygon integration with holes ── "
    "Error in `fmesher.read(prefix, "manifold")`: File"
    This appears to have been caused by unexpected floating point arithmetic
    issues for special test geometries, that have therefore been changed.
    We have also added a 30s timeout on all fmesher calls to prevent long-running
    tests for potential similar problems in the future.
* Checks for new version, 2.4.0:
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
    
  For R 4.0.5 and older, separate INLA package building may be needed to install
  the newer INLA version suggested by inlabru, but for the newer R versions
  direct installation can be done
    
* Spurious doi problem. The 10.1214/17-AOAS1078 doi is listed on
    https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-4/Point-process-models-for-spatio-temporal-distance-sampling-data-from/10.1214/17-AOAS1078.full
  which works, but https://doi.org/10.1214/17-AOAS1078 may give an invalid redirect
  response with a malformed version of the target URL.
  Found the following (possibly) invalid DOIs:
    DOI: 10.1214/17-AOAS1078
      From: inst/CITATION
      Status: libcurl error code 6:
      	Could not resolve host: projecteuclid.org%5Cjournals%5Cannals-of-applied-statistics%5Cvolume-11%5Cissue-4%5CPoint-process-models-for-spatio-temporal-distance-sampling-data-from%5C10.1214
      Message: Error


## Downstream dependencies
inlabru does not have any reverse dependencies
