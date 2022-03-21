## Test environments
* ubuntu 20.04 (local, with INLA), R 4.1.2, R devel
* ubuntu 20.04 (on github, with INLA), R 4.1.2, R devel
* macOS-latest (on github, with INLA), R 4.1.2
* windows-latest (on github, with INLA), R 4.1.2
* win-builder, R devel
* R-hub (Fedora Linux, clang, gfortran) R devel

## Submission notes
* Bugfix and minor feature release 2.5.0
* NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
* CRAN checks for old version, 2.4.0:
  - donttest: "failure: length > 1 in coercion to logical" (see partial reports below)
    These errors were traced to the INLA package, and are fixed there since
    version 22.03.03; The CRAN installations of INLA need to be upgraded
    to fix this issue.  The "stable" INLA version hasn't been changed yet, but
    the inlabru package points to the "testing" version, see
    https://www.r-inla.org/download-install
    For some platforms (notably with old GLIBC versions), inla.binary.install()
    needs to be run after the regular package installation, to install compatible
    program binaries.
* Checks for new version, 2.5.0 (with latest INLA, 22.03.16):
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
  direct installation can be done. Bugs relating to vector-if-statements
  were fixed in INLA version 22.03.03
    
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



## Failure reports

### INLA 22.02.16-2 and INLA 22.02.28-1

 ----------- FAILURE REPORT -------------- 
 --- failure: length > 1 in coercion to logical ---
 --- srcref --- 
: 
 --- package (from environment) --- 
INLA
 --- call from context --- 
FUN(X[[i]], ...)
 --- call from argument --- 
is.null(x) || is.na(x)
 --- R stacktrace ---
where 1: FUN(X[[i]], ...)
where 2: lapply(X = X, FUN = FUN, ...)
where 3: sapply(marginals.random, function(x) (is.null(x) || is.na(x)))
where 4: inla.collect.random(results.dir, debug)

