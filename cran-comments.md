## Test environments
* ubuntu 20.04 (local, with INLA), R 3.6.2, R 4.0.4, R devel
* ubuntu 20.04 (on github, with INLA), R 4.0.4, R devel
* macOS-latest (on github, with INLA), R 4.0.4
* windows-latest (on github, with INLA), R 4.0.4
* win-builder, R 4.0.4, R devel
* R-hub, R devel

## Submission notes
* Bugfix release 2.3.1
* CRAN checks for old version, 2.3.0:
  - ERROR: examples that relied on PROJ6 or later now check if PROJ6 is available
  - NOTE: Rd xrefs; remaining unavaiable xrefs have been removed
  - WARN: Error in readRDS; package-unrelated problem with the CRAN system
          for Check: Rd cross-references
* Checks for new version, 2.3.1:
  - NOTE: Additional_repositories used for non-CRAN Suggested package INLA

## R CMD check results

With the new version, 2.3.1, there were no ERRORs or WARNINGs on the test systems,
except for a spurious "Suggests orphaned package: ggmap" message on R-hub, and
windows on github but not win-builder; see below for details.

Comments:

* The non-CRAN Suggested package INLA has been extensively tested with inlabru
  locally and in github actions for both Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:

  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/testing

* When checking on R-hub, it says that
     Suggests orphaned package: 'ggmap'
  To the best of my knowledge, this is a problem with R-hub itself.
  ggmap _was_ orphaned at some point, but is no longer orphaned.
  rhub appears to be confused about this, based on the discussion here:
  https://community.rstudio.com/t/orphaned-package-on-windows-build/84165/4

## Downstream dependencies
inlabru does not have any reverse dependencies
