## Submission notes

* Minor release 2.7.0, including new features, bug fixes, and
  fixes for package dependencies (spatstat.core removed)
* CRAN checks for old version, 2.6.0:
  NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
* Checks for new version, 2.7.0 (with latest INLA, 22.11.28-1):
  - Spurious error message about potentially invalid doi, see below

## Test environments

* ubuntu 22.04 (local, with INLA), R 4.2.2, R 4.1.3, R devel
* ubuntu 20.04 (on github, with INLA), R 4.2.2, R 4.1.3, R devel
* macOS-latest (on github, with INLA), R 4.2.2
* windows-latest (on github, with INLA), R 4.2.2
* win-builder; R devel, 4.2.2, 4.1.3
* R-hub;
    Windows Server R-devel
    Fedora Linux R-devel
    Ubuntu Linux 20.04 R-release
* For the github platforms, separate tests were also
  done without installing packages in Suggests.

## R CMD check results

Comments:

* The non-CRAN Suggested package INLA has been extensively tested with inlabru
  locally and in github actions for both Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:
```
Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/testing
``` 
* Spurious URL issue noted by `urlchecker` and some other checkers:
```
> urlchecker::url_check()
✖ Error: README.md:37:31 503: Service Unavailable
[doi:10.1111/2041-210X.13168](https://doi.org/10.1111/2041-210X.13168),
                              ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```                              
  Manually accessing https://doi.org/10.1111/2041-210X.13168 correctly leads to
```
  https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13168
```
  Manually accessing https://doi.org/10.1214/17-AOAS1078 correctly leads to
```
  https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/
    issue-4/Point-process-models-for-spatio-temporal-distance-sampling-data-from/10.1214/17-AOAS1078.full
```

## Downstream dependencies

We tested the reverse dependencies (bmstdr, intSDM, PointedSDMs, rSPDE) and
found no problems.
