## Submission notes

* Minor release 2.9.0, including new features, speed improvements, bug fixes,
  and package dependency updates (Moved 'sp' from Depends to Imports to prepare
  for future move to Suggests, and added dependency on the new 'fmesher' package)
* CRAN checks for old version, 2.8.0:
  NOTE: Additional_repositories is used for non-CRAN Suggested package INLA
  NOTE: Documented arguments not in \usage in documentation object
* Checks for new version, 2.9.0 (with latest INLA, 23.08.26):
  - Spurious error message about potentially invalid doi, see below

## R CMD check results

Comments:

* All "Documented arguments not in \usage in documentation object" cases have
  been corrected

* The non-CRAN Suggested package INLA has been extensively tested with inlabru
  locally and in github actions for both Linux, Windows, and macOS.
  The needed repository specification is included in the package DESCRIPTION:
```
Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/testing
``` 
* Spurious URL issue noted by `urlchecker`:
```
> urlchecker::url_check()
✖ Error: README.md:37:31 403: Forbidden
[doi:10.1111/2041-210X.13168](https://doi.org/10.1111/2041-210X.13168),
                              ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```                              
  Manually accessing https://doi.org/10.1111/2041-210X.13168 correctly leads to
```
  https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13168
```

## revdepcheck results

We checked 7 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
