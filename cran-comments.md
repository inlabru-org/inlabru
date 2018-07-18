## Test environments
* ubuntu 18.04 (local), R release (3.5.1)
* ubuntu 18.04 (local), R devel (2018-07-11)
* ubuntu 18.04 (local), R 3.4.4
* ubuntu 14.04 (on travis-ci), R 3.5.0
* TODO: win-builder (devel and release)

On the local installs, --run-donttest was activated to also test
long running examples and tests, with INLA 18.07.12 installed.

## R CMD check results

0 errors | 0 warnings | 0 notes

This is a bugfix release.
  
  * Corrected INLA repository URL
  * Namespace corrections and documentation updates
  * Computational method robustification
