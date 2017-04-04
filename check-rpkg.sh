##/bin/sh -e

echo Checking $2 with $1
$1 --vanilla CMD check --as-cran $2
