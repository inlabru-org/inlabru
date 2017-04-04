##
##
## First:
##   make stage
## This clones the repository from $(PKGPATH) (henceforth called the origin.)
## into a staging clone, called stage.
##
## Keep the staging clone up-to-date by pulling from the origin:
##   make stage-update
##
## Which version to build and check?
##   make branch-origin What it the current origin branch?
##   make branch        What is the current staging branch?
##   make branch.origin Set the staging branch to the current origin branch?
##   make branch.name   Set the staging branch to "name".
## The latter 2 command automatically run 'make stage-update' first
##
## Build a package:
##   make rpkg
##   make rpkg-Rdev
##   make rpkg-source
##
## Build with a different make tool
## (bmake and pmake are the same; install the bmake Ubuntu package):
##   make bmake-rpkg
##   make bmake-rpkg-Rdev:
##
## Check the latest built *source* package (from make rpkg-source):
##   make check
##   make check-Rdev
##
## Install the latest built package:
##   make install
##   make install-source
##
## List config variables:
##   make config

PKG=inlabru
PKGPATH=../inlabru/

STAGE=stage
BUILD=build

SED = sed
CPR = rsync -a --delete
R = R
RDEV = ~/R-dev/bin/R

FILENAME_SRC := $(shell ls -rt $(PKG)_*.tar.gz | grep -v _R_ | tail -1 ) 
FILENAME_BIN := $(shell ls -rt $(PKG)_*.tar.gz | grep _R_ | tail -1 )

BRANCH := $(shell cd $(STAGE); git branch | grep "^\*" | $(SED) "s/^\* //")
BRANCH_ORIGIN := $(shell cd $(PKGPATH); git branch | grep "^\*" | $(SED) "s/^\* //")

default:

branch:
	@echo "Current branch: "$(BRANCH)
branch-origin:
	@echo "Current branch on local origin: "$(BRANCH_ORIGIN)

branch.origin: stage-update
	@cd $(STAGE) ; git checkout -f $(BRANCH_ORIGIN)
	@make branch

branch.%: stage-update
	@cd $(STAGE) ; git checkout -f $*
	@make branch

stage-cleanup:
	@cd $(STAGE) ; git checkout -f $(BRANCH)
stage-update: stage-cleanup
	@cd $(STAGE) ; git pull

stage:
	@git clone $(PKGPATH) $(STAGE)/
	@make branch

prepare: stage-cleanup clean-build

roxy: stage-cleanup
	@cd $(STAGE) ; Rscript -e "library(devtools); document()"

build-dir:
	mkdir -p $(BUILD)
rpkg: prepare build-dir
	$(R) CMD INSTALL -l $(BUILD) --build $(STAGE)
	$(MAKE) clean-build
rpkg-Rdev: prepare build-dir
	$(RDEV) CMD INSTALL -l $(BUILD) --build $(STAGE)
	$(MAKE) clean-build

rpkg-source: prepare
	$(R) CMD build --build $(STAGE)
rpkg-source-Rdev: prepare
	$(RDEV) CMD build --build $(STAGE)

clean-stage:
	rm -rf $(STAGE)
clean-build:
	rm -rf $(BUILD)

clean: clean-stage clean-build

## Check building with a non-GNU-make:
bmake-rpkg: prepare build-dir
	@./bmake-check.sh $(R) $(BUILD) $(STAGE)
	$(MAKE) clean-build
bmake-rpkg-Rdev: prepare build-dir
	@./bmake-check.sh $(RDEV) $(BUILD) $(STAGE)
	$(MAKE) clean-build

## Check the latest built source package:
check:
	@./check-rpkg.sh $(R) $(FILENAME_SRC)
check-Rdev:
	@./check-rpkg.sh $(RDEV) $(FILENAME_SRC)

## Install latest built packages:
install:
	$(R) CMD INSTALL $(FILENAME_BIN)
install-source:
	$(R) CMD INSTALL $(FILENAME_SRC)

config:
	@echo PKG $(PKG)
	@echo PKGPATH $(PKGPATH)
	@echo BRANCH $(BRANCH)
	@echo BRANCH_ORIGIN $(BRANCH_ORIGIN)
	@echo FILENAME_SRC $(FILENAME_SRC)
	@echo FILENAME_BIN $(FILENAME_BIN)

.PHONY: default prepare roxy rpkg rpkg-source config \
	bmake-rpkg bmake-rpkg-Rdev check check-Rdev \
	stage stage-cleanup stage-update \
	clean-stage clean-build clean \
	build-dir
