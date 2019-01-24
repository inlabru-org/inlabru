FROM ubuntu:18.04
MAINTAINER fbachl
## Based on https://rtask.thinkr.fr/blog/installation-of-r-3-5-on-ubuntu-18-04-lts-and-tips-for-spatial-packages/

# CONFIGURE TIMEZONE
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Utilities
RUN apt-get update && apt-get install -y gnupg ca-certificates pandoc

# Enable UBUNTU GIS repository
RUN echo  'deb http://ppa.launchpad.net/ubuntugis/ubuntugis-unstable/ubuntu bionic main' >> /etc/apt/sources.list.d/ubuntugis.sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 314DF160

# Add R 3.5 repository
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' >> /etc/apt/sources.list.d/cran35.sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

RUN apt-get -y update
RUN apt-get -y upgrade

# Needed for add-apt-repository:
RUN apt-get -y install software-properties-common

# Libs needed for devtools:
RUN apt-get -y install libssl-dev

# Install GEOSPATIAL UBUNTU PACKAGES
RUN apt-get -y install libgdal-dev libproj-dev libgeos-dev libudunits2-dev libv8-dev libcairo2-dev libnetcdf-dev libgeos++-dev

# INSTALL R
RUN apt-get -y install r-base r-base-core r-recommended

# Add repository of precompiled R packages
RUN add-apt-repository -y ppa:marutter/c2d4u3.5
RUN apt-get update
RUN apt-get -y install r-cran-rgl


## Roxygen and Devtools and its dependencies

RUN R -e "install.packages(c('roxygen2'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('devtools'), repos='https://cran.rstudio.com/')"

## inlabru requirements (except INLA)
RUN R -e "install.packages(setdiff(c('rgdal', 'rgeos', 'sp', 'testthat', 'ggmap', 'rgl', 'sphereplot', 'raster', 'dplyr', 'maptools', 'mgcv', 'shiny', 'spatstat', 'spatstat.data', 'RColorBrewer', 'graphics', 'knitr', 'rmarkdown'), installed.packages()[,1]))"

## Pkgdown
RUN R -e "install.packages('pkgdown')"

## INLA
RUN R -e "install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)"

