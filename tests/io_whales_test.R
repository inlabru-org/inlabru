##
## This file is for testing the I/O functions for the whale data set 
##
## Author: Fabian Bachl (FEB)
##

source("io_whales.R") 

## Load sightings
  
  sight = io_whales.sightings(years=c(2000:2006),E=1,BfLim=5,PDLim=6)  

## Load effort

  eff = io_whales.effort()

## Load boundaries

  bound = io_whales.boundary()   # which internally uses io_star.boundary()

## Load coast

  coast = io_whales.coast() # which internally uses io_star.coast()

