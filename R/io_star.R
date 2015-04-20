##
## Functions to read and write STAR data
##
## Author: Fabian Bachl (FEB)
##


io_star.getDataDir = function()
{ ## Returns the path the STAR boundary and coast line data is stored in
  return("data/")
}



io_star.boundary = function()
{ ## Read the START boundary file "BoundSTAR.dat"
  boundSTART = read.table(paste(io_star.getDataDir(), "BoundSTAR.dat", sep = .Platform$file.sep), header = F, skip = 1)
  colnames(boundSTART) = c("lat", "lon")
  
  # RETURN
  return(boundSTART)
}

io_star.coast = function()
{ ## Read the START coast line file "AreaSTAR2.dat"
  
  coast = read.table(paste(io_star.getDataDir(), "AreaSTAR2.dat", sep = .Platform$file.sep), sep = ",", 
                     colClasses = "numeric", comment.char = "*", col.names = c("lat", 
                                                                               "lon"))
  ## original dataset is not numberic
  coast$lat = as.numeric(as.character(coast$lat))
  coast$lon = as.numeric(as.character(coast$lon))
  
  # RETURN
  return(coast)
}
