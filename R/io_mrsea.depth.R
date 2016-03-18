# INTERNAL DATA STORAGE
io_mrsea.depth.getDataDir = function() {return(system.file("data",package="iDistance"))}


#' Load \link{mrsea.depth} survey data from raw data sets
#'
#' @aliases io_mrsea.depth.pkgdata.load
#' @export
#' @return \code{mrsea.depth} the \link{mrsea.depth} data set
#' @examples \\dontrun{mrsea.depth = io_mrsea.depth.pkgdata.load();}
#' @author Lindesay Scott-Hayward <\email{lass@st-andrews.ac.uk}>
#'

io_mrsea.depth.pkgdata.load = function() { 
  
  library("MRSea")
  data("predict.data.re")
  preddata = predict.data.re
  preddata = preddata[preddata$season==1 & preddata$impact==0,]
  colnames(preddata)[c(2,3)] = c("x","y")
  
  data(mrsea)
  depth = covdata.import(preddata, "depth", mrsea)
  depth$mesh.p4s = "+proj=utm +zone=32"
  return(depth)
  
}


#' Regenerate \link{mrsea.depth} data and store it to \code{mrsea.depth.RData}
#' 
#' Uses \code{\link{io_mrsea.depth.pkgdata.load}} to load the data and stores
#' the result to mrsea.depth.RData. Thereby the data that is distributed with 
#' our package is generated.
#'
#' @aliases io_mrsea.depth.pkgdata.save
#' @export
#' @return NULL
#' @examples \\dontrun{io_mrsea.depth.pkgdata.save();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_mrsea.depth.pkgdata.save = function(){
  ## save the data we will include in the R package
  mrsea.depth = io_mrsea.depth.pkgdata.load()
  save(mrsea.depth,file=paste0(io_mrsea.depth.getDataDir(),"/mrsea.depth.RData"))
}
