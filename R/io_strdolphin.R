# INTERNAL DATA STORAGE
io_strdolphin.getDataDir = function() {return(system.file("data",package="iDistance"))}


#' Load \link{strdolphin} survey data from raw data sets
#'
#' @aliases io_strdolphin.pkgdata.load
#' @export
#' @return \code{strdolphin} the \link{strdolphin} data set
#' @examples \\dontrun{strdolphin = io_strdolphin.pkgdata.load();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.pkgdata.load = function()
{ 
  effort = io_strdolphin.effort(rem.intermediate=FALSE)
  par3d.args = io_strdolphin.par3d()
  strdolphin = list(sighting=io_strdolphin.sighting(BfLim=5,PDLim=6),
                   transect = as.transect(effort),
                   effort = effort,
                   mesh = io_strdolphin.mesh(),
                   coast.boundary = io_star.coast(),
                   sea.boundary = io_star.boundary(),
                   inner.boundary = io_strdolphin.inner.boundary(),
                   geometry = "geo",
                   mesh.coords = c("lon","lat"),
                   time.coords = "year",
                   par3d.args = par3d.args,
                   ips = io_strdolphin.ips(),
                   ips.yearly = io_strdolphin.ips.yearly()
  )
  class(strdolphin) = c("strdolphin","etpdata","dsdata","list")
  return(strdolphin)
}

#' Regenerate \link{strdolphin} data and store it to \code{strdolphin.RData}
#' 
#' Uses \code{\link{io_strdolphin.pkgdata.load}} to load the data and stores
#' the result to strdolphin.RData. Thereby the data that is distributed with 
#' our package is generated.
#'
#' @aliases io_strdolphin.pkgdata.save
#' @export
#' @return NULL
#' @examples \\dontrun{io_strdolphin.pkgdata.save();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.pkgdata.save = function(){
  ## save the data we will include in the R package
  strdolphin = io_strdolphin.pkgdata.load()
  save(strdolphin,file=paste0(io_strdolphin.getDataDir(),"/strdolphin.RData"))
}

#' Load \link{strdolphin} integration points
#' 
#'
#' @aliases io_strdolphin.ips
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.ips = function(){
  load(file= paste(io_strdolphin.getDataDir(), "etp.ips.RData",sep=.Platform$file.sep))
  return(ips)
}

#' Load \link{strdolphin} integration points (yearly)
#' 
#'
#' @aliases io_strdolphin.ips.yearly
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.ips.yearly = function(){
  load(file= paste(io_strdolphin.getDataDir(), "etp.ips.yearly.RData",sep=.Platform$file.sep))
  return(ips.yearly)
}


#' Load \link{strdolphin} inner boundary
#' 
#'
#' @aliases io_strdolphin.inner.boundary
#' @export
#' @return boundary inner boundary
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.inner.boundary = function(){
  load(file= paste(io_strdolphin.getDataDir(), "etp.inner.boundary.Rdata",sep=.Platform$file.sep))
  return(innerBnd)
}

#' Load \link{strdolphin} par3d parameters
#' 
#'
#' @aliases io_strdolphin.par3d
#' @export
#' @return params pard3d() parameters for rgl plot of strdolphin
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.par3d = function(){
  #par3d(zoom=0.62,FOV=1)
  #par3d.args = par3d(no.readonly=TRUE)
  #save(par3d.args,file="data/strdolphin.par3d.RData")
  load(file= paste(io_strdolphin.getDataDir(), "etp.par3d.RData",sep=.Platform$file.sep))
  params = par3d.args
  return(params)
}


#' Load \link{strdolphin} sightings
#' 
#' Load whale sightings from "Bmus86_06.csv" and apply filters (see parameters).
#'
#' @aliases io_strdolphin.sighting
#' @export
#' @return sightings Whale sighting
#' @examples \\dontrun{sightings = io_strdolphin.sightings(years=c(2000:2006),E=1,BfLim=5);}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.sighting = function(years=NULL,E=NULL,BfLim=NULL,PDLim=NULL)
{ 
  # Load the file
  sightdat = read.csv(file= paste(io_strdolphin.getDataDir(), "strdolphin.sighting.csv",sep=.Platform$file.sep),head=TRUE,sep=",")

  # Remove NA lines
  sightdat = sightdat[!is.na(sightdat$perpdist),]
  
  # rename
  names(sightdat)[names(sightdat) %in% "perpdist"] = "PD"
  
  # Filter out on.eff=0
  sightdat = sightdat[sightdat$on.eff==1, ]
  
  # Year selection
  if (is.numeric(years)) {
    #print("io_strdolphin.sightings(): Filtering years")
    sightdat = sightdat[sightdat$Year %in% years, ]
  }
  
  # Get rid of the non-LT-mode sightings
  if (is.numeric(E)) {
    #print("io_strdolphin.sightings(): Filtering E")
    sightdat = sightdat[sightdat$E==E, ]
  }
  
  # Tim commented that it is not LT mode if Bft>5
  if (is.numeric(BfLim)) {
    #print("io_strdolphin.sightings(): Filtering Bf")
    sightdat = sightdat[sightdat$bft<=BfLim,]
  } 
  
  # Distance limit
  if (is.numeric(PDLim)) {
    #print("io_strdolphin.sightings(): Filtering PD")
    sightdat = sightdat[sightdat$PD<=PDLim,]
  } 
  
  # Rename date fields
  names(sightdat)[which("Year" == names(sightdat))] = "year"
  names(sightdat)[which("Month" == names(sightdat))] = "month"
  names(sightdat)[which("Day" == names(sightdat))] = "day"
  names(sightdat)[which("Time" == names(sightdat))] = "time"
  
  # Rename "Lat" and "Long" to "lat" and "lon"
  names(sightdat)[which(names(sightdat) %in% c("Lat","Long"))] = c("lat","lon")
  
  # Rename "PD" to "distance"
  names(sightdat)[which(names(sightdat) %in% c("PD"))] = c("distance")
  
  # RETURN
  return(sightdat)
} 


#' Generate \link{strdolphin} transect lines from ETP effort data
#'
#' @aliases io_strdolphin.transect
#' @export
#' @return tr Data frame with transect lines
#' @examples \\dontrun{transect = io_strdolphin.transect();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'


io_strdolphin.transect = function(years=FALSE,rem.intermediate=TRUE)
{ 
  eff = io_strdolphin.effort(years=years,rem.intermediate=rem.intermediate)
  
  eff = eff[!is.na(eff[,1]),]
  
  
  tr = cbind(eff[1:(dim(eff)[1]-1),],eff[2:dim(eff)[1],c("lat","lon","bft")])
  colnames(tr)[6] = "start.lat"
  colnames(tr)[7] = "start.lon"
  colnames(tr)[8] = "start.bft"
  colnames(tr)[10] = "end.lat"
  colnames(tr)[11] = "end.lon"
  colnames(tr)[12] = "end.bft"
  
  # remove entries with star.bft=NA (no transects!)
  tr = tr[!is.na(tr$start.bft),]
  
  
  class(tr) = c("transect","data.frame")
  return(tr)
}

#' Load \link{strdolphin} effort
#' 
#' Load effort data from "ETP effort.csv"
#'
#' @aliases io_strdolphin.effort
#' @export
#' @return effort ETP effort related to \link{strdolphin} data
#' @examples \\dontrun{transects = io_strdolphin.transects();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.effort = function(years=FALSE,rem.intermediate=TRUE)
{ 
  # read effort data
  tr = read.csv(file=paste(io_strdolphin.getDataDir(),"ETP effort.csv",sep=.Platform$file.sep))
  
  # Remove intermediate GPS records and NA lines
  na.rows = which(is.na(tr[,1]))
  if (rem.intermediate) {
    end.idx = na.rows-1
    start.idx = c(1,na.rows[1:length(na.rows)-1]+1)
    idx = as.vector(rbind(start.idx,end.idx,na.rows))
    tr = tr[idx,]
  } else {
    #tr = tr[!is.na(tr[,1]),]
  }
  
  # Remove bad transects with lat/lon == 0
  tr = tr[!(tr[,"lat"]==0 & tr[,"lon"]==0) ,]
  
  # Filter years
  if (is.numeric(years)) {
    tr = tr[tr$year %in% years,]
  }
  
  #
  # Return
  # 
  
  class(tr) = c("etpeffort","data.frame")
  return(tr)
} 


#' Load \link{strdolphin} boundary
#' 
#' @aliases io_strdolphin.boundary
#' @export
#' @return boundary 
#' @examples \\dontrun{boundary = io_strdolphin.boundary();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.boundary = function(...) {return(io_star.boundary(...))}

#' Load \link{strdolphin} coast
#' 
#' @aliases io_strdolphin.coast
#' @export
#' @return coast
#' @examples \\dontrun{coast = io_strdolphin.coast();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.coast = function(...) {return(io_star.coast(...))}

#' Load \link{strdolphin} mesh
#' 
#' @aliases io_strdolphin.mesh
#' @export
#' @return mesh An INLA mesh
#' @examples \\dontrun{mesh = io_strdolphin.mesh(); plot(mesh)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_strdolphin.mesh = function(...) 
{
  load(file=paste(io_strdolphin.getDataDir(), "mesh_stitchInterior.Rdata",sep="/"))
  return(mesh)
}
