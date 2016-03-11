# INTERNAL DATA STORAGE
io_sbdolphin.getDataDir = function() {return(system.file("data",package="iDistance"))}


#' Load \link{sbdolphin} survey data from raw data sets
#'
#' @aliases io_sbdolphin.pkgdata.load
#' @export
#' @return \code{sbdolphin} the \link{sbdolphin} data set
#' @examples \\dontrun{sbdolphin = io_sbdolphin.pkgdata.load();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.pkgdata.load = function()
{ 
  effort = io_sbdolphin.effort(rem.intermediate=FALSE)
  par3d.args = io_sbdolphin.par3d()
  dset = list(sighting=io_sbdolphin.sighting(BfLim=5,PDLim=6),
                transect = as.transect(effort),
                effort = effort,
                mesh = io_sbdolphin.mesh(),
                coast.boundary = io_star.coast(),
                sea.boundary = io_star.boundary(),
                inner.boundary = io_sbdolphin.inner.boundary(),
                geometry = "geo",
                mesh.coords = c("lon","lat"),
                time.coords = "year",
                par3d.args = par3d.args,
                ips = io_sbdolphin.ips(),
                ips.yearly = io_sbdolphin.ips.yearly()
  )
  
  # Convert into dsdata format
  cat("Converting ETP data into dsdata. This might take a couple of minutes.")
  dset$effort = as.effort.etpeffort(dset$effort, dset$sighting)
  dset$sighting = NULL
  class(dset) = c("dsdata", "list")
  
  # Sort transects and segments
  srt = sort(segment.id(dset$effort), index.return = TRUE)
  dset$effort = dset$effort[srt$ix,]
  srt = sort(transect.id(dset$effort), index.return = TRUE)
  dset$effort = dset$effort[srt$ix,]
  
  # Discard effort outside the mesh
  msk = is.inside(dset$mesh, dset$effort, mesh.coords = paste0("start.", dset$mesh.coords)) & is.inside(dset$mesh, dset$effort, mesh.coords = paste0("end.", dset$mesh.coords))
  dset$effort = dset$effort[msk,]
  
  # Attach CRS projection strings
  dset$p4s = "+proj=longlat"
  dset$mesh.p4s = "+proj=longlat"
  
  return(dset)
}

#' Regenerate \link{sbdolphin} data and store it to \code{sbdolphin.RData}
#' 
#' Uses \code{\link{io_sbdolphin.pkgdata.load}} to load the data and stores
#' the result to sbdolphin.RData. Thereby the data that is distributed with 
#' our package is generated.
#'
#' @aliases io_sbdolphin.pkgdata.save
#' @export
#' @return NULL
#' @examples \\dontrun{io_sbdolphin.pkgdata.save();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.pkgdata.save = function(){
  ## save the data we will include in the R package
  sbdolphin = io_sbdolphin.pkgdata.load()
  save(sbdolphin,file=paste0(io_sbdolphin.getDataDir(),"/sbdolphin.RData"))
}

#' Load \link{sbdolphin} integration points
#' 
#'
#' @aliases io_sbdolphin.ips
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.ips = function(){
  load(file= paste(io_sbdolphin.getDataDir(), "etp.ips.RData",sep=.Platform$file.sep))
  return(ips)
}

#' Load \link{sbdolphin} integration points (yearly)
#' 
#'
#' @aliases io_sbdolphin.ips.yearly
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.ips.yearly = function(){
  load(file= paste(io_sbdolphin.getDataDir(), "etp.ips.yearly.RData",sep=.Platform$file.sep))
  return(ips.yearly)
}


#' Load \link{sbdolphin} inner boundary
#' 
#'
#' @aliases io_sbdolphin.inner.boundary
#' @export
#' @return boundary inner boundary
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.inner.boundary = function(){
  load(file= paste(io_sbdolphin.getDataDir(), "etp.inner.boundary.Rdata",sep=.Platform$file.sep))
  return(innerBnd)
}

#' Load \link{sbdolphin} par3d parameters
#' 
#'
#' @aliases io_sbdolphin.par3d
#' @export
#' @return params pard3d() parameters for rgl plot of sbdolphin
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.par3d = function(){
  #par3d(zoom=0.62,FOV=1)
  #par3d.args = par3d(no.readonly=TRUE)
  #save(par3d.args,file="data/sbdolphin.par3d.RData")
  load(file= paste(io_sbdolphin.getDataDir(), "etp.par3d.RData",sep=.Platform$file.sep))
  params = par3d.args
  return(params)
}


#' Load \link{sbdolphin} sightings
#' 
#' Load whale sightings from "Bmus86_06.csv" and apply filters (see parameters).
#'
#' @aliases io_sbdolphin.sighting
#' @export
#' @return sightings Whale sighting
#' @examples \\dontrun{sightings = io_sbdolphin.sightings(years=c(2000:2006),E=1,BfLim=5);}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.sighting = function(years=NULL,E=NULL,BfLim=NULL,PDLim=NULL)
{ 
  # Load the file
  sightdat = read.csv(file= paste(io_sbdolphin.getDataDir(), "sbdolphin.sighting.csv",sep=.Platform$file.sep),head=TRUE,sep=",")

  # Remove NA lines
  sightdat = sightdat[!is.na(sightdat$perpdist),]
  
  # rename
  names(sightdat)[names(sightdat) %in% "perpdist"] = "PD"
  
  # Filter out on.eff=0
  sightdat = sightdat[sightdat$on.eff==1, ]
  
  # Year selection
  if (is.numeric(years)) {
    #print("io_sbdolphin.sightings(): Filtering years")
    sightdat = sightdat[sightdat$Year %in% years, ]
  }
  
  # Get rid of the non-LT-mode sightings
  if (is.numeric(E)) {
    #print("io_sbdolphin.sightings(): Filtering E")
    sightdat = sightdat[sightdat$E==E, ]
  }
  
  # Tim commented that it is not LT mode if Bft>5
  if (is.numeric(BfLim)) {
    #print("io_sbdolphin.sightings(): Filtering Bf")
    sightdat = sightdat[sightdat$bft<=BfLim,]
  } 
  
  # Distance limit
  if (is.numeric(PDLim)) {
    #print("io_sbdolphin.sightings(): Filtering PD")
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


#' Generate \link{sbdolphin} transect lines from ETP effort data
#'
#' @aliases io_sbdolphin.transect
#' @export
#' @return tr Data frame with transect lines
#' @examples \\dontrun{transect = io_sbdolphin.transect();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'


io_sbdolphin.transect = function(years=FALSE,rem.intermediate=TRUE)
{ 
  eff = io_sbdolphin.effort(years=years,rem.intermediate=rem.intermediate)
  
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

#' Load \link{sbdolphin} effort
#' 
#' Load effort data from "ETP effort.csv"
#'
#' @aliases io_sbdolphin.effort
#' @export
#' @return effort ETP effort related to \link{sbdolphin} data
#' @examples \\dontrun{transects = io_sbdolphin.transects();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.effort = function(years=FALSE,rem.intermediate=TRUE)
{ 
  # read effort data
  tr = read.csv(file=paste(io_sbdolphin.getDataDir(),"ETP effort.csv",sep=.Platform$file.sep))
  
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


#' Load \link{sbdolphin} boundary
#' 
#' @aliases io_sbdolphin.boundary
#' @export
#' @return boundary 
#' @examples \\dontrun{boundary = io_sbdolphin.boundary();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.boundary = function(...) {return(io_star.boundary(...))}

#' Load \link{sbdolphin} coast
#' 
#' @aliases io_sbdolphin.coast
#' @export
#' @return coast
#' @examples \\dontrun{coast = io_sbdolphin.coast();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.coast = function(...) {return(io_star.coast(...))}

#' Load \link{sbdolphin} mesh
#' 
#' @aliases io_sbdolphin.mesh
#' @export
#' @return mesh An INLA mesh
#' @examples \\dontrun{mesh = io_sbdolphin.mesh(); plot(mesh)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sbdolphin.mesh = function(...) 
{
  load(file=paste(io_sbdolphin.getDataDir(), "mesh_stitchInterior.Rdata",sep="/"))
  return(mesh)
}
