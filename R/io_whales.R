# INTERNAL DATA STORAGE
io_whales.getDataDir = function() {return(system.file("data",package="iDistance"))}


#' Load \link{whales} survey data from raw data sets
#'
#' @aliases io_whales.pkgdata.load
#' @export
#' @return \code{whales} the \link{whales} data set
#' @examples \\dontrun{whales = io_whales.pkgdata.load();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.pkgdata.load = function()
{ 
  effort = io_whales.effort(rem.intermediate=FALSE)
  par3d.args = io_whales.par3d()
  whales = list(sighting=io_whales.sighting(BfLim=5,PDLim=6),
                transect = as.transect(effort),
                effort = effort,
                mesh = io_whales.mesh(),
                coast.boundary = io_star.coast(),
                sea.boundary = io_star.boundary(),
                inner.boundary = io_whales.inner.boundary(),
                geometry = "geo",
                mesh.coords = c("lon","lat"),
                time.coords = "year",
                par3d.args = par3d.args,
                ips = io_whales.ips(),
                ips.yearly = io_whales.ips.yearly()
  )
  class(whales) = c("whales","etpdata","dsdata","list")
  return(whales)
}

#' Regenerate \link{whales} data and store it to \code{whales.RData}
#' 
#' Uses \code{\link{io_whales.pkgdata.load}} to load the data and stores
#' the result to whales.RData. Thereby the data that is distributed with 
#' our package is generated.
#'
#' @aliases io_whales.pkgdata.save
#' @export
#' @return NULL
#' @examples \\dontrun{io_whales.pkgdata.save();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.pkgdata.save = function(){
  ## save the data we will include in the R package
  whales = io_whales.pkgdata.load()
  save(whales,file=paste0(io_whales.getDataDir(),"/whales.RData"))
}

#' Load \link{whales} integration points
#' 
#'
#' @aliases io_whales.ips
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.ips = function(){
  load(file= paste(io_whales.getDataDir(), "bwhales.ips.RData",sep=.Platform$file.sep))
  return(ips)
}

#' Save \link{whales} integration points
#' 
#'
#' @aliases io_whales.ips.save
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.ips.save = function(){
  data(whales)
  ips = projection.integration(data, idx = idx, 
                               mesh = mesh, mesh.split=FALSE, mesh.coords = c("lon","lat"), geometry = "geo", 
                               truncation = 6, n.distance = 5, distance.scheme ="equidistant", fake.distance=TRUE,
                               transect.scheme = "gaussian", n.transect = 1, projection = "linear",
                               group.by = NULL)
  
  save(ips,file= paste(io_whales.getDataDir(), "bwhales.ips.RData",sep=.Platform$file.sep))
  return(ips)
}

#' Load \link{whales} integration points (yearly)
#' 
#'
#' @aliases io_whales.ips.yearly
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.ips.yearly = function(){
  load(file= paste(io_whales.getDataDir(), "bwhales.ips.yearly.RData",sep=.Platform$file.sep))
  return(ips.yearly)
}


#' Save \link{whales} integration points (yearly)
#' 
#'
#' @aliases io_whales.ips.yearly.save
#' @export
#' @return ips integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.ips.yearly.save = function(){
  data(whales)
  ips.yearly = projection.integration(data, idx = idx, 
                               mesh = mesh, mesh.split=FALSE, mesh.coords = c("lon","lat"), geometry = "geo", 
                               truncation = 6, n.distance = 5, distance.scheme ="equidistant", fake.distance=TRUE,
                               transect.scheme = "gaussian", n.transect = 1, projection = "linear",
                               group.by = "year")
  
  save(ips.yearly,file= paste(io_whales.getDataDir(), "bwhales.ips.yearly.RData",sep=.Platform$file.sep))
  return(ips)
}


#' Load \link{whales} inner boundary
#' 
#'
#' @aliases io_whales.inner.boundary
#' @export
#' @return boundary inner boundary
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.inner.boundary = function(){
  load(file= paste(io_whales.getDataDir(), "bwhales.inner.boundary.Rdata",sep=.Platform$file.sep))
  return(innerBnd)
}

#' Load \link{whales} par3d parameters
#' 
#'
#' @aliases io_whales.par3d
#' @export
#' @return params pard3d() parameters for rgl plot of whales
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.par3d = function(){
  #par3d(zoom=0.62,FOV=1)
  #par3d.args = par3d(no.readonly=TRUE)
  #save(par3d.args,file="data/bwhales.par3d.RData")
  load(file= paste(io_whales.getDataDir(), "bwhales.par3d.RData",sep=.Platform$file.sep))
  params = par3d.args
  return(params)
}
  

#' Load \link{whales} sightings
#' 
#' Load whale sightings from "Bmus86_06.csv" and apply filters (see parameters).
#'
#' @aliases io_whales.sighting
#' @export
#' @return sightings Whale sighting
#' @examples \\dontrun{sightings = io_whales.sightings(years=c(2000:2006),E=1,BfLim=5);}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.sighting = function(years=NULL,E=NULL,BfLim=NULL,PDLim=NULL)
{ 
  # Load the file
  sightdat = read.csv(file= paste(io_whales.getDataDir(), "bluewhales.sighting.csv",sep=.Platform$file.sep),head=TRUE,sep=",")
  # Columns prior to 2007 update:
  #[1] "Sp1"      "Sp2"      "Sp3"      "Sp4"      "Cruz"     "Sght"     "E"        "Obs"      "year"    
  #[10] "month"    "day"      "time"     "lat"      "lon"      "GS"       "TotGS"    "distance" "Bf"      
  #[19] "SH"       "SD"       "RF"       "HS"       "VS"       "WSp"      "WDi"      "Cue"      "Me"      
  #[28] "Ph"       "Bi"       "Mx"       "Crs"      "SST"      "Vis"      "Angl"     "Retcl"    "RadD"    
  #[37] "InitID"   "MSp"
  
  # rename
  names(sightdat)[names(sightdat) %in% "perpdist"] = "PD"
  
  # Filter out on.eff=0
  sightdat = sightdat[sightdat$on.eff==1, ]
  
  # Year selection
  if (is.numeric(years)) {
    #print("io_whales.sightings(): Filtering years")
    sightdat = sightdat[sightdat$Year %in% years, ]
  }
  
  # Get rid of the non-LT-mode sightings
  if (is.numeric(E)) {
    #print("io_whales.sightings(): Filtering E")
    sightdat = sightdat[sightdat$E==E, ]
  }
  
  # Tim commented that it is not LT mode if Bft>5
  if (is.numeric(BfLim)) {
    #print("io_whales.sightings(): Filtering Bf")
    sightdat = sightdat[sightdat$bft<=BfLim,]
  } 
  
  # Distance limit
  if (is.numeric(PDLim)) {
    #print("io_whales.sightings(): Filtering PD")
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


#' Generate \link{whales} transect lines from ETP effort data
#'
#' @aliases io_whales.transect
#' @export
#' @return tr Data frame with transect lines
#' @examples \\dontrun{transect = io_whales.transect();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'


io_whales.transect = function(years=FALSE,rem.intermediate=TRUE)
{ 
  eff = io_whales.effort(years=years,rem.intermediate=rem.intermediate)

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

#' Load \link{whales} effort
#' 
#' Load effort data from "ETP effort.csv"
#'
#' @aliases io_whales.effort
#' @export
#' @return effort ETP effort related to \link{whales} data
#' @examples \\dontrun{transects = io_whales.transects();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.effort = function(years=FALSE,rem.intermediate=TRUE)
{ 
  # read effort data
  tr = read.csv(file=paste(io_whales.getDataDir(),"ETP effort.csv",sep=.Platform$file.sep))

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


#' Load \link{whales} boundary
#' 
#' @aliases io_whales.boundary
#' @export
#' @return boundary 
#' @examples \\dontrun{boundary = io_whales.boundary();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.boundary = function(...) {return(io_star.boundary(...))}

#' Load \link{whales} coast
#' 
#' @aliases io_whales.coast
#' @export
#' @return coast
#' @examples \\dontrun{coast = io_whales.coast();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.coast = function(...) {return(io_star.coast(...))}

#' Load \link{whales} mesh
#' 
#' @aliases io_whales.mesh
#' @export
#' @return mesh An INLA mesh
#' @examples \\dontrun{mesh = io_whales.mesh(); plot(mesh)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_whales.mesh = function(...) 
{
  load(file=paste(io_whales.getDataDir(), "mesh_stitchInterior.Rdata",sep="/"))
  return(mesh)
}
