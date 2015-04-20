# INTERNAL DATA STORAGE
io_sst.getDataDir = function() {return(system.file("data",package="iDistance"))}


#' Load sea surface temperature (SST) from SODA data
#'
#' The SODA (Simple Ocean Data Assimilation) model is summarized at http://apdrc.soest.hawaii.edu/datadoc/soda_2.2.4.php
#' the format of the URL query is described at http://coastwatch.pfeg.noaa.gov/erddap/griddap/documentation.html
#' The SST data may also be downloaded using the interface at http://coastwatch.pfeg.noaa.gov/erddap/griddap/hawaii_d90f_20ee_c4cb.html.
#' Specify time (year and months 01 to 12), latitude (-30.25 to 40.25), longitude (190.25 to 290.75), variable (sst) and file type (NetCDF),
#' then click ?Submit?.  This procedure should produce the same URL query and same NetCDF file as the code below.

#' @aliases io_sst.load
#' @export
#' @param year Vector of years to return SST for
#' @param month Vector of months to select from the SST data sets
#' @param source.dir Direcory to load SST NetCDF files from (default: system temporary directory)
#' @return \code{sst} sea surface temperature
#' @examples \\dontrun{}
#' @author Tim Gerrodette Oct 2014 <\email{tim.gerrodette@@noaa.gov}>
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

io_sst.load = function(year=2006,month=1,source.dir=tempdir(),lat.stride=1,lon.stride=1,lat.range=c(-30.25,40.25),lon.range=c(-169.25,-69.25)){
  require(RNetCDF)
  # We have to convert the given lon.range to the [0,360]-system used for SODA data
  converted.lon.range = lon.range
  converted.lon.range[lon.range<0] = converted.lon.range[lon.range<0]+360
  sst.list = list()
  for (k in 1:length(year)){
    yr = year[k]
    filename = paste("soda-sst-",as.character(yr),"_lo",as.character(converted.lon.range[1]),"-",as.character(converted.lon.range[2]),"_la",as.character(lat.range[1]),"-",as.character(lat.range[2]),"_str",as.character(lat.stride),"-",as.character(lon.stride),".nc",sep="")
    full.filename = file.path(source.dir,"..",filename)
    if (!file.exists(full.filename)){
      io_sst.download(target.dir=source.dir,year=yr,lat.range=c(-30.25,40.25),lon.range=converted.lon.range,lat.stride=lat.stride,lon.stride=lon.stride,filename=full.filename)
    }
    sst.nc <- open.nc(full.filename)
    month.idx <- utcal.nc("seconds since 1970-01-01",var.get.nc(sst.nc,0))[,2] # convert the NetCDF object into more conventional R objects
    
    sst = read.nc(sst.nc)
    close.nc(sst.nc)
    
    # select months
    if (is.numeric(month)) { msel = month }
    else { msel = month[[k]] }
    sst$temp = sst$temp[,,msel] 
    sst$month <- as.numeric(month.idx[msel])
    
    # annotate with year
    sst$year = yr
    
    # change coordinates system
    sst$longitude[sst$longitude>180] = sst$longitude[sst$longitude>180] - 360
    lola.idx = which(names(sst) %in% c("latitude","longitude"))
    names(sst)[lola.idx] = c("lat","lon")
    
    #append to list
    sst.list[[k]] = sst
  }
  #sst <- array(abind::abind(lapply(sst.list,function(x) {x$temp}),along=4),
  #             dim=c(length(sst.list[[1]]$longitude),length(sst.list[[1]]$latitude),length(sst.list[[1]]$time),length(year)),
  #             dimnames=list(lon=sst.list[[1]]$longitude,lat=sst.list[[1]]$latitude,month=month.idx,year=year))
  #sst = sst[,,month,]
  #class(sst) = c("sst","array")  
  return(sst.list)
}

#' Download sea surface temperature (SST) data in NetCDF format
#' 
#' Create a formatted URL query to download SODA monthly mean SST data from Univ of Hawaii Asia-Pacific Data Research Center server, one year at a time
#' in code below, specify the years of the data desired, then execute the other R commands.
#' the SODA model produces estimates at 0.5 degree spatial resolution for both longitude and latitude
#' data are written to files named sst-YYYY.nc in the provided directory in NetCDF format, where YYYY is the year requested
#'
#' 
#' @aliases io_sst.download(dir)
#' @export
#' @param target.dir target directory for NetCDF SST files
#' @param year Vector of years to load SST data for
#' @author Tim Gerrodette Oct 2014 <\email{tim.gerrodette@@noaa.gov}>
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_sst.download = function(target.dir=tempdir(),year=2006,lat.range=c(-30.25,40.25),lon.range=c(190.25,290.75),lon.stride=1,lat.stride,filename=NA) {
  #' year <- c(1986:1990,1992,1993,1998:2000,2003,2006)
  for (yr in year) {
    url.query <- paste("http://coastwatch.pfeg.noaa.gov/erddap/griddap/hawaii_d90f_20ee_c4cb.nc?temp[(",as.character(yr),"-01-15T00:00:00Z):1:(",as.character(yr),"-12-15T00:00:00Z)][(5.01):1:(5.01)][(",as.character(lat.range[1]),"):",as.character(lat.stride),":(",as.character(lat.range[2]),")][(",as.character(lon.range[1]),"):",as.character(lon.stride),":(",as.character(lon.range[2]),")]",sep="")
    if (is.na(filename)) {
      filename = paste("soda-sst-",as.character(yr),".nc",sep="")
      ncfile <- file.path(target.dir,"..",filename)
    }
    else {
      ncfile = filename
    }
    download.file(url.query,ncfile,mode="wb")
  }
}

#' Save SST data linearized on mesh triangles
#'
#' @aliases io_sst.save
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

io_sst.save = function(year=1986:2007,month=6:12){
  load(file=paste(io_sst.getDataDir(), "mesh_stitchInterior.Rdata",sep="/")) # load mesh
  isst = list()
  k=1
  for (y in year){
    sst = io_sst.load(year=year,month=month)
    msst = sst
    msst[[1]]$temp = apply(sst[[1]]$temp,MARGIN=c(1,2),function(x) {return(mean(x))})
    
    loc = data.frame(data$mesh$loc[,c(1,2)])
    colnames(loc) = data$mesh.coords
    isst[[k]] = interpolate.sst.internal(msst[[1]],loc,extrapolate=TRUE)
    k = k+1
  }
  sst = list()
  sst$values = do.call(cbind,isst)
  colnames(sst$values) = year
  sst$mesh = mesh
  sst$time.coords = "year"
  sst$mesh.coords = c("lon","lat")

  save(sst,file=paste0(io_sst.getDataDir(),"/sst.RData"))
}
