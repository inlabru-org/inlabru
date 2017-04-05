# Sea surface temperature (SST) from SODA data
#
# The SODA (Simple Ocean Data Assimilation) model is summarized at http://apdrc.soest.hawaii.edu/datadoc/soda_2.2.4.php
# The data can be accessed via \link{io_sst.load}, which internally calls \link{io_sst.load} if the SODA NetCDF files can not be found on the local machine.
#
# SST data can be plotted via \link{plot.sst} (see also the example)
# Other methods: \link{mean.sst},\link{interpolate.sst}
#
# @aliases plot.sst
# @examples \\dontrun{sst = io_sst.load(year=c(2006,2007)); plot(sst)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
# @name sst
NULL


# Plot SST data
#
# @aliases plot.sst
# @export
# @examples \\dontrun{sst = io_sst.load(year=c(2006,2007),month=c(5,6)); plot(sst)}
# @author Tim Gerrodette Oct 2014 <\email{tim.gerrodette@@noaa.gov}>
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
# 

plot.sst = function(sst){
  warning("This plot code is merely a placeholder.")
  dims = dim(sst)
  par(mfrow=c(dims[4],dims[3]),mai=c(0.1,0.1,0.1,0.1))
  for (i in 1:dims[4]){
    for (k in 1:dims[3]) {
      image(sst[,,k,i],col=gray(0:20/20),xaxt='n', yaxt='n',main=month.name[as.numeric(dimnames(sst)$month[k])]) # ]
    }
  }
}


# Mean of sea surface temperature (SST)
# 
# @aliases mean.sst
# @export
# @param sst sea surface temperature
# @param \code{mweights} Monthwise averaging weights (must equal the number of months \code{sst} incorporates)
# @return \code{sst} weighted yearly mean sea surface temperature
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#

mean.sst = function(sst,mweights=rep(1,dim(sst)[3])) {
  wfun = function(x,w) {weighted.mean(x,mweights)}
  sst = apply(sst,MARGIN=c(1,2,4),FUN=wfun)
  class(sst) = c("sst","array")
  return(sst)
}


# Interpolate sea surface temperature (SST)
# 
# @aliases interpolate.sst
# @export
# @param sst Spatial SST estimate or list of such
# @param loc Locations to interpolate SST at
# @return sst.at.loc SST at the provided locations
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

interpolate.sst = function(sst,loc) {
  if (length(dim(sst$temp))>2){ stop(print("This method is not implemented to work with temporal SST slabs.")) }
  if (length(names(sst))==0){
    # We got a list of SST fields
    sst.at.loc = lapply(sst,function(x) interpolate.sst.internal(x,loc))
  }
  else {
    sst.at.loc = interpolate.sst.internal(sst,loc)
  }
  
return(sst.at.loc)
}

interpolate.sst.internal = function(sst,loc,extrapolate=TRUE) {
  
  # Define lattice and mesh
  lattice = inla.mesh.lattice(x=sst$lon,y=sst$lat)
  mesh = inla.mesh.create(lattice=lattice,extend=list(n=5),boundary=lattice$segm)
  
  # extrapolate SST at locations with NA values
  if (extrapolate) { sst = extrapolate.sst(sst,mesh) }
  
  # SST as vector (grid indexing)
  sst.on.grid = as.vector(sst$temp)
  
  # Project to values at the locations defined above
  A = inla.spde.make.A(mesh,loc=as.matrix(loc[,c("lon","lat")]))
  inv.idx <- rep(NA, mesh$n)
  inv.idx[mesh$idx$lattice] = seq_len(mesh$n)
  sst.at.loc = as.vector(A %*% sst.on.grid[inv.idx])
  
  return(sst.at.loc)
}


# Extract sea surface temperature (SST) for given locations (and time)
# 
# @aliases get.sst
# @export
# @param loc Locations to interpolate SST at
# @return sst.at.loc SST at the provided locations
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.sst = function(loc,method=yearly.mean.sst,...) {
  if (is.character(method)){
    if (method=="ymean"){ method = yearly.mean.sst }
  }
  return(method(loc,...))
}

yearly.mean.sst = function(loc,month=6:12,...){
  year = unique(loc$year)
  sst = io_sst.load(year=year,month=month,...)
  msst = lapply(sst,function(x){m = apply(x$temp,MARGIN=c(1,2),mean);x$temp=m;return(x)})
  isst = rep(NA,dim(loc)[1])
  for (ny in 1:length(msst)){
    y = year[ny]
    idx = loc$year==y
    if (any(idx)){
      isst[idx] = interpolate.sst(msst[[ny]],loc[idx,])
    }  
  }
  return(isst)
}

yearly.dev.sst = function(loc,month=6:12,...){
  year = unique(loc$year)
  sst = io_sst.load(year=year,month=month,...)
  msst = lapply(sst,function(x){m = apply(x$temp,MARGIN=c(1,2),mean);x$temp=m;return(x)})
  temp = lapply(msst,function(x){return(x$temp)})
  arr = array(unlist(temp),dim=c(dim(temp[[1]]),length(temp)))
  tmean = apply(arr,MARGIN=c(1,2),mean)
  msst = lapply(msst,function(x){m = x$temp - tmean;x$temp=m;return(x)})
  
  isst = rep(NA,dim(loc)[1])
  for (ny in 1:length(msst)){
    y = year[ny]
    idx = loc$year==y
    if (any(idx)){
      isst[idx] = interpolate.sst(msst[[ny]],loc[idx,])
    }  
  }
  return(isst)
}

yearly.sdev.sst = function(loc,month=6:12,...){
  year = unique(loc$year)
  sst = io_sst.load(year=year,month=month,...)
  msst = lapply(sst,function(x){
    m = apply(x$temp,MARGIN=c(1,2),mean)
    x$temp=m -mean(m,na.rm=TRUE)
    return(x)})
  
  isst = rep(NA,dim(loc)[1])
  for (ny in 1:length(msst)){
    y = year[ny]
    idx = loc$year==y
    if (any(idx)){
      isst[idx] = interpolate.sst(msst[[ny]],loc[idx,])
    }  
  }
  return(isst)
}



# Yearly weighted average SST for given locations
#
# @aliases yearly.wmean.sst
# @export
# @param loc Geographic locations annotated with year, day and month
# @param wloc If defined, locations to determine weights from
# @return sst SST
# @examples \\dontrun{}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
# 
yearly.wmean.sst = function(loc,wloc=FALSE,...){
  year = unique(loc$year)
  if (is.data.frame(wloc)){ weights = weights.sst(wloc) }
  else { weights = weights.sst(loc) }
  
  # interpolated SST
  isst = rep(NA,dim(loc)[1])
  
  for (ny in 1:length(year)){
    y = year[ny]
    sst = io_sst.load(year=y,month=as.numeric(colnames(weights)),...)
    w = weights[as.character(y),]
    wmsst = apply(sst[[1]]$temp,MARGIN=c(1,2),function(x) {return(weighted.mean(x,w))})
    sst[[1]]$temp = wmsst
    idx = loc$year==y
    if (any(idx)){
      isst[idx] = interpolate.sst(sst[[1]],loc[idx,])
    }  
  }
  return(isst)
}


# Weights for yearly averages of SST data
#
# @aliases weights.sst
# @export
# @param effdat effort data
# @return sstweights Weights for each month and each year of the given effort data
# @examples \\dontrun{}
# @author Yuan (Joyce) Yuan <\email{yy84@@st-andrews.ac.uk}>
# 

weights.sst = function(effdat)
{
  yr <- unique(effdat$year[!is.na(effdat$year)]) ## all the yrs with effort
  month.min <- sapply(yr,function(i)min(effdat$month[effdat$year==i],na.rm=T)) # months of first effort by year
  day.min <- sapply(yr,function(i) min(effdat$day[effdat$year==i&effdat$month==month.min[which(yr==i)]],na.rm=T))
  month.max <- sapply(yr,function(i)max(effdat$month[effdat$year==i],na.rm=T)) # months of last effort by year
  day.max <- sapply(yr,function(i) max(effdat$day[effdat$year==i&effdat$month==month.max[which(yr==i)]],na.rm=T))
  effortdaycount = data.frame(year = yr, monthstart = month.min, daystart =day.min, monthend = month.max, dayend = day.max)
  #for the year with effort
  monthall = c(min(month.min):max(month.max))
  days.effortyr = matrix(data=NA, nrow=length(yr), ncol=length(monthall), dimnames= list(yr, monthall))
  ## days for each month, note that the start and end month differ by year
  days.effortyr[effortdaycount$monthstart!=7, "7"]= 0
  days.effortyr[effortdaycount$monthstart!=7, "8"]= 31-effortdaycount[effortdaycount$monthstart==8,]$daystart +1
  days.effortyr[effortdaycount$monthstart==7, "7"]=31 - effortdaycount[effortdaycount$monthstart==7,]$daystart +1
  days.effortyr[effortdaycount$monthstart==7, "8"]=31
  days.effortyr[ , "9"] =30
  days.effortyr[ , "10"] =31
  days.effortyr[effortdaycount$monthend!=12, "12"]= 0
  days.effortyr[effortdaycount$monthend!=12, "11"]= effortdaycount[effortdaycount$monthend!=12,]$dayend
  days.effortyr[effortdaycount$monthend==12, "11"]= 30
  days.effortyr[effortdaycount$monthend==12, "12"]= effortdaycount[effortdaycount$monthend==12,]$dayend
  # transform days count into weights
  total = apply(days.effortyr,1, sum)
  weights.effortyr = t(apply(matrix(c(1:nrow(days.effortyr)), ncol=1),1, function(X)days.effortyr[X,]/total[X]))
  rownames(weights.effortyr) = yr
  ## weights for the gap year
  ## for those gap yrs with no effort, Tim has agreed with start day in July is 30 and end day in Dec is 6
  yrall = c(range(effdat$year)[1]:range(effdat$year)[2])
  yrgap = setdiff(yrall, yr)
  days.gapyr = matrix(data=NA, nrow = length(yrgap), ncol=length(monthall), dimnames = list(yrgap,monthall))
  days.gapyr[, as.character(c(8:11))] = matrix(rep(c(31,30,31,30), each=nrow(days.gapyr)), nrow= nrow(days.gapyr), byrow=F)
  days.gapyr[,as.character(7)] = rep(2, nrow(days.gapyr))
  days.gapyr[,as.character(12)] = rep(6,nrow(days.gapyr))
  days.gapyr ## for gap years, days are all the same for each month across years
  weights.gapyr = days.gapyr/sum(days.gapyr[1,])
  ## combine yreffort and yrgap together
  sstweights = rbind(weights.gapyr, weights.effortyr )
  ## make the rows in order with all the years
  sstweights = sstweights[as.character(yrall),]
  return(sstweights)
}


# Extrapolate sea surface temperature (SST) at regions where the SODA defines the sst to be NA
# 
# @aliases extrapolate.sst
# @export
# @param sst Spatial SST estimates with NA values
# @param mesh inla.mesh to extrapolate values on
# @return sst Spatial SST with filled in conditional expectations for NA values
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>


extrapolate.sst = function(sst,mesh) {
  # SST as vector
  sst.on.grid = as.vector(sst$temp)  

  inv.idx <- rep(NA, mesh$n)
  inv.idx[mesh$idx$lattice] = seq_len(mesh$n)
  
  sst.on.mesh = sst.on.grid[inv.idx]  
  na.vec = is.na(sst.on.mesh)
  
  u.msk = na.vec # mask of variable to compute conditional expectation for
  c.msk = !(u.msk) # mask of variables to condition on

  # precision matrices
  Q = inla.spde.precision(inla.spde2.matern(mesh),theta=c(0, log(1e-6)))
  Q_uu = Q[u.msk,u.msk] 
  Q_cc = Q[c.msk,c.msk]
  Q_uc = Q[u.msk,c.msk]  
  
  # Given values and marginal mean
  x_c = sst.on.mesh[!na.vec] # values of variables to condition on
  mu_c = rep(0,length(x_c)) # marginal mean 
  mu_u = rep(0,sum(u.msk)) # marginal mean
  
  # compute conditional expectation mu_uc via Cholesky decomposition
  z = Q_uu%*%mu_u - Q_uc %*% (x_c-mu_c)
  L_uu = chol(Q_uu) # upper triangular!
  mu_uc = as.vector(solve(L_uu,solve(t(L_uu),as.vector(z))))
  
  # fill in values
  sst.on.mesh[na.vec] = mu_uc
  
  # transform back to grid indices
  sst$temp = matrix(sst.on.mesh[mesh$idx$lattice],nrow = nrow(sst$temp),ncol = ncol(sst$temp))
  
  # return
  return(sst)
}