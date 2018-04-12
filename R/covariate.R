# Create covariate data set (covdata)
#
# @aliases make.covdata
# @name make.covdata
# @export
# @param mesh
# @param values
# @param mesh.coords
# @param time.coords
# @examples \\dontrun{}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

make.covdata = function(mesh, values, mesh.coords, time.coords){
  covdata = list(mesh=mesh, values=values, mesh.coords=mesh.coords, time.coords=time.coords)
  class(covdata) = c("covdata","list")
  return(covdata)
}

# Create covariate data set (covdata)
#
# Note: documentation does not match code
#
# @aliases covdata.import
# @name covdata.import
# @export
# @param mesh
# @param values
# @param mesh.coords
# @param time.coords
# @examples \\dontrun{}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

covdata.import = function(dframe, colname, data){
  
  covloc = as.matrix(dframe[,data$mesh.coords])
  scale = max(abs(dframe[,colname]))
  dframe[,colname] = dframe[,colname]/scale
  
  spde.mdl = INLA::inla.spde2.matern(mesh = data$mesh, alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.1)
  A = INLA::inla.spde.make.A(data$mesh, loc = covloc)
  stack = INLA::inla.stack(data=list(y = dframe[,colname]), A=list(A,1), 
                     effects = list(spde=1:spde.mdl$n.spde, m = rep(1,nrow(dframe))))
  
  result = INLA::inla( formula = y ~ f(spde, model = spde.mdl) -1, 
                       family = "gaussian",
                       data = INLA::inla.stack.data(stack),
                       control.predictor = list(A = INLA::inla.stack.A(stack)),
                       verbose = FALSE)
        
  value = scale * result$summary.random$spde[,"mode", drop = FALSE]
  
  depth = make.covdata(mesh = data$mesh, values = value, mesh.coords = data$mesh.coords, time.coords = NULL)
}


# Plot covariate data set
#
# Note: documentation does not match code
#
# @aliases plot.covdata
# @export
# @param covdata A covariate data set
# @examples \\dontrun{data(sst); plot.covariate(sst, time = 1)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

plot.covdata = function(covdata, time = 1, fun=NULL, ...){
  if (is.null(covdata$geometry)) { covdata$geometry = "euc" }
  col = covdata$values[, time]
  pixelplot.mesh(data = covdata, col = col, ...)
}


# Make covariate function
#
# Note: documentation does not match code
#
# @aliases make.covariate
# @name make.covariate
# @export
# @param covdata A covariate data set
# @return A function intended to act as a covariate for a distance sampling model
# @examples \\dontrun{data(sst); plot.covariate(sst, time = 1)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

make.covariate = function(cdata, method = NULL, ...){

    if ( class(cdata)[1] == "covdata") { 
      if ( is.null(method) ) { method = get.value }
      return( function(loc){ method(cdata, loc) } )
    }  

    else if (class(cdata)[1] == "SpatialPolygonsDataFrame") {
      if ( is.null(method) ) { method = shapefile.to.covariate }
      return( method(cdata, ...) )
    }
  
}



# Get covariate value at given location/time
#
# @aliases get.value
# @export
# @param covariate A covariate data set
# @param loc Locations (spatial/temporal) to obtain the covariate value at
# @return \code{values} Covariate values
# @examples \\dontrun{data(sst); get.value(sst,loc=data.frame(lon=-110,lat=0,year=1986))}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.value = function(covariate,loc){
  times = unique(loc[,covariate$time.coords])
  values = numeric(length=nrow(loc))
  if (length(times)>0) {
    for (t in times){
      msk = loc[,covariate$time.coords]==t
      A = INLA::inla.spde.make.A(covariate$mesh,loc=as.matrix(loc[msk,covariate$mesh.coords]))
      inside = is.inside(covariate$mesh, loc = loc[,covariate$mesh.coords], mesh.coords = covariate$mesh.coords)
      values[msk] = as.vector(A%*%covariate$values[,as.character(t)])
      values[msk & !inside] = NA
    }
  } else {
    A = INLA::inla.spde.make.A(covariate$mesh,loc=as.matrix(loc[,covariate$mesh.coords]))
    inside = is.inside(covariate$mesh, loc = loc[,covariate$mesh.coords], mesh.coords = covariate$mesh.coords)
    values = as.vector(A%*%covariate$values[,1])
    values[!inside] = NA
  }
  return(values)
}

# Get maximum value of covariate
#
# @aliases get.max
# @export
# @param covariate A covariate data set
# @param loc Locations (spatial/temporal) to obtain the maximal covariate value at. If not provided, return overall maximum.
# @return \code{values} Maximal covariate values
# @examples \\dontrun{data(sst); get.max(sst)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.max = function(covariate,loc=NULL){
  if (is.null(loc)){ return(max(as.vector(covariate$values))) }
  else {
    vals = get.value(covariate,loc)
    return(max(vals))
  }
}

# Get minimum value of covariate
#
# @aliases get.min
# @export
# @param covariate A covariate data set
# @param loc Locations (spatial/temporal) to obtain the maximal covariate value at. If not provided, return overall maximum.
# @return \code{values} Maximal covariate values
# @examples \\dontrun{data(sst); get.min(sst)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.min = function(covariate,loc=NULL){
  if (is.null(loc)){ return(max(as.vector(covariate$values))) }
  else {
    vals = get.value(covariate,loc)
    return(min(vals))
  }
}


# Get temporal mean
#
# @aliases get.tmean
# @export
# @param covariate A covariate data set
# @param loc Locations (spatial/temporal) to obtain temporal mean covariate value at. If not provided, return overall maximum.
# @return \code{values} Temporal mean
# @examples \\dontrun{data(sst); get.min(sst)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.tmean = function(covariate,loc=NULL){
  tm = apply(covariate$values,MARGIN=1,mean)
  loc2 = loc[,covariate$mesh.coords,drop=FALSE]
  A = INLA::inla.spde.make.A(covariate$mesh,loc=as.matrix(loc2))
  values = as.vector(A%*%tm)
  return(values)
}


# Get spatial mean of covariate for given time
#
# Note: If this function seems to run forever use
# refine = list(max.edge=K) with some integer K>1
#
# WARNING: This is an approximation!
#
#
# @aliases get.smean
# @export
# @param covariate A covariate data set
# @param loc Locations (spatial/temporal) to obtain the spatial mean covariate value at.
# @param weights Averaging weights, default all ones
# @return \code{values} Spatial mean or spatial weighted mean
# @examples \\dontrun{data(sst); get.smean(sst)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.smean = function(covariate,loc=NULL,weights=NULL){
  sm <- numeric(length=ifelse(is.null(loc), 0, dim(loc)[1]))
  
  if (is.null(weights)) {
    weights <- rep(1, covariate$mesh$n)
  }
  weights <- weights * Matrix::diag(INLA::inla.mesh.fem(covariate$mesh, order=1)$c0)
  weights <- weights / sum(weights)
  
  loc2 <- data.frame(covariate$mesh$loc[, 1:length(covariate$mesh.coords),
                                        drop=FALSE])
  colnames(loc2) <- covariate$mesh.coords
  
  if (is.null(covariate$time.coords)){
    values <- get.value(covariate, loc2)
    sm <- sum(values * weights)
  } else {
    for (t in unique(loc[,covariate$time.coords])){
      loc2[[covariate$time.coords]] <- t
      values <- get.value(covariate, loc2)
      sm[loc[covariate$time.coords]==t] <- sum(values * weights)
    }
  }
  return(sm)
}



# Get mean
#
# @aliases get.mean
# @export
# @param covariate A covariate data set
# @param loc Locations (spatial/temporal) to obtain mean covariate value at. If not provided, return overall maximum.
# @param timepoints ????
# @param weights ????
# @return \code{values} Temporal mean
# @examples \\dontrun{data(sst); get.min(sst)}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.mean = function(covariate,loc=NULL, timepoints=NULL, weights=NULL){
  if (is.null(timepoints)) {
    timepoints <- as.numeric(unique(colnames(covariate$values)))
  }
  loc <- data.frame(time = timepoints)
  colnames(loc) <- covariate$time.coords
  values <- get.smean(covariate, loc=loc, weights=weights)
  return(mean(values))
}
