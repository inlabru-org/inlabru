#' Plot covariate data set
#'
#' @aliases plot.covariate
#' @export
#' @param covariate A covariate data set
#' @examples \\dontrun{data(sst); plot.covariate(sst)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

plot.covariate = function(covariate,time=NULL,fun=NULL){

  # The grid
  xlim = range(covariate$mesh$loc[,1])
  ylim = range(covariate$mesh$loc[,2])
  x = seq(xlim[1],xlim[2],length.out=500)
  y = seq(ylim[1],ylim[2],length.out=500)
  grid = expand.grid(x=x,y=y)

  # Locations/time
  if (is.null(time)) { time = colnames(covariate$values)[[1]] }
  loc = data.frame(lon = grid$x, lat = grid$y)
  loc[[covariate$time.coords]] = rep(time,nrow(loc))

  # Covariate values
  if (is.null(fun)) { fun = get.value }
  val = do.call(fun,list(covariate=covariate,loc=loc))

  # Plot
  require(lattice)
  levelplot(val ~ grid$x + grid$y,
            col.regions = topo.colors(100),
            panel=function(...){ panel.levelplot(...) },
            xlab = covariate$mesh.coords[[1]],
            ylab = covariate$mesh.coords[[2]])

}


#' Get covariate value at given location/time
#'
#' @aliases get.value
#' @export
#' @param covariate A covariate data set
#' @param loc Locations (spatial/temporal) to obtain the covariate value at
#' @return \code{values} Covariate values
#' @examples \\dontrun{data(sst); get.value(sst,loc=data.frame(lon=-110,lat=0,year=1986))}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.value = function(covariate,loc){
  times = unique(loc[,covariate$time.coords])
  values = numeric(length=nrow(loc))
  for (t in times){
    msk = loc[,covariate$time.coords]==t
    A = inla.spde.make.A(covariate$mesh,loc=as.matrix(loc[msk,covariate$mesh.coords]))
    inside = is.inside(covariate$mesh, loc = loc[msk,covariate$mesh.coords], mesh.coords = covariate$mesh.coords)
    values[msk] = as.vector(A%*%covariate$values[,as.character(t)])
    values[msk & !inside] = NA
  }
  return(values)
}

#' Get maximum value of covariate
#'
#' @aliases get.max
#' @export
#' @param covariate A covariate data set
#' @param loc Locations (spatial/temporal) to obtain the maximal covariate value at. If not provided, return overall maximum.
#' @return \code{values} Maximal covariate values
#' @examples \\dontrun{data(sst); get.max(sst)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.max = function(covariate,loc=NULL){
  if (is.null(loc)){ return(max(as.vector(covariate$values))) }
  else {
    vals = get.value(covariate,loc)
    return(max(vals))
  }
}

#' Get minimum value of covariate
#'
#' @aliases get.min
#' @export
#' @param covariate A covariate data set
#' @param loc Locations (spatial/temporal) to obtain the maximal covariate value at. If not provided, return overall maximum.
#' @return \code{values} Maximal covariate values
#' @examples \\dontrun{data(sst); get.min(sst)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.min = function(covariate,loc=NULL){
  if (is.null(loc)){ return(max(as.vector(covariate$values))) }
  else {
    vals = get.value(covariate,loc)
    return(min(vals))
  }
}


#' Get temporal mean
#'
#' @aliases get.tmean
#' @export
#' @param covariate A covariate data set
#' @param loc Locations (spatial/temporal) to obtain temporal mean covariate value at. If not provided, return overall maximum.
#' @return \code{values} Temporal mean
#' @examples \\dontrun{data(sst); get.min(sst)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.tmean = function(covariate,loc=NULL){
  tm = apply(covariate$values,MARGIN=1,mean)
  loc2 = loc[,covariate$mesh.coords,drop=FALSE]
  A = inla.spde.make.A(covariate$mesh,loc=as.matrix(loc2))
  values = as.vector(A%*%tm)
  return(values)
}


#' Get spatial mean of covariate for given time
#'
#' Note: If this function seems to run forever use
#' refine = list(max.edge=K) with some integer K>1
#'
#' WARNING: This is an approximation!
#'
#'
#' @aliases get.smean
#' @export
#' @param covariate A covariate data set
#' @param loc Locations (spatial/temporal) to obtain the spatial mean covariate value at.
#' @param weights Averaging weights, default all ones
#' @return \code{values} Spatial mean or spatial weighted mean
#' @examples \\dontrun{data(sst); get.smean(sst)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.smean = function(covariate,loc=NULL,weights=NULL){
  sm <- numeric(length=ifelse(is.null(loc), 0, dim(loc)[1]))

  if (is.null(weights)) {
    weights <- rep(1, covariate$mesh$n)
  }
  weights <- weights * Matrix::diag(inla.mesh.fem(covariate$mesh, order=1)$c0)
  weights <- weights / sum(weights)

  loc2 <- data.frame(covariate$mesh$loc[, 1:length(covariate$mesh.coords),
                                        drop=FALSE])
  colnames(loc2) <- covariate$mesh.coords
  if (!(covariate$time.coords %in% colnames(loc))) {
    stop("Time coordinates missing")
  }
  for (t in unique(loc[,covariate$time.coords])){
    loc2[[covariate$time.coords]] <- t
    values <- get.value(covariate, loc2)
    sm[loc[covariate$time.coords]==t] <- sum(values * weights)
  }
  return(sm)
}



#' Get mean
#'
#' @aliases get.mean
#' @export
#' @param covariate A covariate data set
#' @param loc Locations (spatial/temporal) to obtain mean covariate value at. If not provided, return overall maximum.
#' @param timepoints ????
#' @param weights ????
#' @return \code{values} Temporal mean
#' @examples \\dontrun{data(sst); get.min(sst)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

get.mean = function(covariate,loc=NULL, timepoints=NULL, weights=NULL){
  if (is.null(timepoints)) {
    timepoints <- as.numeric(unique(colnames(covariate$values)))
  }
  loc <- data.frame(time = timepoints)
  colnames(loc) <- covariate$time.coords
  values <- get.smean(covariate, loc=loc, weights=weights)
  return(mean(values))
}
