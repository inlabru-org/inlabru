
sample.lgcp = function(mesh, weights, geometry = "euc", strategy = NULL, R = 6371) {
  
if (is.null(strategy)) {
  if ( geometry == "euc" ) { strategy = "bounding-rectangle" }
  if ( geometry == "geo" ) { strategy = "sliced-spherical"}
}

  
if ( strategy == "bounding-rectangle") {
    
    ###########################################################################
    ## Code for simulating LGCPs on triangulated domains (planar or spherical)
    ## used throughout Simpson et al. (2011)
    ## Author: Daniel Simpson
    ## This code comes with no warranty or guarantee of any kind.
    ###########################################################################
    
    # Construct bounding rectangle
    loc <- mesh$loc
    xmin = min(loc[,1])
    xmax = max(loc[,1])
    ymin = min(loc[,2])
    ymax = max(loc[,2])
    area =(xmax- xmin)*(ymax- ymin)
    
    #Simulate number of points
    lambda_max <- max(weights)
    Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
    
    #Simulate uniform points on the bounding rectangle
    x <- runif(n=Npoints, min=xmin, max=xmax)
    y <- runif(n=Npoints, min=ymin, max=ymax)
    
    points <-cbind(x,y)
    
    #Do some thinning
    A <- inla.mesh.project(mesh,points)$A
    weights =exp( weights-lambda_max)
    pointValues = as.vector(A%*%weights)
    keep = which(runif(Npoints) < pointValues)
    
    return(points[keep,])

} else if ( strategy == "spherical" ) {
  # Simulate number of points
  area =   4*pi*R^2
  lambda_max <- max(weights)
  Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
  if (Npoints > 5000000) { stop(paste0("Too many points!: ",Npoints)) }

  # Choose z uniformly distributed in [-1,1].
  z <- runif(n=Npoints, min=-1, max=1)
  # Choose t uniformly distributed on [0, 2*pi).
  t <- runif(n=Npoints, min=0, max=2*pi)
  r = sqrt(1-z^2)
  x = r * cos(t)
  y = r * sin(t)
  
  loc = euc.to.geo(data.frame(x=x,y=y,z=z), R = 1)
  points = as.matrix(loc[,c("lon","lat")])
  A = inla.mesh.project(mesh,points)$A
  weights = exp(weights - lambda_max)
  pointValues = as.vector(A%*%weights)
  points = points[runif(Npoints) < pointValues,]
  return(points)

} else if ( strategy == "sliced-spherical" ) {
  
  # Simulate number of points
  lon.range = range(mesh$loc[,1])
  lon.rel = abs((lon.range[2]-lon.range[1])/360)
  area = (4*pi*R^2)*lon.rel
  
  lambda_max <- max(weights)
  weights = exp(weights - lambda_max)
  Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
  # cat(Npoints)
  sampled.points = list()
  nblocks = max(1,Npoints/1000000)
  sp = ceiling(seq(1,Npoints+1, length.out = max(2,nblocks)))
  n.sampled = 0
  
  for (k in 1:(length(sp)-1)) {
    n.points = sp[k+1] - sp[k]
    
    # Choose z uniformly distributed in [-1,1].
    z <- runif(n = n.points, min=-1, max=1) # transforms into longitude
    # Choose t uniformly distributed on [0, 2*pi).
    t <- runif(n = n.points, min = 2*pi*(lon.range[1]+180)/360, max = 2*pi*(lon.range[2]+180)/360)
    r = sqrt(1-z^2)
    x = r * cos(t)
    y = r * sin(t)
    
    loc = euc.to.geo(data.frame(x=x,y=y,z=z), R = 1)
    loc[,2] = loc[,2] - 180
    
    points = as.matrix(loc[,c("lon","lat")])
    A <- inla.mesh.project(mesh,points)$A
    pointValues = as.vector(A%*%weights)
    keep = runif(n.points) < pointValues
    sampled.points[[k]] = points[keep,]
    n.sampled = n.sampled + sum(pointValues == 0)
  }
#   cat(paste0("sp:",sp))
#   cat(paste0("n.sampled:",n.sampled))
#   print(paste0("Fraction of points kept: ", n.sampled/Npoints) )
  return(do.call(rbind,sampled.points))
}



  
}
