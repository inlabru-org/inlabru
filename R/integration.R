# Split lines at mesh edges
#
# @aliases split.lines
# @export
# @param mesh An inla.mesh object
# @param sp Start points of lines
# @param ep End points of lines
# @param filter.zero.length Filter out segments with zero length? (Bool)
# @param ... argments to int.quadrature
# @return List of start and end points resulting from splitting the given lines
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
# 
split.lines = function(mesh, sp, ep, filter.zero.length = TRUE) {

  # locations for splitting
  loc = as.matrix(rbind(sp,ep))
  idx = 1:dim(sp)[1]

  # Filter out segments not on the mesh
  t1 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0))]
  loc = as.matrix(rbind(sp,ep))

  # Split them segments into parts
  if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
  np = dim(sp)[1]
  sp.idx = t(rbind(1:np,np+1:np))
  splt = inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
  #plot(data$mesh)
  #points(loc)
  #points(splt$split.loc,col="blue)

  sp = splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]] # Start point of new segments
  ep = splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]] # End points of new segments
  idx = idx[splt$split.idx[,1]]
  origin = splt$split.origin

  # Filter out zero length segments
  if ( filter.zero.length ) {
    sl = apply((ep-sp)^2,MARGIN=1,sum)
    sp = sp[!(sl==0),]
    ep = ep[!(sl==0),]
    origin = origin[!(sl==0)]
    idx = idx[!(sl==0)]
  }

  return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))

}


# Gaussian quadrature and other integration point constructors
# 
# Contruct integration points for each of lines defined by the start and end points provided.
# The following schemes are available: 
# "equidistant" : Equidistant integration points without boundary. All weights are identical and sum uf to the length of a line.
# "gaussian": Points and weight according to the Gaussian quadrature rule. Currently only n=1 and n=2 are supported (Exact integration for linear and quadratic functions).
# "twosided-gaussian": Experimental
#    
# @aliases int.quadrature
# @export
# @param sp Start points of lines
# @param ep End points of lines
# @param scheme Integration scheme (gaussian or equdistant)
# @param n Number of integration points
# @return List with integration poins (ips), weights (w) and weights including line length (wl)
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.quadrature = function(sp=NULL,ep=NULL,scheme="gaussian",n.points=1,geometry="euc",coords = NULL){

  if ( is.null(colnames(sp)) & !is.null(coords)) { sp = data.frame(sp) ; colnames(sp) = coords }
  if ( is.null(colnames(ep)) & !is.null(coords) & !(length(ep)==0) ) { ep = data.frame(ep) ; colnames(ep) = coords  }
  
  if (is.vector(sp)) { 
    n.lines = 1
    len = function(a) {abs(a)}
  }
  else { 
    n.lines = dim(sp)[1] 
    len = function(x) {apply(x^2,MARGIN=1,function(y){return(sqrt(sum(y)))})}
  }
  if (scheme == "gaussian"){
    if (n.points==1){
      # points
      if (geometry == "euc"){
        ips = sp + (ep-sp)/2
        # weights
        w = rep(1,n.lines)
        wl = w*len(ep-sp)
      }  
      else if (geometry == "geo"){
        colnames(sp) = coords
        sp.euc = geo.to.euc(data.frame(sp))
        colnames(ep) = coords
        ep.euc = geo.to.euc(data.frame(ep))
        ips = sp.euc + (ep.euc-sp.euc)/2
        ips = euc.to.geo(data.frame(ips))[,coords,drop=FALSE]
        # weights
        w = rep(1,n.lines)
        wl = w*dist.geo(data.frame(sp),data.frame(ep))
      }

      # PD vector
      # pdv = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
      # Index of line a point comes from
      line.idx = 1:n.lines
      
    } else if (n.points==2){
      if (geometry == "geo") {stop("Geometry geo not supported")}
      # points
      ips1 = sp + (-0.5*sqrt(1/3)+1/2) * (ep-sp)
      ips2 = sp + (0.5*sqrt(1/3)+1/2) * (ep-sp)
      ips = rbind(ips1,ips2)
      
      # weights
      w = rep(1,dim(ips)[1])/2
      wl = w*len(ep-sp)
      
      # PD vector
      #dvec = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
      #dvec = rbind(dvec,dvec)
      
      # Index of line a point comes from
      line.idx = c(1:n.lines,1:n.lines)
      
    } else if (n.points==3){
      if (geometry == "geo") {stop("Geometry geo not supported")}
      # points
      ips1 = sp + (-0.5*sqrt(3/5)+1/2) * (ep-sp)
      ips2 = sp + 0.5*(ep-sp)
      ips3 = sp + (0.5*sqrt(3/5)+1/2) * (ep-sp)
      ips = rbind(ips1,ips2,ips3)
      
      # weights
      w = c(rep(5/9,n.lines),rep(8/9,n.lines),rep(5/9),n.lines)/2
      wl = w*len(ep-sp)
      
      # Index of line a point comes from
      line.idx = c(1:n.lines,1:n.lines)
    }
    
    else { stop("Gaussian quadrature with n>3 not implemented") }
  }
  else if (scheme == "twosided-gaussian"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    ips1 = int.quadrature(sp,sp+0.5*(ep-sp),scheme="gaussian",n.points=n.points, geometry, coords)
    ips2 = int.quadrature(sp+0.5*(ep-sp),ep,scheme="gaussian",n.points=n.points, geometry, coords)
    ips = rbind(ips1$ips,ips2$ips)
    w = 0.5*rbind(ips1$w,ips2$w)
    wl = 0.5*rbind(ips1$wl,ips2$wl)
    line.idx = rbind(ips1$line.idx,ips2$line.idx)
  }
  else if (scheme == "equidistant"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    # points
    mult = seq(0-1/(2*n.points),1+1/(2*n.points),length.out=n.points+2)[2:(n.points+1)]
    nips = list()
    for (k in 1:n.points){
      nips[[k]] = sp + mult[k]*(ep-sp)
    }
    ips = do.call(rbind,nips)
    # weights
    w = rep(1,dim(ips)[1])/n.points
    wl = w*(len(ep-sp)[rep(1:n.lines,n.points)])
    # Index of line a point comes from
    line.idx = rep(1:n.lines,n.points)
  }
  else if (scheme == "trapezoid"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    # points
    mult = seq(0,1,length.out=n.points)
    nips = list()
    for (k in 1:n.points){
      nips[[k]] = sp + mult[k]*(ep-sp)
    }
    ips = do.call(rbind,nips)
    # weights
    w = rep(1,dim(ips)[1])/(n.points-1)
    w[1] = 0.5 * w[1]
    w[length(w)] = 0.5 * w[length(w)] 
    wl = w * (len(ep-sp)[rep(1:n.lines,n.points)])
    # Index of line a point comes from
    line.idx = rep(1:n.lines,n.points)
  }
  else if (scheme == "fixed"){ 
    ips = data.frame(tmp = sp)
    #colnames(ips) = coords
    w = rep(1,length(sp))
    wl = rep(1,length(sp))
    line.idx = rep(NaN,length(sp))
  }
  return(list(ips=ips,w=w,wl=wl,line.idx=line.idx))
}


# Expand integration points over additional dimensions
# 
# See int.quadrature on how to formulate the further arguments.
#
# @aliases int.expand
# @export
# @param ips Integration points to expand
# @param ... lists with arguments for multiple calls of int.quadrature
# @return Integration points
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.expand = function(ips, ...) {
  dots = list(...)
  x = list()
  weight = list()
  for (k in 1: length(dots)){
    if (is.numeric(dots[[k]])) {
      weight[[k]] = rep(1,length(dots[[k]]))
    }
    else {
      quad = do.call(int.quadrature,dots[[k]])
      x[[k]] = quad$ips
      weight[[k]] = quad$wl
      names(x)[[k]] = names(dots)[[k]]
    }
  }
  if (length(dots)>1){
    grd = do.call(expand.grid, x)
    w = apply(do.call(expand.grid, weight), MARGIN = 1, prod)
  }
  else {
    grd = data.frame(x[[1]])
    names(grd) = names(x)[[1]]
    w = weight[[1]]
  }
  rips = ips[as.vector(t(matrix(1:nrow(ips),ncol = nrow(grd),nrow = nrow(ips)))) , ]
  rgrd = grd[rep(seq_len(nrow(grd)), nrow(ips)), , drop = FALSE]
  rips$weight = rips$weight * w[rep(seq_len(nrow(grd)), nrow(ips))]
  return(cbind(rips[,!(colnames(rips)=="weight"),drop=FALSE],rgrd,weight=rips$weight))
}
