#' Automatic integration scheme selection
#'
#'
#' @aliases select.integration
#' @export
#' @param data A distance sampling data set
#' @return scheme An integration scheme 
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

select.integration = function(data){
  if (data$geometry=="euc"){ scheme = grid.integration }
  else if (data$geometry=="geo"){ scheme = cylindric.integration }
  else {stop("Unable to select integration scheme")}
  return(scheme)
}

#' INLA inference in log Gaussian Cox processes for distance sampling data.
#'
#'
#' @aliases pproc_int.scheme.default
#' @export
#' @import pracma
#' @param data A distance sampling data set, e.g. \link{whales}.
#' @return \code{int.points}, Integration points with respective weights.
#' @author Yuan (Joyce) Yuan <\email{yy84@@st-andrews.ac.uk}>
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

pproc_int.scheme.default = function(data){
  
  eff = data$effort
  
  #
  # Joyce's code
  #
  
  
  ####################
  ## NAs separate transect segments
  allID.NAeff = which(is.na(eff$lat), arr.ind=T)
  allID.NAeff =allID.NAeff [-1]
  #head(allID.NAeff)
  nLTsegment = length(allID.NAeff) -1
  
  w=6
  dLT=12
  dPD=3
  
  
  ### for the 1st LT ######
  p0lonlat = c(eff[1, c("lon", "lat")])
  p1lonlat = c(eff[2, c("lon", "lat")])
  
  IntMatrixLonLatList.1stLT=designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)$intgridlonlat
  
  ### the following lapply() starts from the 2nd transect line ###
  ####### it is much faster to use lapply(); tried loop, cannot deal with
  IntMatrixLonLatList.exclude1stOne = lapply(matrix(data=c(2:nLTsegment), ncol=1),
                                             function(idx.LTsegment){
                                               p0lonlat = c(eff[(allID.NAeff[idx.LTsegment-1]+1), c("lon", "lat")])
                                               p1lonlat = c(eff[(allID.NAeff[idx.LTsegment]-1), c("lon", "lat")])
                                               designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)$intgridlonlat})
  ## append the 1st LT sement to the rest obtained from lapply() above
  
  IntMatrixLonLatList = append(list(IntMatrixLonLatList.1stLT),IntMatrixLonLatList.exclude1stOne)
  #save(IntMatrixLonLatList,  file = paste(resultdir, "IntMatrixLonLatList", version, ".Rdata",sep=""))
  ##############################
  # debug
  # p0lonlat = c(eff[2, c("lon", "lat")])
  # p1lonlat = c(eff[3, c("lon", "lat")])
  # IntMatrixLonLatList.1stLT=designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)$intgridlonlat
  #
  # for(idx.LTsegment in c(2:nLTsegment))
  # {
  #   p0lonlat = c(eff[(allID.NAeff[idx.LTsegment-1]+1), c("lon", "lat")])
  #   p1lonlat = c(eff[(allID.NAeff[idx.LTsegment]-1), c("lon", "lat")])
  #   designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)$intgridlonlat
  # }
  
  #########################################################################################
  ## unlist the above list.result using do.call("rbind", listobject)
  IntPointMatrixLonLat = do.call("rbind", IntMatrixLonLatList )
  #save(IntPointMatrixLonLat, file=paste(resultdir, "IntPointMatrixLonLat", version, ".Rdata",sep=""))
  #nrow(IntPointMatrixLonLat)
  ## check NAs in the output
  idNA = which(is.na(IntPointMatrixLonLat[,"lon"]), arr.ind=T)
  #length(idNA)
  
  
  # #####################################################
  # ## add NA to separate between LT segments
  # load(paste(resultdir, "IntMatrixLonLatList", version, ".Rdata",sep=""))
  # length(IntMatrixLonLatList)
  # head(IntMatrixLonLatList[[1]][[2]])
  #
  # IntMatrixLonLatListAddNA = NULL
  # IntMatrixLonLatListAddNA = lapply(c(1:length(IntMatrixLonLatList)),
  #                                   function(X){
  #                                     IntMatrixLonLatListAddNA[[X]] = rbind(IntMatrixLonLatList[[X]], NA)})
  # IntPointMatrixLonLatAddNA = do.call("rbind",IntMatrixLonLatListAddNA)
  # save(IntPointMatrixLonLatAddNA, file=paste(resultdir, "IntPointMatrixLonLatAddNA", version, ".Rdata",sep=""))
  ###################################################
  
  #####################################################
  #######  for the design matrix on tangent plane
  ### for the 1st LT ######
  p0lonlat = c(eff[1, c("lon", "lat")])
  p1lonlat = c(eff[2, c("lon", "lat")])
  IntMatrixv1v2v3List.1stLT=designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)$intgridv1v2v3
  
  IntMatrixv1v2v3List.exclude1stOne = lapply(matrix(data=c(2:nLTsegment), ncol=1),
                                             function(idx.LTsegment){
                                               p0lonlat = c(eff[(allID.NAeff[idx.LTsegment-1]+1), c("lon", "lat")])
                                               p1lonlat = c(eff[(allID.NAeff[idx.LTsegment]-1), c("lon", "lat")])
                                               designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)$intgridv1v2v3})
  ## append the 1st LT sement to the rest obtained from lapply() above
  IntMatrixv1v2v3List = append(list(IntMatrixv1v2v3List.1stLT),IntMatrixv1v2v3List.exclude1stOne)
  #save(IntMatrixv1v2v3List,  file = paste(resultdir, "IntMatrixv1v2v3List", version, ".Rdata",sep=""))
  #length(IntMatrixv1v2v3List)
  ############################
  
  ipd = rep(unique(IntMatrixv1v2v3List[[1]]$dPD), nrow(IntPointMatrixLonLat))
  int.points = data.frame(IntPointMatrixLonLat,PD=ipd)
  
  #
  # Return
  #
  
  return(int.points)
  
}

pproc_int.scheme.default.tidy = function(data,w=6,dLT=12*8,dPD=3){
  
  tr = data$transect
  IntMatrixLonLatList = lapply(matrix(data=c(1:dim(tr)[1]), ncol=1),
                               function(i){
                                 p0lonlat = list(lon=tr[i,"start.lon"], lat=tr[i,"start.lat"])
                                 p1lonlat = list(lon=tr[i,"end.lon"], lat=tr[i,"end.lat"])
                                 designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)})
  
  LL = do.call("rbind",lapply(IntMatrixLonLatList,function(x){x$intgridlonlat}))
  VV = do.call("rbind",lapply(IntMatrixLonLatList,function(x){x$intgridv1v2v3}))
  
  int.points = cbind(LL,VV)
  colnames(int.points)[8] = "PD"
  #
  # Return
  #
  
  return(int.points)
  
}


pproc_int.scheme.reduce.daily= function(data,w=6,dLT=12*8,dPD=3){
  
  tr = transect.reduce.by.day(data$transect)
  
  IntMatrixLonLatList = lapply(matrix(data=c(1:dim(tr)[1]), ncol=1),
                               function(i){
                                 p0lonlat = list(lon=tr[i,"start.lon"], lat=tr[i,"start.lat"])
                                 p1lonlat = list(lon=tr[i,"end.lon"], lat=tr[i,"end.lat"])
                                 designintegration.oneLTsegement(p0lonlat, p1lonlat, w, dLT, dPD)})
  
  LL = do.call("rbind",lapply(IntMatrixLonLatList,function(x){x$intgridlonlat}))
  VV = do.call("rbind",lapply(IntMatrixLonLatList,function(x){x$intgridv1v2v3}))
  
  int.points = cbind(LL,VV)
  colnames(int.points)[8] = "PD"
  
  return(int.points)
  
}

#' Midpoint integration
#'
#'
#' @aliases cylindric.integration
#' @export
#' @return int.points integration points
#'   $lat : Latitude
#'   $lon : Longitude
#'   $PD : distance from observer
#'   $weight: integration weight
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
midpoint.integration = function(data,group,distance){
  tr = data$transect
  spoint = startpoint(tr,data)
  epoint = endpoint(tr,data)
  weight = 2*dist(epoint,spoint)/length(distance)
  ipoint = ssum(spoint,epoint)/2
  ips = list()
  for (k in 1:length(distance)) {
    ips[[k]] = data.frame(ipoint,weight=weight,distance=distance[k])
  }
  return(do.call(rbind,ips))
}

#' Construct integration points in cylindric space
#'
#'
#' @aliases cylindric.integration
#' @export
#' @return int.points integration points
#'   $lat : Latitude
#'   $lon : Longitude
#'   $PD : distance from observer
#'   $weight: integration weight
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
cylindric.integration = function(data,group=NULL,...){
  tr = as.transect(data)
  if (is.null(group)){group = group.by.length(tr,data,max.ratio=500,max.arclen=5000)}
  else if (is.function(group)) { group = group(tr,data) }

  # transform to cylindric space
  cyl.list = lapply(unique(group),
                    function(g) {
                      a = to.cyl.transect(tr[g==group,],data)
                      b = linedata.etptransect(tr[g==group,],data,"year")[1]
                      return(c(a,list(year=b)))
                    }
  )
  
  # construct integration points
  geo.ips = lapply(cyl.list,function(x) {
    ips = ischeme(x$startpoints,x$endpoints,...)
    geo = cyl.to.geo(ips,B=x$B)
    geo = cbind(geo,ips,year=x$year)
    return(geo)})
  
  geo.ips = do.call(rbind, geo.ips)
  return(geo.ips)
}


#' Create integration points from cylindric coordinates of a group of transect lines.
#'
#' @aliases ischeme
#' @export
#' @param spoints Start points of transects
#' @param epoints End points of transects
#' @return Integration points ($lat,$lon) with integration weights ($weight) and perpendicular distance ($PD)
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

ischeme = function(spoints,epoints,n.dist=5,truncation=6,R=6371,project=FALSE,...){
  if (project){
    spoints$z = 0
    epoints$z = 0
  }
  # Finns integration scheme
  len = dist(spoints,epoints)
  tim = tdist(spoints,epoints)
  w.hat.s = sum(len)
  w.hat.t = sum(tim)
  s.hat = (1/w.hat.s) * apply(cbind(len,len) * ssum(spoints,epoints)/2,MARGIN=2,sum)
  t.hat = (1/w.hat.t) * apply(cbind(tim,tim) * tsum(spoints,epoints)/2,MARGIN=2,sum)

  # distances
  PD = seq(0,truncation,length.out=n.dist+1)+(truncation/(2*n.dist))
  PD = PD[1:n.dist]
  
  # weights
  w =  R*w.hat.s * truncation/n.dist
  
  ret = data.frame(data.frame(t(s.hat)),
                   distance = PD,
                   weight = w)
  class(ret) = class(spoints)
  return(ret)
}

#' Temporary implementation of integration point construction
#'
#' @aliases temporary.integration.scheme
#' @export
#' @param effort ETP effort data with NA-lines removed
#' @return Integration points ($lat,$lon) with integration weights ($weight), perpendicular distance ($PD) and year ($year)
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

temporary.integration.scheme = function(effort,max.ratio=300,...){
  transect = as.transect.etpeffort(effort)
  tmpdata = list(effort=effort,transect=transect)
  class(tmpdata) = c("etpdata","dsdata","data.frame")
  gfun = function(x,y) {return(group.by.length.transect(x,y,max.ratio=max.ratio,max.arclen=5000))}
  group = nested.group(transect,tmpdata,"year",gfun)
  ips = cylindric.integration(tmpdata,group,...)
  return(ips)
}


#' Integration on a grid (for Euclidean data)
#'
#' @aliases grid.integration
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>


grid.integration = function(data,group,n.distance=3,n.length=3,truncation=3){
  distance = seq(-truncation,truncation,length.out=(n.distance+2))
  distance = distance[2:(length(distance)-1)]
  lens = seq(0,1,length.out=(n.length+2))
  lens = lens[2:(length(lens)-1)]
  tr = as.transect(data)
  spoints = startpoint(data,tr)
  epoints = endpoint(data,tr)
  ips = list()
  i = 1
  for (k in 1:nrow(spoints)){
    for (d in distance){
      for (l in lens) {
        spoint = spoints[k,] ; class(spoint) = class(spoints)
        epoint = epoints[k,] ; class(epoint) = class(epoints)
        locvec = ssum(epoint,-spoint)
        orthvec = normalize.euc(data.frame(x=locvec$y,y=-locvec$x));
        pt = spoint + orthvec*d + l*locvec
        pt$weight = 2*dist(epoint,spoint)*truncation/(length(distance)*length(lens))
        pt$distance = abs(d)
        ips[[i]] = pt
        i=i+1
      }
    }
  }
  ips = do.call(rbind,ips)
  return(ips)
}


#' Triangle indices of points given a mesh
#'
#' @aliases triangle
#' @export
#' @param mesh A inla.mesh
#' @param loc Locations using the coordinate system of the mesh 
#' @return tri Triangle indices
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

triangle = function(mesh,loc){
  mcross = function(a,b) { return( (a[,1])*(b[2]) - (a[,2])*(b[1]) )}
  tri = numeric(length=dim(loc)[1])
  tv = mesh$graph$tv  
  for (j in 1:nrow(tv)) {
    v = mesh$loc[tv[j,],c(1,2)]
    
    a = v[1,]
    b = v[2,]
    c = v[3,]
    
    ap =  loc - cbind(rep(a[1],nrow(loc)),rep(a[2],nrow(loc)))
    bp =  loc - cbind(rep(b[1],nrow(loc)),rep(b[2],nrow(loc)))
    cp =  loc - cbind(rep(c[1],nrow(loc)),rep(c[2],nrow(loc)))
    
    ab = a-b
    bc = b-c
    ca = c-a
    
    c1 = mcross(ap,ab)
    c2 = mcross(bp,bc)
    c3 = mcross(cp,ca)
    
    # AP × AB, BP × BC, and CP × CA must have same sign
    inside = which( ((sign(c1) == -1) & (sign(c2)==-1) & (sign(c3)==-1)) | ((sign(c1) == 1) & (sign(c2)==1) & (sign(c3)==1)))
    
    tri[inside] = j
  }
  return(tri)
}


ip.project = function(projection, mesh, ips, w, mesh.coords = c("x","y")){

  #
  # Project integration points
  #
  if (projection == "single") {
    # 1) By weighted average
    tips = data.frame(ips)
    colnames(tips) = data$mesh.coords
    tips = data.frame(tips,tri=tri,weight=w)
    fun = function(x) { merge.points(x,data$mesh.coords) }
    ips = recurse.rbind(fun,tips,c("tri"))
    
  } else if ( projection == "bezier" ) {
    
    # Obtain baryocentric coordinates of integration points:
    res = inla.fmesher.smorg(mesh$loc, mesh$graph$tv,points2mesh=ips)
    tri = res$p2m.t # re-calculate triangles
    ips.bc = res$p2m.b # barycentric coordinates
    ips.B = bezier.basis(ips.bc)
    ctl.points = diag(3) # for now, the control points are the triangle vertices
    ctl.B = bezier.basis(ctl.points)
    
    ips.wB = ips.B * matrix(rep(w,size(ips.B)[2]),ncol=size(ips.B)[2]) # point-wise product of ips.B and w
    colsum = function(x) { return( apply(x,MARGIN=2,sum) ) }
    rhs = simplify2array(by(ips.wB,tri,colsum,simplify=TRUE)) # for each triangle, sum up wB
    
    bmp = barycentric.midpoints(mesh,sort(unique(tri))) # edge midpoint barycentric coordinates of triangles
    bmp.bez = lapply(bmp,bezier.basis) # basis evaluated at those midpoints
    stri = sort(unique(tri))
    
    wc = matrix(NA,nrow=length(unique(tri)),ncol=6)
    neg = vector(length=length(unique(tri)))
    for (k in 1:length(stri)){
      t = stri[k]
      y = rhs[,k]
      b1 = bmp.bez$b1[k,]
      b2 = bmp.bez$b2[k,]
      b3 = bmp.bez$b3[k,]
      
      cns.euc = mesh$loc[mesh$graph$tv[t,],]
      cns.bez = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv[t,,drop=FALSE],points2mesh=cns.euc)$p2m.b
      cns.basis = bezier.basis(cns.bez)
      
      
      if (sum(b1)==0 | sum(b2)==0 | sum(b3)==0){
        warning("BUG!!!!!!!!!")
      } else { B = t(rbind(cns.basis,rbind(b1,b2,b3))) }
      
      wc[k,] = solve(B,y)
      
    }
    
    
    ips = rbind(bmp$p1,bmp$p2,bmp$p3,bmp$mp1,bmp$mp2,bmp$mp3)[,1:2]
    colnames(ips) = data$mesh.coords
    weight = as.vector(wc)
    ips = data.frame(ips,weight)
    
    # Accumulate
    wsum = function(x) { w = sum(x$weight) ; x=x[1,]; x$weight = w; return( x ) }
    tips = by(ips,round(ips[1:2],digits=2),wsum)
    #tips = tips[!simplify2array(lapply(tips,is.null))]
    ips = do.call(rbind,tips)
    
  }  else if ( projection == "linear" ) {
    
    # Obtain baryocentric coordinates of integration points:
    res = inla.fmesher.smorg(mesh$loc, mesh$graph$tv,points2mesh=ips[,mesh.coords])
    tri = res$p2m.t # re-calculate triangles
    ips.bc = res$p2m.b # barycentric coordinates
    ips.B = ips.bc
    # ips.B = bezier.basis(ips.bc)
    #ctl.points = diag(3) # for now, the control points are the triangle vertices
    #ctl.B = bezier.basis(ctl.points)
    
    ips.wB = 0.5 * ips.B * matrix(rep(w,size(ips.B)[2]),ncol=size(ips.B)[2]) # point-wise product of ips.B and w
    colsum = function(x) { return( apply(x,MARGIN=2,sum) ) }
    rhs = simplify2array(by(ips.wB,tri,colsum,simplify=TRUE)) # for each triangle, sum up wB
    
    bmp = barycentric.vertices(mesh,sort(unique(tri))) # edge midpoint barycentric coordinates of triangles
    bmp.bez = bmp
    #bmp.bez = lapply(bmp,bezier.basis) # basis evaluated at those midpoints
    stri = sort(unique(tri))
    
    wc = matrix(NA,nrow=length(unique(tri)),ncol=6)
    neg = vector(length=length(unique(tri)))
    for (k in 1:length(stri)){
      t = stri[k]
      y = rhs[,k]
      b1 = bmp.bez$b1[k,]
      b2 = bmp.bez$b2[k,]
      b3 = bmp.bez$b3[k,]
      
      #cns.euc = mesh$loc[mesh$graph$tv[t,],]
      #cns.bez = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv[t,,drop=FALSE],points2mesh=cns.euc)$p2m.b
      #cns.basis = bezier.basis(cns.bez)
      
      
      if (sum(b1)==0 | sum(b2)==0 | sum(b3)==0){
        warning("BUG!!!!!!!!!")
      } else { B = t(rbind(rbind(b1,b2,b3))) }
      
      wc[k,] = solve(B,y)
      
    }
    
    
    ips = rbind(bmp$p1,bmp$p2,bmp$p3,bmp$mp1,bmp$mp2,bmp$mp3)[,1:2]
    colnames(ips) = data$mesh.coords
    weight = as.vector(wc)
    ips = data.frame(ips,weight)
    
    # Accumulate
    wsum = function(x) { w = sum(x$weight) ; x=x[1,]; x$weight = w; return( x ) }
    tips = by(ips,round(ips[1:2],digits=2),wsum)
    #tips = tips[!simplify2array(lapply(tips,is.null))]
    ips = do.call(rbind,tips)
    
  } else { # 3) project to midpoints of triangle  
  }

ret = ips
return(ret)
  
}

bezier.basis = function(b){
  return( cbind(b[,1]^2 , b[,2]^2, b[,3]^2 , b[,1]*b[,2], b[,1]*b[,3], b[,2]*b[,3]  ))
}

barycentric.midpoints = function(mesh,tri){
  mp1 = mesh$loc[mesh$graph$tv[tri,1],] + 0.5*(mesh$loc[mesh$graph$tv[tri,2],]-mesh$loc[mesh$graph$tv[tri,1],])
  mp2 = mesh$loc[mesh$graph$tv[tri,1],] + 0.5*(mesh$loc[mesh$graph$tv[tri,3],]-mesh$loc[mesh$graph$tv[tri,1],])
  mp3 = mesh$loc[mesh$graph$tv[tri,2],] + 0.5*(mesh$loc[mesh$graph$tv[tri,3],]-mesh$loc[mesh$graph$tv[tri,2],])
  
  orth = mesh$loc[mesh$graph$tv[tri,2],]-mesh$loc[mesh$graph$tv[tri,1],]
  mp1 = mp1 + 50*eps(1)*cbind(-orth[,2],orth[,1],orth[,3])
  
  orth = mesh$loc[mesh$graph$tv[tri,3],]-mesh$loc[mesh$graph$tv[tri,1],]
  mp2 = mp2 + 50*eps(1)*cbind(-orth[,2],orth[,1],orth[,3])
  
  orth = mesh$loc[mesh$graph$tv[tri,3],]-mesh$loc[mesh$graph$tv[tri,2],]
  mp3 = mp3 + 50*eps(1)*cbind(-orth[,2],orth[,1],orth[,3])
  
  b1 = matrix(NA,nrow=size(mp1)[1],ncol=3)
  b2 = matrix(NA,nrow=size(mp2)[1],ncol=3)
  b3 = matrix(NA,nrow=size(mp3)[1],ncol=3)
  for (k in 1:length(tri)){
      tr = tri[k]
      tv = mesh$graph$tv[c(tr,tr),]
      b1[k,] = c(0.5,0.5,0)# inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp1[k,]))$p2m.b
      b2[k,] = c(0.5,0,0.5)# inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp2[k,]))$p2m.b
      b3[k,] = c(0,0.5,0.5)# inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp3[k,]))$p2m.b
    }
  #orth = mesh$loc[mesh$graph$tv[tri,3],]-mesh$loc[mesh$graph$tv[tri,2],]
  #mp3 = mp3 + 10*eps(1)*cbind(-orth[,2],orth[,1],orth[,3])
  #inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp3[k,]))$p2m.b
  #plot(data$mesh)
  #points(mp1); points(mp2); points(mp3)
  b = list(b1=b1,b2=b2,b3=b3,mp1=mp1,mp2=mp2,mp3=mp3,
           p1=mesh$loc[mesh$graph$tv[tri,1],],
           p2=mesh$loc[mesh$graph$tv[tri,2],],
           p3=mesh$loc[mesh$graph$tv[tri,3],])

  return(b)
}




barycentric.vertices = function(mesh,tri){
  mp1 = mesh$loc[mesh$graph$tv[tri,1],] 
  mp2 = mesh$loc[mesh$graph$tv[tri,2],] 
  mp3 = mesh$loc[mesh$graph$tv[tri,3],] 
  
  b1 = matrix(NA,nrow=size(mp1)[1],ncol=3)
  b2 = matrix(NA,nrow=size(mp2)[1],ncol=3)
  b3 = matrix(NA,nrow=size(mp3)[1],ncol=3)
  for (k in 1:length(tri)){
    tr = tri[k]
    tv = mesh$graph$tv[c(tr,tr),]
    b1[k,] =  inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp1[k,]))$p2m.b
    b2[k,] =  inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp2[k,]))$p2m.b
    b3[k,] =  inla.fmesher.smorg(loc=mesh$loc,tv=tv,points2mesh=t(mp3[k,]))$p2m.b
  }

  b = list(b1=b1,b2=b2,b3=b3,mp1=mp1,mp2=mp2,mp3=mp3,
           p1=mesh$loc[mesh$graph$tv[tri,1],],
           p2=mesh$loc[mesh$graph$tv[tri,2],],
           p3=mesh$loc[mesh$graph$tv[tri,3],])
  
  return(b)
}



#' Split lines at mesh edges
#'
#' @aliases split.lines
#' @export
#' @param mesh An inla.mesh object
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param filter.zero.length Filter out segments with zero length? (Bool)
#' @param ... argments to int.quadrature
#' @return List of start and end points resulting from splitting the given lines
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

split.lines = function(mesh, sp, ep, filter.zero.length = TRUE) {
  
  # locations for splitting
  loc = as.matrix(rbind(sp,ep))
  idx = 1:dim(sp)[1]
  
  # Filter out segments not on the mesh
  t1 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
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
  
  return(list(sp=sp,ep=ep,split.origin=origin,idx=idx))
  
}

#' Replicate lines with different distances to base line
#'
#' @aliases replicate.lines
#' @export
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param truncation distance at which we truncate sightings
#' @param ... argments to int.quadrature
#' @return List of 1) Start and end points of replicated lines 2) their distances to the original line
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

replicate.lines = function(sp,ep,truncation,scheme="equidistant",n=3,fake.distance=FALSE,geometry="euc") {
  
  if (fake.distance) {
    quad = int.quadrature(0,truncation,scheme,n)
  } else {
    quad = int.quadrature(-truncation,truncation,scheme,n)
  }
  
  dst = quad$ips
  idx = 1:dim(sp)[1]
  
  if (geometry == "euc"){
    # Normal vectors
    dvec = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
    nrm = apply(dvec^2,MARGIN=1,function(x) {return(sqrt(sum(x)))})
    dvec = dvec/rep(nrm,dim(dvec)[2])
    
    tmp.sp = list()
    tmp.ep = list()
    distance = list()
    for (k in 1:length(dst)){
      if (fake.distance){
        tmp.sp[[k]] = sp
        tmp.ep[[k]] = ep
      }
      else {
        tmp.sp[[k]] = sp + dvec*dst[k]
        tmp.ep[[k]] = ep + dvec*dst[k]
      }
      distance[[k]] = rep(dst[k],dim(sp)[1])
    }
    sp  = do.call(rbind,tmp.sp)
    ep  = do.call(rbind,tmp.ep)
    distance = abs(do.call(c,distance))
    idx = rep(idx,length(dst))
    return(list(sp=sp,ep=ep,distance=distance,idx=idx))
  
  } else if (geometry == "geo"){
    stop("not supported: geometry geo")
  }
}

#' Gaussian quadrature and other integration point constructors
#' 
#' Contruct integration points for each of lines defined by the start and end points provided.
#' The following schemes are available: 
#' "equidistant" : Equidistant integration points without boundary. All weights are identical and sum uf to the length of a line.
#' "gaussian": Points and weight according to the Gaussian quadrature rule. Currently only n=1 and n=2 are supported (Exact integration for linear and quadratic functions).
#' "twosided-gaussian": Experimental
#'    
#' @aliases int.quadrature
#' @export
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param scheme Integration scheme (gaussian or equdistant)
#' @param n Number of integration points
#' @return List with integration poins (ips), weights (w) and weights including line length (wl)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.quadrature = function(sp,ep,scheme="gaussian",n.points=1,geometry="euc",coords=c("x","y")){
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
    
    else { stop("Gaussian quadrature with n>2 not implemented") }
  }
  else if (scheme == "twosided-gaussian"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    ips1 = int.quadrature(sp,sp+0.5*(ep-sp),scheme="gaussian",n.points=n.points)
    ips2 = int.quadrature(sp+0.5*(ep-sp),ep,scheme="gaussian",n.points=n.points)
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
  return(list(ips=ips,w=w,wl=wl,line.idx=line.idx))
}


#' Wrapper for \link{int.quadrature}
#' 
#' Returns a data frame instead of a list. See documentation of \link{int.quadrature}.
#' 
#' @aliases int.1d
#' @param ... see 
#' @export \link{int.quadrature}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.1d = function(...){
  quad = int.quadrature(...)
  ips = data.frame(distance = quad$ips, weight = quad$wl)
  return(ips)
}

#' Line transect integration
#' 
#' Creates a set of integration points with weights from a \link{dsdata} structure. 
#' The integration points can be based on the transect lines or the line segments of the survey (parameter "on").
#' By default, the integration points are arranged on a grid and their number in direction of the transect and perpendicular
#' to the transect are given by the parameters n.length and n.distance. The boundary of the integration in perpendicular
#' direction is given by distance.truncation. Since transect lines can be long compared to the mesh that is used to model
#' the data it can be useful to set the parameter line.split to TRUE. This means that transect lines / segments are split
#' into parts at the edges of the mesh before the integration points are constructed (for each of these parts). The
#' mesh that is used for this procedure is (by default) the mesh of the data set (data$mesh). However, it can by replaced
#' by a mesh provided as a parameter. Additionaly, the given mesh can automatically be refined (mesh.refine) if a denser set 
#' of integration points is required. The latter also holds for integration points constructed using the option 
#' "projection = TRUE". Hereby, after their contruction, integration points are projected onto points at the mesh vertices.
#'
#' @aliases int.points
#' @export
#' @param data Either a dsdata/etpdata data set (e.g. whales) or a data.frame describing effort data
#' @param on Either "transect" or "segment". This determines on which of these the integration is based on. Alternatively a two column index matrix, first column: column index of transect start points in effort data, second column: column index of transect end points in effort data
#' @param line.split TRUE or FALSE, determines if lines that cross mesh edged should be splitted
#' @param mesh Mesh used to construct the integration points. By default the mesh of the given data set.
#' @param mesh.split Split mesh triangles into four sub-triangles for refined integration.
#' @param mesh.coords Character description of the mesh coordinates, e.g. c("lon","lat")
#' @param geometry Either "geo" (geographic)  or "euc" (Euclidean)
#' @param length.scheme Integration scheme along the line (transect/segment)
#' @param n.length Number of integration points along the line
#' @param distance.scheme Integration scheme along perpendicular distance
#' @param n.distance Number of integration points perpendicular distance
#' @param distance.truncation Truncation for perpendicular distance (i.e. integration limit)
#' @param fake.distance Wether or not integration points stay on the transect line and distances are faked.
#' @param projection Type of projection. Currently only "linear" works.
#' @param group.by Create independent integration points for sub-groups of the data, e.g. group.by=c("year","month")
#' @return List with integration poins (ips), weights (w) and weights including line length (wl)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>


int.points = function(data, 
                  on = "transect",
                  line.split = FALSE,
                  mesh = data$mesh, 
                  mesh.split = FALSE,
                  mesh.coords = data$mesh.coords,
                  geometry = data$geometry, 
                  length.scheme = "gaussian", 
                  n.length = 1, 
                  distance.scheme = "equidistant", 
                  n.distance = 5, 
                  distance.truncation = 1, 
                  fake.distance = FALSE,
                  projection = NULL,
                  group.by = NULL,
                  filter.zero.length = TRUE
){
  
  if (mesh.split) { mesh = mesh.split(mesh) }
  
  if (is.data.frame(data)){
    covariates = data
    # Segment/transect start and end points 
    sp = strip.coords(data[, paste0("start.",mesh.coords)])
    ep = strip.coords(data[, paste0("end.",mesh.coords)])
    idx = data.frame(1:nrow(sp), 1:nrow(sp))
    
  } else {
    
    # Get segments/transect indices
    if (on == "transect" ){ idx = as.transect(data) }
    else if (on == "segment" ) { idx = as.segment(data) }
    else ( idx = on )
    
    # Covariates
    covariates = data$effort
    # Segment/transect start and end points 
    sp = startpoint(idx,data)[,mesh.coords]
    ep = endpoint(idx,data)[,mesh.coords]
  }
  
  # Filter out bogus stuff
  loc = as.matrix(rbind(sp,ep))
  t1 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  if (any(t1==0) | any(t2==0)) { warning("filtering...")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0)),]
  
  
  # Replicate segments for different distances (out:sp,ep,distance)
  line.rep = replicate.lines(sp,ep,distance.truncation,scheme=distance.scheme,n=n.distance,fake.distance=fake.distance)
  sp = line.rep$sp
  ep = line.rep$ep
  distance = line.rep$distance
  #split.origin=splt$split.origin,idx=idx
  idx = idx[line.rep$idx,]
  
  # Split lines
  if ( line.split ) {
    line.spl = split.lines(mesh, sp, ep, filter.zero.length)
    sp = line.spl$sp
    ep = line.spl$ep
    distance = distance[line.spl$split.origin]
    idx = idx[line.spl$split.origin,]
  }
  
  # Determine integration points along lines
  quad = int.quadrature(sp, ep, scheme = length.scheme, n.points = n.length, geometry=geometry,coords=mesh.coords)
  ips = quad$ips 
  w = 2*distance.truncation*quad$wl/n.distance
  distance = distance[quad$line.idx]
  idx = idx[quad$line.idx,]
  
  # Wrap everything up and perform projection according to distance and given group.by argument
  ips = data.frame(ips)
  colnames(ips) = mesh.coords
  ips= cbind(ips, distance = distance, weight = w)
  
  #
  # OLD SCHEME
  #
  
  #   if (!is.null(group.by)) { ips = cbind(ips,covariates[idx[,1],group.by,drop=FALSE])}
  #   
  #   projfun = function(dframe) { data.frame(ip.project(projection = projection, mesh, as.matrix(dframe[,mesh.coords]), dframe$weight,mesh.coords=mesh.coords),
  #                                           dframe[1,c("distance",group.by),drop=FALSE]) }
  #   ips = recurse.rbind(projfun,ips,c("distance",group.by))
  if ( !is.null(projection)){
    fu = function(x) data.frame(cbind(project.weights(x, mesh = data$mesh, mesh.coords = data$mesh.coords),distance = x$distance[1]))
    ips = recurse.rbind(fu, ips, c(list("distance"),group.by))
  }
  if (length(idx) == dim(ips)[1]) {ips = data.frame(ips, idx = idx[,1])} # only return index if it makes sense
  return(ips)
}

#' Prjection integration done FAST. Not implemented completely, has to deal with distances.

project.weights = function(ips, mesh, mesh.coords = NULL){
  
  w = ips[,"weight"]
  
  if ( is.null(mesh.coords) ) {
    ips.loc = ips
  } else {
    ips.loc = as.matrix(ips[,mesh.coords])
  }
  
  res = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = ips.loc)
  tri = res$p2m.t 
  
  # Filter integration points outside mesh
  if ( any(tri==0) ) {
    ips = ips[!tri==0,]
    w = w[!tri==0]
    res$p2m.b = res$p2m.b[!tri==0,]
    warning("Integration points outside mesh. Omitting.")
  }
  
  nw = w * res$p2m.b
  w.by = by(as.vector(nw), as.vector(mesh$graph$tv[tri,]), sum, simplify = TRUE)
  nw = as.vector(w.by)
  
  if ( is.null(mesh.coords)) {
    ret = data.frame(mesh$loc[as.numeric(names(w.by))], w = nw)
  } else {
    nips = mesh$loc[as.numeric(names(w.by)),1:length(mesh.coords)]
    colnames(nips) = mesh.coords
    ret = cbind(nips, weight = nw)
  }
  
  return(ret)
  
}


#' Expand integration points over additional dimensions
#' 
#' See int.quadrature on how to formulate the further arguments.
#'
#' @aliases int.expand
#' @export
#' @param ips Integration points to expand
#' @param ... lists with arguments for multiple calls of int.quadrature
#' @return Integration points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

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
  return(cbind(rgrd,rips))
}

