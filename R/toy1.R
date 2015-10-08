#' A toy data set
#' 
#' @name toy1
NULL

io_toy1.getDataDir = function() {return(system.file("data",package="iDistance"))}

io_toy1.pkgdata.save = function(){
  ## save the data we will include in the R package
  toy1 = generate.toy1()
  save(toy1,file=paste0(io_toy1.getDataDir(),"/toy1.RData"))
}

generate.toy1 = function(lambda=10,dfun=function(x){exp(-(0.5*(x)^2)*2)},truncation = 3){
  tr1s1 = data.frame(strat="1",trans="1.1",
                     seg="1.1.1",start.x=-1,start.y=1,end.x=7,end.y=6,det=NA)
  
  toy = list()
  toy$geometry = "euc"
  effort = rbind(tr1s1); 
  class(effort) = c("effort","data.frame")
  toy$effort = effort
  detections = make.detections(as.segment(effort),toy,lambda=lambda,dfun=dfun,truncation=truncation) #plot(detections[,c("x","y")],xlim=c(0,7),ylim=c(-2,5))
  
  # join effort with detections
  effort = join(effort,detections)
  
  # sort columns of effort
  other.idx = !(names(effort) %in% c("strat","trans","seg","det"))
  effort = cbind(effort[,c("strat","trans","seg","det")],effort[,other.idx])
  
  # set class
  class(effort) = c("effort","data.frame")
  toy$effort = effort
  
  # Mesh
  lattice = inla.mesh.lattice(x = seq(-3, 9, length.out=7),
                              y = seq(-3, 9, length.out=5))
  toy$mesh = do.call(inla.mesh.create,list(lattice = lattice, refine = list(max.edge = 3)))
  
  toy$mesh.coords = c("x","y")
  
  class(toy) = c("toy1","dsdata","list")
  return(toy)
}



make.detections = function(seg,data,lambda=1,dfun=function(x){1},truncation=1.5){
  spoints = startpoint(seg,data)
  epoints = endpoint(seg,data)
  detections = list()
  ids = id(seg,data)
  for (k in 1:nrow(spoints)){
    spoint = spoints[k,] ; class(spoint) = class(spoints)
    epoint = epoints[k,] ; class(epoint) = class(epoints)
    trlen = dist(spoint,epoint)
    lambda.total = 2*lambda*truncation*trlen
    
    # Sample the number of potential detections
    n.smp = rpois(1,lambda.total)
    # Sample observation distance for each potential detection (signed for side of transect line)
    distance = truncation*2*(runif(n.smp)-0.5)
    # Sample relative location along transect line
    loc = sort(runif(n.smp))
    # Filter out detections not acutally observed 
    dprob = dfun(abs(distance))
    do.observe = dprob>runif(n.smp)
    loc = loc[do.observe]
    distance = distance[do.observe]
    # Calculate detection locations
    locvec = ssum(epoint,-spoint)
    orthvec = normalize.euc(data.frame(x=locvec$y,y=-locvec$x));
    orthvec = orthvec[rep(1,length(distance)),]
    pts = spoint[rep(1,length(distance)),] + orthvec*distance + loc*locvec[rep(1,length(distance)),]
    obs.loc = spoint[rep(1,length(distance)),] + loc*locvec[rep(1,length(distance)),]
    colnames(obs.loc) = c("obs.x","obs.y")
    seg.id = ids[k]
    det.id = paste0(seg.id,".",1:length(distance))
    detections[[k]] = cbind(seg=seg.id,det=det.id,pts,distance=abs(distance), obs.loc)
  }
  return(do.call(rbind,detections))  
}
