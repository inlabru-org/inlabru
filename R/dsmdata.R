import.dsmdata = function(dsmdata){
  
  # Extract data
  segdata <- dsmdata$segdata
  distdata <- dsmdata$distdata
  obsdata <- dsmdata$obsdata
  preddata <- dsmdata$preddata
  
  # 1. DEFINE BLOCKS
  # Take account of depth in defining blocks - Block.Label is identifier
  # More than one covariate can be listed
  # If no covariates, then covar.col=NA and transects are used (to be done)
  segdata <- define.blocks.f(seg=segdata,covar.col=6,geometry="euc")
  #segdata[1:2, ]
  
  # 2. GET THE ANGLE OF DIRECTION FOR EACH SEGMENT
  segdata <- get.direction.segment.f(data=segdata,geometry="euc")
  #segdata[1:2, ]
  
  # 3. GET THE START AND END POINTS OF ALL SEGMENTS
  # NOTE: x and y get renamed to mid.x and mid.y to avoid confusion with detection coordinates
  segdata <- start.end.points.segments.f(seg=segdata,use.tran=FALSE,geometry="euc")
  #segdata[1:2, ]
  
  # 4. AMALGAMATE SEGMENTS INTO NEW BLOCKS
  blocks <- get.blocks.f(seg=segdata,geometry="euc")
  #blocks[1:2, ]
  
  # 5. ADD SEGMENT AND BLOCK LABELS TO DETECTIONS
  # so that detections and blocks data can be combined
  distdata <- add.labels.to.obs.f(dists=distdata,obs=obsdata,seg=segdata)
  #distdata[1:2, ]
  
  
  newdata <- combine.dsmdata.f(blocks=blocks,dists=distdata)
  # Only one strata for this data so need to add this information
  newdata$strat <- 1
  #newdata[1:3, ]
  
  # Generate location for detections if location missing
  # Use do.plot=T to see the location generated for each detection
  det.new.coords <- generate.obs.location.f(seg=segdata,dists=distdata,geometry="euc",do.plot=F)
  #det.new.coords[1:3, ]
  
  
  # Prediction data to mesh
  if (require(splancs)) {
    loc = as.matrix(preddata[,c("x","y")])
    seg = inla.nonconvex.hull(loc, convex = -0.01)
    mesh = inla.mesh.create(interior = seg, refine = list(max.edge = (max(as.vector(loc)) - min(as.vector(loc)))/30))
  } else {
    mesh = NULL
  }
  
  
  dset = make.dsdata(effort = newdata, mesh = mesh)
  
  return(dset)
}





