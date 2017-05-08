
int.slines = function(data, mesh, group = NULL, project = TRUE) {
  
  # Extract start and end coordinates
  qq = coordinates(data)
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  
  idx = do.call(rbind, lapply(1:length(qq), function(k) do.call(cbind, lapply(qq[[k]], function(x) rep(k, nrow(x)-1) ))))
  idx = cbind(idx, idx)
  
  # Filter out points outside the mesh...
  loc = as.matrix(rbind(sp,ep))
  t1 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  if (any(t1==0) | any(t2==0)) { warning("Found spatial lines with start or end point ouside of the mesh. Omitting.")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0)),]
  
  # Split at mesh edges
  line.spl = split.lines(mesh, sp, ep, TRUE)
  sp = line.spl$sp
  ep = line.spl$ep
  idx = idx[line.spl$split.origin,]
  
  # Determine integration points along lines
  
  sp3d = within(data.frame(sp), Z <- 0)
  coordinates(sp3d) = c("X1","X2","Z")
  proj4string(sp3d) = proj4string(data)
  ep3d = within(data.frame(ep), Z <- 0)
  coordinates(ep3d) = c("X1","X2","Z")
  proj4string(ep3d) = proj4string(data)
  
  sp3d = spTransform(sp3d, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  ep3d = spTransform(ep3d, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  mp3d = SpatialPoints((coordinates(sp3d) + coordinates(ep3d))/2, proj4string = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  
  ips = coordinates(spTransform(mp3d, CRS(proj4string(data))))[,c(1,2)]
  w = spDists(coordinates(spTransform(sp3d, CRSobj = CRS("+proj=longlat")))[,c(1,2)],
              coordinates(spTransform(ep3d, CRSobj = CRS("+proj=longlat")))[,c(1,2)],
              diagonal = TRUE, longlat = TRUE)
  
  # Wrap everything up and perform projection according to distance and given group argument
  ips = data.frame(ips)
  colnames(ips) = c("x","y")
  
  # Weights
  ips = cbind(ips, weight = w)  
  if ( "weight" %in% names(data) ) { ips$weight = ips$weight * data$weight[idx[,1]] }
  
  
  # if ( !is.null(group) ) { ips = cbind(ips, as.data.frame(data)[idx[,1],group,drop=FALSE]) }
  
  coordinates(ips) = c("x","y")
  if (!is.null(coordnames(data))) coordnames(ips) = coordnames(data)
  proj4string(ips) = proj4string(data)
  
  # Project to mesh vertices
  if ( project ) ips = vertex.projection(ips, mesh, columns = "weight", group = group)
  
  ips
}

# Determine integration points 
# 
# @aliases polygon.integration
# @export
# @param data A \link{dsdata} object
# @param poly A function turning segments into polygons, e.g. \link{swath}
# @param project If TRUE, project integration points to mesh vertices
# @param group A grouping for the polygons, e.g. "year" for temporal integration schemes. This can also be a function mapping segments to a group.
# @param ... Additional arguments passed on to the poly function.
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

polygon.integration = function(data, poly, project = TRUE, group = NULL, ...){
    poly.args = list(...)
    poly = do.call(poly, c(list(data), poly.args))
    loc = poly[,data$mesh.coords]
    ips = int.polygon(data$mesh, loc, group = poly$seg)
    colnames(ips)[match("group", colnames(ips))] = "seg"

    if ( project ) {
      # Was a grouping supplied? If so, calculate the group index for each integration point
      if ( is.null(group) ) {
        ips.grp = rep(1, nrow(ips))
        group.name = "int.group"
      } else if ( is.character(group) ) {
        ips.grp = segdata(data)[match(ips$seg, segdata(data)$seg), group]
        group.name = group
      } else {
        ips.grp = group(segdata(data)[match(ips$group, segdata(data)$seg),])
        group.name = colnames(ips.grp)[1]
      }
      # Project the integration points for each group independently
      ips = lapply(unique(ips.grp), function(g){  
                  gips = data.frame(project.weights(ips[ips.grp==g,], data$mesh, mesh.coords = data$mesh.coords))
                  gips[[group.name]] = g
                  return(gips)
                  })
      ips = do.call(rbind,ips)
    }
}
  



# Integration points for polygons inside an inla.mesh
# 
# @aliases int.polygon
# @export
# @param mesh An inla.mesh object
# @param loc Locations defining the polygons
# @param group If loc defines multiple polygons then this is the ID of the group for each location in loc
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.polygon = function(mesh, loc, group = NULL){
  
if ( is.null(group) ) { group = rep(1, nrow(loc)) }
ipsl = list()
# print(paste0("Number of polygons to integrate over: ", length(unique(group)) ))
for ( g in unique(group) ) {
  gloc = loc[group==g, ]
  # Check where the polygon intersects with the mesh edges
  sp = gloc[1:(nrow(gloc)-1),]
  ep = gloc[2:nrow(gloc),]
  sloc = split.lines(mesh, sp, ep, filter.zero.length = FALSE)$split.loc[,1:2]
  if (!is.null(sloc)){ colnames(sloc) = colnames(loc) }
  
  # plot(mesh) ; points(sp) ; points(ep) ; points(sloc)
  bloc = rbind(gloc)
  bnd = inla.mesh.segment(loc = bloc)
  imesh = inla.mesh.create(boundary = bnd, loc = rbind(mesh$loc[,1:2], sloc[,1:2]))
  # plot(imesh) ; points(sp) ; points(ep) ; points(gloc)
  # plot(imesh) ; points(sp) ; points(ep) ; points(gloc) ; plot(mesh, add = TRUE)
  
  ips = data.frame(imesh$loc[,1:2])
  colnames(ips) = c("x","y")
  ips$weight = diag(as.matrix(inla.mesh.fem(imesh)$c0))
  # ips = as.data.frame(project.weights(ips, mesh, mesh.coords = c("x","y")))
  ips$group = g
  ipsl = c(ipsl, list(ips))
}
return(do.call(rbind,ipsl))
}



