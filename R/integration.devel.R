#' Determine integration points 
#' 
#' @aliases polygon.integration
#' @export
#' @param data A \link{dsdata} object
#' @param poly A function turning segments into polygons, e.g. \link{swath}
#' @param project If TRUE, project integration points to mesh vertices
#' @param group A grouping for the polygons, e.g. "year" for temporal integration schemes. This can also be a function mapping segments to a group.
#' @param ... Additional arguments passed on to the poly function.
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

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
  



#' Integration points for polygons inside an inla.mesh
#' 
#' @aliases int.polygon
#' @export
#' @param mesh An inla.mesh object
#' @param loc Locations defining the polygons
#' @param group If loc defines multiple polygons then this is the ID of the group for each location in loc
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

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

  imesh = inla.mesh.2d(loc.domain = sloc, 
                       max.edge =  sqrt(2*max(diff(range(sloc[,1])), diff(range(sloc[,2])))^2))
  
  ii = is.inside(imesh, mesh$loc)

  imesh = inla.mesh.2d(loc.domain = sloc,
                       mesh$loc[ii,],
                       max.edge =  sqrt(2*max(diff(range(sloc[,1])), diff(range(sloc[,2])))^2))
  
  
  ips = data.frame(imesh$loc[,1:2])
  colnames(ips) = c("x","y")
  ips$weight = diag(as.matrix(inla.mesh.fem(imesh)$c0))
  ips = as.data.frame(project.weights(ips, mesh, mesh.coords = c("x","y")))
  ips$group = g
  ipsl = c(ipsl, list(ips))
}
return(do.call(rbind,ipsl))
}


circle = function(data, radius, n = 6) {
  
  mp = segdata(data)[,paste0("start.",dset$mesh.coords)]
  seg = segdata(data)[,"seg"]
  
  df = list()
  for (k in 1:length(seg)){
    df[[k]] = data.frame(x = radius*sin(seq(0,2*pi, length.out=n)),
                         y = radius*cos(seq(0,2*pi, length.out=n)),
                         seg = seg[k])
  }
  df = do.call(rbind, df)
  colnames(df) = c(data$mesh.coords, "seg")
  return(df)
}


swath = function(data, width = NULL) {
  
  sp = segdata(data)[,paste0("start.",dset$mesh.coords)]
  ep = segdata(data)[,paste0("end.",dset$mesh.coords)]
  seg = segdata(data)[,"seg"]
  
  if ( is.null(width) ) { width = segdata(data)$width }
  
  colnames(sp) = data$mesh.coords
  colnames(ep) = data$mesh.coords
  
  v = ep-sp
  vo = data.frame(x = v[,2],y = -v[,1])
  v1 = data.frame(sp + width/2 * normalize.euc(vo), seg)
  v2 = data.frame(sp - width/2 * normalize.euc(vo), seg)
  v3 = data.frame(ep - width/2 * normalize.euc(vo), seg)
  v4 = data.frame(ep + width/2 * normalize.euc(vo), seg)
  df = rbind(v1, v4, v3, v2)
  colnames(df) = c(data$mesh.coords, "seg")
  return(df)
}



