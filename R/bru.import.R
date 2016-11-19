
import.mexdolphin = function() {
  
  library(dsm)
  data(mexdolphins)
  data.p4s = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  dset = import.dsmdata(mexdolphins, covar.col = 8)
  mexdolphin = as.spatial.dsdata(dset, cnames = c("x","y"), crs = CRS(data.p4s))
  
  # Target CRS
  target.p4s = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs +towgs84=0,0,0"
  
  # Units to km
  mexdolphin$points = spTransform(mexdolphin$points, CRS = CRS(target.p4s))
  mexdolphin$samplers = spTransform(mexdolphin$samplers, CRS = CRS(target.p4s))  
  mexdolphin$points$distance = mexdolphin$points$distance / 1000
  mexdolphin$mesh$loc = mexdolphin$mesh$loc/1000
  
  # Set mesh crs
  mexdolphin$mesh$crs = inla.CRS(target.p4s)
  
  # Coordnames is set to NULL by spTransofrm ??
  coordnames(mexdolphin$samplers) = c("x","y")
  
  # Set weights (transect width)
  mexdolphin$samplers$weight = 16
  
  ##### Prediction polygon
  polyloc = as.data.frame(mexdolphin$mesh$loc[mexdolphin$mesh$segm$int$idx,c(1,2)])
  colnames(polyloc) = c("x","y")
  po = Polygon(polyloc, hole=FALSE)
  pos = Polygons(list(po), ID = "c")
  predpoly = SpatialPolygons(list(pos))
  df = data.frame(weight = 1)
  rownames(df) = "c"
  predpolyd = SpatialPolygonsDataFrame(predpoly, data = df)
  # plot(predpolyd)
  mexdolphin$ppoly = predpolyd
  proj4string(mexdolphin$ppoly) = proj4string(mexdolphin$points)
  
  ##### Simulate a whole population #####
  distance = seq(0, 8, length.out = 20)
  fml = coordinates + distance ~ df.lsigma + g(spat, model = inla.spde2.matern(mexdolphin$mesh), mesh = mexdolphin$mesh)
  pred = expression(log(1-exp(-(distance/(exp(df.lsigma)))^-1)) + spat + Intercept)
  init.tutorial()
  r = lgcp(points = mexdolphin$points, samplers = mexdolphin$samplers, model = fml, predictor = pred, mesh = mexdolphin$mesh)
  llambda = log(8) + r$summary.random$spat$mean + r$summary.fixed["Intercept","mean"]
  # sum( exp(llambda) * diag(inla.mesh.fem(mexdolphin$mesh)$c0))
  smexdolphin = mexdolphin
  smexdolphin$llambda = llambda
  pts = data.frame(sample.lgcp(mexdolphin$mesh, llambda))
  
  ### remove samples outside prediction polygon. The sp package over() method somehow does not work.
  loc = mexdolphin$ppoly@polygons[[1]]@Polygons[[1]]@coords
  seg = inla.mesh.segment(loc[rev(1:nrow(loc)),])
  msh = inla.mesh.2d(boundary=seg, max.edge = 100)
  # plot(msh)
  # points(pts, col = 1+is.inside(msh, coordinates(pts)))
  pts = pts[is.inside(msh, coordinates(pts)),]
  coordinates(pts) = c("x","y")
  proj4string(pts) = proj4string(mexdolphin$points)
  mexdolphin$simulated = SpatialPointsDataFrame(pts, data.frame(size = rep(1, length(pts))))
  mexdolphin$lambda = exp(llambda)
  
  #### Depth covariate #####
  depth = mexdolphins$preddata[,c("x","y","depth")]
  coordinates(depth) = c("x","y")
  proj4string(depth) = data.p4s
  # gg.map(depth) + gg.point(depth, CRS("+proj=longlat")) + 
  # gg.polygon(mexdolphin$ppoly, CRS("+proj=longlat"))
  mexdolphin$depth = spTransform(depth, CRS(target.p4s))
 
  # return
  mexdolphin 
}
