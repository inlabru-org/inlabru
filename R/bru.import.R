
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


save.mexdolphin = function() {
  mexdolphin = import.mexdolphin()
  save("mexdolphin", file=paste0(system.file("data",package="inlabru"),"/mexdolphin.RData"))
}


import.seals = function(sealfile = "WestIce2012.csv", icefile = "reflectance_0.0025deg_grid_modis_20120328_1310.tif") {
  
  #' Load seal data  
  
  seals = read.csv(sealfile)
  
  #' Turn seal data into spatial object
  
  coordinates(seals) = c("longitude","latitude")
  proj4string(seals) = CRS("+proj=longlat")
  
  #' To avoid confusion I remove the x-y coordinates created by Martin
  
  seals = seals[,setdiff(names(seals), c("x","y"))]
  
  #' Change CRS
  
  target.p4s = "+proj=utm +zone=27 +units=km"
  seals = spTransform(seals, CRS(target.p4s))
  coordnames(seals) = c("x","y")
  
  #' Select a strip
  
  seals = seals[(coordinates(seals)[,2]>7947) & (coordinates(seals)[,2]<7954), ] # strip 9
  # seals = seals[(coordinates(seals)[,2]>7931) & (coordinates(seals)[,2]<7938.5), ] # strip 12
  # seals = seals[(coordinates(seals)[,2]>7925) & (coordinates(seals)[,2]<7930), ] # strip 13
  #seals = seals[(coordinates(seals)[,2]>7920) & (coordinates(seals)[,2]<7924), ] # strip 14
  #seals = seals[(coordinates(seals)[,2]>7910) & (coordinates(seals)[,2]<7920), ] # strip 15
  
  #' Add total number of seals
  
  seals$all = seals$harps + seals$hooded
  # plot(seals)
  
  #' Build a mesh. This mesh will be fine at the photo locations but coarse elsewhere
  
  bnd = inla.nonconvex.hull(coordinates(seals), resolution = 170, convex = 0.5)
  bnd2 = inla.nonconvex.hull(coordinates(seals), resolution = 100, convex = 0.7)
  mesh = inla.mesh.2d(boundary = list(bnd, bnd2), max.edge = c(0.2,3))
  mesh$crs = inla.CRS(projargs = CRS(target.p4s))
  # ggplot() + gg(mesh) + gg(seals) + coord_equal()
  # mesh$n
  
  #' Let's plot the observed counts
  
  # ggplot() + gg(mesh) + gg(seals, mapping = aes(longitude, latitude, color = log(all/area)), size = 3) + coord_equal()
  
  #' Ice Covariate
  
  ice = rgdal::readGDAL(icefile)
  ice = spTransform(ice, CRS(target.p4s))
  ii = is.inside(mesh, coordinates(ice))
  ice = ice[as.vector(ii),]
  
  # ggplot() +  gg(ice, mapping = aes(x, y, color = band1), size = 1) + gg(mesh) + coord_equal() + scale_color_gradientn(colours = topo.colors(100)) 
  
  #' Interpolate ice covariate
  
  icecv = covariate(ice, predictor = band1, mesh = mesh)
  plot(icecv)
  
  #' Add band1 covariate to seals data frame
  
  ice.band1 = evaluator(icecv)
  seals$ice = ice.band1(x=coordinates(seals)[,1], y = coordinates(seals)[,2])
  
  #' Plot seal count and ice together
  # 
  #   plot.spatial(mesh = mesh, col = data.frame(mode = icecv$values), property = "mode", nx = 2000) +
  #     scale_fill_gradientn(colours = topo.colors(100)) +
  #     gg(seals, mapping = aes(x, y, color = log(all/area)), size = 3) +
  #     scale_color_gradientn(colours = heat.colors(100))
  
  
  #' Create a data set
  seals = list(mesh = mesh, points = seals, ice.data = ice, ice.cv = icecv)
}

save.seals = function(...) {
  # save.seals("/home/fbachl/devel/r/seals/WestIce2012.csv", "/home/fbachl/devel/r/seals/reflectance_0.0025deg_grid_modis_20120328_1310.tif")
  seals = import.seals(...)
  save("seals", file=paste0(system.file("data",package="inlabru"),"/seals.RData"))
}


import.gorillas = function() {
  
  # Load Gorilla data from spatstat
  data(gorillas,package="spatstat")
  
  # Create SpatialPoints representing nest locations
  nests = as.data.frame(gorillas)
  coordinates(nests) = c("x","y")
  proj4string(nests) = CRS("+proj=utm +zone=32N +datum=WGS84") # from the Gorillas help file
  
  #' Turn the observation window into spatial polygon
  boundary = spoly(points = gorillas$window$bdry, crs = CRS("+proj=utm +zone=32N +datum=WGS84"))
  
  #' Build mesh
  bnd = inla.mesh.segment(loc = data.frame(gorillas$window$bdry[[1]]))
  mesh = inla.mesh.2d(interior = bnd, max.edge = 250)
  mesh$crs = inla.CRS(proj4string(nests))
  
  #' Turn covariates int SpatialGridDataFrame
  gcov = list()
  for ( nm in names(gorillas.extra) ) { 
    gcov[[nm]] = as(gorillas.extra[[nm]], "SpatialGridDataFrame") 
    proj4string(gcov[[nm]]) = proj4string(nests)
    coordnames(gcov[[nm]]) = c("x","y")
    names(gcov[[nm]]) = nm
  }
  
  # Create a plot sampling data set
  set.seed(121)
  plotpts = plotsample(nests, boundary, x.ppn=0.6, y.ppn=0.6, nx=5.4, ny=5.4)
  counts = point2count(plotpts$plots,plotpts$dets)
  x = coordinates(counts)[,1]
  y = coordinates(counts)[,2]
  count = counts@data$n
  
  # Make gam data frame
  gnestcount_9x9_60pc = data.frame(x = x, y = y, count = count, exposure = counts$area)
  gnestplots_9x9_60pc = plotpts$plots
  gnestdets_9x9_60pc = plotpts$dets
  sample_9x9_60pc = list( counts = gnestcount_9x9_60pc, 
                          plots = gnestplots_9x9_60pc, 
                          nests = gnestdets_9x9_60pc)
  
  # plot to check:
  # ggplot() +gg(boundary) +gg(gnestplots_9x9_60pc) +  gg(gnestdets_9x9_60pc,pch="+",cex=4) +
  #   geom_text(aes(label=count, x=x, y=y)) + coord_fixed()
  # 
  
  ### Make final gorilla data set
  gorillas = list(nests = nests,
                  mesh = mesh,
                  boundary = boundary,
                  plotsample = sample_9x9_60pc,
                  gcov = gcov)
  
  return(gorillas)
}


save.gorillas = function() {
  gorillas = import.gorillas()
  save("gorillas", file=paste0(system.file("data",package="inlabru"),"/gorillas.RData"))
}

