#
# Load the data points (whale sightings, transects and mesh designed by Joyce)
#
  source("../R/io_whales.R")
  whales = io_whales.pkgdata.load()

#
# run LGCP 
#
  require(INLA)
  source("../R/pproc_geometry.R")
  source("../R/pproc_lgcp.R")
  result = pproc_lgcp(whales,int.subsample=TRUE)

#
# Visualize result
#
  plot(whales$mesh)

  require(lattice)
  proj <- inla.mesh.projector(whales$mesh,dims = c(300,300))
  x1.mean <- inla.mesh.project(proj, field=mdl$INLA$result$summary.ran$i$mean)
  par(new=FALSE)
  levelplot(x1.mean,xlab='', ylab='',col.regions=topo.colors(100), scale=list(draw=FALSE))
