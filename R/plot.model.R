
# Play (animate) spatial field
#
# Animates a spatial field using RGL. 
# 
# @aliases play.spatial
# @export
# @param group Example: group = list(year = c(1,2)) animates the field for years 1 and 2
# @param ... Parameters passed on to \link{plot.spatial}
#

play.spatial = function(group = list(), rgl, ...){
  if ( !rgl ) { rgl = TRUE}
  globe()
  sargs = list(...)
  myanim = function(time, ...) {
    rgl::par3d(skipRedraw = TRUE)
    grp = list()
    grp[[names(group)[[1]]]] = group[[1]][(floor(time) %% 2)+1]
    do.call(pixelplot.mesh, c(sargs, list(group = grp, add = TRUE, rgl = TRUE))) ;
    rgl::par3d(skipRedraw = FALSE)
    return("")
  }
  
  rgl::play3d(myanim, duration = 10, startTime = 0, fps = 1)
  
}