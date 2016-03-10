# This filie contains several helpers to plot distance sampling data using ggplot2
# gg.mesh() : A geom for plotting the mesh
# gg.seg() : A geom for plotting segment lines
# gg.det() : A geom for plotting detections
# gg.swath() : A geom for plotting the segment area


#' Plot inla.mesh using ggplot2
#'
#' 
#' @aliases gg.mesh
#' @name gg.mesh
#' @export
#' @param data a \link{dsdata} object
#' @return A ggplot2 object
#' 

gg.mesh = function(data, color = rgb(0,0,0,0.1), ...) {
  if ( class(data)[1] == "inla.mesh" ) { mesh = data } else { mesh = data$mesh }
  df = rbind(data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,2],c(1,2)]),
               data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
               data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]))
  gg = geom_segment(data = df, aes(x=a.1,y=a.2,xend=b.1,yend=b.2), color = color, ...)
  return(gg)
}


#' Plot segments using ggplot2
#'
#' @aliases plot gg.seg
#' @export
#' @param data a \link{dsdata} object
#' @return A ggplot2 object
#' 

gg.seg = function(data, ...) {
  df = data.frame(startpoint(as.segment(data),data)[,data$mesh.coords], 
                  endpoint(as.segment(data),data)[,data$mesh.coords])
  colnames(df) = c("x","y","xend","yend")
  gg = geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend), ...)
  return(gg)
}


#' Plot detections using ggplot2
#'
#' 
#' @aliases gg.det
#' @export
#' @param data a \link{dsdata} object
#' @param color Color used to plot detections
#' @return A ggplot2 object
#' 

gg.det = function(data, color = "red", ...) {
  df = detdata(data)
  gg = geom_point(data = df, aes_string(x = data$mesh.coords[1], 
                                        y = data$mesh.coords[2]), 
                  color = color, ...)
  return(gg)
}

#' Plot segment swath
#'
#' 
#' @aliases gg.swath
#' @export
#' @param data a \link{dsdata} object
#' @param width The width of the swath
#' @param fill Filling color of the swath
#' @return A ggplot2 object
#' 

gg.swath = function(data, width = 1, fill = rgb(0,0,0,0.1), ...) {
  sp = startpoint(as.segment(data),data)
  ep = endpoint(as.segment(data),data)
  v = ep-sp
  vo = data.frame(x = v[,2],y = -v[,1])
  v1 = data.frame(sp + width * normalize.euc(vo), seg = 1:dim(vo)[1])
  v2 = data.frame(sp - width * normalize.euc(vo), seg = 1:dim(vo)[1])
  v3 = data.frame(ep - width * normalize.euc(vo), seg = 1:dim(vo)[1])
  v4 = data.frame(ep + width * normalize.euc(vo), seg = 1:dim(vo)[1])
  df = rbind(v1, v2, v3, v4)
  gg = geom_polygon(data = df, aes(group = seg), fill = fill)

  return(gg)
}



