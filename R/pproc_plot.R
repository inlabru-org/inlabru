#' Plot a LGCP model.
#'
#' 
#' @aliases plot plot.pproc_lgcp
#' @export
#' @param mdl A \link{pproc_lgcp} model
#' 

plot.pproc_lgcp = function(mdl,rgl=TRUE,field="spde",property="mean",add=FALSE,colorbar=TRUE,int.points=FALSE,logscale=TRUE,col=NULL,...){
  if (any("geometry" %in% names(mdl$data))){ geometry = mdl$data$geometry }
  
  sgh.points = detdata(mdl$data)
  if (rgl) {
    require(rgl)
    #if (is.null(field)) { col = logintensity.pproc(mdl,points=loc,covariates=list(sst = sstfun, clon=clon,clat=clat),property=property)}
    if (field %in% names(mdl$INLA$result$summary.ran)){
      col = mdl$INLA$result$summary.ran[[field]][[property]]
    } else {
      ind = inla.stack.index(mdl$stack, field)$data
      col = mdl$INLA$result$summary.fitted.values[ind,property]
    }
    
    if (geometry=="geo"){
      if (!add){ rgl.open() }
      bg3d(color="white")
      rgl.earth()
      rgl.sphmesh(mdl$data$mesh,add=TRUE,radius=1.001,col=col)
      rgl.sphpoints(long=sgh.points$lon+360,lat=sgh.points$lat,radius=1.01,col="red",size=5)
      if ( int.points ) { rgl.sphpoints(long=mdl$int.points$lon+360,lat=mdl$int.points$lat,radius=1.01,col=rgb(0,0,1),size=3) }
      if (colorbar){
        # This is a temporary workaround, rgl people are working on something like
        cb.col = cm.colors(100)
        rgl.linestrips(x=-0.85,y=0.85,z=seq(-1,1,length.out=length(cb.col)),col=cb.col,lwd=50)
        rgl.texts(x=-0.79,y=0.79,z=c(-0.97,0.97),c(format(min(col),digits=3),format(max(col),digits=3)),col="black",cex=1)
      }
    }
    else {
      plot(mdl$data$mesh,rgl=rgl,col=col,...)
      rgl.points(x=sgh.points$x,y=sgh.points$y,z=0.02,col=rgb(1,0,0),size=8)
      rgl.points(x=mdl$int.points$x,y=mdl$int.points$y,z=0.01,col=rgb(0,0,0.6),size=4)
    }
    if (any("par3d.args" %in% names(mdl$data))) { do.call(par3d,mdl$data$par3d.args) }
  }
  else {
    require(lattice)
    proj <- inla.mesh.projector(mdl$data$mesh,dims = c(300,300))
    xlim = range(mdl$data$mesh$loc[,1])
    ylim = range(mdl$data$mesh$loc[,2])
    x = seq(xlim[1],xlim[2],length.out=300)
    y = seq(ylim[1],ylim[2],length.out=300)
    grid = expand.grid(x=x,y=y)
    if (!is.null(col)) {
      A = inla.spde.make.A(mdl$data$mesh,loc=cbind(grid$x,grid$y))
      col = A%*%as.vector(col)
      msk = apply(abs(A),MARGIN=1,sum)>0
      col[!msk] = NA
    }
    else if (field %in% names(mdl$INLA$result$summary.ran)){
      col = inla.mesh.project(proj, field=mdl$INLA$result$summary.ran[[field]][[property]])
    } else {
      ind = inla.stack.index(mdl$stack, field)$data
      if (property == "qfrac") {
        col1 = inla.mesh.project(proj, field=mdl$INLA$result$summary.fitted.values[ind,"0.975quant"])
        col2 = inla.mesh.project(proj, field=mdl$INLA$result$summary.fitted.values[ind,"0.025quant"])
        col = col2-col1
        
      } else {
        
        col = inla.mesh.project(proj, field=mdl$INLA$result$summary.fitted.values[ind,property])
      }
    }
      
    if (!logscale) { col = exp(col) }
    
    print(paste0("min:", min(col,na.rm=TRUE),", max:",max(col,na.rm=TRUE), ", mean:",mean(col,na.rm=TRUE)))
    
    sgh.points = detdata(mdl$data)
    co1 = mdl$data$mesh.coords[1]
    co2 = mdl$data$mesh.coords[2]
    levelplot(col ~ grid$x + grid$y,col.regions=topo.colors(100),
              panel=function(...){
                panel.levelplot(...)
                if ( int.points ) { panel.points(x=mdl$int.points[,co1],y=mdl$int.points[,co2], pch=".",col="blue",cex=2) }
                panel.points(x=sgh.points[,co1],y=sgh.points[,co2], col=rgb(0,0,0),pch=".",cex=2)
              },
              xlab = co1, ylab = co2,...
    ) 
  }
}