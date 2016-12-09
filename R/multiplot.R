#' @title Multiple ggplots on a page.
#'
#' @description
#' Renders multiple ggplots on a single page. 
#'
#' @param ... Comma-separated \code{ggplot} objects.
#' @param plotlist A list of \code{ggplot} objects - an alternative to the comma-separated argument above.
#' @param cols Number of columns of plots on the page.
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored. If the layout is 
#' something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then plot 1 will go in the upper left, 
#' 2 will go in the upper right, and 3 will go all the way across the bottom.
#'
#' @source 
#' \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'  
#' @examples 
#' library(inlabru)
#' library(raster)
#' data(gorillanests)
#' plotpts = plotsample(gnests,gnestboundary,x.ppn=0.5,y.ppn=0.5,nx=5,ny=5)
#' p1 = ggplot() +gg(plotpts$plots) +gg(plotpts$dets) +gg(gnestboundary)
#' countdata = point2count(plotpts$plots,plotpts$dets)
#' x=coordinates(countdata)[,1]
#' y=coordinates(countdata)[,2]
#' count=countdata@data$n
#' p2 = ggplot() +gg(gnestboundary) +gg(plotpts$plots) +  geom_text(aes(label=count, x=x, y=y))
#' multiplot(p1,p2,cols=2)
#' 
#' @export
#
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}