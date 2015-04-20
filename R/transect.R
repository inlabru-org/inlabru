#' Plot transects
#'
#' @aliases plot.transect
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
plot.transect = function(tr,data=FALSE,rgl=TRUE,col=NULL,add=FALSE,radius=1.1,...){
  if ("geometry" %in% names(data)) { geometry=data$geometry } else {geometry="euc"}
  if (geometry=="geo"){
    if (is.null(col)) {col="black"}
    if (rgl){
      if (!add){
        #dev = open3d()
        rgl.earth()
      } else {
        #dev = NULL
      }
      # draw transects
      
      rgl.sphlines2(startpoint(tr,data),endpoint(tr,data),add=add,radius=radius,col=col,lwd=2)
    }
    else {
      xl = c(-13,32)
      yl = c(-160,-75)
      
      if (is.numeric(col)){
        grpcol = col
      } else {
        spoint = startpoint(tr)
        epoint = endpoint(tr)
        lgrp = length.group.transect(tr,max.arclen=Inf,max.ratio=1.1,n.lookahead=Inf)
        grpcol = 1+mod(lgrp,254)
      }   
      plot.new()
      for (k in 1:500){
        plot(c(spoint[k,"lat"],epoint[k,"lat"]),c(spoint[k,"lon"],epoint[k,"lon"]),xlim=xl,ylim=yl,type="l",xlab="",ylab="",col=grpcol[k],lwd=3)
        par(new=TRUE)
      }
    }
  }
  else if (geometry=="euc"){
    epoint = endpoint(tr,data)
    spoint = startpoint(tr,data)
    xlim = range(rbind(epoint[,1],spoint[,1]))
    ylim = range(rbind(epoint[,2],spoint[,2]))
    for (k in 1:dim(spoint)[1]){
      plot(c(spoint[k,1],epoint[k,1]),c(spoint[k,2],epoint[k,2]),type="l",xlab="",ylab="",lwd=3,...)# xlim=xlim,ylim=ylim
      par(new=TRUE)
    }
  }
  return(invisible())
}


# 
# GENERICS
#

startpoint = function(...){UseMethod("startpoint")}
endpoint = function(...){UseMethod("endpoint")}
len = function(...){UseMethod("len")}
tim = function(...){UseMethod("tim")}
gaplen = function(...){UseMethod("gaplen")}
numel = function(...){UseMethod("numel")}
linedata = function(...){UseMethod("linedata")}


#' Length of a transect line, from start point to end points.
#' 
#' ! Ignores gaps between segments
#' 
#' @aliases len.transect
#' @export
#' @param tr Transects
#' @param dsdata Distance sampling data set
#' @example \dontrun{ toy=toy1() ; plot(toy) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

len.transect = function(tr,data) {
  if (data$geometry == "euc") {
    return(dist.euc(startpoint(tr,data),endpoint(tr,data))) 
  } else if ((data$geometry == "geo")) {
    return(dist.geo(startpoint(tr,data),endpoint(tr,data)))
  } else {stop("Unknown geometry")}
}


#' Length of gaps between transects
#' 
#' @aliases gaplen.transect
#' @export
#' @param tr Transects
#' @param dsdata Distance sampling data set
#' @example \dontrun{ toy=toy1() ; plot(toy) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

gaplen.transect = function(tr,data,...) {
  spoints = startpoint(tr,data,...)
  epoints = endpoint(tr,data,...)
  if (data$geometry == "euc") {
    return(dist.euc(epoints[1:numel(tr)-1,],spoints[2:numel(tr),])) 
  } else if ((data$geometry == "geo")) {
    return(dist.geo(epoints[1:numel(tr)-1,],spoints[2:numel(tr),])) 
  } else {stop("Unknown geometry")}
}

#' Length of gaps between segments
#' 
#' @aliases gaplen.segment
#' @export
#' @param seg Transect segments
#' @param dsdata Distance sampling data set
#' @example \dontrun{ toy=toy1() ; plot(toy) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

gaplen.segment = function(seg,data,...) {
  spoints = startpoint(seg,data,...)
  epoints = endpoint(seg,data,...)
  if (data$geometry == "euc") {
    return(dist.euc(epoints[1:numel(tr)-1,],spoints[2:numel(tr),])) 
  } else if ((data$geometry == "geo")) {
    return(dist.geo(epoints[1:numel(tr)-1,],spoints[2:numel(tr),])) 
  } else {stop("Unknown geometry")}
}

tim.transect = function(tr,...) {
  return(rep(0,numel(tr))) 
}


#' Group transects up to a maximal arc length
#' 
#'
#' @aliases group.by.length.transect
#' @export 
#' @param tr Transects
#' @param data Data set providing information on transects
#' @param max.arclength Maximal arc length between start and end point of a group
#' @param max.ratio A factor describing the maximal deviation of the gruoped transect from the arc between the group's start and end point
#' @return \code{group} A vector of group indicators (integer)
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
group.by.length = function(...){UseMethod("group.by.length")}
group.by.length.transect = function(tr,data,max.arclen=1000,max.ratio=1.1,n.lookahead = 100) {
  transect.length = len(tr,data)
  gap.length = c(0,gaplen(tr,data))
  group = rep(0,numel(tr))
  n.tr = numel(tr)
  idx1 = 1
  g = 1
  
  while (TRUE) {
    endidx = min(n.tr,idx1+n.lookahead)
    travel.length = cumsum(transect.length[idx1:endidx]+c(0,gap.length[(idx1+1):endidx]))
    arc.length = dist.geo(startpoint(tr,data)[idx1,],endpoint(tr,data)[idx1:endidx,])
    arc.length.start = dist.geo(startpoint(tr,data)[idx1,],startpoint(tr,data)[idx1:endidx,])
    
    #is.smooth = (travel.length<=max.ratio*arc.length)
    is.smooth = (travel.length<=sqrt(max.ratio*arc.length))
    is.short = (arc.length < max.arclen)
    is.proper = c(1,arc.length[1:(length(arc.length)-1)]<arc.length[2:(length(arc.length))])
    is.proper2 = c(1,arc.length.start[1:(length(arc.length.start)-1)]<arc.length.start[2:(length(arc.length.start))])
    
    idx2 = idx1 + max(0,match(FALSE,is.smooth & is.short & is.proper & is.proper2,nomatch=2)-2)
    
    if (idx2>=n.tr){
      group[idx1:length(group)] = g
      break
    } 
    else {
      group[idx1:idx2] = g
      g = g+1
      idx1 = idx2+1
    } 
  }
  return(group)
}

#' Transform geographic transect coordinates to cylindric coordinates. 
#' 
#' The start point of the first transect end the end point of the last transect serve to define the equator of
#' the cylindric transformation
#' 
#'
#' @aliases to.cyl.transect
#' @export 
#' @param tr Transects
#' @param data Data set providing information on transects
#' @return A list with cylindric coordinates of the start ($startpoints) and end ($endpoints) points of the transects. 
#' $B is the basis used to rotate euclidian vectors such that the cylindric equator is the x-y plane in the new coordinate system.
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

to.cyl.transect = function(tr,data){
  spoints = startpoint(tr,data,keep=TRUE)
  epoints = endpoint(tr,data,keep=TRUE)
  B = basis.geo(spoints[1,],epoints[dim(epoints)[1],])
  
  cyl.spoints = geo.to.cyl(spoints,B=B)
  class.spoints = "cyl"
  cyl.epoints = geo.to.cyl(epoints,B=B)
  class.epoints = "cyl"

  return(list(startpoints=cyl.spoints,endpoints=cyl.epoints,B=B))
}



#' Nested grouping of transects
#'
#' @aliases nested.group
#' @export 
#' @param tr Transects
#' @param data Data
#' @param ... criteria to group by. Either strings (identifying a transect property column) or grouping functions
#' @return group group index
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
nested.group = function(...){UseMethod("nested.group")}
nested.group.transect = function(tr,data,crit,crit2=FALSE,...){
  
  if (is.character(crit)) { grp = as.numeric(as.factor(linedata(tr,data,crit))) } 
  else { grp = crit(tr,data)  }
  
  # nested ?
  if (nargs()==3) { return(grp) }
  else {
    nested.grp = list()
    new.grp = rep(0,length(grp))
    last.group = 0
    for (g in unique(grp)) {
      idx = grp==g
      nested.grp[[g]] = nested.group.transect(tr[idx,],data,crit2,...) + last.group
      last.group = max(nested.grp[[g]])
      new.grp[idx] = nested.grp[[g]]
    }
    return(new.grp)
  }
}

internal.nested.group = function(tr,data,crit){

}
