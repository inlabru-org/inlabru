
#' ETP distance sampling data format
#' 
#' Fields:
#' - effort
#' - sighting
#' - geometry = "euc"
#' - mesh.coords = c("lon","lat")
#' 
#' @name dsdata
#' @example whales
NULL

detdata.etpdata = function(data,detection=NULL,...){ 
  if (is.null(detection)) {
    return(data$sighting)
  } 
  else {
    return(data$sighting[detection,])
  }
}


as.transect.etpeffort = function(effort){
  #na.rows = which(is.na(effort$cruise))
  #start.idx = c(1,na.rows[1:length(na.rows)-1]+1)
  #end.idx = na.rows-1
  end.idx = which(effort$effort=="E")
  start.idx = c(1,end.idx[1:length(end.idx)-1]+1)
  #
  tr = data.frame(start=start.idx,end=end.idx)
  class(tr) = c("etptransect","transect","data.frame")
  return(tr)
}


startpoint.etptransect = function(tr,data=data,keep=FALSE) {
  if (keep){ return(data$effort[tr$start,]) }
  else { return(data$effort[tr$start,c("lat","lon")])}
}

endpoint.etptransect = function(tr,data=data,keep=FALSE) {
  if (keep){ return(data$effort[tr$end,]) } 
  else { return(data$effort[tr$end,c("lat","lon")]) } 
}

numel.etptransect = function(tr,data=data) { return(dim(tr)[1]) }

linedata.etptransect = function(tr,data,fields){
  return(data$effort[tr$start,fields])
}
