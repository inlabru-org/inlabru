#' @title Summarise DIC and WIC from \code{lgcp} objects.
#'
#' @description
#' Calculates DIC and WAIC differences and produces an ordered summary. 
#'
#' @param ... Comma-separated objects of class \code{lgcp}
#' @param order.by If 'DIC', orders output by DIC differences; If 'DIC', orders output 
#' by DIC differences.
#'
#' @return A data frame with each row containing the model name, DIC, WAIC, deltaDIC, and
#' deltaWAIC.
#'  
#' @export
deltaIC = function(...,order.by="DIC"){
  if(order.by != "DIC" & order.by != "WAIC") {
    warning("Invalid order.by argument: ordering by DIC")
    order.by = DIC
  }
  names <- as.character(substitute(list(...)))[-1L]
  model = eval(list(...))
  nmod = length(model)
  dic = waic = rep(NA,nmod)
  for(i in 1:nmod) {
    mod = model[[i]]
    if(!inherits(mod,"lgcp")) stop("Object ",i," is not of class lgcp.")
    if(is.null(mod$dic$dic)) stop("Object ",i," does not have a DIC.")
    if(is.null(mod$waic$waic)) stop("Object ",i," does not have a WAIC.")
    dic[i] = mod$dic$dic
    waic[i] = mod$waic$waic
  }
  if(order.by == "DIC") ord = order(dic)
  if(order.by == "WAIC") ord = order(waic)
  dic = dic[ord]
  waic = waic[ord]
  names = names[ord]
  ddic = dic - min(dic)
  dwaic = waic - min(waic)
  outdat = data.frame(Model=names,DIC=dic,WAIC=waic,Delta.DIC=ddic,Delta.WAIC=dwaic)
  return(outdat)
}
