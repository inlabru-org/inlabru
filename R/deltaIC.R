#' @title Summarise DIC and WIC from \code{lgcp} objects.
#'
#' @description
#' Calculates DIC and WAIC differences and produces an ordered summary. 
#'
#' @param ... Comma-separated objects of class \code{lgcp}
#' @param criterion If 'DIC', plots DIC differences; If 'DIC', plots DIC differences.
#'
#' @return A data frame with each row containing the model name, DIC, WAIC, deltaDIC, and
#' deltaWAIC.
#'  
#' @export
deltaIC = function(...,criterion="DIC"){
  if(criterion != "DIC" & criterion != "WAIC") {
    warning("Invalid criterion argument: using DIC")
    criterion = "DIC"
  }
  names <- as.character(substitute(list(...)))[-1L]
  model = eval(list(...))
  nmod = length(model)
  dic = waic = rep(NA,nmod)
  for(i in 1:nmod) {
    mod = model[[i]]
    if(criterion=="DIC" & is.null(mod$dic$dic)) stop("Object ",i," does not have a DIC.")
    if(criterion=="WAIC" & is.null(mod$waic$waic)) stop("Object ",i," does not have a WAIC.")
    dic[i] = mod$dic$dic
    waic[i] = mod$waic$waic
  }
  if(criterion == "DIC") ord = order(dic)
  if(criterion == "WAIC") ord = order(waic)
  dic = dic[ord]
  waic = waic[ord]
  names = names[ord]
  ddic = dic - min(dic)
  dwaic = waic - min(waic)
  if(criterion=="DIC") 
    outdat = data.frame(Model=names,DIC=dic,Delta.DIC=ddic)
  if(criterion=="WAIC") 
    outdat = data.frame(Model=names,DIC=waic,Delta.DIC=dwaic)
  return(outdat)
}
