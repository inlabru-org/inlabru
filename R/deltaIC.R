#' @title Summarise DIC and WAIC from \code{lgcp} objects.
#'
#' @description
#' Calculates DIC and WAIC differences and produces an ordered summary. 
#'
#' @param ... Comma-separated objects inheriting from class \code{inla} and obtained from a run of \link[INLA]{inla}, \link{bru} or \link{lgcp}
#' @param criterion If 'DIC', plots DIC differences; If 'WAIC', plots WAIC differences.
#'
#' @return A data frame with each row containing the model name, DIC, WAIC, deltaDIC, and
#' deltaWAIC.
#'  
#' @export
#' 
#' @examples
#' 
#' \donttest{
#' # Generate some data
#' input.df <- data.frame(x=cos(1:10))
#' input.df <- within(input.df, y <- 5 + 2*cos(1:10) + rnorm(10, mean=0, sd=0.1))
#' 
#' # Fit two models
#' fit <- bru(y ~ x, "gaussian", input.df)
#' fit2 <- bru(y ~ x, "Poisson", input.df)
#' 
#' # Compare DIC
#' 
#' deltaIC(fit, fit2)
#' }

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
