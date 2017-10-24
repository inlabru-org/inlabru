iinla.env = new.env()
iinla.env$log = sprintf("inlabru @ %s", date())



#' @title Global setting for tutorial sessions
#' 
#' @description Increases verbosity and sets the inference strategy to empirical Bayes.
#'
#' @aliases init.tutorial
#' @export
#' 
#' @return NULL
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' # Determine current bru default:
#' bo = bru.options()
#' 
#' # By default, INLA's integration strategy is set to the INLA default 'auto':
#' bo$inla.options$control.inla
#' 
#' # Now, let's run init.tutorial() to make empirical Bayes the default 
#' integration method when \code{bru} calls \code{inla}
#' 
#' init.tutorial()
#' 
#' # Check if it worked:
#' bru.options()$inla.options$control.inla
#' 
#' }
#' 
init.tutorial = function() {
  cat("Setting defaults for tutorial session. \n")
  iinla.setOption("iinla.verbose", list(TRUE))
  iinla.setOption("control.inla", list(list(int.strategy = "eb")))
  iinla.setOption("control.compute", list(list(config = TRUE, dic = TRUE, waic = TRUE)))
}

# @title Print inlabru log
# 
# @description Print inlabru log
#
# @aliases bru.log
# @export
# 
# @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
# 

bru.log = function() {
  sprintf(iinla.env$log)
}

msg = function(txt) {
  if ( iinla.getOption("iinla.verbose") ) { cat(paste0(txt,"\n")) }
  logentry(txt)
}

logentry = function(txt) { iinla.env$log = c(iinla.env$log, paste0(Sys.time(),": ", txt)) }

iinla.getOption = function (option = c("control.compute", "control.inla","iinla.verbose")) 
{
  if (missing(option)) 
    stop("argument is required.")
  envir = iinla.env
  option = match.arg(option, several.ok = TRUE)
  if (exists("iinla.options", envir = envir)) 
    opt = get("iinla.options", envir = envir)
  else opt = list()

  default.opt = list(iinla.verbose = FALSE,
                     control.compute = list(config = TRUE),
                     control.inla = list(int.strategy = "auto"))
  res = c()
  for (i in 1:length(option)) {
    if (option[i] %in% names(opt)) {
      res[[option[i]]] = opt[[option[i]]]
    }
    else {
      res[[option[i]]] = default.opt[[option[i]]]
    }
  }
  if (length(res) == 1) { res = res[[1]] }
  return(res)
}



iinla.setOption = function (...) 
{
  iinla.setOption.core = function(option = c("control.compute","control.inla","iinla.verbose"), value) {
    envir = iinla.env
    option = match.arg(option, several.ok = FALSE)
    if (!exists("iinla.options", envir = envir)) 
      assign("iinla.options", list(), envir = envir)
    if (is.character(value)) {
      eval(parse(text = paste("iinla.options$", option, 
                              "=", shQuote(value), sep = "")), envir = envir)
    }
    else {
      eval(parse(text = paste("iinla.options$", option, 
                              "=", ifelse(is.null(value), "NULL", value), 
                              sep = "")), envir = envir)
    }
    return(invisible())
  }
  called = list(...)
  len = length(names(called))
  if (len > 0L) {
    for (i in 1L:len) {
      do.call(iinla.setOption.core, args = list(names(called)[i], 
                                               called[[i]]))
    }
  }
  else {
    iinla.setOption.core(...)
  }
  return(invisible())
}
