#' Prediction from fitted inla model
#' 
#' Takes a fitted inla object produced by the function inla() and produces predictions given a 
#' new set of values for the model covariates or the original values used for the model fit. 
#' The predictions can be based on any R expression that is valid given these values/covariates 
#' and the posterior of the estimated effects.
#' 
#' IMPORTANT: The inla object provided has to have an additional field called "formula". This
#' is the formula used to call the inla() function and will be used to convert the inla object
#' to a bru object. Thereafter, preict.bru() is called to perform the prediction.
#'
#' @aliases predict.inla
#' @param object An object obtained by calling \link{bru} or \link{lgcp}.
#' @param ... Arguments passed on to predict.bru().
#' @return A \code{prediction} object.
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples 
#' \dontrun{
#' 
#' # Generate some data
#' 
#' input.df <- data.frame(x=cos(1:10))
#' input.df <- within(input.df, y <- 5 + 2*cos(1:10) + rnorm(10, mean=0, sd=0.1))
#'
#' # Fit a Gaussian likelihood model
#' 
#' formula = y ~ x
#' fit <- inla(formula, "gaussian", data = input.df, control.compute=list(config = TRUE))
#' 
#' # Add formula to inla object
#' 
#' fit$formula = y ~ x
#'
#' Estimate posterior statistics of exp(x), where x is the random effect.
#' 
#' xpost = predict(fit, NULL, ~ exp(x))
#' xpost
#' plot(xpost)
#' 
#' }

predict.inla = function(object, ...) {
  fit$sppa$model = make.model(object$formula)
  class(fit) = c("bru","inla")
  predict(fit, ...)
}



extract.summary = function(result, property) {
  ret = list()
  
  for ( label in rownames(result$summary.fixed) ) {
    ret[[label]] = result$summary.fixed[label, property]
  }
  
  for ( label in names(result$summary.random) ) {
    ret[[label]] = result$summary.random[[label]][, property]
  }
  
  for ( label in rownames(result$summary.hyperpar) ) {
    new.label = gsub(" ","_", x = label, fixed = TRUE)
    ret[[new.label]] = result$summary.hyperpar[label, property]
  }
  
  
  # For factors we add a data.frame with column names equivalent to the factor levels
  fac.names = names(effect(result$model))[unlist(lapply(result$model$effects, function(e) {e$model == "factor"}) )]
  for (name in fac.names) {
    tmp = unlist(ret[startsWith(names(ret), name)])
    names(tmp) = lapply(names(tmp), function(nm) {substring(nm, nchar(name)+1)})
    ret[[name]] = tmp
  }
  
  
  ret
}



##
## A wrapper for inla.posterior.sample()
##
## Converts each sample into a list of sub-samples representing the
## latent variables
##
## Example: samples[[1]]$somefield is the value of the field "somefield"
##          in the first sample
##

inla.posterior.sample.structured = function(result,n){

  # Workaround for older versions of INLA
  if ("hyper.user.scale" %in% formalArgs(inla.posterior.sample)) {
    samples = INLA::inla.posterior.sample(n, result)
  } else { samples = INLA::inla.posterior.sample(n, result, FALSE) }
  
  ssmpl = list()
  for (i in 1:length(samples)) {
    smpl.latent = samples[[i]]$latent
    smpl.hyperpar = samples[[i]]$hyperpar
    vals = list()

    # Extract simulated predictor and fixed effects
    for (name in unique(c("Predictor", result$names.fixed))) { 
      vals[[name]] = extract.entries(name,smpl.latent) 
    }
    
    # For fixed effects that were modeled via factors we attach an extra vector holding the samples
    fac.names = names(effect(result$model))[unlist(lapply(result$model$effects, function(e) {e$model == "factor"}) )]
    for (name in fac.names) {
      vals[[name]] = smpl.latent[startsWith(rownames(smpl.latent), name),]
      names(vals[[name]]) = lapply(names(vals[[name]]), function(nm) {substring(nm, nchar(name)+1)})
    }
    
      
    # Extract simulated latent variables. If the model is "clinear", however, extract the realisations
    # from the hyperpar field.
    if (length(result$summary.random) > 0) {
      for (k in 1:length(result$summary.random)) {
        name = unlist(names(result$summary.random[k]))
        model = result$model.random[k]
        if (!(model=="Constrained linear")) { vals[[name]] = extract.entries(name,smpl.latent) }
        else { vals[[name]] = smpl.hyperpar[paste0(paste0("Beta_intern for ",name)," -- in user scale")] }
      }
    }
    if(length(smpl.hyperpar)>0){
      names(smpl.hyperpar) = sapply(names(smpl.hyperpar), function(nm) gsub(" ","_", x = nm, fixed = TRUE))
    }
    ssmpl[[i]] = c(vals, smpl.hyperpar)
  }

#
# Return
#
return(ssmpl)
}

extract.entries = function(name,smpl){
  ename = gsub("\\.", "\\\\.", name)
  ename = gsub("\\(", "\\\\(", ename)
  ename = gsub("\\)", "\\\\)", ename)
  ptn = paste("^", ename, "[\\:]*[\\.]*[0-9]*$", sep="")
  return(smpl[grep(ptn,rownames(smpl))])
}

# Stack multiple observations
#
#
# @aliases inla.stack.y
# @export
# @param ... observation vectors
# @return y observation vector
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#

inla.stack.y = function(...) {
  all.y = lapply(list(...),function(x) return(x$data$data$y))
  nrow.y = lapply(list(...),function(x) return(nrow(x$data$data)))
  y = list()
  for (j in 1:nargs()){
    ny = nrow.y[[j]]
    y[[j]] = cbind( matrix(NA,nrow=ny,ncol=j-1) , all.y[[j]], matrix(NA,nrow=ny,ncol=nargs()-j) )
  }
  return(do.call(rbind,y))
}

# Stack multiple exposures
#
#
# @aliases inla.stack.e
# @export
# @param ... observation vectors
# @return e observation vector
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

inla.stack.e = function(...) {
  all.y = lapply(list(...),function(x) return(x$data$data$e))
  nrow.y = lapply(list(...),function(x) return(nrow(x$data$data)))
  y = list()
  for (j in 1:nargs()){
    ny = length(all.y[[j]])
    if (ny==0) { 
      e.tmp = rep(NA,nrow.y[[j]])
    } else { 
      e.tmp = all.y[[j]] 
    }
    y[[j]] = cbind( matrix(NA,nrow=nrow.y[[j]],ncol=j-1) , e.tmp , matrix(NA,nrow=nrow.y[[j]],ncol=nargs()-j) )
  }
  return(do.call(rbind,y))
}


# Join stacks intended to be run with different likelihoods
#
# @aliases inla.stack.mjoin
# 

inla.stack.mjoin = function(..., compress = TRUE, remove.unused = TRUE){
  y = inla.stack.y(...)
  e = inla.stack.e(...)
  mstack = inla.stack.join(...,compress = compress, remove.unused = remove.unused)
  mstack$data$y = y
  mstack$data$e = e
  return(mstack)
}


# Retrieve data from stack. Other than inla.stack.data this will give
# an observation vector y with multiple columns
#
# @aliases inla.stack.mdata
# 

inla.stack.mdata = function(stack){
  mdata = inla.stack.data(stack)
  if (!is.null(stack$data$y)) {
    mdata$y.inla = stack$data$y
  } else {
    mdata$y.inla = stack$data$data$y
  }
  
  return(mdata)
}

# Combine stacks by adding up predictors
#
# @aliases inla.stack.add 
# 

inla.stack.add = function(...) {
  stacks = list(...)
  stk3 = inla.stack.sum(stacks[[1]]$data$data, 
                        A = lapply(stacks,function(x) {return(inla.stack.A(x))}), 
                        effects = lapply(stacks,function(x) {return(x$effects$data)}))
}



plotmarginal.inla = function(result,varname="Intercept", link = function(x){x}, add = FALSE, ggp = TRUE, lwd=3,...){
  vars = variables.inla(result)
  ovarname = varname
  if (paste0("Beta for ", varname) %in% rownames(vars)) { varname = paste0("Beta for ", varname)}
  
  if (varname %in% c(result$names.fixed, rownames(result$summary.hyperpar) )) {
    
    if (vars[varname,"type"] == "fixed"){
      marg = inla.tmarginal(link, result$marginals.fixed[[varname]])
    }
    else if (vars[varname,"type"] == "random"){
      marg = inla.tmarginal(link, result$marginals.random[[varname]])
    }
    else if (vars[varname,"type"] == "hyperpar"){
      marg = inla.tmarginal(link, result$marginals.hyperpar[[varname]])
    }
    uq = inla.qmarginal(0.975, marg)
    uqy = inla.dmarginal(uq, marg)
    lq = inla.qmarginal(0.025, marg)
    lqy = inla.dmarginal(lq, marg)
    inner.x = seq(lq, uq, length.out = 100)
    inner.marg = data.frame(x = inner.x, y = inla.dmarginal(inner.x, marg))

    df = data.frame(marg)
    ggplot(data = df, aes_string(x="x",y="y")) + geom_path() + geom_ribbon(ymin = 0,aes_string(ymax = "y"), alpha = 0.1) +
      geom_segment(x = lq, y = 0, xend = lq, yend = lqy) +
      geom_segment(x = uq, y = 0, xend = uq, yend = uqy) +
      geom_ribbon(data = inner.marg, ymin = 0, aes_string(ymax = "y"), alpha = 0.1) +
      xlab(ovarname) + ylab("pdf")

  } else {
    
    df = result$summary.random[[varname]]
    colnames(df) = c("ID","mean","sd","lower","mid","upper","mode","kld")
    p <- ggplot(df, aes_string("ID", "mode"))
    p + geom_crossbar(aes_string(ymin = "lower", ymax = "upper")) + ylab("mod and quantiles") + xlab(paste0(varname," ID"))

  }
}

variables.inla = function(result, include.random=TRUE){
  handle.missing <- function(col.names) {
    cbind(data.frame(type = character(0),
                     model = character(0),
                     as.data.frame(
                         matrix(NA, 0, length(col.names),
                                dimnames=list(c(), col.names)))))
  }

  handle.missing.columns <- function(data, col.names) {
    missing.names <- setdiff(col.names, colnames(data))
    if (length(missing.names) > 0) {
      df <- as.data.frame(matrix(NA, nrow(data), length(missing.names),
                                 dimnames=list(NULL, missing.names)))
      return(cbind(data, df))
    } else {
      return(data)
    }
}

  handle.data.frame <- function(data, type, model, col.names) {
    ##    rownames(data) ???
    cbind(data.frame(type = rep(type, nrow(data)),
                     model = rep(model, nrow(data))),
          handle.missing.columns(data, col.names))
  }

  ## Get column names, handling possibly empty output:
  fixed.missing <- (is.null(result$summary.fixed) ||
                      (ncol(result$summary.fixed) == 0))
  hyperpar.missing <- (is.null(result$summary.hyperpar) ||
                         (ncol(result$summary.hyperpar) == 0))
  random.missing <- (is.null(result$summary.random) ||
                       (length(result$summary.random) == 0) ||
                         (ncol(result$summary.random[[1]]) == 0))
  col.names <- c()
  if (!fixed.missing) {
    col.names <- union(col.names, colnames(result$summary.fixed))
  }
  if (!hyperpar.missing) {
    col.names <- union(col.names, colnames(result$summary.hyperpar))
  }
  if (!random.missing) {
    col.names <- union(col.names, colnames(result$summary.random[[1]]))
  }
  if (fixed.missing) {
    fixed <- handle.missing(col.names)
  } else {
    fixed <- handle.data.frame(result$summary.fixed,
                               "fixed", "fixed", col.names)
  }
  if (hyperpar.missing) {
    hyperpar <- handle.missing(col.names)
  } else {
    hyperpar <- handle.data.frame(result$summary.hyperpar,
                                  "hyperpar", NA, col.names)
  }
  if (random.missing) {
    random <- list(handle.missing(col.names))
  } else {
    if (!include.random) {
      random <-
        lapply(seq_along(result$summary.random),
               function(x) {
                 dat <- data.frame(mean = mean(result$summary.random[[x]]$mean),
                                   sd = mean(result$summary.random[[x]]$sd))
                 output <- handle.data.frame(dat,
                                             "random",
                                             result$model.random[x],
                                             col.names)
                 rownames(output) <-
                   names(result$summary.random)[x]
               output
             })
    } else {
      random <-
        lapply(seq_along(result$summary.random),
               function(x) {
                 output <- handle.data.frame(result$summary.random[[x]],
                                             "random",
                                             result$model.random[x],
                                             col.names)
                 rownames(output) <-
                   paste(names(result$summary.random)[x],
                         seq_len(nrow(result$summary.random[[x]])))
                 output
               })
    }
  }

  variables <- do.call(rbind, c(list(fixed, hyperpar), random))
  return(variables)
}

