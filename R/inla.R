#' Prediction from fitted inla model
#' 
#' Takes a fitted inla object produced by the function inla() and produces predictions given a 
#' new set of values for the model covariates or the original values used for the model fit. 
#' The predictions can be based on any R expression that is valid given these values/covariates 
#' and the posterior of the estimated effects.
#' 
#' @aliases predict.inla
#' @export
#' @param object A \code{bru} object obtained by calling \link{bru} or \link{lgcp}.
#' @param ... Arguments passed on to \link{predict.bru}.
#' @return A \code{prediction} object.
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples
#' \donttest{
#' # Some features use the INLA package.
#' if (require("INLA", quietly = TRUE)) {
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
#' # Estimate posterior statistics of exp(x), where x is the fixed effect.
#' 
#' xpost = predict(fit, NULL, ~ exp(x))
#' xpost
#' plot(xpost)
#' 
#' }
#' }

predict.inla <- function(object, ...) {
  object$sppa$model <- make.model(object$.args$formula)
  class(object) <- c("bru", class(object))
  predict(object, ...)
}



#' Sampling based on bru posteriors
#' 
#' @description 
#' Takes a fitted \code{inla} object produced by \code{INLA::inla()} and produces samples given 
#' a new set of values for the model covariates or the original values used for the model fit. The
#' samples can be based on any R expression that is valid given these values/covariates and the joint
#' posterior of the estimated random effects.
#'  
#' @aliases generate.inla
#' @export
#' @family sample generators
#' @param object An \code{inla} object obtained by calling \code{INLA::inla()}.
#' @param ... additional arguments passed on to\link{generate.bru}.
#' 
#' @return List of generated samples
#' @seealso \link{predict.inla}
#' 
#' @author Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @examples
#' \donttest{
#' # Some features use the INLA package.
#' if (require("INLA", quietly = TRUE)) {
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
#' # Generate samples from the posterior distribution of exp(x), where x is the fixed effect.
#' 
#' xpost = generate(fit, NULL, ~ exp(x), n.samples = 2)
#' xpost
#' plot(xpost[[1]])
#' 
#' }
#' }

generate.inla = function(object,
                         ...)
{
  object$sppa$model <- make.model(object$.args$formula)
  class(object) <- c("bru", class(object))
  generate(object, ...)
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
  if ("hyper.user.scale" %in% formalArgs(INLA::inla.posterior.sample)) {
    samples = INLA::inla.posterior.sample(n = n, result = result)
  } else {
    samples = INLA::inla.posterior.sample(n = n,
                                          result = result,
                                          intern = FALSE)
  }
  
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
      names(smpl.hyperpar) <- vapply(names(smpl.hyperpar),
                                     function(nm) gsub(" ","_", x = nm, fixed = TRUE),
                                     "name")
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

# Expand observation vectors/matrices in stacks into to a multicolumn matrix for multiple likelihoods
#
# @aliases inla.stack.mexpand
# @name inla.stack.mexpand
# @export
# @param ... List of stacks that contain vector observations
#            (existing multilikelihood observation matrices are also permitted)
# @param old.names A vector of strings with the names of the observation vector/matrix for each stack.
#        If a single string, this is assumed for all the stacks. (default "BRU.response")
# @param new.name The name to be used for the expanded observation matrix,
#        possibly the same as an old name. (default "BRU.response")
# @return a list of modified stacks with multicolumn observations
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}> and Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#

inla.stack.mexpand <- function(...,
                               old.names = "BRU.response",
                               new.name = "BRU.response") {
  stacks <- list(...)
  if (length(old.names) == 1) {
    old.names <- rep(old.names, length(stacks))
  }
  y.cols <- unlist(lapply(seq_along(stacks),
                          function(x, stacks, old.names) {
                            LHS <- INLA::inla.stack.LHS(stacks[[x]])[[old.names[x]]]
                            ifelse(is.vector(LHS), 1, ncol(LHS))
                          },
                          stacks = stacks, old.names = old.names))
  y.offset <- c(0, cumsum(y.cols))
  y.cols.total <- sum(y.cols)
  for (j in 1:length(stacks)){
    LHS <- INLA::inla.stack.LHS(stacks[[j]])
    RHS <- INLA::inla.stack.RHS(stacks[[j]])
    A <- INLA::inla.stack.A(stacks[[j]])
    # Access the raw tag indexing information
    tags <- list(data = stacks[[j]]$data$index,
                 effects = stacks[[j]]$effects$index)
    
    # Expand the observation vector/matrix into a multilikelihood observation matrix:
    y.rows <- ifelse(is.vector(LHS[[old.names[j]]]),
                     length(LHS[[old.names[j]]]),
                     nrow(LHS[[old.names[j]]]))
    LHS[[new.name]] <-
      cbind( matrix(NA, nrow = y.rows, ncol=y.offset[j]),
             LHS[[old.names[j]]],
             matrix(NA, nrow = y.rows, ncol=y.cols.total - y.offset[j+1]) )
    
    # Create the modified stack, with model compression disabled to prevent modifications:
    stacks[[j]] <-
      INLA::inla.stack.sum(data = LHS, A = A, effects = RHS, compress = FALSE, remove.unused = FALSE)
    # Since the row indexing is unchanged, copy the tag index information:
    stacks[[j]]$data$index <- tags$data
    stacks[[j]]$effects$index <- tags$effects
  }
  stacks
}

# Stack multiple exposures
# Obsolete. Do not use. Exposure vector stacking is handled automatically by inla.stack
#
# @aliases inla.stack.e
# @param ... observation vectors
# @return e observation vector
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

inla.stack.e = function(...) {
  stop("This function is obsolete and should not be used.")
}


# Join stacks intended to be run with different likelihoods
#
# @param ... List of stacks that contain vector observations
#            (existing multilikelihood observation matrices are also permitted)
# @param old.names A vector of strings with the names of the observation vector/matrix for each stack.
#        If a single string, this is assumed for all the stacks. (default "BRU.response")
# @param new.name The name to be used for the expanded observation matrix,
#        possibly the same as an old name. (default "BRU.response")
# @aliases inla.stack.mjoin
# 

inla.stack.mjoin = function(..., compress = TRUE, remove.unused = TRUE,
                            old.names = "BRU.response", new.name = "BRU.response"){
  stacks = inla.stack.mexpand(..., old.names = old.names, new.name = new.name)
  do.call(INLA::inla.stack.join,
          c(stacks, list(compress = compress, remove.unused = remove.unused)))
}


# Retrieve data from stack.
# Obsolete. Do not use. Exposure vector stacking is handled automatically by inla.stack
#
# The special "BRU.response" object should have been constructed before, via inla.stack.mjoin.
# Since all the work to setup the stack properly was already done by inla.stack.mjoin,
# this function just passes its arguments through to INLA::inla.stack.data
#
# @aliases inla.stack.mdata
# @param stack The stack to extract data from
# @param ... Additional named objects to add to the list of data objects.
#            Typically used for spde model objects.
# @return A list of data and effect vectors/matrices, and optional extra objects
# 

inla.stack.mdata = function(stack, ...){
  warning("inla.stack.mdata: This function is obsolete and should not be used.")
  INLA::inla.stack.data(stack, ...)
}

# Combine stacks by adding up predictors "horizontally".
# Only the data section from the first stack is preserved.
# Only the tag information from the data section first stack is preserved.
# TODO: The effects tag information is not computed.
# NOTE: This function appears to be unused.
#
# @aliases inla.stack.add 
# 

inla.stack.add = function(..., compress = TRUE, remove.unused = TRUE) {
  stacks <- list(...)
  stack <- INLA::inla.stack.sum(data = INLA::inla.stack.LHS(stacks[[1]]), 
                                A = lapply(stacks, function(x) { INLA::inla.stack.A(x) }), 
                                effects = lapply(stacks, function(x) { INLA::inla.stack.RHS(x) }),
                                compress = compress, remove.unused = remove.unused)
  stack$data$index <- stacks[[1]]$data$index
  stack$effects$index <- NULL
  stack
}



plotmarginal.inla = function(result,varname="Intercept", link = function(x){x}, add = FALSE, ggp = TRUE, lwd=3,...){
  vars = variables.inla(result)
  ovarname = varname
  if (paste0("Beta for ", varname) %in% rownames(vars)) { varname = paste0("Beta for ", varname)}
  
  if (varname %in% c(result$names.fixed, rownames(result$summary.hyperpar) )) {
    
    if (vars[varname,"type"] == "fixed"){
      marg = INLA::inla.tmarginal(link, result$marginals.fixed[[varname]])
    }
    else if (vars[varname,"type"] == "random"){
      marg = INLA::inla.tmarginal(link, result$marginals.random[[varname]])
    }
    else if (vars[varname,"type"] == "hyperpar"){
      marg = INLA::inla.tmarginal(link, result$marginals.hyperpar[[varname]])
    }
    uq = INLA::inla.qmarginal(0.975, marg)
    uqy = INLA::inla.dmarginal(uq, marg)
    lq = INLA::inla.qmarginal(0.025, marg)
    lqy = INLA::inla.dmarginal(lq, marg)
    inner.x = seq(lq, uq, length.out = 100)
    inner.marg = data.frame(x = inner.x, y = INLA::inla.dmarginal(inner.x, marg))

    df = data.frame(marg)
    ggplot(data = df, aes_string(x="x",y="y")) + geom_path() +
      geom_ribbon(ymin = 0,aes_string(ymax = "y"), alpha = 0.1) +
      geom_segment(x = lq, y = 0, xend = lq, yend = lqy) +
      geom_segment(x = uq, y = 0, xend = uq, yend = uqy) +
      geom_ribbon(data = inner.marg, ymin = 0, aes_string(ymax = "y"), alpha = 0.1) +
      xlab(ovarname) + ylab("pdf")

  } else {
    
    df = result$summary.random[[varname]]
    colnames(df) = c("ID","mean","sd","lower","mid","upper","mode","kld")
    p <- ggplot(df, aes_string("ID", "mode"))
    p + geom_crossbar(aes_string(ymin = "lower", ymax = "upper")) +
      ylab("mod and quantiles") + xlab(paste0(varname," ID"))

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

