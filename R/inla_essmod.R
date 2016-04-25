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

  samples = inla.posterior.sample(n,result,hyper.user.scale = TRUE)

  ssmpl = list()
  for (i in 1:length(samples)) {
    smpl.latent = samples[[i]]$latent
    smpl.hyperpar = samples[[i]]$hyperpar
    vals = list()

    # Extract simulated predictor and fixed effects
    for (name in c("Predictor",result$names.fixed)) { vals[[name]] = extract.entries(name,smpl.latent) }

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
    ssmpl[[i]] = c(vals, list(hyper=smpl.hyperpar))
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
  ptn = paste("^", ename, "\\.[0-9]*$", sep="")
  return(smpl[grep(ptn,rownames(smpl))])
}

#' Stack multiple observations
#'
#'
#' @aliases inla.stack.y
#' @export
#' @param ... observation vectors
#' @return y observation vector
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

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

#' Stack multiple exposures
#'
#'
#' @aliases inla.stack.e
#' @export
#' @param ... observation vectors
#' @return e observation vector
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

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


#' Join stacks intended to be run with different likelihoods
#'
#' @aliases inla.stack.mjoin
#' 

inla.stack.mjoin = function(..., compress = TRUE, remove.unused = TRUE){
  y = inla.stack.y(...)
  e = inla.stack.e(...)
  mstack = inla.stack.join(...,compress = compress, remove.unused = remove.unused)
  mstack$data$y = y
  mstack$data$e = e
  return(mstack)
}


#' Retrieve data from stack. Other than inla.stack.data this will give
#' an observation vector y with multiple columns
#'
#' @aliases inla.stack.mdata
#' 

inla.stack.mdata = function(stack){
  mdata = inla.stack.data(stack)
  mdata$y.inla = stack$data$y
  return(mdata)
}

#' Combine stacks by adding up predictors
#'
#' @aliases inla.stack.add 
#' 

inla.stack.add = function(...) {
  stacks = list(...)
  stk3 = inla.stack.sum(stacks[[1]]$data$data, 
                        A = lapply(stacks,function(x) {return(inla.stack.A(x))}), 
                        effects = lapply(stacks,function(x) {return(x$effects$data)}))
}


#' Shortcut to refine an inla.mesh object
#'
#'
#' @aliases inla.mesh.refine
#' @export
#' @param mesh an inla.mesh object
#' @param int.points integration points
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

mesh.refine = function(mesh,refine=list(max.edge=1)){
rmesh = inla.mesh.create(loc=mesh$loc,interior=inla.mesh.interior(mesh),boundary=inla.mesh.boundary(mesh),refine=refine)
return(rmesh)
}

#' Split triangles of a mesh into four triangles
#'
#' Warning: does not reconstruct interior boundary
#' Warning2: Works in euclidean coordinates. Not suitable for sphere.
#'
#' @aliases mesh.split
#' @export
#' @param mesh an inla.mesh object
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

mesh.split = function(mesh,n=1){

  p1 = mesh$loc[mesh$graph$tv[,1],]
  p2 = mesh$loc[mesh$graph$tv[,2],]
  p3 = mesh$loc[mesh$graph$tv[,3],]

  m1 = p1 + 0.5*(p2-p1)
  m2 = p1 + 0.5*(p3-p1)
  m3 = p2 + 0.5*(p3-p2)
  all.loc = rbind(mesh$loc,m1,m2,m3)

  bnd.mid = mesh$loc[mesh$segm$bnd$idx[,1],] + 0.5 * ( mesh$loc[mesh$segm$bnd$idx[,2],] - mesh$loc[mesh$segm$bnd$idx[,1],]  )
  all.bnd = rbind(mesh$segm$bnd$loc,bnd.mid)

  #   int.mid = mesh$loc[mesh$segm$int$idx[,1],] + 0.5 * ( mesh$loc[mesh$segm$int$idx[,2],] - mesh$loc[mesh$segm$int$idx[,1],]  )
  #   all.int = rbind(int.mid)
  #
  #   plot(mesh)
  #   points(mesh$loc[mesh$segm$int$idx[1,1],])
  #   points(mesh$loc[mesh$segm$int$idx[,2],])
  #
  #   int = rbind(mesh$loc[mesh$segm$int$idx[,1],],mesh$loc[mesh$segm$int$idx[,2],])
  #   points(int)
  #   points(all.int)

  mesh2 = inla.mesh.create(loc = all.loc, boundary = all.bnd )

  if (n == 1) { return(mesh2) }
  else { return(mesh.split(mesh2,n-1))}
}

# plot.marginal = function(...){UseMethod("plot.marginal")}
plot.marginal.inla = function(result,varname="Intercept", link = function(x){x}, add = FALSE, ggp = TRUE, lwd=3,...){
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
    if ( ggp ) {
      require( ggplot2 )
      df = data.frame(marg)
      ggplot(data = df, aes(x=x,y=y)) + geom_path() + geom_ribbon(ymin = 0,aes(ymax = y), alpha = 0.1) +
        geom_segment(x = lq, y = 0, xend = lq, yend = lqy) +
        geom_segment(x = uq, y = 0, xend = uq, yend = uqy) +
        geom_ribbon(data = inner.marg, ymin = 0, aes(ymax = y), alpha = 0.1) +
        xlab(ovarname) + ylab("pdf")
    } else {
      if (!add) {
        plot(marg,type='l',xlab=varname,ylab="Posterior density",...)
      } else {
        points(marg, type='l', xlab="",ylab="", ...)
      }
      lheight = max(marg[,"y"])
      lines(x=c(vars[varname,"mode"],vars[varname,"mode"]),y=c(0,lheight),col="blue",lwd=lwd,...)
      lines(x=c(vars[varname,"mean"],vars[varname,"mean"]),y=c(0,lheight),col="red",lwd=lwd,...)
      lines(x=c(vars[varname,"0.025quant"],vars[varname,"0.025quant"]),y=c(0,lheight),col=rgb(0,0.6,0),lwd=lwd,...)
      lines(x=c(vars[varname,"0.975quant"],vars[varname,"0.975quant"]),y=c(0,lheight),col=rgb(0,0.6,0),lwd=lwd,...)
    }
  } else {
    if ( require(ggplot2) ){
      df = result$summary.random[[varname]]
      colnames(df) = c("ID","mean","sd","lower","mid","upper","mode","kld")
      p <- ggplot(df, aes(ID, mode))
      p + geom_crossbar(aes(ymin = lower, ymax = upper)) + ylab("mod and quantiles") + xlab(paste0(varname," ID"))
    } else {
      uq = result$summary.random[[varname]][,"0.975quant"]
      lq = result$summary.random[[varname]][,"0.025quant"]
      md = result$summary.random[[varname]][,"mode"]
      plot(link(md), ylim = c(min(lq),max(uq)))
      points(link(uq), col = "red")
      points(link(lq), col = "red")
    }
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

