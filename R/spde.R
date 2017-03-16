#' @title Posteriors of SPDE hyper parameters and Matern correlation or covariance function.
#' 
#' @description
#' Plots the posterior distribution of the range, log(range), variance, or log(variance) 
#' parameter of a model's SPDE component. Can also plot Matern correlation or covariance function.
#' \code{inla.spde.result}. 
#'
#' @param result An object inheriting from \code{inla}.
#' @param name Name of the SPDE, see \code{names(result$summary.random)}
#' @param what One of "range", "log.range", "variance", "log.variance", "matern.correlation" or "matern.covariance".
#' @return A posterior
#'  
#' @export

spde.posterior = function(result, name, what = "range") {
  spdespec = result$sppa$model$effects[[name]]$inla.spde
  spderesult <- inla.spde.result(result, name, spdespec)
  
  if ( what == "matern.correlation" || what == "matern.covariance") {
    
    if ( what == "matern.correlation" ) { 
      corr = TRUE 
      ylab = "Matern Correlation"
    } else { 
      corr = FALSE
      ylab = "Matern Covariance"
      }

    if ( class(result$sppa$model$effects[[name]]$mesh) == "inla.mesh.1d") { d = 1 } else { d = 2 }
    
    materncov = function(dist, kappa,d,corr){
      inla.matern.cov(nu=1, kappa, x=dist,d=d, corr = corr)
    }
    
    kappaQ = inla.qmarginal(p=c(0.025,0.5,0.975), marginal=spderesult$marginals.kappa[[1]])
    
    xmax = exp(spderesult$summary.log.range.nominal[["0.975quant"]]) * 1.2
    x = seq(0,xmax,length=200)
    
    # median
    med = apply(matrix(data=x, ncol=1), 1, 
                              function(X)materncov(dist=X, kappa =kappaQ[2],d=d,corr=corr ))
    # lower band
    uq = apply(matrix(data=x, ncol=1), 1, 
                              function(X)materncov(dist=X, kappa =kappaQ[1],d=d,corr=corr))
    # upper band
    lq = apply(matrix(data=x, ncol=1), 1,
                              function(X)materncov(dist=X, kappa =kappaQ[3],d=d,corr=corr))

    df = data.frame(x = x, median = med, lq = lq, uq = uq)
    attr(df, "type") = "1d"
    attr(df, "misc") = list(dims = "x", predictor = c("distance", ylab))
    class(df) = list("prediction","data.frame")
    df
    
  } else {
  
    marg = switch(what,
                  range = spderesult$marginals.range.nominal[[1]],
                  log.range = spderesult$marginals.log.range.nominal[[1]],
                  variance = spderesult$marginals.variance.nominal[[1]],
                  log.variance = spderesult$marginals.log.variance.nominal[[1]]
                  )
    if(is.null(marg)) stop("Invalid varname: ",what,". must be one of 'range', 
                           'log.range',  'variance',  'log.variance', 
                           'matern.correlation', matern.covariance")
    
    med = inla.qmarginal(0.5, marg)
    uq = inla.qmarginal(0.975, marg)
    lq = inla.qmarginal(0.025, marg)
    inner.x = seq(lq, uq, length.out = 100)
    inner.marg = data.frame(x = inner.x, y = inla.dmarginal(inner.x, marg))
    colnames(inner.marg) = c(what, "pdf")
    df = data.frame(marg)
    colnames(df) = c(what, "pdf")
    attr(df, "type") = "0d"
    attr(df, "summary") = list(uq = uq, lq = lq, median = med, inner.marg = inner.marg)
    class(df) = list("prediction","data.frame")
    df
  }
}