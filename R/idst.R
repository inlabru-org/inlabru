#' iDistance wrapper for INLA
#' 
#'
#' @aliases idst
#' @export
#' @param data A \link{dsdata} object
#' @param model A \link{model} object
#' @param ips A data.frame of integration points
#' @param stack A stack configuration (Not implemented yet, use NULL)
#' @param n Number of \link{inla} iterations
#' @param idst.verbose If TRUE, be verbose (use verbose=TRUE) to make INLA verbose
#' @param ... Arguments passed on to \link{inla}
#' @return A \link{inla} object


idst = function(data, model, ips = NULL, stack = NULL, predict = NULL, n = 1, idst.verbose = FALSE, ...){
  
  model = join(model)
  update.model(model, result = NULL) # This will initialize the history
  
  if ( is.null(ips) ) {
    if ( "ips" %in% names(data) ) { ips = data$ips }
    else { stop("Parameter 'ips' not provided and data set data has no field 'ips'" ) }
  }
  
  if ( is.null(stack) ) {
    det.stack <- detection.stack(data, model = model)
    int.stack <- integration.stack(data, scheme = ips, model = model)
    
    if ( !is.null(predict) ) {
      if ("model" %in% names(predict)) { predict = list(predict) } # Old syntax
      pstacks = list()
      for (pr in predict) {
        if ( !("data" %in% names(pr)) ) { pr$data = data }
        if ( !("loc" %in% names(pr)) ) { pr$loc = data$mesh$loc[,c(1,2)] ; colnames(pr$loc) = data$mesh.coords }
        pstacks = c(pstacks, list(do.call(prediction.stack, pr)))
      }
      stk <- do.call(inla.stack, c(list(det.stack, int.stack), pstacks))
    } else { 
      stk <- inla.stack(det.stack, int.stack)
    }
    
  }

  for ( k in 1:n ) {
    
    result <- tryCatch( inla(formula = model$formula, 
                     family = "poisson",
                     data = c(inla.stack.data(stk), list.data(model)),
                     control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                     E = inla.stack.data(stk)$e, ...), 
                   error = function(e) { 
                     if (k == 1) { stop(e) }
                     else { stop(paste0("INLA crashed during iteration ",k,". It is likely that there is a convergence problem.")) }
                     }
                   )
    if ( idst.verbose ) { cat(".") }
    # Update model
    update.model(model, result)
    
    if ( k > 1) {
      # Update stacks
      det.stack <- detection.stack(data, model = model)
      int.stack <- integration.stack(data, scheme = ips, model = model)
      

      if ( !is.null(predict) ) {
        if ("model" %in% names(predict)) { predict = list(predict) } # Old syntax
        pstacks = list()
        for (pr in predict) {
          if ( !("data" %in% names(pr)) ) { pr$data = data }
          if ( !("loc" %in% names(pr)) ) { pr$loc = data$mesh$loc[,c(1,2)] ; colnames(pr$loc) = data$mesh.coords }
          pstacks = c(pstacks, list(do.call(prediction.stack, pr)))
        }
        stk <- do.call(inla.stack, c(list(det.stack, int.stack), pstacks))
      } else { 
        stk <- inla.stack(det.stack, int.stack)
      }
      
    }
  }
  result$stack = stk
  return(result)
}
