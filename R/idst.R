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
#' @param ... Arguments passed on to \link{inla}
#' @return A \link{inla} object


idst = function(data, model, ips = NULL, stack = NULL, n = 1, ...){
  
  model = join(model)
  
  if ( is.null(ips) ) {
    if ( "ips" %in% names(data) ) { ips = data$ips }
    else { stop("Parameter 'ips' not provided and data set data has no field 'ips'" ) }
  }
  
  if ( is.null(stack) ) {
    det.stack <- detection.stack(data, model = model)
    int.stack <- integration.stack(data, scheme = ips, model = model)
    stk <- inla.stack(det.stack, int.stack)  
  }

  for ( k in 1:n ) {
    result <- inla(formula = model$formula, 
                   family = "poisson",
                   data = c(inla.stack.data(stk), list.data(model)),
                   control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                   E = inla.stack.data(stk)$e,
                   ...)
    if ( k > 1) {
      update.model(model, result)
      det.stack <- detection.stack(data, model = model)
      int.stack <- integration.stack(data, scheme = ips, model = model)
      stk <- inla.stack(det.stack, int.stack)  
    }
  }
  return(result)
}
