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


idst = function(data, model, ips = NULL, stack = NULL, predict = NULL, n = 1, ...){
  
  model = join(model)
  
  if ( is.null(ips) ) {
    if ( "ips" %in% names(data) ) { ips = data$ips }
    else { stop("Parameter 'ips' not provided and data set data has no field 'ips'" ) }
  }
  
  if ( is.null(stack) ) {
    det.stack <- detection.stack(data, model = model)
    int.stack <- integration.stack(data, scheme = ips, model = model)
    
    if ( !is.null(predict) ) {
      if ( is.model(predict) ) {
        loc = data$mesh$loc[,c(1,2)] ; colnames(loc) = data$mesh.coords
        predict = list(data = data, model = predict, loc = loc, tag = "prediction")
      }
      if ( !("data" %in% names(predict)) ) { predict$data = data }
      pred.stack = do.call(prediction.stack,predict)
      stk <- inla.stack(det.stack, int.stack, pred.stack)
    } else { 
      stk <- inla.stack(det.stack, int.stack)
    }
    
  }

  for ( k in 1:n ) {
    result <- inla(formula = model$formula, 
                   family = "poisson",
                   data = c(inla.stack.data(stk), list.data(model)),
                   control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                   E = inla.stack.data(stk)$e,
                   ...)
    # Update model
    update.model(model, result)
    
    if ( k > 1) {
      # Update stacks
      det.stack <- detection.stack(data, model = model)
      int.stack <- integration.stack(data, scheme = ips, model = model)
      
      if ( !is.null(predict) ) {
        if ( is.model(predict) ) {
          loc = data$mesh$loc[,c(1,2)] ; colnames(loc) = data$mesh.coords
          predict = list(data = data, model = predict, loc = loc, tag = "prediction")
        }
        if ( !("data" %in% names(predict)) ) { predict$data = data }
        pred.stack = do.call(prediction.stack,predict)
        stk <- inla.stack(det.stack, int.stack, pred.stack)
      } else { 
        stk <- inla.stack(det.stack, int.stack)
      }
    }
  }
  result$stack = stk
  return(result)
}
