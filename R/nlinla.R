
nlinla.taylor = function(expr, data) {
  effects = colnames(data)
  grd = matrix(NA, ncol = length(effects), nrow = nrow(data))
  const = numeric()
  const2 = numeric()
  for (k in 1:nrow(data)) {
    myenv = new.env()
    invisible(lapply(effects, function(x) myenv[[x]] = data[k,x]))
    tmp = numericDeriv(expr[[1]], effects, myenv)
    grd[k,] = attr(tmp, "gradient")
    const[k] = tmp[1] - sum(grd[k,] * t(as.vector(data[k,effects])))
  }
  grd =  data.frame(grd)
  colnames(grd) = effects
  return(list(gradient = grd, const = const))
}

nlinla.epunkt = function(model, data, result = NULL) {
  if ( is.null(result) ){
    df = data.frame(matrix(0, nrow = nrow(data), ncol = length(model$effects)))
    colnames(df) = model$effects
    df
  } else {
    evaluate.model(model, result, data, link = identity, do.sum = FALSE)  
  }
}

nlinla.reweight = function(A, weights, model, data){
  expr = model$expr
  epkt = nlinla.epunkt(model, data, result = model$result)
  ae = nlinla.taylor(expr, epkt)
  for ( k in 1:length(A) ) {
    nm = names(A)[k]
    if (!(nm=="")){ A[[k]] = A[[k]] * ae$gradient[[nm]] }
  }
  for ( k in 1:length(weights) ){
    for (j in 1:ncol(weights[[k]])) {
      nm = names(weights[[k]])[j]
      if ( nm %in% names(ae$gradient) ) {
        weights[[k]][[nm]] = weights[[k]][[nm]] * ae$gradient[[nm]]
      }
    }
  }
  
  return(list(A = A, weights = weights, const = ae$const))
}

