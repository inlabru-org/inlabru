
nlinla.taylor = function(expr, epunkt, data, env) {
  
  if ( "offset" %in% names(epunkt) ) {
    stop("One of your model components is an offset. 
          However, you are using a non-linear predictor (formula),
          which would set the respective term to zero. 
          Please remove the offset component and add its value
          to the predictor formula.")
  }
  
  if (nrow(epunkt)<=1000) {
    effects = colnames(epunkt)
    df = as.data.frame(data)
    df = df[,setdiff(names(df),names(epunkt)),drop=FALSE]
    wh = cbind(df, epunkt)
    myenv = new.env()
    invisible(lapply(colnames(wh), function(x) myenv[[x]] = wh[,x]))
    invisible(lapply(setdiff(names(env), names(myenv)), function(x) myenv[[x]] = env[[x]]))
    tmp = numericDeriv(expr[[1]], effects,  rho = myenv)
    gra = attr(tmp, "gradient")
    nr = nrow(gra)
    ngrd = matrix(NA,nrow=nr,ncol=length(effects))
    for (k in 1:length(effects)){
      # as.matrix required since diag() will not work if gra has single entry
      ngrd[,k] = diag(as.matrix(gra[, ((k-1)*nr+1):(k*nr), drop=FALSE]))
    }
    nconst = as.vector(tmp) - rowSums(ngrd * epunkt)
    ngrd =  data.frame(ngrd)
    colnames(ngrd) = effects
    return(list(gradient = ngrd, const = nconst))
  } else {
    blk = round(1:nrow(epunkt)/1000)
    qq = by(1:nrow(epunkt), blk, function(idx) {nlinla.taylor(expr,epunkt[idx,],data[idx,], env)})
    nconst = do.call(c, lapply(qq, function(x) { x$const }))
    ngrd = do.call(rbind, lapply(qq, function(x) { x$grad }))
    return(list(gradient = ngrd, const = nconst))
  }
}

nlinla.epunkt = function(model, data, result = NULL) {
  # This function determines the current point around which
  # to perform the taylo approximation
  # (1) If result is NULL set all all effects to 0
  # (2) If result is a data.frame, use the entries as to where to approximate
  # (3) if result is an inla object, use these estimates as to where to approximate
  
  dfdata = data.frame(data) # data as data.frame (may have been supplied as Spatial* object)
  if ( is.null(result) ){
    df = data.frame(matrix(0, nrow = nrow(dfdata), ncol = length(model$effects)))
    colnames(df) = elabels(model)
    df
  } else if ( !inherits(result, "inla") & is.data.frame(result) ) {
    # If result contains only a single row data frame repeat it to match the data
    if ( (nrow(result) == 1) & (nrow(dfdata)>1) ) {
        result = result[rep(1,nrow(dfdata)),,drop = FALSE]
    }
    # Check if all variables have been supplied. Those that aren't are set to 0
    for (eff in setdiff(names(model$effects), names(result))) {
      result[[eff]] = 0
    }
    return(result)
  } else { 
    evaluate.model(model, result, data, property = "mean") 
  }
}

nlinla.reweight = function(A, model, data, expr, result){
  epkt = nlinla.epunkt(model, data, result = result)
  ae = nlinla.taylor(expr, epkt, data, environment(model$formula))
  for ( k in 1:length(A) ) {
    nm = names(A)[k]
    if ( !(is.null(nm) || nm == "")){ A[[k]] = A[[k]] * ae$gradient[[nm]] }
  }
  return(list(A = A, const = ae$const))
}

