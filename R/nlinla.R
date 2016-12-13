
nlinla.taylor = function(expr, epunkt, data) {
  if (nrow(epunkt)<=1000) {
    effects = colnames(epunkt)
    df = as.data.frame(data)
    df = df[,setdiff(names(df),names(epunkt)),drop=FALSE]
    wh = cbind(df, epunkt)
    myenv = new.env()
    invisible(lapply(colnames(wh), function(x) myenv[[x]] = wh[,x]))
    tmp = numericDeriv(expr[[1]], effects,  rho = myenv)
    gra = attr(tmp, "gradient")
    nr = nrow(gra)
    ngrd = matrix(NA,nrow=nr,ncol=length(effects))
    for (k in 1:length(effects)){ ngrd[,k] = diag(gra[,((k-1)*nr+1):(k*nr)]) }
    nconst = as.vector(tmp) - rowSums(ngrd * epunkt)
    ngrd =  data.frame(ngrd)
    colnames(ngrd) = effects
    return(list(gradient = ngrd, const = nconst))
  } else {
    blk = round(1:nrow(epunkt)/1000)
    qq = by(1:nrow(epunkt), blk, function(idx) {nlinla.taylor(expr,epunkt[idx,],data[idx,])})
    nconst = do.call(c, lapply(qq, function(x) { x$const }))
    ngrd = do.call(rbind, lapply(qq, function(x) { x$grad }))
    return(list(gradient = ngrd, const = nconst))
  }
}

nlinla.epunkt = function(model, data, result = NULL) {
  if ( is.null(result) ){
    df = data.frame(matrix(0, nrow = nrow(data.frame(data)), ncol = length(model$effects)))
    colnames(df) = elabels(model)
    df
  } else { evaluate.model(model, result, data, property = "mode") }
}

nlinla.reweight = function(A, model, data){
  expr = model$expr
  epkt = nlinla.epunkt(model, data, result = model$result)
  ae = nlinla.taylor(expr, epkt, data)
  for ( k in 1:length(A) ) {
    nm = names(A)[k]
    if ( !(is.null(nm) || nm == "")){ A[[k]] = A[[k]] * ae$gradient[[nm]] }
  }
  return(list(A = A, const = ae$const))
}

