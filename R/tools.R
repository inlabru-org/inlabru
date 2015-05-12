recurse.rbind = function(fun,dframe,cols){
  if (length(cols) == 1) {
    cols = cols[[1]]
    tmp = list()
    k = 1
    for (cl in unique(dframe[,cols])){
      tmp[[k]] = fun(dframe[dframe[,cols]==cl,])
      k = k+1
    }
    ret = do.call(rbind,tmp)
  } else {
    col = cols[[1]]
    tmp = list()
    k = 1
    for (cl in unique(dframe[,col])){
      tmp[[k]] = recurse.rbind(fun,dframe[dframe[,col]==cl,],cols[2:length(cols)])
      k = k+1
    }
    ret = do.call(rbind,tmp)
  }
  return(ret)
}


