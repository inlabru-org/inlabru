#' iDistance models for INLA
#' 
#' Facilitates the spatial modeling approaches using INLA.
#' A \code{model} has a formula that describes one or more of INLA's \code{f} objects.
#' It describes how given a data frame of locations or other coordinates translates into
#' evaluating the predictors of the formula. For manually setting up a \code{model} see the
#' constructor \link{make.model}. Useful operators on \code{model} objects are:
#' \itemize{
#'  \item{\link{make.model}: }{Create a model}
#'  \item{\link{summary}: }{Summarize a model}
#'  \item{\link{evaluate}: }{ Determine the linear predictor's value at some location }
#'  \item{\link{list.data}: }{ List environmental variables needed to call INLA}
#' }
#' 
#' 
#' @name model
NULL


#####################################
# GENERICS
#####################################

evaluate = function(...){UseMethod("evaluate")}
list.covariates = function(...){UseMethod("list.covariates")}
list.A = function(...){UseMethod("list.A")}
list.indices = function(...){UseMethod("list.indices")}
list.data = function(...){UseMethod("list.data")}


#####################################
# Constructor
#####################################

#' Create an inlabru \link{model}
#' 
#'
#' @export
#' @param fml A formula
#' @param points a data.frame of data points
#' @return A \link{model} object
make.model = function(fml, data) {
  submodel = list()
  covariates = list()
  
  # Define function for shifting the effect name into the g-call
  shift.names = function(fml) {
    tms = terms(fml)
    labels = attr(tms, "term.labels")
    effects = character()
    for (k in 1:length(labels)){
      lb = labels[[k]]
      gpd = getParseData(parse(text=lb))
      # Determine the name of the effect
      ename = getParseText(gpd, id=1)
      is.fixed = (gpd[1,"token"] == "SYMBOL")
      if ( !(ename %in% c("g","f","offset")) & !is.fixed ) {
        effects[[k]] = ename
        # Replace effect name by g(ename
        labels[[k]] = gsub(paste0(ename,"("), paste0("g(",ename,", "), lb, fixed = TRUE)
      }
    }
    labels
  }
  
  lbl = shift.names(fml)
  
  for ( k in 1:length(lbl) ) {
    
    lb = lbl[[k]]
    
    #
    # We expect to see three types of terms
    # g
    # f 
    # offset
    # fixed
    
    if (substr(lb,1,2) == "g(") { ########################## g() term
      
      # Extract mesh and spde model
      ge = eval(parse(text = lb), envir = environment(fml))
      
      # Replace function name by INLA f function
      lb = gsub("g\\(","f(",lb)
      
      # Remove extra mesh argument
      lb = gsub("[,][ ]*mesh[ ]*=[^),]*", "", lb)
      
      # Remove extra map argument
      pat = paste0("", deparse(ge$map))
      lb = gsub(deparse(ge$map), "", lb, fixed =TRUE)
      lb = gsub("[,][ ]*map[ ]*=[^),]*", "", lb)
      
      # Remove extra A.msk argument
      pat = paste0("", deparse(ge$A.msk))
      lb = gsub(deparse(ge$map), "", lb, fixed =TRUE)
      lb = gsub("[,][ ]*A.msk[ ]*=[^),]*", "", lb)
      
      lbl[[k]] = lb
      
      submodel[[ge$f$label]] = c(ge$f, list(mesh = ge$mesh, 
                                            A.msk = ge$A.msk,
                                            mesh.coords = ge$f$term,
                                            map = ge$map,
                                            inla.spde = ge$model))
      
    } else if (substr(lb,1,2) == "f(") { ##################### f() term
      
    } else if (substr(lb,1,7) == "offset(") { ################ offset() term
      
    } else { ################################################# fixed effect term
      
      gpd = getParseData(parse(text=lb))
      
      if (gpd[1,"token"] == "SYMBOL") {
        if (!(gpd[1,"text"] %in% c(names(environment(fml)), names(data)))) { environment(fml)[[gpd[1,"text"]]] = 0 }
        if ( gpd[1,"text"] %in% names(data) ) {
          
        } else {
          covariates[[lbl[k]]] = function(x) {
            v = rep(1, nrow(as.data.frame(x)))
            ret = data.frame(v)
            colnames(ret) = effect
            return(ret)
          }
          environment(covariates[[lbl[k]]]) = new.env()
          assign("effect", lbl[k], envir = environment(covariates[[lbl[k]]]))
        }
        
      } else if (gpd[1,"token"] == "expr") {
        old.label = lb
        lb = paste0(gpd[2,"text"],".effect")
        effects[k] = paste0(gpd[2,"text"],".effect")
        covariates[[lb]] = function(...) {do.call(function(...) {eval(parse(text=old.label), envir = list(...))}, as.list(...))}
      }
      
      
      # New label
      lbl[[k]] = paste0("f(",lb,", model='linear')")
      ge = eval(parse(text = lbl[[k]]))
      
      submodel[[lb]] = c(eval(parse(text = lbl[[k]])), 
                         list(mesh = ge$mesh, 
                              A.msk = ge$A.msk,
                              mesh.coords = ge$f$term,
                              map = as.name(lb)))
      
    }
  }
  
  # Did the old formula have an intercept?
  if ( attr(terms(fml), "intercept") == 0 ) {icpt = "-1"} else { icpt = ""}
  
  # Combine labels into new formula
  new.rhs = as.formula(paste0(". ~ ", paste0(lbl, collapse = " + "), icpt))
  new.fml = update.formula(fml, new.rhs)
  environment(new.fml) = environment(fml)
  
  # Return
  ret = list(effects = submodel, formula = new.fml, in.formula = fml)
  class(ret) = c("model","list")
  ret
}


#####################################
# OPERATORS ON MODELS
#####################################

labels = function(mdl) { names(mdl$effects) }
predictor = function(mdl) { parse(text=paste0(labels(m),collapse="+"))}
effect = function(mdl, name = NULL) {
  if (is.null(name)) { 
    mdl$effects } 
  else {
    mdl$effects[[name]]
  }
   
}


#' Model summary
#'
#' @aliases summary.model
#' @export
#' @param mdl An \link{inlabru} \link{model}
#' 
summary.model = function(mdl) {
  for (label in labels(mdl)) {
    eff = effect(mdl,label)
    cat(paste0(label, ":\n"))
    cat(paste0("\t model: ",eff$model, "\n"))
    cat(paste0("\t map: ",eff$map,", class = ", class(eff$map),"\n"))
    cat(paste0("\t covariate function: ",eff$covariate, "\n"))
    cat("\n")
    
  }
  cat(paste0("--- FORMULA ---\n\n"))
  print(mdl$formula)
}



#' List data needed to run INLA
#'
#' @aliases list.data.model
#' @export
#' 

list.data.model = function(model){
  # THIS IS A WORKAROUND. If a formula is constructed and contains the stack of 
  # a previous run INLA will confuse the data from the new stack with the old stack
  # due to the data = ... as.list(environment(jmdl$formula))) ... parameter.
  # INLA reports the following error: 
  # Error in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 2379, 2739
  
  assign("stack", NULL, envir = environment(model$formula))
  
  # Formula environment as list
  elist = as.list(environment(model$formula))
  
  # Remove previous inla results. For some reason these slow down the next INLA call.
  elist = elist[unlist(lapply(elist, function(x) !inherits(x, "inla")))]
}

#' List of covariates effects needed to run INLA
#'
#' @aliases list.covariates.model
#' 

list.covariates.model = function(mdl, pts){
    
  all.varnames = c(all.vars(mdl$formula), all.vars(mdl$expr)) 
  covar.data = list()
    
    # Points may be annotated with covariates
    
    for (vname in all.varnames){
      if (vname %in% names(data.frame(pts))) {
        covar.data[[vname]] = data.frame(pts)[,vname]
      }
    }
    covar.data = data.frame(do.call(cbind,covar.data))
}

#' List of A matrices needed to run INLA
#'
#' @aliases list.A.model
#' 

list.A.model = function(mdl, points){
  A.lst = list()

  # For each effect create indices
  for ( eff in effect(mdl) ) {
    if ( is.null(eff$n) ) {
      A.lst[[eff$label]] = Matrix::Diagonal(nrow(data.frame(points)))
    } else if ( is.null(eff$mesh) ) {
      idx = mapper(eff$map, points)
      A = matrix(0, nrow = nrow(data.frame(points)), ncol=eff$n)
      A[cbind(1:nrow(A), idx)] = 1
      A.lst[[eff$label]] = A
    } else {
      
      loc = mapper(eff$map, points)
      loc = as.matrix(loc) # inla.spde.make.A requires matrix format as input

      if (!("n.group" %in% names(eff$inla.spde))) { ng = 1 } else { ng = eff$inla.spde$n.group }
      if (ng > 1) {
        group = as.matrix(points[,eff$time.coords])
      } else { group = NULL }
      A = inla.spde.make.A(eff$mesh, loc = loc, group = group)
      # Mask columns of A
      if (!is.null(eff$A.msk)) { A = A[, as.logical(eff$A.msk), drop=FALSE]}
      # Weights for models with A-matrix are realized in the follwoing way:
      if (!is.null(eff$weights)) {
        w = eff$weights
        A = as.matrix(A)*as.vector(w)
      }
      A.lst[[eff$label]] = A

    }
  }
  
  
  
  
  
  
  return(A.lst)
}


#' List of spde indexing effects needed to run INLA
#'
#' @aliases list.indices.model
#' 

list.indices.model = function(mdl, points){
  
  # The list of indices to be returned
  idx = list()

  # For each effect create indices
  for ( eff in effect(mdl) ) {
    name = eff$label
    
    # Only create indices if the number n of effect variables is known
    if (!is.null(eff$n)){
      
      # Effects with meshes need a special treatment ...
      if ( is.null(eff$mesh) ) { 
        idx[[name]] = seq_len(eff$n) 
      } else {
        if ( "m" %in% names(eff$mesh) ) {
          idx[[name]] = 1:eff$mesh$m # support inla.mesh.1d models
          # If a is masked, correct number of indices
          if (!is.null(eff$A.msk)) { lst[[name]] = 1:sum(eff$A.msk)}
        } else {
          ng = eff$inla.spde[[name]]$n.group
          if (is.null(ng)) { ng = 1 }
          idx[[name]] = inla.spde.make.index(name, n.spde = eff$inla.spde$n.spde, n.group = ng)
        }
      }
    } else {
      idx[[name]] = mapper(eff$map, points)
    }
    
  }
  idx
}

#' Evaluate model at given locations
#' 
#' Compute an approximation to the linear predictor at given locations and gicen coordinates.
#'
#' @aliases evaluate.model evaluate
#' @export
#' @param model An iDistance \link{model}
#' @param result The result of an \link{inla} run or a sample obtained from \link{inla.posterior.sample.structured}
#' @param loc Locations and covariates needed to evaluate the model. If \code{NULL}, SPDE models will be evaluated at the mesh coordinates.
#' @param property Property of the model compnents to obtain value from. Default: "mode". Other options are "mean", "0.025quant", "0.975quant" and "sd".
#' @param use.covariate DEPRECATED (will be irgnored)

evaluate.model = function(model, result, loc, property = "mode", do.sum = TRUE, link = identity, n = 1, predictor = model$expr, use.covariate = NULL) {
  
  # Obtain A-matrices
  A = list.A.model(model, loc)
  
  # Do we otain our values from sampling or from a property of a summary?
  if ( property == "sample") {
    if ( inherits(result, "inla") ) {
      smp = inla.posterior.sample.structured(result, n = n) 
    } else {
      smp = result
      n = length(smp)
    }
  } 
 
  posts = list()
  
  for ( eff in effect(model) ){
    name = eff$label
    
    # Values from samples
    if ( property == "sample") {

      post = lapply(smp, function(s) {
        mult = as.vector(s[[name]])
        if (length(mult) == 1) { as.vector(A[[name]] %*% rep(mult, nrow(A[[name]]))) }
        else { rowSums(t(t(as.matrix(A[[name]])) * mult)) }
      })
    
    # Values from inla result
    } else {
      
      # Fixed effects
      if (name %in% rownames(result$summary.fixed)) {
        post = result$summary.fixed[name, property]
        post = rep(post, nrow(data.frame(loc)))
        
      # Hyper parameter
      } else if (paste0("Beta for ",name) %in% rownames(result$summary.hyperpar)) {
        
        post = result$summary.hyperpar[paste0("Beta for ",name), property]
        post = rep(post, nrow(data.frame(loc)))
        
        # Random effects
      } else if (name %in% names(result$summary.random)) {
        
        post = A[[name]] %*% result$summary.random[[name]][,property]
        post = as.vector(post)
        
      } else {
        stop(sprintf("Effect with name %s not found in inla result.", name))
      }
      
    }
    posts[[name]] = post
  }
  
  if ( property == "sample") {
    if (!is.null(predictor) && do.sum) { 
      ret = do.call(Map, c(list(function(...){apply(cbind(...),MARGIN=1,identity)}), posts))
      ret = lapply(ret, function(r) {eval(predictor, envir = cbind(as.data.frame(t(r)), as.data.frame(loc)))})
    } else {
      ret = do.call(Map, c(list(function(...){apply(cbind(...),MARGIN=1,sum)}), posts))
    }
    ret = lapply(ret, link)
  } else {
    ret = do.call(cbind, posts)
    if (!is.null(predictor) && do.sum) { ret = eval(predictor, c(posts, as.list(data.frame(loc)))) } 
    else if ( do.sum ) { ret = apply(ret, MARGIN = 1, sum) }
    ret = link(ret)
  }
  
  return(ret)
}




#' A wrapper for the \link{inla} \link{f} function
#' 
#' g makes the mesh an explicit argument and catches the provided model for later usage.
#' See \link{f} for details on model specifications.
#'
#' @aliases g
#' @export
#' @param mesh An inla.mesh
#' @param model See \link{f} model specifications
#' @param ... Passed on to inla \link{f}
#' @return A list with mesh, model and the return value of the f-call

g = function(..., map, A.msk, mesh = NULL, model = "linear") {
  if (as.character(substitute(map))[[1]] == "" ) { map = NULL } else { map = substitute(map) }
  if (as.character(substitute(A.msk))[[1]] == "" ) { A.msk = NULL } else { A.msk = A.msk }
  
  fvals = f(..., model = model)
  
  if (fvals$model == "spde2" & is.null(map)) { map = expression(coordinates) }
  
  ret = list(mesh = mesh, 
             map = map,
             A.msk = A.msk,
             model = model, 
             f = fvals)
}


mapper = function(map, points) {
  # If map is a function, return map(points)
  # If map is not a function but a column in points, return points$map
  # Otherwise assume map is an expression that can be evaluated with points as an environment
  cat("mapping")
  if (is.function(map)) { loc = map(points) } 
  else if (class(map) == "call") { loc = eval(map, data.frame(points)) } 
  else {
    fetcher = get0(as.character(map))
    if (is.function(fetcher)) { loc = fetcher(points) } 
    else {
      if(as.character(map) %in% names(as.data.frame(points))) {
        loc = as.data.frame(points)[,as.character(map)] 
      } else {
        loc = rep(1, nrow(data.frame(points))) 
      }
    }
  }
}
