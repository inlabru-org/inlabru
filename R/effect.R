####################################################################################################
# GENERICS
####################################################################################################

#' @export
#' @rdname value.component
value = function(...){UseMethod("value")}
#' @export
#' @rdname amatrix.component
amatrix = function(...){UseMethod("amatrix")}
#' @export
#' @rdname map.component
map = function(...){UseMethod("map")}
#' @export
#' @rdname index.component
index = function(...){UseMethod("index")}


####################################################################################################
# CONSTRUCTORS
####################################################################################################
#' inlabru latent model component construction
#' 
#' @description 
#' 
#' Similar to glm(), gam() and inla() \link{bru} uses formula objects to describe response data and latent 
#' (unknonw) components of the model to be fitted. However, in addition to the syntax compatible with 
#' \link[INLA]{inla}, bru components offer addtitional functionality which facilitates modeling.
#' 
#' In inlabru, latent components can be constructed using R formulae or explicit parameter. For background 
#' information on the formulae inlabru accepts please see \link{component.formula}. For more details on the
#' model parameters and how inlabru employs INLA's \code{f} function please see \link{component.character}
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' @family Component constructor
#' @param object A formula or a label (character)
#' @param ... Arguments passed on to \link{component.formula} or \link{component.character}
#' @export 

component = function(object, ...){UseMethod("component")}


#' inlabru latent model component construction using formulae
#' 
#' @description 
#' 
#' Similar to glm(), gam() and inla() \link{bru} uses formula objects to describe response data and latent 
#' (unknonw) components of the model to be fitted. However, in addition to the syntax compatible with 
#' \link[INLA]{inla}, bru components offer addtitional functionality which facilitates modeling.
#' 
#' @details 
#' 
#' \link{bru} will understand formulae describing fixed effect models just like the other methods. For instance, the
#' formula \code{y ~ x} will fit the linear combination of an effect named \code{x} and an intercept to
#' the response \code{y} with respect to the likelihood family stated when calling \link{bru}. Mathematically,
#' the linear predictor \eqn{\eta} would be written down as 
#' 
#' \deqn{\eta = \beta * x + c,} 
#' 
#' where:
#' 
#' \itemize{
#' \item{\eqn{c} }{is the \emph{intercept}}
#' \item{\eqn{x }}{is a \emph{covariate}}
#' \item{\eqn{\beta} }{is a \emph{random variable} associated with \eqn{x} and} 
#' \item{\eqn{\psi = \beta * x }}{ is called the \emph{random effect} of \eqn{x}}
#' }
#' 
#' A problem that arises when using this kind of R formula is that it does not clearly relect the mathematical
#' formula. For instance, when providing the formula to inla, the resulting object will refer to the random
#' effect \eqn{\psi = \beta * x } as \code{x}. Hence, it is not clear if \code{x} refers to the covariate 
#' or the effect of the covariate.
#' 
#' @section Naming random effects:
#' 
#' In INLA, a simple random effect model would be expressed as 
#' 
#' \itemize{\item{\code{formula = y ~ f(x, model = "linear")},}}
#' 
#' where \link[INLA]{f} is the inla specific function to set up random effects of all kinds. The underlying 
#' predictor would again be \eqn{\eta = \beta * x + c} but the result of fitting the model would state 
#' \code{x} as the random effect's name. bru allows to rewrite this formula in order to explicitly state 
#' the name of the random effect and the name of the associated. This is achived by replacing \code{f}
#' with an arbitrary name that we wish to assign to the effect, e.g.
#' 
#' \itemize{\item{\code{components = y ~ psi(x, model = "linear")}.}}
#'
#' Being able to discriminate between \eqn{x} and \eqn{\psi} is relevant because of two functionalities
#' bru offers. The formula parameters of both, \link{bru} and the prediction method \link{predict.bru}
#' are interpreted in the mathematical sense. For instance, \code{predict} may be used to analyze the
#' an analytical combination of the covariate \eqn{x} and the intercept using
#' 
#' \itemize{\item{\code{predict(fit, data.frame(x=1)), ~ exp(x + Intercept)}.}}
#' 
#' On the other hand, predict may be used to only look at a transformation of the random effect \eqn{\psi}
#' 
#' \itemize{\item{\code{predict(fit, NULL, ~ exp(psi))}.}}
#'  
#' @aliases component.formula
#' @export
#' @family Component constructor
#' @param object A formula describing latent model components.
#' @param ... Ignored arguments (S3 compatibility)
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples
#' # As an example, let us create a linear component. Here, the component is 
#' # called "myLinearEffectOfX" while the covariate the component acts on is 
#' # called "x". Note that a list of components is returned because the 
#' # formula may define multiple components
#' 
#' eff = component(~ myLinearEffectOfX(map = x, model = "linear"))
#' summary(eff[[1]])
#' 


component.formula = function(object, ...) {
  code = code.components(object)
  parsed = lapply(code, function(x) parse(text=x))
  components = lapply(parsed, function(component.expression) eval(component.expression, envir = environment(object)))
  names(components) = lapply(components, function(x) x$label)
  components
}

#' inlabru latent model component construction using parameters
#' 
#' This function is inlabru's equivalent to INLA's \code{f} function but adds functionality that
#' is unique to inlabru. 
#'  
#' @aliases component.character
#' @family Component constructor
#' @export
#' @method component character
#' @param object A string giving the component its name
#' @param data EXPERIMENTAL
#' @param model Either one of "offset", "factor", "linear" or a model accepted by INLA's \code{f} function
#' @param map EXPERIMENTAL
#' @param n EXPERIMENTAL
#' @param season.length EXPERIMENTAL
#' @param group EXPERIMENTAL
#' @param values EXPERIMENTAL
#' @param A.msk Boolean vector for masking (deactivating) columns of the A-matrix
#' @param ... EXPERIMENTAL
#' @return An component object
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples
#' \donttest{
#' if (require("INLA", quietly = TRUE)) {
#'     
#' # As an example, let us create a linear component. Here ,the component is 
#' # called "myEffectOfX" while the covariate the component acts on is called "x":
#' 
#' eff = component("myEffectOfX", model = "linear", map = x)
#' summary(eff)
#' 
#' # A more complicated component:
#' eff = component("myEffectOfX", model = inla.spde2.matern(inla.mesh.1d(1:10)), map = x)
#' }
#' }

component.character = function(object,
                             data,
                             model,
                             map,
                             n = NULL,
                             season.length = NULL,
                             group = NULL,
                             replicate = NULL,
                             values = NULL,
                             A.msk = NULL,
                             ...){
  
  # INLA models:
  # itypes = c(linear, iid, mec, meb, rgeneric, rw1, rw2, crw2, seasonal, besag, besag2, bym, bym2, besagproper, 
  #            besagproper2, fgn, fgn2, ar1, ar1c, ar, ou, generic, generic0, generic1, generic2, generic3, spde, 
  #            spde2, spde3, iid1d, iid2d, iid3d, iid4d, iid5d, 2diid, z, rw2d, rw2diid, slm, matern2d, copy,
  #            clinear, sigm, revsigm, log1exp, logdist) 
  #
  # Supported:
  # btypes = c("offset", "factor", "linear", "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")
  
  # The label
  label = object
  
  # Model type as character
  model.type = model
  if ( inherits(model, "inla.spde") ) { model.type = "spde" }
  
  # Default component (to be filled)
  component = list(label = label,
                inla.formula = NA,
                type = model.type,
                map = substitute(map),
                mesh = NA,
                group.char = as.character(substitute(group)), # Name of the data column holding the group index
                replicate.char = as.character(substitute(replicate)), # Name of the data column holding the replicate index
                values = values,
                A.msk = A.msk,
                model = model,
                env = parent.frame())
  
  # Main bit
  if ( model.type %in% c("offset") ) {
    component$inla.formula = as.formula(paste0("~ . + offset(offset)"))
  } 
  else if ( model.type %in% c("factor") ) {
    component$inla.formula = as.formula(paste0("~ . + ", label))
  }
  else {
    
    # Construct a call to the f function from the parameters provided
    # Ultimately, this call will be converted to the model formula presented to INLA
    fcall = sys.call()
    fcall[[1]] = "f"
    fcall[[2]] = as.symbol(label)
    
    # Remove parameters inlabru supports but INLA doesn't
    fcall = fcall[!(names(fcall) %in% c("map","A.msk", "mesh"))]
    
#    # For SPDE models we need a little nasty trick
#   if ( model.type %in% c("spde") ) {
#     #tmp = NA
#     #fcall$group = as.symbol("tmp") 
#   }
    
    # Arguments to f as list
    f.args = as.list(fcall[2:length(fcall)])
    component$f.args = f.args
    
    # A trick for "Copy" models
    if ("copy" %in% names(f.args)) { 
      f.args[["copy"]] = NULL
      fcall$model = NULL
    }
    
    # Protect arguments that may need access to actual data
    # TODO: make a more general system, that also handles ngroup etc.
    if ("group" %in% names(f.args)) {
      f.args[["group"]] <- substitute(group)
    }
    if ("replicate" %in% names(f.args)) {
      f.args[["replicate"]] <- substitute(replicate)
    }
    
    # Call f and store the results
    # TODO: Investigate what outputs from INLA::f require data, and access that
    # information with some other method. ::f should only be called by INLA,
    # and never by inlabru, since it sets up temporary files that will not be
    # removed, and requires data not available at this point!
    # Only a small handful of component$f elements are actually accessed
    # by inlabru code.
##    fvals <- do.call(INLA::f, f.args, envir = parent.frame())
    component$f = NULL
    
    # Second part of the SPDE model trick above
    # TODO: generalise group and replicate handling for all models
    if ( model.type %in% c("spde") ) { 
      fcall$group = as.symbol(paste0(label, ".group"))
    }
    
    # Generate the formula that will be presented to INLA
    component$inla.formula = as.formula(paste0("~ . + ", paste0(deparse(fcall), collapse = "")))
    
    # Set the default mesh used for interpolation
    component$mesh = make.default.mesh(component, model, model.type, fvals)

    # Set ngroup and nrep gefaults
    if (is.null(component$f$ngroup)) { component$ngroup = 1 } else { component$ngroup = component$f$ngroup }
    if (is.null(component$f$nrep)) { component$nrep = 1 } else { component$nrep = component$f$nrep }
    
  }

  class(component) = c("component","list")
  component
}

# Picks a default mesh given the type of model and parameters provided
make.default.mesh = function(component, model, model.type, fvals = NULL){
  
  miss.msg = paste0("component '%s' (type '%s') requires argument '%s'. ",
                    "Check out f() for additional information on this argument.")
  
  if ( model.type %in% c("linear", "clinear") ) {
    mesh = NA
  }
  else if ( model.type %in% c("iid") ) {
    if ( is.null(fvals$n) ) {
      stop(sprintf(miss.msg, component$label,
                   model.type, "n"))
      }
    mesh = INLA::inla.mesh.1d(1:fvals$n)
  }
  else if ( model.type %in% c("seasonal") ) {
    if ( is.null(fvals$season.length) ) {
      stop(sprintf(miss.msg, component$label,
                   model.type, "season.length"))
    }
    mesh = INLA::inla.mesh.1d(1:fvals$season.length)
  }
  else if ( model.type %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
    # TODO: Check that equal/unequal spacing in these models works
    #       in the desired way, and document what this means.
    if ( is.null(fvals$values) ) {
      stop(sprintf(miss.msg, component$label,
                   model.type, "values"))
    }
    mesh = INLA::inla.mesh.1d(sort(unique(fvals$values)))
  }
  else if ( model.type %in% c("spde") ) {
    mesh = model$mesh
  } 
  else {
    stop(paste0("component type '", model.type, "' not implemented."))
  }
  
  mesh
}

#' Convert components to R code
#'  
#' @aliases code.components
#' @keywords internal
#' @param components A \link{formula} describing latent model components.
#' @param fname Chracter setting the name of the function that will interpret the components.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

code.components = function(components) {
  fname = "component"
  tms = terms(components)
  codes = attr(tms, "term.labels")
  
  # Check for offset()
  isoff = as.vector(unlist(lapply(rownames(attr(tms, "factors")), function(s) substr(s,1,6)=="offset")))
  if (any(isoff)) {
    codes[[length(codes)+1]] = rownames(attr(tms, "factors"))[isoff]
  }
  
  for (k in 1:length(codes)){
    code = codes[[k]]
    
    # Function syntax or fixed component?
    ix = regexpr("(", text = code, fixed = TRUE)
    is.offset = FALSE
    if (ix > 0) {
      label = substr(code, 1, ix-1)
      is.fixed = FALSE
      if (label == "offset") {is.offset = TRUE} 
    } else {
      label = code
      is.fixed = TRUE
    }
    
    # Make code
    if ( is.fixed ) {
      codes[[k]] = sprintf("%s(\"%s\", map = %s, model = 'linear')", fname, label, label)
    }
    else if ( is.offset ) {
      codes[[k]] = gsub(paste0(label,"("), paste0(fname,"(\"",label,"\", model = \"offset\", map = "), code, fixed = TRUE)
    }
    else {
      codes[[k]] = gsub(paste0(label,"("), paste0(fname,"(\"",label,"\", "), code, fixed = TRUE)
    }

    
  }
  
  codes
}


####################################################################################################
# OPERATORS
####################################################################################################


#' Summarize an component
#'  
#' @aliases summary.component
#' @export
#' @method summary component
#' @keywords internal
#' @param object An component.
#' @param ... ignored.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

summary.component = function(object, ...) {
  
  eff <- list('Label' = object$label,
              'Type' = object$type,
              'Map' = sprintf("%s [class: %s]",
                              deparse(object$map),
                              class(object$map)),
              'INLA formula' = as.character(object$inla.formula))
  class(eff) <- c('summary.component', 'list')
  eff
}

#' Print the summary of an component
#'  
#' @aliases print.summary.component
#' @export
#' @method print summary.component
#' @keywords internal
#' @param x A 'summary.component' object.
#' @author Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @rdname summary.component

print.summary.component = function(x, ...) {
  for (name in names(x)) {
    # Split TAB character to attempt proper printing in RStudio,
    # but even though this makes a difference on the command line,
    # there's no component when inside the function. /FL
    cat(name, ":", "\t", x[[name]], "\n", sep="")
  }
  invisible(x)
}



#' Evaluate an component
#' 
#' Calculates a latent component given some data and the state of the component's internal random variables.
#' 
#' TODO: Improve speed for iterated calls by making 'mapped' a parameter 
#' 
#' @aliases value.component
#' @export
#' @method value component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param state Either a numeric vector or a list with a numeric entry whose name is equal to the name parameter.
#' @param A A matrix overriding the default projection matrix.
#' @param ... Unused.
#' @return A numeric vector of component values
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'


value.component = function(component, data, state, A = NULL, ...) {
  
  # Convenience: extract state if a list of states was provided
  if ( is.list(state) &!is.data.frame(state)) { state = state[[component$label]] }
  
  # Obtain covariates
  mapped = mapper(component$map, data, component, component$env)
  
  if (is.data.frame(mapped)) { 
    if (component$label %in% names(mapped)) {
      mapped = mapped[,component$label,drop=TRUE]
    } else {
      mapped = mapped[,1,drop=TRUE] 
    }
  }
    
  # Make A-matrix (if not provided)
  if ( is.null(A) ) { A = amatrix(component, data) }
  
  # Determine component depending on the type of latent model
  if ( component$type %in% c("linear", "clinear") ) {

    values = A %*% (state * mapped)
  }
  else if ( component$type %in% c("offset") ) {

    values = A %*% mapped
  }
  else if ( component$type %in% c("factor") ) {

    values = A %*% state[mapped]
  } 
  else if ( component$type %in% c("iid", "seasonal") ) {
    
    values = A %*% state[mapped]
    
  }
  else if ( component$type %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
    
    values = A %*% state
    
  }
  else if ( component$type %in% c("spde") ) {

    values = A %*% state
    
  } else {
    stop(paste0("Evaluation of ", component$type, " not implemented."))
  }
  
  as.vector(values)
}




#' Construct A-matrix
#' 
#' Constructs the A-matrix for a given component and and some data
#'  
#' @aliases amatrix.component
#' @export
#' @method amatrix component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

amatrix.component = function(component, data, ...) {
  
  if ( component$type %in% c("spde") ) {
    if ( component$map == "coordinates" ) {
      if ( is.na(proj4string(data)) | is.null(component$mesh$crs) ) {
        loc = data
      } else {
        loc = stransform(data, crs = component$mesh$crs)
      }
      
    } else {
      loc = mapper(component$map, data, component)
      if ( !is.matrix(loc) & !inherits(loc,"Spatial") ) loc = as.matrix(loc)
    }
    
    
    if (component$ngroup > 1) {
      group = data[[component$group.char]]
    } else { group = NULL }
    
    A = INLA::inla.spde.make.A(component$mesh, 
                               loc = loc, 
                               group = group, 
                               n.group = component$ngroup, 
                               n.repl = component$nrep)
  } 
  else if ( component$type %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
    if ( is.null(component$values) ) { stop("Parameter 'values' not set for component '", component$label, "'.") }
    mesh = INLA::inla.mesh.1d(component$values)
    loc = mapper(component$map, data, component)
    A = INLA::inla.spde.make.A(mesh, loc)
  }
  else {
    if ( is.null(data) ) { A = 1 }
    else { A = Matrix::Diagonal(nrow(data.frame(data))) }
  }

  # Mask columns of A
  if (!is.null(component$A.msk)) { 
    A = A[, as.logical(component$A.msk), drop=FALSE]
  }
  
  # Weight rows of A
  # if (!is.null(component$weights)) {
  #   A = as.matrix(A)*as.vector(component$weights)
  # }
  
  A
}

#' Obtain covariate
#' 
#' @section Simple covariates and the map parameter:
#' 
#' It is not unusual for a random effect act on a transformation of a covariate. In other frameworks this
#' would mean that the transformed covariate would have to be calculated in advance and added to the 
#' data frame that is usually provided via the \code{data} parameter. inlabru provides the option to do
#' this transformation automatically. For instance, one might be interested in the effect of a covariate
#' \eqn{x^2}. In inla and other frameworks this would require to add a column \code{xsquared} to the
#' input data frame and use the formula
#' 
#' \itemize{\item{\code{formula = y ~ f(xsquared, model = "linear")},}}
#' 
#' In inlabru this can be achived using two ways of using the \code{map} parameter. 
#' 
#' \itemize{
#' \item{\code{components = y ~ psi(map = x^2, model = "linear")}}
#' \item{\code{components = y ~ psi(map = mySquareFun(x), model = "linear")},}
#' \item{\code{components = y ~ psi(map = myOtherSquareFun, model = "linear")},}
#' 
#' }
#' 
#' In the first example inlabru will interpret the map parameter as an expression to be evaluated within
#' the data provided. Since \eqn{x} is a knonwn covariate it will know how to calculate it. The second
#' example is an expression as well but it uses a function alled \code{mySquareFun}. This function is
#' defined by user but has wo be accessible within the work space when setting up the compoonents. 
#' The third example provides the function \code{myOtherSquareFun} directly and not within an expression.
#' In this case, inlabru will call the function using the data provided via the  \code{data} parameter. 
#' inlabru expects that the output of this function is a data.frame with "psi" being the name of the 
#' single existing column. For instance, 
#' 
#' \code{myOtherSquareFun = function(data) {
#'                             data = data[,"x", drop = FALSE] ; 
#'                             colnames(data) = "psi" ;
#'                             return(data)}}
#' 
#' @section Spatial Covariates:
#' 
#' When fitting spatial models it is common to work with covariates that depend on space, e.g. sea
#' surface temperature or elevation. Although it is straight forward to add this data to the input
#' data frame or write a covariate function like in the previous section there is an even more 
#' convenient way in inlabru. Spatial covariates are often stored as \code{SpatialPixelDataFrame}, 
#' \code{SpatialPixelDataFrame} or \code{RasterLayer} objects. These can be provided directly via 
#' the map parameter if the input data is a \code{SpatialPointsDataFrame}. inlabru will automatically
#' evaluate and/or interpolate the coariate at your data locations when using code like
#' 
#' \itemize{\item{\code{components = y ~ psi(mySpatialPixels, model = "linear")}.}}
#'
#' @section Coordinates:
#' 
#' A common spatial modelling component when using inla are SPDE models. An important feature of
#' inlabru is that it will automatically calculate the so called A-matrix which maps SPDE 
#' values at the mesh vertices to values at the data locations. For this purpose, the map parameter
#' can be se to \code{coordinates}, which is the \code{sp} package function that extracts point 
#' coordinates from the SpatialPointsDataFrame that was provided as input to bru. The code for
#' this would look as follows:
#' 
#' \itemize{\item{\code{components = y ~ mySPDE(map = coordinates, model = inla.spde2.matern(...))}.}}
#'  
#' @aliases map.component
#' @export
#' @method map component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the component's map.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

map.component = function(component, data, ...) {
  
  mp = mapper(component$map, data, component, component$env)

}



#' Obtain indices
#' 
#' Idexes into to the components
#'  
#' @aliases index.component
#' @export
#' @method index component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the component's map.
#' @param ... Unused.
#' @return a data.frame of indices or list of indices into the components latent variables
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

index.component = function(component, data, ...) {
  
  if ( component$type %in% c("spde") ) {
    if ( "m" %in% names(component$mesh) ) {
      # idx = 1:component$mesh$m # support inla.mesh.1d models
      idx = INLA::inla.spde.make.index(name = component$label, n.spde = component$mesh$m)
      # If a is masked, correct number of indices
      if ( !is.null(component$A.msk) ) { idx = 1:sum(component$A.msk) }
    } else {
      
      idx = INLA::inla.spde.make.index(component$label, 
                                       n.spde = component$model$n.spde, 
                                       n.group = component$ngroup,
                                       n.repl = component$nrep)
    }
  }
  else if ( component$type %in% c("factor") ) {
    idx = map.component(component, data)
    if (!is.data.frame(idx)) { idx = data.frame(idx) ; colnames(idx) = component$label }
    idx[,1] = as.factor(paste0(component$label, idx[,1]))
    #idx[,1] = as.factor(idx[,1])
  }
  else {
    idx = map.component(component, data)
    if (!is.data.frame(idx)) {
      idx = data.frame(idx)
      colnames(idx) = component$label
    }
  }

  idx
  
}



mapper = function(map, points, eff, env = NULL) {
  
  # Evaluate the map with the points as an environment
  emap = tryCatch(eval(map, c(data.frame(points), as.list(env))), error = function(e) {})
  
  # 0) Eval failed. map everything to 1. This happens for automatically added Intercept
  # 1) If we obtain a function, apply the function to the points
  # 2) If we obtain a SpatialGridDataFrame extract the values of that data frame at the point locations using the over() function
  # 3) Else we obtain a vector and return as-is. This happens when map states a column of the data points
  
  if ( is.null(emap) ) { 
    loc = data.frame(x = rep(1, max(1,nrow(as.data.frame(points)))))
    colnames(loc) = deparse(map)
  }
  else if ( is.function(emap) ) { loc = emap(points) }
  else if ( inherits(emap, "SpatialGridDataFrame") | inherits(emap, "SpatialPixelsDataFrame")) {
    if ( length(eff$group.char) == 0 ) {
      loc = over(points, emap)[,1,drop=FALSE]
      colnames(loc) = eff$label
    } else {
      layr = points[[eff$group.char]]
      loc = vector()
      for (l in unique(layr)) { loc[layr == l] = over(points[layr == l,], emap)[,l] }
      loc = data.frame(loc = loc)
      colnames(loc) = eff$label
    }
    
  }
  else if ( eff$label == "offset" && is.numeric(emap) && length(emap)==1 ) { loc = data.frame(offset = rep(emap, nrow(points)))}
  else { loc = emap }
  
  # Always return as many rows as points has
  if (is.vector(loc) && (length(loc) == 1) && (nrow(as.data.frame(points))>1)) {
    loc = rep(loc, nrow(as.data.frame(points)))
  }
  
  # Check if any of the locations are NA. If we are dealing with SpatialGridDataFrame try
  # to fix that by filling in nearest neighbor values.
  if ( any(is.na(loc)) & ( inherits(emap, "SpatialGridDataFrame") | inherits(emap, "SpatialPixelsDataFrame") )) {
    
    warning(sprintf("Map '%s' has returned NA values. As it is a SpatialGridDataFrame I will try to fix this by filling in values of spatially nearest neighbors. In the future, please design your 'map=' argument as to return non-NA for all points in your model domain/mesh. Note that this can also significantly increase time needed for inference/prediction!",
                    deparse(map), eff$label))
    
    BADpoints = points[as.vector(is.na(loc)),]
    GOODpoints = points[as.vector(!is.na(loc)),]
    dst = rgeos::gDistance(SpatialPoints(GOODpoints),SpatialPoints(BADpoints), byid=T)
    nn = apply(dst, MARGIN = 1, function(row) which.min(row)[[1]])
    loc[is.na(loc)] = loc[!is.na(loc)][nn]
    colnames(loc) = eff$label
  }
  
  # Check for NA values.    
  if ( any(is.na(loc)) ) {
    stop(sprintf("Map '%s' of component '%s' has returned NA values. Please design your 'map='
                 argument as to return non-NA for all points in your model domain/mesh.",
                 as.character(map)[[1]], eff$label))
  }
  
  loc
}



