# GENERICS ----

#' @export
#' @rdname value.component
value <- function(...) {
  UseMethod("value")
}
#' @export
#' @rdname amatrix_eval.component
amatrix_eval <- function(...) {
  UseMethod("amatrix_eval")
}
#' @export
#' @rdname input_eval.component
input_eval <- function(...) {
  UseMethod("input_eval")
}
#' @export
#' @rdname index_eval.component
index_eval <- function(...) {
  UseMethod("index_eval")
}



# CONSTRUCTORS ----

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

component <- function(object, ...) {
  UseMethod("component")
}


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
#' eff <- component(~ myLinearEffectOfX(map = x, model = "linear"))
#' summary(eff[[1]])
#' # Equivalent shortcuts:
#' eff <- component(~ myLinearEffectOfX(x, model = "linear"))
#' eff <- component(~ myLinearEffectOfX(x))
component.formula <- function(object, ...) {
  code <- code.components(object)
  parsed <- lapply(code, function(x) parse(text = x))
  components <- lapply(parsed,
                       function(component.expression) {
                         eval(component.expression,
                              envir = environment(object))
                       })
  names(components) <- lapply(components, function(x) x$label)
  class(components) <- c("component_list", "list")
  components
}

bru_input <- function(input, label = NULL, layer = NULL, selector = NULL) {
  inp <- list(input = input,
              label = label,
              layer = layer,
              selector = selector)
  class(inp) <- c("bru_input", "list")
  inp
}

bru_subcomponent <- function(input = NULL,
                             mapper = NULL,
                             model = NULL,
                             n = NULL,
                             values = NULL,
                             season.length = NULL,
                             weights = NULL) {
  type <- model
  if (inherits(model, "inla.spde")) {
    type <- "spde"
  } else if (inherits(model, "inla.rgeneric")) {
    type <- "rgeneric"
  } else if (inherits(model,
                      c("clinear", "sigm", "revsigm",
                        "log1exp", "logdist"))) {
    type <- "specialnonlinear"
  } else if (is.character(model)) {
    type <- model
  } else {
    type <- "unknown"
  }
  subcomponent <-
    list(input = input,
         mapper = mapper,
         model = model,
         type = type,
         n = n,
         values = values,
         season.length = season.length,
         weights = weights)
  class(subcomponent) <- c("bru_subcomponent", "list")
  subcomponent
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
#' @param replicate EXPERIMENTAL
#' @param values EXPERIMENTAL
#' @param A.msk Boolean vector for masking (deactivating) columns of the A-matrix
#' @param ... EXPERIMENTAL
#' @return An component object
#' 
#' @details 
#' Deprecated: \code{map} (now \code{main}), \code{mesh} (now \code{mapper})
#'
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'
#' @examples
#' \donttest{
#' if (require("INLA", quietly = TRUE)) {
#'
#'   # As an example, let us create a linear component. Here ,the component is
#'   # called "myEffectOfX" while the covariate the component acts on is called "x":
#'
#'   eff <- component("myEffectOfX", main = x, model = "linear")
#'   summary(eff)
#'
#'   # A more complicated component:
#'   eff <- component("myEffectOfX",
#'                    main = x,
#'                    model = inla.spde2.matern(inla.mesh.1d(1:10)))
#' }
#' }

component.character <- function(object,
                                # Main model parameters
                                main = NULL, # This must be kept as 1st arg.
                                model = NULL,
                                mapper = NULL,
                                main_layer = NULL,
                                main_selector = NULL,
                                n = NULL,
                                values = NULL,
                                season.length = NULL,
                                # Weights
                                weights = NULL,
                                weights_layer = NULL,
                                # Group model parameters
                                group = NULL,
                                group_layer = NULL,
                                control.group = NULL,
                                # Replicate model parameters
                                replicate = NULL,
                                replicate_layer = NULL,
                                A.msk = NULL,
                                # Deprecated parameters
                                # map -> main
                                map = NULL,
                                # mesh -> mapper
                                mesh = NULL,
                                ...) {

  # INLA models:
  # itypes = c(linear, iid, mec, meb, rgeneric, rw1, rw2, crw2, seasonal, besag, besag2, bym, bym2, besagproper,
  #            besagproper2, fgn, fgn2, ar1, ar1c, ar, ou, generic, generic0, generic1, generic2, generic3, spde,
  #            spde2, spde3, iid1d, iid2d, iid3d, iid4d, iid5d, 2diid, z, rw2d, rw2diid, slm, matern2d, copy,
  #            clinear, sigm, revsigm, log1exp, logdist)
  #
  # Supported:
  # btypes = c("offset", "factor", "linear", "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")

  # The label
  label <- object
  
  if (is.null(model)) {
    # TODO: may need a special marker to autodetect factor variables,
    #       which can only be autodetected in a data-aware pass.
    model <- "linear"
  }

  if (!is.null(substitute(group))) {
    if (is.null(control.group)) {
      control.group <- INLA::inla.set.control.group.default()
    }
    group_model <- control.group$model
  } else {
    group_model <- NULL
  }

  if (is.null(main_layer)) {
    main_layer <- label
  }

  if (!is.null(substitute(map))) {
    if (is.null(substitute(main))) {
      main <- substitute(map)
      warning("Use of 'map' is deprecated and may be disabled; use 'main' instead.")
    } else {
      warning("Deprecated 'map' overridden by 'main'.")
    }
  }

  if (!is.null(mesh)) {
    if (is.null(mapper)) {
      mapper <- mesh
      warning("Use of 'mesh' is deprecated and may be disabled; use 'mapper' instead.")
    } else {
      warning("Deprecated 'mesh' overridden by 'mapper'.")
    }
  }

  # Default component (to be filled)
  component <- list(
    label = label,
    inla.formula = NULL,
    main = bru_subcomponent(
      input = bru_input(substitute(main),
                        label = label,
                        layer = main_layer,
                        selector = main_selector),
      mapper = mapper,
      model = model,
      n = n,
      values = substitute(values),
      season.length = season.length,
      weights =
        if (is.null(substitute(weights))) {
          NULL
        } else {
          bru_input(substitute(weights),
                    label = paste0(label, ".weights"),
                    layer = weights_layer)
        }),
    group = bru_subcomponent(
      input = bru_input(substitute(group),
                        label = paste0(label, ".group"),
                        layer = group_layer),
      mapper = NULL,
      n = NULL,
      model = group_model),
    replicate = bru_subcomponent(
      input = bru_input(substitute(replicate),
                        label = paste0(label, ".repl"),
                        layer = replicate_layer),
      mapper = NULL,
      n = NULL,
      model = "iid"),
    A.msk = A.msk,
    env = parent.frame()
  )
  
  # Main bit
  if (component$main$type %in% c("offset")) {
    component$inla.formula <- as.formula(paste0("~ . + offset(offset)"))
  } else if (component$main$type %in% c("factor")) {
    component$inla.formula <- as.formula(paste0("~ . + ", label))
  }
  else {
    # Construct a call to the f function from the parameters provided
    # Ultimately, this call will be converted to the model formula presented to INLA
    fcall <- sys.call()
    fcall[[1]] <- "f"
    fcall[[2]] <- as.symbol(label)
    fcall[["model"]] <- as.symbol(model)
    
    # Remove parameters inlabru supports but INLA doesn't,
    # and substitute parameters that inlabru will transform
    fcall <- fcall[!(names(fcall) %in% c("main",
                                         "mapper",
                                         "group_mapper",
                                         "replicate_mapper",
                                         "A.msk",
                                         "map",
                                         "mesh"))]

    # A trick for "Copy" models
    if ("copy" %in% names(fcall)) {
      ### TODO:
      ### Check this. Where is the copy information stored and used?
      ### Should postprocess components to link them together.
      fcall[["copy"]] <- NULL
      fcall$model <- NULL
    }

    # Replace arguments that will be evaluated by a mapper
    # TODO: make a more general system, that also handles ngroup etc.
    suffixes <- list("group" = "group",
                     "replicate" = "repl",
                     "values" = "values",
                     "weights" = "weights")
    for (arg in names(suffixes)) {
      if (arg %in% names(fcall)) {
        fcall[[arg]] <- as.symbol(paste0(label, ".", suffixes[[arg]]))
      }
    }

    # These values were previously initialised by an INLA::f call, but ::f
    # should only be called by INLA, and never by inlabru, since it sets up
    # temporary files that will not be removed, and also requires data not
    # available at this point!  Until multi-stage model initialisation is
    # implemented, require the user to explicitly provide these values.
    if (!is.null(component$main$n)) {
      fcall[["n"]] <- component$main$n
    }
    if (!is.null(season.length)) {
      fcall[["season.length"]] <- season.length
    }
#    # TODO: This is a temporary compatibility:
#    fcall[["values"]] <- as.symbol(paste0(label, ".values"))
    
    component$inla.formula <-
      as.formula(paste0("~ . + ",
                        as.character(parse(text = deparse(fcall)))))
    
    ## TODO: All the below needs to be in the second pass.
    
    # Set the default mesh used for interpolation
    # TODO: Adjust this for a first and second pass of initalisation
    component$main <-
      make_mapper(component$main,
                  label,
                  input_values = NULL,
                  strict = FALSE)
    component$group <-
      make_mapper(component$group,
                  paste(label, "group"),
                  input_values = NULL,
                  strict = FALSE)
    component$replicate <-
      make_mapper(component$replicate,
                  paste(label, "replicate"),
                  input_values = NULL,
                  strict = FALSE)
    
    # Set ngroup and nrep defaults
    if (is.null(component$group$n)) {
      fcall[["ngroup"]] <- 1
    } else {
      fcall[["ngroup"]] <- component$group$n
    }
    if (is.null(component$replicate$n)) {
      fcall[["nrep"]] <- 1
    } else {
      fcall[["nrep"]] <- component$replicate$n
    }

    # Generate the formula that will be presented to INLA
    component$inla.formula <-
      as.formula(paste0("~ . + ",
                        paste0(deparse(fcall),
                               collapse = "\n")))
    
    component$fcall <- fcall
    
  }

  class(component) <- c("component", "list")
  component
}

# Defines a default mapper given the type of model and parameters provided
# Checks mapping$mapper, model$mesh (for spde models), mapping$values,
# input_values or mapping$n, in that order.
make_mapper <- function(subcomp,
                        label,
                        input_values = NULL,
                        strict = TRUE) {
  miss.msg <- paste0(
    "component '%s' (type '%s') requires argument '%s'. ",
    "Check out f() for additional information on this argument."
  )

  if (!is.null(subcomp[["mapper"]])) {
    if (inherits(subcomp[["mapper"]], "inla.mesh.1d")) {
      subcomp$n <- subcomp[["mapper"]]$m
    } else if (inherits(subcomp[["mapper"]], "inla.mesh")) {
      subcomp$n <- subcomp[["mapper"]]$n
    }
  } else if (subcomp[["type"]] %in% c("spde")) {
    subcomp$mapper <- subcomp[["model"]]$mesh
    subcomp$n <- subcomp[["model"]]$n.spde
  } else if (subcomp[["type"]] %in% "linear") {
    subcomp$n <- 1
    subcomp$values <- 1
  } else if (subcomp[["type"]] %in% c("iid", "exchangeable")) {
    if (is.null(subcomp[["n"]])) {
      if (strict) {
        stop(sprintf(
          miss.msg, label,
          subcomp[["type"]], "n"
        ))
      }
    } else {
      subcomp$mapper <- INLA::inla.mesh.1d(seq_len(subcomp[["n"]]))
    }
  }
  else if (subcomp[["type"]] %in% c("seasonal")) {
    if (is.null(subcomp[["season.length"]])) {
      if (strict) {
        stop(sprintf(
          miss.msg, label,
          subcomp[["type"]], "season.length"
        ))
      }
    } else {
      subcomp$mapper <- INLA::inla.mesh.1d(seq_len(subcomp[["season.length"]]))
    }
  } else {
    if (!is.null(subcomp[["values"]])) {
      subcomp$mapper <- INLA::inla.mesh.1d(sort(unique(subcomp[["values"]])))
      subcomp$n <- subcomp$mapper$n
      subcomp$values <- subcomp$mapper$loc
    } else if (!is.null(input_values)) {
      subcomp$mapper <- INLA::inla.mesh.1d(sort(unique(input_values),
                                                na.last = NA))
      subcomp$n <- subcomp$mapper$n
      subcomp$values <- subcomp$mapper$loc
    } else if (strict) {
      stop(sprintf(
        miss.msg, label,
        subcomp[["type"]], "values"
      ))
    }
  }

  subcomp
}

#' Convert components to R code
#'
#' @aliases code.components
#' @keywords internal
#' @param components A \link{formula} describing latent model components.
#' @param fname Chracter setting the name of the function that will interpret the components.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

code.components <- function(components) {
  fname <- "component"
  tms <- terms(components)
  codes <- attr(tms, "term.labels")

  # Check for offset()
  isoff <- as.vector(unlist(lapply(rownames(attr(tms, "factors")), function(s) substr(s, 1, 6) == "offset")))
  if (any(isoff)) {
    codes[[length(codes) + 1]] <- rownames(attr(tms, "factors"))[isoff]
  }

  for (k in 1:length(codes)) {
    code <- codes[[k]]

    # Function syntax or fixed component?
    ix <- regexpr("(", text = code, fixed = TRUE)
    is.offset <- FALSE
    if (ix > 0) {
      label <- substr(code, 1, ix - 1)
      is.fixed <- FALSE
      if (label == "offset") {
        is.offset <- TRUE
      }
    } else {
      label <- code
      is.fixed <- TRUE
    }

    # Make code
    if (is.fixed) {
      codes[[k]] <- sprintf("%s(\"%s\", main = %s, model = 'linear')",
                            fname, label, label)
    }
    else if (is.offset) {
      codes[[k]] <-
        gsub(paste0(label, "("),
             paste0(fname, "(\"", label, "\", model = \"offset\", main = "),
             code, fixed = TRUE)
    }
    else {
      codes[[k]] <- gsub(paste0(label, "("),
                         paste0(fname, "(\"", label, "\", "),
                         code, fixed = TRUE)
    }
  }

  codes
}


# OPERATORS ----


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

summary.component <- function(object, ...) {
  eff <- list(
    "Label" = object$label,
    "Type" = object$type,
    "Main" = sprintf(
      "%s [class: %s]",
      deparse(object$main$input),
      class(object$main$input)
    ),
    "INLA formula" = as.character(object$inla.formula)
  )
  class(eff) <- c("summary.component", "list")
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

print.summary.component <- function(x, ...) {
  for (name in names(x)) {
    # Split TAB character to attempt proper printing in RStudio,
    # but even though this makes a difference on the command line,
    # there's no component when inside the function. /FL
    cat(name, ":", "\t", x[[name]], "\n", sep = "")
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


value.component <- function(component, data, state, A = NULL, ...) {

  # Convenience: extract state if a list of states was provided
  if (is.list(state) & !is.data.frame(state)) {
    state <- state[[component$label]]
  }

  # Obtain covariates
  main_input <- input_eval(component, part = "main", data = data, env = component$env)
  group_input <- input_eval(component, part = "group", data = data, env = component$env)
  repl_input <- input_eval(component, part = "repl", data = data, env = component$env)
  weights_input <- input_eval(component, part = "weights", data = data, env = component$env)

  if (is.data.frame(main_input)) {
    message("'main_input' is a data.frame! This code should be moved to input_eval.bru_input().")
    if (component$label %in% names(main_input)) {
      main_input <- main_input[, component$label, drop = TRUE]
    } else {
      main_input <- main_input[, 1, drop = TRUE]
    }
  }

  # Make A-matrix (if not provided)
  if (is.null(A)) {
    A <- amatrix(component, data)
  }

  # Determine component depending on the type of latent model
  if (component$type %in% c("linear", "clinear")) {
    values <- A %*% (state * mapped)
  }
  else if (component$type %in% c("offset")) {
    values <- A %*% mapped
  }
  else if (component$type %in% c("factor")) {
    values <- A %*% state[mapped]
  }
  else if (component$type %in% c("iid", "seasonal")) {
    values <- A %*% state[mapped]
  }
  else if (component$type %in% c("rw1", "rw2", "ar", "ar1", "ou")) {
    values <- A %*% state
  }
  else if (component$type %in% c("spde")) {
    values <- A %*% state
  } else {
    stop(paste0("Evaluation of ", component$type, " not implemented."))
  }

  as.vector(values)
}





#' Construct A-matrix
#'
#' Constructs the A-matrix for a given component and and some data
#'
#' @aliases amatrix_eval.component
#' @export
#' @method amatrix_eval component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

amatrix_eval.bru_subcomponent <- function(subcomp, data, env = NULL, ...) {
  ## TODO: Check for offset component

  if (!is.null(subcomp$weights)) {
    weights <- input_eval(subcomp$weights, data, env = env)
  } else {
    weights <- 1.0
  }
  if (inherits(subcomp$mapper, c("inla.mesh", "inla.mesh.1d"))) {
    val <- input_eval(subcomp$input, data, env = env)
    if (!is.matrix(val) & !inherits(val, "Spatial")) {
      val <- as.matrix(val)
    }
    A <- INLA::inla.spde.make.A(subcomp$mapper, loc = val, weights = weights)
  } else if (subcomp$model %in% c("linear", "clinear")) {
    val <- input_eval(subcomp$input, data, env = env)
    A <- Matrix::sparseMatrix(i = seq_len(nrow(data)),
                              j = rep(1, nrow(data)),
                              x = val * weights,
                              dims = c(nrow(data), subcomp$n))
  } else if (is.null(subcomp$mapper)) {
    if (subcomp$n > 1) {
      stop(paste0("Missing mapper (NULL) for subcomponent '", subcomp$label, "'"))
    }
    A <- Matrix::Diagonal(1, 1.0)
  } else {
    stop(paste0("Unsupported mapper of class '",
                paste0(class(component$main$mapper), collapse = "', '"), "'",
                " for subcomponent '", subcomp$label, "'"))
  }
  
  A
}



#' Construct A-matrix
#'
#' Constructs the A-matrix for a given component and and some data
#'
#' @aliases amatrix_eval.component
#' @export
#' @method amatrix_eval component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

amatrix_eval.component <- function(component, data, ...) {
  ## TODO: Check for offset component
  
  A.main <- amatrix_eval(component$main, data, component$env)
  A.group <- amatrix_eval(component$group, data, component$env)
  A.repl <- amatrix_eval(component$replicate, data, component$env)
  
  A <- INLA::inla.row.kron(A.repl,
                           INLA::inla.row.kron(A.group, A.main))
  
  # Mask columns of A
  if (!is.null(component$A.msk)) {
    A <- A[, as.logical(component$A.msk), drop = FALSE]
  }

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
#' @aliases input_eval.component
#' @export
#' @method input_eval component
#' @keywords internal
#' @param component An component.
#' @param part The input part to evaluate; \code{'main'}, \code{'group'},
#'   \code{'replicate'}, \code{'weights'}
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the component's map.
#' @param ... Unused.
#' @return An vector or a coordinate matrix
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

input_eval.component <- function(component,
                                 part = c("main", "group", "replicate", "weights"),
                                 data,
                                 ...) {
  part <- match.arg(part)
  if (part %in% "weights") {
    val <- input_eval(component[["main"]]$weights, data, component$env,
                      label = paste(component$label, "(main, weights)"), ...)
  } else {
    val <- input_eval(component[[part]]$input, data, component$env,
                      label = paste0(component$label, " (", part, ")"), ...)
  }
  val
}

#' @export

input_eval.component_list <-
  function(components,
           part = c("main", "group", "replicate", "weights"),
           data,
           ...)
  {
    part <- match.arg(part)
    result <- lapply(components, function(x) input_eval(x, part = part, data = data, ...))
    names(result) <- names(components)
    result
  }





input_eval.bru_input <- function(input, data, env = NULL, label = NULL, ...) {

  # Evaluate the map with the data as an environment
  if (is.null(env)) {
    emap <- tryCatch(eval(input$input, envir = data.frame(data), enclos = parent.frame()),
                     error = function(e) {})
  } else {
    emap <- tryCatch(eval(input$input, envir = data.frame(data), enclos = env),
                     error = function(e) {})
  }
    
  # 0) Eval failed. map everything to 1. This happens for automatically
  #    added Intercept, and for some components that cannot be evaluated
  #    based on the data, e.g. missing columns for a certain likelihood in a
  #    multilikelihood model.
  # 1) If we obtain a function, apply the function to the data
  # 2) If we obtain a SpatialGridDataFrame extract the values of that data
  #    frame at the point locations using the over() function
  # 3) Else we obtain a vector and return as-is. This happens when input
  #    references a column of the data points
  
  n <- nrow(as.data.frame(data))

  if (is.null(emap)) {
    val <- rep(1, n)
  } else if (is.function(emap)) {
    # Allow but detect failures:
    val <- tryCatch(emap(data),
                    error = function(e) {})
    if (is.null(val)) {
      # TODO: figure out if we need to do something else in this case.
      # Returning NULL seems like a good way to keep track of these failures
      # that are likely to happen for multilikelihood models; A component only
      # needs to be evaluable for at least one of the likelihoods.
    } else if (identical(as.character(input$input), "coordinates")) {
      # Return SpatialPoints instead of a matrix
      val <- as.data.frame(val)
      coordinates(val) <- seq_len(ncol(val))
      # Allow proj4string failures:
      p4s <- tryCatch(proj4string(data),
                      error = function(e) {})
      if (!is.null(p4s)) {
        proj4string(val) <- CRS(p4s)
      }
    }
  }
  else if (inherits(emap, "SpatialGridDataFrame") |
           inherits(emap, "SpatialPixelsDataFrame")) {
    if (is.null(input[["selector"]])) {
      layer <-
        if (is.null(input[["layer"]])) {
          1
        } else {
          input[["layer"]]
        }
      val <- over(data, emap)[, layer, drop = TRUE]
    } else {
      layer <- data[[input$selector]]
      val <- numeric(n)
      for (l in unique(layer)) {
        val[layer == l] <- over(data[layer == l, , drop = FALSE],
                                emap)[, l, drop = TRUE]
      }
    }
  } else if ((input$label == "offset") &&
           is.numeric(emap) &&
           (length(emap) == 1)) {
    val <- rep(emap, n)
  } else {
    val <- emap
  }

  # Always return as many rows as data has
  if (is.vector(val) && (length(val) == 1) && (n > 1)) {
    val <- rep(val, n)
  }

  # Check if any of the locations are NA. If we are dealing with SpatialGridDataFrame try
  # to fix that by filling in nearest neighbor values.
  # # TODO: Check how to deal with this fully in the case of multilikelihood models
  if (any(is.na(as.data.frame(val))) &&
      (inherits(emap, "SpatialGridDataFrame") ||
       inherits(emap, "SpatialPixelsDataFrame"))) {
    warning(sprintf(
      "Map '%s' has returned NA values. As it is a SpatialGridDataFrame I will try to fix this by filling in values of spatially nearest neighbors. In the future, please design your 'map=' argument as to return non-NA for all points in your model domain/mesh. Note that this can also significantly increase time needed for inference/prediction!",
      deparse(input$input), label
    ))

    if (is.null(ncol(val))) {
      ok <- !is.na(val)
    } else {
      nok <- is.na(val)
      ok <- (rowSums(nok) == 0)
    }
    dst <- rgeos::gDistance(SpatialPoints(data[ok, , drop = FALSE]),
                            SpatialPoints(data[!ok, , drop = FALSE]),
                            byid = TRUE)
    nn <- apply(dst, MARGIN = 1, function(row) which.min(row)[[1]])
    if (is.null(ncol(val))) {
      val[!ok] <- val[ok][nn]
    } else {
      val[nok] <- val[ok, , drop = FALSE][nn, , drop = FALSE][nok[!ok,]]
    }
  }

  # Check for NA values.
  if (any(is.na(as.data.frame(val)))) {
    # TODO: remove this check and make sure NAs are handled properly elsewhere
    stop(sprintf(
      "Input '%s' of component '%s' has returned NA values. Please design your 
                 argument as to return non-NA for all points in your model domain/mesh.",
      as.character(input$input)[[1]], label
    ))
  }

  val
}




#' Obtain indices
#'
#' Idexes into to the components
#'
#' @aliases index_eval.component
#' @export
#' @method index_eval component
#' @keywords internal
#' @param component An component.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the component's map.
#' @param ... Unused.
#' @return a data.frame of indices or list of indices into the components latent variables
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>,
#'   Finn Lindgren \email{finn.lindgren@@gmail.com}
#'

index_eval.component <- function(component, ...) {
  g.n <-
    if (is.null(component$group$n)) {
      1
    } else {
      component$group$n
    }
  r.n <-
    if (is.null(component$replicate$n)) {
      1
    } else {
      component$replicate$n
    }
  main <- component$main$values
  group <- seq_len(g.n)
  repl <- seq_len(r.n)
  idx <- data.frame(
    main = rep(main, times = g.n * r.n),
    group = rep(rep(group, each = length(main)),
                times = r.n),
    repl = rep(repl, each = length(main) * g.n)
  )
  colnames(idx) <- paste0(component$label, c("", ".group", ".repl"))
  idx
}

