# GENERICS ----

#' Construct A-matrix
#'
#' Constructs the A-matrix for components and data
#' @export
#' @rdname amatrix_eval
amatrix_eval <- function(...) {
  UseMethod("amatrix_eval")
}
#' Obtain covariate values
#'
#' @export
#' @rdname input_eval
input_eval <- function(...) {
  UseMethod("input_eval")
}
#' Obtain indices
#'
#' Indexes into to the components
#'
#' @export
#' @rdname index_eval
index_eval <- function(...) {
  UseMethod("index_eval")
}
#' @rdname add_mappers.component
add_mappers <- function(...) {
  UseMethod("add_mappers")
}

#' Methods for `bru_mapper` objects
#'
#' @details
#' * `bru_mapper` Generic mapper S3 constructor. Default constructor only
#' adds `bru_mapper` to the class definition, unless the `mapper` input already
#' inherits from `bru_mapper`.
#' @param \dots Arguments passed on to other methods
#' @param mapper A mapper S3 object, normally inheriting from `bru_mapper`
#' @export
#' @seealso [bru_mapper_methods] for specific method implementations.
#' @rdname bru_mapper
bru_mapper <- function(...) {
  UseMethod("bru_mapper")
}
#' @details
#' * `ibm_n` Generic.
#' @export
#' @rdname bru_mapper
ibm_n <- function(mapper, ...) {
  UseMethod("ibm_n")
}
#' @details
#' * `ibm_values` Generic.
#' @export
#' @rdname bru_mapper
ibm_values <- function(mapper, ...) {
  UseMethod("ibm_values")
}
#' @details
#' * `ibm_amatrix` Generic. Implementations must return a (sparse) matrix of size `NROW(input)`
#' by `ibm_n(mapper)`
#' @param input The values for which to produce a mapping matrix
#' @export
#' @rdname bru_mapper
ibm_amatrix <- function(mapper, input, ...) {
  UseMethod("ibm_amatrix")
}



# CONSTRUCTORS ----

#' inlabru latent model component construction
#'
#' @description
#'
#' Similar to `glm()`, `gam()` and `inla()`, [bru] models can be constructed via
#' a formula-like syntax, where each latent effect is specified. However, in
#' addition to the parts of the syntax compatible with `INLA::inla`, `bru`
#' components offer additional functionality which facilitates modelling, and
#' the predictor expression can be specified separately, allowing more complex
#' and non-linear predictors to be defined. The formula syntax is just a way to
#' allow all model components to be defined in a single line of code, but the
#' definitions can optionally be split up into separate component definitions.
#' See Details for more information.
#'
#' The `component` methods all rely on the [component.character()] method, that
#' defines a model component with a given label/name. This method is normally
#' the only method that user code might directly access, since the
#' `component.formula` and other methods are called inside [bru()] instead.
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{Finn.Lindgren@@gmail.com}
#' @family component constructors
#' @param object A formula (for `component.formula`), a character label
#' (for `component.character`), or a list of component objects
#' (for `component.list`; not yet supported as `bru` input (TODO)).
#' @param \dots Arguments passed on to the methods
#' @export

component <- function(object, ...) {
  UseMethod("component")
}


#' @family component constructors
#' @export
#' @param main = NULL
#' main takes an R expression that evaluates to where the latent variables should be evaluated (coordinates, indices, continuous scalar (for rw2 etc)).
#'   Arguments starting with weights, group, replicate behave similarly to main, but for the corresponding features of `INLA::f()`.
#' @param model Either one of "offset", "factor_full", "factor_contrast", "linear" or a model name or
#' object accepted by INLA's `f` function. If set to NULL, then "linear" is used.
#' @param mapper
#' Information about how to do the mapping from the values evaluated in `main`,
#' and to the latent variables. Auto-detects spde model objects in model and
#' extracts the mesh object to use as the mapper, and auto-generates mappers
#' for indexed models. (Default: NULL, for auto-determination)
#' @param main_layer,main_selector
#' The `_layer` is a numeric index or character name of which
#' layer/variable to extract from a covariate data object given in `main`
#' (Default: The effect component name, if it exists i the covariate object,
#' otherwise the first column of the covariate data frame)
#'
#' The `_selector` is character name of a variable whose contents determines which layer to
#' extract from a covariate for each data point. Overrides the `layer`.
#' (Default: NULL)
#' @param n The number of latent variables in the model. Should be auto-detected
#' for most or all models (Default: NULL, for auto-detection).
#' An error is given if it can't figure it out by itself.
#' @param values Specifies for what covariate/index values INLA should build
#' the latent model. Normally generated internally based on the mapping details.
#' (Default: NULL, for auto-determination)
#' @param season.length Passed on to `INLA::f()` for model `"seasonal"`
#' (TODO: check if this parameter is still fully handled)
# Weights
#' @param weights,weights_layer,weights_selector
#' Optional specification of effect scaling weights.
#' Same syntax as for `main`.
# Group model parameters
#' @param group,group_mapper,group_layer,group_selector
#' Optional specification of kronecker/group model indexing.
#' @param control.group `list` of kronecker/group model parameters, currently
#' passed directly on to `INLA:f`
# Replicate model parameters
#' @param replicate,replicate_mapper,replicate_layer,replicate_selector
#' Optional specification of indices for an independent
#' replication model. Same syntax as for `main`
#' @param A.msk TODO: check/fix/deprecate this parameter.
#' Likely doesn't work at the moment, and I've found no examples that use it.
# Deprecated parameters
#' @param map Deprecated. Use `main` instead.
#' @param mesh Deprecated. Use `mapper` instead.
#' @param envir_extra TODO: check/fix this parameter.
#'
#' @details The `component.character` method is inlabru's equivalent to INLA's
#' `f` function but adds functionality that is unique to inlabru.
#'
#' @family component constructors
#' @rdname component
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
#'     main = x,
#'     model = inla.spde2.matern(inla.mesh.1d(1:10))
#'   )
#' }
#' }
#'
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
                                weights_selector = NULL,
                                # Group model parameters
                                group = NULL,
                                group_mapper = NULL,
                                group_layer = NULL,
                                group_selector = NULL,
                                control.group = NULL,
                                # Replicate model parameters
                                replicate = NULL,
                                replicate_mapper = NULL,
                                replicate_layer = NULL,
                                replicate_selector = NULL,
                                A.msk = NULL,
                                # Deprecated parameters
                                # map -> main
                                map = NULL,
                                # mesh -> mapper
                                mesh = NULL,
                                envir_extra = NULL,
                                ...) {

  # INLA models:
  # itypes = c(linear, iid, mec, meb, rgeneric, rw1, rw2, crw2, seasonal, besag, besag2, bym, bym2, besagproper,
  #            besagproper2, fgn, fgn2, ar1, ar1c, ar, ou, generic, generic0, generic1, generic2, generic3, spde,
  #            spde2, spde3, iid1d, iid2d, iid3d, iid4d, iid5d, 2diid, z, rw2d, rw2diid, slm, matern2d, copy,
  #            clinear, sigm, revsigm, log1exp, logdist)
  #
  # Supported:
  # btypes = c("offset", "factor_full", "factor_contrast", "linear", "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")

  # The label
  label <- object

  if (is.null(model)) {
    # TODO: may need a special marker to autodetect factor variables,
    #       which can only be autodetected in a data-aware pass.
    model <- "linear"
  }

  # Force evaluation of explicit inputs
  force(values)

  if (!is.null(substitute(group))) {
    if (is.null(control.group)) {
      control.group <- INLA::inla.set.control.group.default()
    }
    group_model <- control.group$model
  } else {
    group_model <- "exchangeable"
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

  envir <- parent.frame()
  if (is.null(envir_extra)) {
    envir_extra <- new.env(envir)
  }

  # Default component (to be filled)
  component <- list(
    label = label,
    inla.formula = NULL,
    main = bru_subcomponent(
      input = bru_input(substitute(main),
        label = label,
        layer = main_layer,
        selector = main_selector
      ),
      mapper = mapper,
      model = model,
      n = n,
      values = values,
      season.length = season.length,
      weights =
        if (is.null(substitute(weights))) {
          NULL
        } else {
          bru_input(substitute(weights),
            label = paste0(label, ".weights"),
            layer = weights_layer,
            selector = weights_selector
          )
        }
    ),
    group = bru_subcomponent(
      input = bru_input(substitute(group),
        label = paste0(label, ".group"),
        layer = group_layer,
        selector = group_selector
      ),
      mapper = group_mapper,
      n = NULL,
      model = group_model
    ),
    replicate = bru_subcomponent(
      input = bru_input(substitute(replicate),
        label = paste0(label, ".repl"),
        layer = replicate_layer,
        selector = replicate_selector
      ),
      mapper = replicate_mapper,
      n = NULL,
      model = "iid"
    ),
    A.msk = A.msk,
    env = envir,
    env_extra = envir_extra
  )

  # Main bit
  if (component$main$type %in% c("offset")) {
    component$inla.formula <- as.formula(paste0("~ . + offset(offset)"),
      env = envir
    )
  } else {
    # Construct a call to the f function from the parameters provided
    # Ultimately, this call will be converted to the model formula presented to INLA
    fcall <- sys.call()
    fcall[[1]] <- "f"
    fcall[[2]] <- as.symbol(label)
    names(fcall)[2] <- ""
    if (is.null(names(fcall)) || (names(fcall)[3] == "")) {
      # 'main' is the only parameter allowed to be nameless,
      # and only if it's the first parameter (position 3 in this fcall)
      names(fcall)[3] <- "main"
    }

    # Store model name or object in the environment
    model_name <- paste0("BRU_", label, "_main_model")
    fcall[["model"]] <- as.symbol(model_name)
    assign(model_name, component$main$model, envir = component$env_extra)

    # Remove parameters inlabru supports but INLA doesn't,
    # and substitute parameters that inlabru will transform
    fcall <- fcall[!(names(fcall) %in% c(
      "main",
      "mapper",
      "group_mapper",
      "replicate_mapper",
      "A.msk",
      "map",
      "mesh",
      "main_layer",
      "group_layer",
      "replicate_layer",
      "weights_layer",
      "main_selector",
      "group_selector",
      "replicate_selector",
      "weights_selector"
    ))]

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
    suffixes <- list(
      "group" = "group",
      "replicate" = "repl",
      "weights" = "weights"
    )
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

    # Make sure 'values' is setup properly.
    if (is.null(component$main$values)) {
      fcall <- fcall[!("values" %in% names(fcall))]
    } else {
      values_name <- paste0("BRU_", label, "_values")
      fcall[["values"]] <- as.symbol(values_name)
      assign(values_name, component$main$values, envir = component$env_extra)
    }

    # Setup factor precision parameter
    if (identical(component$main$type, "factor")) {
      if (is.null(fcall[["hyper"]])) {
        # TODO: allow configuration of the precision via prec_linear
        factor_hyper_name <- paste0("BRU_", label, "_main_factor_hyper")
        fcall[["hyper"]] <- as.symbol(factor_hyper_name)
        assign(factor_hyper_name,
          list(prec = list(initial = log(1e-6), fixed = TRUE)),
          envir = component$env_extra
        )
      }
    }

    component$inla.formula <-
      as.formula(paste0(
        "~ . + ",
        as.character(parse(text = deparse(fcall)))
      ),
      env = envir
      )

    component$fcall <- fcall
  }

  class(component) <- c("component", "list")
  component
}





#' inlabru latent model component construction using formulae
#'
#' @description
#'
#' Similar to glm(), gam() and inla() [bru] uses formula objects to describe response data and latent
#' (unknonw) components of the model to be fitted. However, in addition to the syntax compatible with
#' [inla][INLA::inla], bru components offer addtitional functionality which facilitates modeling.
#'
#' @details
#'
#' [bru] will understand formulae describing fixed effect models just like the other methods. For instance, the
#' formula `y ~ x` will fit the linear combination of an effect named `x` and an intercept to
#' the response `y` with respect to the likelihood family stated when calling [bru]. Mathematically,
#' the linear predictor \eqn{\eta} would be written down as
#'
#' \deqn{\eta = \beta * x + c,}
#'
#' where:
#'
#' \itemize{
#' \item{\eqn{c} }{is the *intercept*}
#' \item{\eqn{x }}{is a *covariate*}
#' \item{\eqn{\beta} }{is a *random variable* associated with \eqn{x} and}
#' \item{\eqn{\psi = \beta * x }}{ is called the *random effect* of \eqn{x}}
#' }
#'
#' A problem that arises when using this kind of R formula is that it does not clearly relect the mathematical
#' formula. For instance, when providing the formula to inla, the resulting object will refer to the random
#' effect \eqn{\psi = \beta * x } as `x`. Hence, it is not clear if `x` refers to the covariate
#' or the effect of the covariate.
#'
#' @section Naming random effects:
#'
#' In INLA, a simple random effect model would be expressed as
#'
#' \itemize{\item{`formula = y ~ f(x, model = "linear")`,}}
#'
#' where [f][INLA::f] is the inla specific function to set up random effects of all kinds. The underlying
#' predictor would again be \eqn{\eta = \beta * x + c} but the result of fitting the model would state
#' `x` as the random effect's name. bru allows to rewrite this formula in order to explicitly state
#' the name of the random effect and the name of the associated. This is achived by replacing `f`
#' with an arbitrary name that we wish to assign to the effect, e.g.
#'
#' \itemize{\item{`components = y ~ psi(x, model = "linear")`.}}
#'
#' Being able to discriminate between \eqn{x} and \eqn{\psi} is relevant because of two functionalities
#' bru offers. The formula parameters of both, [bru] and the prediction method [predict.bru]
#' are interpreted in the mathematical sense. For instance, `predict` may be used to analyze the
#' an analytical combination of the covariate \eqn{x} and the intercept using
#'
#' \itemize{\item{`predict(fit, data.frame(x=2)), ~ exp(psi + Intercept)`.}}
#'
#' which corresponds to the mathematical expression \eqn{\exp(x \beta + c)}.
#'
#' On the other hand, predict may be used to only look at a transformation of
#' the latent variable \eqn{\beta_\psi}
#'
#' \itemize{\item{`predict(fit, NULL, ~ exp(psi_latent))`.}}
#'
#' #' which corresponds to the mathematical expression \eqn{\exp(\beta)}.
#'
#' @param lhoods TODO: document
#' @param envir TODO: document
#' @export
#' @family component constructors
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname component
#'
#' @examples
#' # As an example, let us create a linear component. Here, the component is
#' # called "myLinearEffectOfX" while the covariate the component acts on is
#' # called "x". Note that a list of components is returned because the
#' # formula may define multiple components
#'
#' eff <- component(~ myLinearEffectOfX(main = x, model = "linear"))
#' summary(eff[[1]])
#' # Equivalent shortcuts:
#' eff <- component(~ myLinearEffectOfX(x, model = "linear"))
#' eff <- component(~ myLinearEffectOfX(x))
component.formula <- function(object, lhoods = NULL, envir = NULL, ...) {
  if (is.null(envir)) {
    envir <- environment(object)
  }
  code <- code.components(object)
  parsed <- lapply(code, function(x) parse(text = x))
  components <- lapply(
    parsed,
    function(component.expression) {
      eval(component.expression,
        envir = envir
      )
    }
  )
  environment(object) <- envir
  component.list(components, lhoods = lhoods, envir = envir)
}

#' @export
#' @rdname component
component.list <- function(object, lhoods = NULL, envir = NULL, ...) {
  stopifnot(all(vapply(object, function(x) inherits(x, "component"), TRUE)))
  if (is.null(envir)) {
    envir <- environment(object)
  }
  names(object) <- lapply(object, function(x) x$label)
  class(object) <- c("component_list", "list")
  environment(object) <- envir
  if (!is.null(lhoods)) {
    object <- add_mappers(object, lhoods = lhoods)
  }
  object
}


#' @export
#' @rdname component
`[.component_list` <- function(x, i) {
  env <- environment(x)
  object <- NextMethod()
  class(object) <- c("component_list", "list")
  environment(object) <- env
  object
}



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param component PARAM_DESCRIPTION
#' @param lhoods PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname add_mappers.component
#' @keywords internal
#' @export

add_mappers.component <- function(component, lhoods) {
  # Filter out lhoods that don't use/support the component
  keep_lh <-
    vapply(lhoods,
      function(lh, label) {
        label %in% parse_inclusion(
          label,
          lh[["include_components"]],
          lh[["exclude_components"]]
        )
      },
      TRUE,
      label = component$label
    )
  lh <- lhoods[keep_lh]
  component$label

  component$main <- add_mapper(component$main,
    label = component$label,
    lhoods = lh,
    env = component$env
  )
  component$group <- add_mapper(component$group,
    label = component$label,
    lhoods = lh,
    env = component$env
  )
  component$replicate <- add_mapper(component$replicate,
    label = component$label,
    lhoods = lh,
    env = component$env
  )

  if (component$main$type %in% c("offset")) {
  } else {
    fcall <- component$fcall

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

    if (is.null(component$main$values)) {
      fcall <- fcall[!("values" %in% names(fcall))]
    } else {
      values_name <- paste0("BRU_", component$label, "_values")
      fcall[["values"]] <- as.symbol(values_name)
      assign(values_name, component$main$values, envir = component$env_extra)
    }

    # Generate the formula that will be presented to INLA
    component$inla.formula <-
      as.formula(paste0(
        "~ . + ",
        paste0(deparse(fcall),
          collapse = "\n"
        )
      ),
      env = component$env
      )

    component$fcall <- fcall
  }

  component
}
add_mappers.component_list <- function(components, lhoods) {
  for (k in seq_along(components)) {
    components[[k]] <- add_mappers(components[[k]], lhoods)
  }
  components
}

bru_input <- function(input, label = NULL, layer = NULL, selector = NULL) {
  inp <- list(
    input = input,
    label = label,
    layer = layer,
    selector = selector
  )
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
  factor_mapping <- NULL
  if (inherits(model, "inla.spde")) {
    type <- "spde"
  } else if (inherits(model, "inla.rgeneric")) {
    type <- "rgeneric"
  } else if (inherits(
    model,
    c(
      "clinear", "sigm", "revsigm",
      "log1exp", "logdist"
    )
  )) {
    type <- "specialnonlinear"
  } else if (is.character(model)) {
    if (identical(model, "factor")) {
      model <- "factor_contrast"
      warning(paste0(
        "Deprecated model 'factor'. Please use 'factor_full' or ",
        "'factor_contrast' instead.\n",
        "Defaulting to 'factor_contrast' that matches the old 'factor' model."
      ),
      immediate. = TRUE
      )
    }
    if (identical(model, "factor_full")) {
      model <- "iid"
      type <- "factor"
      factor_mapping <- "full"
    } else if (identical(model, "factor_contrast")) {
      model <- "iid"
      type <- "factor"
      factor_mapping <- "contrast"
    } else {
      type <- model
    }
  } else {
    type <- "unknown"
  }
  subcomponent <-
    list(
      input = input,
      mapper = mapper,
      model = model,
      type = type,
      n = n,
      values = values,
      season.length = season.length,
      weights = weights,
      factor_mapping = factor_mapping
    )
  class(subcomponent) <- c("bru_subcomponent", "list")

  subcomponent
}

add_mapper <- function(subcomp, label, lhoods = NULL, env = NULL) {
  if (!is.null(subcomp[["mapper"]])) {
    if (!inherits(
      subcomp[["mapper"]],
      c("bru_mapper", "inla.mesh", "inla.mesh.1d")
    )) {
      stop(paste0(
        "Unknown mapper of type '",
        paste0(class(subcomp[["mapper"]]), collapse = ", "),
        "' for ", label
      ))
    }
    subcomp[["mapper"]] <- bru_mapper(subcomp[["mapper"]])
  } else {
    if (!is.null(lhoods)) {
      inp <- lapply(
        lhoods,
        function(lh) {
          input_eval(subcomp$input,
            data = lh$data,
            env = env,
            label = subcomp$input$label,
            null.on.fail = TRUE
          )
        }
      )
      # Check for
      # 1) All NULL; OK, set to 1
      # 2) Some NULL; exclude NULL results
      # TODO: Check for vector/matrix/coordinate inconsistency
      null.results <- vapply(inp, function(x) is.null(x), TRUE)
      if (all(null.results)) {
        inp_values <- 1
        n_values <- 1
      } else {
        if (any(null.results)) {
          inp_ <- inp[!null.results]
        } else {
          inp_ <- inp
        }
        is_spatial <- vapply(inp_, function(x) inherits(x, "Spatial"), TRUE)
        is_matrix <- vapply(inp_, function(x) is.matrix(x), TRUE)
        is_factor <- vapply(inp_, function(x) is.factor(x), TRUE)
        if (any(is_spatial)) {
          if (!all(is_spatial)) {
            stop("Inconsistent input types; spatial and non-spatial")
          }
          message("TODO: Ensure spatial input objects use the same CRS")
          inp_values <- unique(do.call(
            rbind,
            lapply(
              inp_,
              function(x) coordinates(x)
            )
          ))
          n_values <- nrow(inp_values)
        } else if (any(is_matrix)) {
          if (!all(is_matrix)) {
            stop("Inconsistent input types; matrix and non-matrix")
          }
          inp_values <- unique(do.call(rbind, inp_))
          n_values <- nrow(inp_values)
        } else if (any(is_factor)) {
          if (!all(is_factor)) {
            stop("Inconsistent input types; factor and non-factor")
          }
          inp_values <- sort(unique(do.call(c, lapply(inp_, as.character))),
            na.last = NA
          )
          n_values <- length(inp_values)
        } else {
          inp_values <- sort(unique(unlist(inp_)), na.last = NA)
          n_values <- length(inp_values)
        }
      }
      if (n_values < 1) {
        subcomp$n <- 1
        subcomp$values <- NULL
        inp_values <- NULL
      }
      subcomp <- make_mapper(subcomp,
        label,
        input_values = inp_values,
        strict = TRUE
      )
    }
  }
  if (!is.null(subcomp[["mapper"]])) {
    # Check internal consistency of user specified n and the mapper:
    mapper_n <- ibm_n(subcomp[["mapper"]])
    if (!is.null(subcomp[["n"]]) &&
        subcomp[["n"]] != mapper_n) {
      stop(paste0(
        "Size mismatch, n=", subcomp[["n"]], " != ibm_n()=",
        mapper_n, " mapper for label ", label
      ))
    }
    subcomp[["n"]] <- mapper_n
    subcomp[["values"]] <- ibm_values(subcomp[["mapper"]])
  }
  subcomp
}



# Defines a default mapper given the type of model and parameters provided
# Checks subcomp$mapper, subcomp$model$mesh (for spde models), subcomp$n,
# subcomp$values, or input_values, in that order.
make_mapper <- function(subcomp,
                        label,
                        input_values = NULL,
                        strict = TRUE) {
  if (is.null(subcomp[["mapper"]])) {
    if (subcomp[["type"]] %in% c("spde")) {
      subcomp[["mapper"]] <- subcomp[["model"]]$mesh
    }
  }
  if (!is.null(subcomp[["mapper"]])) {
    if (!inherits(
      subcomp[["mapper"]],
      c("bru_mapper", "inla.mesh", "inla.mesh.1d")
    )) {
      stop(paste0(
        "Unknown mapper of type '",
        paste0(class(subcomp[["mapper"]]), collapse = ", "),
        "' for ", label
      ))
    }
    subcomp[["mapper"]] <- bru_mapper(subcomp[["mapper"]])
  } else if (subcomp[["type"]] %in% c("linear", "clinear")) {
    subcomp[["mapper"]] <- bru_mapper_linear()
  } else {
    # No mapper; construct based on input values
    if (!is.null(subcomp[["n"]])) {
      if (!is.null(subcomp[["values"]])) {
        warning(
          "Both 'n' and 'values' provided. Ignoring 'values'",
          immediate. = TRUE
        )
      }
      values <- seq_len(subcomp[["n"]])
    } else if (!is.null(subcomp[["values"]])) {
      values <- subcomp[["values"]]
    } else if (!is.null(input_values)) {
      values <- input_values
    } else {
      stop(paste0("No mapper, no n, and no values given for ", label))
    }

    if (is.factor(values) || is.character(values)) {
      subcomp[["mapper"]] <- bru_mapper_factor(
        values,
        factor_mapping = subcomp[["factor_mapping"]]
      )
    } else {
      values <- sort(unique(values), na.last = NA)
      if (length(values) > 1) {
        subcomp[["mapper"]] <- bru_mapper(INLA::inla.mesh.1d(values))
      } else {
        if (all(values == 1)) {
          subcomp[["mapper"]] <- bru_mapper_const()
        } else {
          subcomp[["mapper"]] <- bru_mapper_linear()
        }
      }
    }
  }
  subcomp
}

#' Convert components to R code
#'
#' @aliases code.components
#' @keywords internal
#' @param components A [formula] describing latent model components.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'

code.components <- function(components, add = "") {
  fname <- "inlabru:::component.character"
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
      codes[[k]] <- paste0(
        fname, '("', label, '", main = ', label,
        ', model = "linear"',
        add, ")"
      )
    }
    else {
      # Add extra code before final bracket
      ix <- max(gregexpr(")", text = code, fixed = TRUE)[[1]])
      code <-
        paste0(
          substr(code, 1, ix - 1),
          add,
          substr(code, ix, nchar(code))
        )

      if (is.offset) {
        codes[[k]] <-
          gsub(paste0(label, "("),
            paste0(fname, "(\"", label, "\", model = \"offset\", main = "),
            code,
            fixed = TRUE
          )
      }
      else {
        codes[[k]] <- gsub(paste0(label, "("),
          paste0(fname, "(\"", label, "\", "),
          code,
          fixed = TRUE
        )
      }
    }
  }

  codes
}

# MAPPERS ----

#' @export
#' @rdname bru_mapper
bru_mapper.default <- function(mapper, ...) {
  if (!inherits(mapper, "bru_mapper")) {
    class(mapper) <- c("bru_mapper", class(mapper))
  }
  mapper
}

#' Implementation methods for mapper objects
#' @param \dots Arguments passed on to other methods
#' @param mapper A mapper S3 object, normally inheriting from `bru_mapper`
#' @seealso [bru_mapper()]
#' @name bru_mapper_methods
#' @rdname bru_mapper_methods
ibm_n.inla.mesh <- function(mapper, ...) {
  mapper$n
}
#' @rdname bru_mapper_methods
ibm_values.inla.mesh <- function(mapper, ...) {
  seq_len(mapper$n)
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_mapper_methods
ibm_amatrix.inla.mesh <- function(mapper, input, ...) {
  if (!is.matrix(input) && !inherits(input, "Spatial")) {
    val <- as.matrix(val)
  }
  INLA::inla.spde.make.A(mapper, loc = input)
}

#' @rdname bru_mapper_methods
ibm_n.inla.mesh.1d <- function(mapper, ...) {
  mapper$m
}
#' @rdname bru_mapper_methods
ibm_values.inla.mesh.1d <- function(mapper, ...) {
  mapper$loc
}
#' @rdname bru_mapper_methods
ibm_amatrix.inla.mesh.1d <- function(mapper, input, ...) {
  INLA::inla.spde.make.A(mapper, loc = input)
}

#' @export
#' @rdname bru_mapper
bru_mapper_linear <- function(...) {
  mapper <- list()
  class(mapper) <- c("bru_mapper_linear", "bru_mapper", "list")
  mapper
}

#' @rdname bru_mapper_methods
ibm_n.bru_mapper_linear <- function(mapper, ...) {
  1L
}
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_linear <- function(mapper, ...) {
  1.0
}

#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_linear <- function(mapper, input, ...) {
  ok <- !is.na(input)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = rep(1, sum(ok)),
    x = input[ok],
    dims = c(NROW(input), ibm_n(mapper))
  )
  A
}

#' @export
#' @rdname bru_mapper
bru_mapper_const <- function(...) {
  mapper <- list()
  class(mapper) <- c("bru_mapper_linear", "bru_mapper", "list")
  mapper
}

#' @rdname bru_mapper_methods
ibm_n.bru_mapper_const <- function(mapper, ...) {
  1L
}
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_const <- function(mapper, ...) {
  1.0
}

#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_const <- function(mapper, input, ...) {
  ok <- !is.na(input)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = rep(1, sum(ok)),
    x = rep(1, sum(ok)),
    dims = c(NROW(input), ibm_n(mapper))
  )
  A
}

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param factor_mapping character; selects the type of factor mapping.
#' * `'contrast'` for leaving out the first factor level.
#' * `'full'` for keeping all levels.
#' @export
#' @rdname bru_mapper
bru_mapper_factor <- function(values, factor_mapping, ...) {
  stopifnot(
    is.factor(values) || is.character(values)
  )
  factor_mapping <- match.arg(factor_mapping, c("full", "contrast"))
  if (is.factor(values)) {
    mapper <- list(
      levels = levels(values),
      factor_mapping = factor_mapping
    )
  } else {
    mapper <- list(
      levels = unique(values),
      factor_mapping = factor_mapping
    )
  }
  class(mapper) <- c("bru_mapper_factor", "bru_mapper", "list")
  mapper
}

#' @rdname bru_mapper_methods
ibm_n.bru_mapper_factor <- function(mapper, ...) {
  length(ibm_values(mapper))
}
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_factor <- function(mapper, ...) {
  if (identical(mapper$factor_mapping, "contrast")) {
    mapper$levels[-1L]
  } else {
    mapper$levels
  }
}

#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_factor <- function(mapper, input, ...) {
  if (is.factor(input)) {
    if (!identical(levels(input), mapper$levels)) {
      input <- factor(as.character(input), levels = mapper$levels)
    }
  } else {
    input <- factor(input, levels = mapper$levels)
  }
  if (mapper$factor_mapping == "full") {
    input <- as.numeric(input)
  } else if (mapper$factor_mapping == "contrast") {
    input <- as.numeric(input) - 1L
  } else {
    stop("Unknown factor mapping '", mapper$factor_mapping, "'.")
  }
  ok <- !is.na(input)
  ok[which(ok)] <- (input[ok] > 0L)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = input[ok],
    x = rep(1.0, sum(ok)),
    dims = c(NROW(input), ibm_n(mapper))
  )
  A
}

# OPERATORS ----


#' Summarize an component
#'
#' @export
#' @method summary component
#' @keywords internal
#' @param object An component.
#' @param ... ignored.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
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
#' @export
#' @keywords internal
#' @param x A 'summary.component' object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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





#' @export
#' @keywords internal
#' @param component An component.
#' @param data A `data.frame` or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname amatrix_eval

amatrix_eval.bru_subcomponent <- function(subcomp, data, env = NULL, ...) {
  ## TODO: Check for offset component

  val <- input_eval(subcomp$input, data, env = env)

  if (!is.null(subcomp$weights)) {
    weights <- input_eval(subcomp$weights, data, env = env)
  } else {
    weights <- 1.0
  }
  if (inherits(subcomp[["mapper"]], "bru_mapper")) {
    A <- weights * ibm_amatrix(subcomp[["mapper"]], input = val)
  } else {
    stop(paste0(
      "Unsupported mapper of class '",
      paste0(class(subcomp[["mapper"]]), collapse = "', '"), "'",
      " for subcomponent '", subcomp$label, "'"
    ))
  }

  A
}



#' @export
#' @keywords internal
#' @param component An component.
#' @param data A `data.frame` or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname amatrix_eval

amatrix_eval.component <- function(component, data, ...) {
  ## TODO: Check for offset component

  A.main <- amatrix_eval(component$main, data, component$env)
  A.group <- amatrix_eval(component$group, data, component$env)
  A.repl <- amatrix_eval(component$replicate, data, component$env)

  A <- INLA::inla.row.kron(
    A.repl,
    INLA::inla.row.kron(A.group, A.main)
  )

  # Mask columns of A
  # TODO: check what this feature is intended for!
  if (!is.null(component$A.msk)) {
    A[, as.logical(component$A.msk)] <- 0.0
  }

  A
}

#' @export
#' @rdname amatrix_eval

amatrix_eval.component_list <- function(components, data, ...) {
  lapply(components, function(x) amatrix_eval(x, data = data, ...))
}

#' @section Simple covariates and the map parameter:
#'
#' It is not unusual for a random effect act on a transformation of a covariate. In other frameworks this
#' would mean that the transformed covariate would have to be calculated in advance and added to the
#' data frame that is usually provided via the `data` parameter. inlabru provides the option to do
#' this transformation automatically. For instance, one might be interested in the effect of a covariate
#' \eqn{x^2}. In inla and other frameworks this would require to add a column `xsquared` to the
#' input data frame and use the formula
#'
#' \itemize{\item{`formula = y ~ f(xsquared, model = "linear")`,}}
#'
#' In inlabru this can be achived using two ways of using the `main` parameter
#' (`map` in version 2.1.13 and earlier).
#'
#' \itemize{
#' \item{`components = y ~ psi(main = x^2, model = "linear")`}
#' \item{`components = y ~ psi(main = mySquareFun(x), model = "linear")`,}
#' \item{`components = y ~ psi(main = myOtherSquareFun, model = "linear")`,}
#'
#' }
#'
#' In the first example inlabru will interpret the map parameter as an expression to be evaluated within
#' the data provided. Since \eqn{x} is a knonwn covariate it will know how to calculate it. The second
#' example is an expression as well but it uses a function alled `mySquareFun`. This function is
#' defined by user but has wo be accessible within the work space when setting up the compoonents.
#' The third example provides the function `myOtherSquareFun` directly and not within an expression.
#' In this case, inlabru will call the function using the data provided via the  `data` parameter.
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
#' convenient way in inlabru. Spatial covariates are often stored as `SpatialPixelDataFrame`,
#' `SpatialPixelDataFrame` or `RasterLayer` objects. These can be provided directly via
#' the map parameter if the input data is a `SpatialPointsDataFrame`. inlabru will automatically
#' evaluate and/or interpolate the coariate at your data locations when using code like
#'
#' \itemize{\item{`components = y ~ psi(mySpatialPixels, model = "linear")`.}}
#'
#' @section Coordinates:
#'
#' A common spatial modelling component when using inla are SPDE models. An important feature of
#' inlabru is that it will automatically calculate the so called A-matrix which maps SPDE
#' values at the mesh vertices to values at the data locations. For this purpose, the map parameter
#' can be se to `coordinates`, which is the `sp` package function that extracts point
#' coordinates from the SpatialPointsDataFrame that was provided as input to bru. The code for
#' this would look as follows:
#'
#' \itemize{\item{`components = y ~ mySPDE(main = coordinates, model = inla.spde2.matern(...))`.}}
#'
#' @export
#' @keywords internal
#' @param component An component.
#' @param part The input part to evaluate; `'main'`, `'group'`,
#'   `'replicate'`, `'weights'`
#' @param data A `data.frame` or Spatial* object of covariates and/or point locations. If null, return the component's map.
#' @param ... Unused.
#' @return An vector or a coordinate matrix
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname input_eval

input_eval.component <- function(component,
                                 part = c("main", "group", "replicate", "weights"),
                                 data,
                                 ...) {
  part <- match.arg(part)
  if (part %in% "weights") {
    val <- input_eval(component[["main"]]$weights, data, component$env,
      label = paste(component$label, "(main, weights)"), ...
    )
  } else {
    val <- input_eval(component[[part]]$input, data, component$env,
      label = paste0(component$label, " (", part, ")"), ...
    )
  }
  val
}

#' @export
#' @rdname input_eval

input_eval.component_list <-
  function(components,
           part = c("main", "group", "replicate", "weights"),
           data,
           ...) {
    part <- match.arg(part)
    result <- lapply(components, function(x) input_eval(x, part = part, data = data, ...))
    names(result) <- names(components)
    result
  }



#' @export
#' @rdname input_eval

input_eval.bru_input <- function(input, data, env = NULL, label = NULL,
                                 null.on.fail = FALSE, ...) {

  # Evaluate the map with the data as an environment
  if (is.null(env)) {
    emap <- tryCatch(eval(input$input, envir = data.frame(data), enclos = parent.frame()),
      error = function(e) {}
    )
  } else {
    emap <- tryCatch(eval(input$input, envir = data.frame(data), enclos = env),
      error = function(e) {}
    )
  }

  # 0) Eval failed. map everything to 1. This happens for automatically
  #    added Intercept, and for some components that cannot be evaluated
  #    based on the data, e.g. missing columns for a certain likelihood in a
  #    multilikelihood model.
  # 1) If we obtain a function, apply the function to the data
  # 2) If we obtain a SpatialGridDataFrame extract the values of that data
  #    frame at the point locations using the over() function
  # 3) Else we obtain a vector and return as-is. This happens when input
  #    references a column of the data points, or some other complete expression

  n <- nrow(as.data.frame(data))

  if (is.null(emap)) {
    if (null.on.fail) {
      return(NULL)
    }
    val <- rep(1, n)
  } else if (is.function(emap)) {
    # Allow but detect failures:
    val <- tryCatch(emap(data),
      error = function(e) {}
    )
    if (is.null(val)) {
      # TODO: figure out if we need to do something else in this case.
      # Returning NULL seems like a good way to keep track of these failures
      # that are likely to happen for multilikelihood models; A component only
      # needs to be evaluable for at least one of the likelihoods.
    } else if (identical(as.character(input$input), "coordinates")) {
      tryCatch(
        {
          # Return SpatialPoints instead of a matrix
          val <- as.data.frame(val)
          coordinates(val) <- seq_len(ncol(val))
          # Allow proj4string failures:
          data_crs <- tryCatch(fm_sp_get_crs(data),
                               error = function(e) {}
          )
          if (!fm_crs_is_null(data_crs)) {
            proj4string(val) <- data_crs
          }
          val
        },
        error = function(e) {
          NULL
        }
      )
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
    }
    val <- eval_SpatialDF(
      emap,
      data,
      layer = layer,
      selector = input[["selector"]]
    )
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
  # to fix that by filling in nearest neighbour values.
  # # TODO: Check how to deal with this fully in the case of multilikelihood models
  # Answer: should respect the lhood "include/exclude" info for the component list
  if (any(is.na(as.data.frame(val))) &&
    (inherits(emap, "SpatialGridDataFrame") ||
      inherits(emap, "SpatialPixelsDataFrame"))) {
    warning(
      paste0(
        "Model input '",
        deparse(input$input),
        "' for '", label, "' returned some NA values.\n",
        "Attempting to fill in spatially by nearest available value.\n",
        "To avoid this basic covariate imputation, supply complete data."
      )
    )
    
    val <- bru_fill_missing(data = emap, where = data, values = val,
                            layer = layer, selector = input[["selector"]])
  }

  # Check for NA values.
  if (any(is.na(as.data.frame(val)))) {
    # TODO: remove this check and make sure NAs are handled properly elsewhere,
    # if possible. Problem: treating NA as "no effect" can have negative side
    # effects. For spatial covariates, can be handled by infill, but a general
    # solution doesn't exist, so not sure how to deal with this.
    stop(sprintf(
      "Input '%s' of component '%s' has returned NA values. Please design your 
                 argument as to return non-NA for all points in your model domain/mesh.",
      as.character(input$input)[[1]], label
    ))
  }

  val
}




#' @export
#' @keywords internal
#' @param component An component.
#' @param ... Unused.
#' @return a list of indices into the components latent variables
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com},
#'   Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname index_eval

index_eval.component <- function(component, ...) {
  main <- component$main$values
  group <- component$group$values
  replicate <- component$replicate$values
  m.n <- length(main)
  g.n <- length(group)
  r.n <- length(replicate)
  idx <- list(
    main = rep(main, times = g.n * r.n),
    group = rep(rep(group, each = m.n), times = r.n),
    repl = rep(replicate, each = m.n * g.n)
  )
  names(idx) <- paste0(component$label, c("", ".group", ".repl"))
  idx
}


#' @export
#' @rdname index_eval

index_eval.component_list <- function(components, ...) {
  lapply(components, function(x) index_eval(x, ...))
}
