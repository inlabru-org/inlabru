# GENERICS ----

#' Construct component linearisations
#'
#' Constructs the linearisation mapper for each component
#' @export
#' @rdname comp_lin_eval
comp_lin_eval <- function(...) {
  UseMethod("comp_lin_eval")
}
#' Obtain component inputs
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
#' Obtain inla index subset information
#'
#' Subsets for `INLA::f()` compatible indexing
#'
#' @export
#' @rdname inla_subset_eval
inla_subset_eval <- function(...) {
  UseMethod("inla_subset_eval")
}
#' @title Add component input/latent mappers
#' @description Add missing mappers between input data and latent variables,
#' based on likelihood data
#' @param \dots Parameters passed on to other methods
#' @export
#' @rdname add_mappers
add_mappers <- function(...) {
  UseMethod("add_mappers")
}




# CONSTRUCTORS ----

#' Latent model component construction
#'
#' @description
#'
#' Similar to `glm()`, `gam()` and `inla()`, [bru()] models can be constructed via
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
#' defines a model component with a given label/name. The user usually
#' doesn't need to call these methods directly, but can instead supply a
#' formula expression that can be interpreted by the [component_list.formula()]
#' method, called inside [bru()].
#'
#' @details
#'
#' As shorthand, [bru()] will understand basic additive formulae describing fixed effect
#' models. For instance, the
#' components specification `y ~ x` will define the linear combination of an
#' effect named `x` and an intercept to
#' the response `y` with respect to the likelihood family stated when calling [bru()]. Mathematically,
#' the linear predictor \eqn{\eta} would be written down as
#'
#' \deqn{\eta = \beta * x + c,}
#'
#' where:
#'
#' \itemize{
#' \item{\eqn{c} }{is the *intercept*}
#' \item{\eqn{x }}{is a *covariate*}
#' \item{\eqn{\beta} }{is a *latent variable* associated with \eqn{x} and}
#' \item{\eqn{\psi = \beta * x }}{ is called the *effect* of \eqn{x}}
#' }
#'
#' A problem that arises when using this kind of R formula is that it does not
#' clearly reflect the mathematical
#' formula. For instance, when providing the formula to inla, the resulting
#' object will refer to the random
#' effect \eqn{\psi = \beta * x } as `x`.
#' Hence, it is not clear when `x` refers to the covariate
#' or the effect of the covariate.
#'
#' @section Naming random effects:
#'
#' In INLA, the `f()` notation is used to define more complex models, but
#' a simple linear effect model can also be expressed as
#'
#' \itemize{\item{`formula = y ~ f(x, model = "linear")`,}}
#'
#' where `f()` is the inla specific function to set up random effects of all kinds. The underlying
#' predictor would again be \eqn{\eta = \beta * x + c} but the result of fitting the model would state
#' `x` as the random effect's name. bru allows rewriting this formula in order to explicitly state
#' the name of the random effect and the name of the associated covariate. This is achieved by replacing `f`
#' with an arbitrary name that we wish to assign to the effect, e.g.
#'
#' \itemize{\item{`components = y ~ psi(x, model = "linear")`.}}
#'
#' Being able to discriminate between \eqn{x} and \eqn{\psi} is relevant because of two functionalities
#' bru offers. The formula parameters of both [bru()] and the prediction method [predict.bru]
#' are interpreted in the mathematical sense. For instance, `predict` may be used to analyze the
#' analytical combination of the covariate \eqn{x} and the intercept using
#'
#' \itemize{\item{`predict(fit, data.frame(x=2)), ~ exp(psi + Intercept)`.}}
#'
#' which corresponds to the mathematical expression \ifelse{html}{\out{e <sup>&#946; + c</sup>}}{\eqn{e^{x \beta + c}}}.
#'
#' On the other hand, predict may be used to only look at a transformation of
#' the latent variable \eqn{\beta_\psi}
#'
#' \itemize{\item{`predict(fit, NULL, ~ exp(psi_latent))`.}}
#'
#' which corresponds to the mathematical expression \ifelse{html}{\out{e <sup>&#946;</sup>}}{\eqn{e^{\beta}}}.
#'
#' @param \dots Parameters passed on to other methods
#'
#' @rdname component
#' @aliases bru_component
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{Finn.Lindgren@@gmail.com}
#' @family component constructors
#' @export
#'
#' @examples
#' # As an example, let us create a linear component. Here, the component is
#' # called "myLinearEffectOfX" while the covariate the component acts on is
#' # called "x". Note that a list of components is returned because the
#' # formula may define multiple components
#'
#' cmp <- component_list(~ myLinearEffectOfX(main = x, model = "linear"))
#' summary(cmp)
#' # Equivalent shortcuts:
#' cmp <- component_list(~ myLinearEffectOfX(x, model = "linear"))
#' cmp <- component_list(~ myLinearEffectOfX(x))
#' # Individual component
#' cmp <- component("myLinearEffectOfX", main = x, model = "linear")
#' summary(cmp)
component <- function(...) {
  UseMethod("component")
}


#' @export
#' @param object A character label for the component
#' @param main
#' `main` takes an R expression that evaluates to where the latent variables
#' should be evaluated (coordinates, indices, continuous scalar (for rw2 etc)).
#' Arguments starting with weights, group, replicate behave similarly to main,
#' but for the corresponding features of `INLA::f()`.
#' @param model Either one of "const" (same as "offset"), "factor_full",
#' "factor_contrast", "linear",
#' "fixed", or a model name or
#' object accepted by INLA's `f` function. If set to NULL, then "linear" is used
#' for vector inputs, and "fixed" for matrix input (converted internally to
#' an iid model with fixed precision)
#' @param mapper
#' Information about how to do the mapping from the values evaluated in `main`,
#' and to the latent variables. Auto-detects spde model objects in model and
#' extracts the mesh object to use as the mapper, and auto-generates mappers
#' for indexed models. (Default: NULL, for auto-determination)
#' @param main_layer,main_selector
#' The `_layer` input should evaluate to a numeric index or character name or
#' vector of which
#' layer/variable to extract from a covariate data object given in `main`.
#' (Default: NULL if `_selector` is given. Otherwise the effect component name,
#'  if it exists in the covariate object, and otherwise the first column of
#'  the covariate data frame)
#'
#' The `_selector` value should be a character name of a variable
#' whose contents determines which layer to extract from a covariate for each
#' data point. (Default: NULL)
#' @param n The number of latent variables in the model. Should be auto-detected
#' for most or all models (Default: NULL, for auto-detection).
#' An error is given if it can't figure it out by itself.
#' @param values Specifies for what covariate/index values INLA should build
#' the latent model. Normally generated internally based on the mapping details.
#' (Default: NULL, for auto-determination)
#' @param season.length Passed on to `INLA::f()` for model `"seasonal"`
#' (TODO: check if this parameter is still fully handled)
# Copy feature
#' @param copy character; label of other component that this component should
#' be a copy of. If the `fixed = FALSE`, a scaling constant is estimated, via a
#' hyperparameter. If `fixed = TRUE`, the component scaling is fixed, by
#' default to 1; for fixed scaling, it's more efficient to express the scaling
#' in the predictor expression instead of making a copy component.
# Weights
#' @param weights,weights_layer,weights_selector
#' Optional specification of effect scaling weights.
#' Same syntax as for `main`.
# Group model parameters
#' @param group,group_mapper,group_layer,group_selector,ngroup
#' Optional specification of kronecker/group model indexing.
#' @param control.group `list` of kronecker/group model parameters, currently
#' passed directly on to `INLA::f`
# Replicate model parameters
#' @param replicate,replicate_mapper,replicate_layer,replicate_selector,nrep
#' Optional specification of indices for an independent
#' replication model. Same syntax as for `main`
#' @param A.msk TODO: check/fix/deprecate this parameter.
#' Likely doesn't work at the moment, and I've found no examples that use it.
#' @param .envir Evaluation environment
#' @param envir_extra TODO: check/fix this parameter.
#'
#' @details The `component.character` method is inlabru's equivalent to INLA's
#' `f` function but adds functionality that is unique to inlabru.
#'
#' Deprecated parameters:
#' * map: Use `main` instead.
#' * mesh: Use `mapper` instead.
#'
#' @rdname component
#' @aliases bru_component
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla(quietly = TRUE)) {
#'   # As an example, let us create a linear component. Here, the component is
#'   # called "myEffectOfX" while the covariate the component acts on is called "x":
#'
#'   cmp <- component("myEffectOfX", main = x, model = "linear")
#'   summary(cmp)
#'
#'   # A more complicated component:
#'   cmp <- component("myEffectOfX",
#'     main = x,
#'     model = INLA::inla.spde2.matern(fm_mesh_1d(1:10))
#'   )
#'
#'   # Compound fixed effect component, where x and z are in the input data.
#'   # The formula will be passed on to MatrixModels::model.Matrix:
#'   cmp <- component("eff", ~ -1 + x:z, model = "fixed")
#'   summary(cmp)
#' }
#' }
#'
component.character <- function(object,
                                # Main model parameters
                                main = NULL, # This must be kept as 1st arg.
                                weights = NULL, # This must be kept as 2nd arg.
                                ..., # Prevent partial matching
                                model = NULL,
                                mapper = NULL,
                                main_layer = NULL,
                                main_selector = NULL,
                                n = NULL,
                                values = NULL,
                                season.length = NULL,
                                # Copy feature
                                copy = NULL,
                                # Weights
                                weights_layer = NULL,
                                weights_selector = NULL,
                                # Group model parameters
                                group = 1L,
                                group_mapper = NULL,
                                group_layer = NULL,
                                group_selector = NULL,
                                ngroup = NULL,
                                control.group = NULL,
                                # Replicate model parameters
                                replicate = 1L,
                                replicate_mapper = NULL,
                                replicate_layer = NULL,
                                replicate_selector = NULL,
                                nrep = NULL,
                                A.msk = NULL,
                                .envir = parent.frame(),
                                envir_extra = NULL) {
  # INLA models:
  # itypes = c(linear, iid, mec, meb, rgeneric, rw1, rw2, crw2, seasonal, besag, besag2, bym, bym2, besagproper,
  #            besagproper2, fgn, fgn2, ar1, ar1c, ar, ou, generic, generic0, generic1, generic2, generic3, spde,
  #            spde2, spde3, iid1d, iid2d, iid3d, iid4d, iid5d, 2diid, z, rw2d, rw2diid, slm, matern2d, copy,
  #            clinear, sigm, revsigm, log1exp, logdist)
  #
  # Supported:
  # btypes = c("const", "offset", "factor_full", "factor_contrast", "linear", "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")

  # The label
  label <- object

  if (is.null(model) && is.null(copy)) {
    # TODO: may need a special marker to autodetect factor and matrix variables,
    #       which can only be autodetected in a data-aware pass.
    model <- "linear"
  }

  # Force evaluation of explicit inputs
  force(values)

  if (!is.null(substitute(group)) &&
    !identical(deparse(substitute(group)), "1L")) {
    if (is.null(control.group)) {
      control.group <- INLA::inla.set.control.group.default()
    }
    group_model <- control.group$model
  } else {
    group_model <- "exchangeable"
  }

  if ("map" %in% names(sys.call())) {
    #    if (!is.null(substitute(map))) {
    if (is.null(substitute(main))) {
      main <- sys.call()[["map"]]
      warning("Use of 'map' is deprecated and may be disabled; use 'main' instead.",
        immediate. = TRUE
      )
    } else {
      warning("Deprecated 'map' overridden by 'main'.",
        immediate. = TRUE
      )
    }
  }

  if ("mesh" %in% names(sys.call())) {
    if (is.null(mapper)) {
      mapper <- list(...)[["mesh"]]
      warning("Use of 'mesh' is deprecated and may be disabled; use 'mapper' instead.")
    } else {
      warning("Deprecated 'mesh' overridden by 'mapper'.")
    }
  }

  if (is.null(envir_extra)) {
    envir_extra <- new.env(parent = .envir)
  }

  # Convert ngroup and nrep to bru_mapper info
  if (!is.null(ngroup)) {
    if (!is.null(group_mapper)) {
      stop("At most one of 'ngroup' and 'group_mapper' should be supplied.")
    }
    group_mapper <- bru_mapper_index(ngroup)
    ngroup <- NULL
  }
  if (!is.null(nrep)) {
    if (!is.null(replicate_mapper)) {
      stop("At most one of 'nrep' and 'replicate_mapper' should be supplied.")
    }
    replicate_mapper <- bru_mapper_index(nrep)
    nrep <- NULL
  }

  # Default component (to be filled)
  component <- list(
    label = label,
    inla.formula = NULL,
    main = bru_subcomponent(
      input = bru_input(
        substitute(main),
        label = label,
        layer = substitute(main_layer),
        selector = main_selector
      ),
      mapper = mapper,
      model = model,
      n = n,
      values = values,
      season.length = season.length
    ),
    group = bru_subcomponent(
      input = bru_input(
        substitute(group),
        label = paste0(label, ".group"),
        layer = substitute(group_layer),
        selector = group_selector
      ),
      mapper = group_mapper,
      n = NULL,
      model = group_model
    ),
    replicate = bru_subcomponent(
      input = bru_input(
        substitute(replicate),
        label = paste0(label, ".repl"),
        layer = substitute(replicate_layer),
        selector = replicate_selector
      ),
      mapper = replicate_mapper,
      n = NULL,
      model = "iid"
    ),
    weights =
      if (is.null(substitute(weights))) {
        NULL
      } else {
        bru_input(
          substitute(weights),
          label = paste0(label, ".weights"),
          layer = substitute(weights_layer),
          selector = weights_selector
        )
      },
    copy = copy,
    A.msk = A.msk,
    env = .envir,
    env_extra = envir_extra
  )

  # Main bit
  # Construct a call to the f function from the parameters provided
  # Ultimately, this call will be converted to the model formula presented to INLA
  fcall <- sys.call()
  fcall[[1]] <- "f"
  fcall[[2]] <- as.symbol(label)
  names(fcall)[2] <- ""
  # 'main' and 'weights' are the only regular parameter allowed to be nameless,
  # and only if they are the first parameters (position 3 and 4 in fcall)
  if (is.null(names(fcall)) || identical(names(fcall)[3], "")) {
    names(fcall)[3] <- "main"
  }
  if (is.null(names(fcall)) || identical(names(fcall)[4], "")) {
    names(fcall)[4] <- "weights"
  }
  unnamed_arguments <- which(names(fcall[-c(1, 2)]) %in% "")
  if (length(unnamed_arguments) > 0) {
    # Without this check, R gives the error
    #   'In str2lang(s) : parsing result not of length one, but 0'
    # in the INLA call instead, which isn't very informative.
    stop(paste0(
      "Unnamed arguments detected in component '", label, "'.\n",
      "  Only 'main' and 'weights' parameters may be unnamed.\n",
      "  Unnamed arguments at position(s) ",
      paste0(unnamed_arguments, collapse = ", ")
    ))
  }

  # Special and general cases:
  if (component$main$type %in% c("offset", "const")) {
    # The offset is included either automatically for ~ . linear models,
    # or explicitly by name in the predictor expression, so no INLA formula
    # component is needed.
    component$inla.formula <- as.formula(paste0("~ ."),
      env = .envir
    )
    component$main$mapper <- bru_mapper_const()
    component$group$mapper <- bru_mapper_index(1L)
    component$replicate$mapper <- bru_mapper_index(1L)
    # Add scalable multi-mapper
    component[["mapper"]] <-
      bru_mapper_pipe(
        list(
          mapper = bru_mapper_multi(list(
            main = component$main$mapper,
            group = component$group$mapper,
            replicate = component$replicate$mapper
          )),
          scale = bru_mapper_scale()
        )
      )
  } else {
    if (!is.null(copy)) {
      # Store copy-model name or object in the environment
      #    model_name <- paste0("BRU_", label, "_copy_model")
      #    fcall[["copy"]] <- component$copy # as.symbol(model_name)
      #    assign(model_name, component$copy, envir = component$env_extra)
    } else {
      # Store model name or object in the environment
      model_name <- paste0("BRU_", label, "_main_model")
      fcall[["model"]] <- as.symbol(model_name)
      assign(model_name, component$main$model, envir = component$env_extra)
    }

    # Remove parameters inlabru supports but INLA doesn't,
    # and substitute parameters that inlabru will transform
    fcall <- fcall[!(names(fcall) %in% c(
      "main",
      "weights",
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

    # Replace arguments that will be evaluated by a mapper
    # TODO: make a more general system
    suffixes <- list(
      "group" = "group",
      "replicate" = "repl"
    )
    for (arg in names(suffixes)) {
      if (arg %in% names(fcall)) {
        fcall[[arg]] <- as.symbol(paste0(label, ".", suffixes[[arg]]))
      }
    }

    if (is.null(copy)) {
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
        assign(
          values_name,
          component$main$values,
          envir = component$env_extra
        )
      }

      # Setup factor precision parameter
      if (component$main$type %in% c("factor", "fixed")) {
        if (is.null(fcall[["hyper"]])) {
          # TODO: allow configuration of the precision via prec.linear
          fixed_hyper_name <- paste0("BRU_", label, "_main_fixed_hyper")
          fcall[["hyper"]] <- as.symbol(fixed_hyper_name)
          assign(
            fixed_hyper_name,
            list(prec = list(
              initial = log(INLA::inla.set.control.fixed.default()$prec),
              fixed = TRUE
            )),
            envir = component$env_extra
          )
        }
      }
    }

    component$inla.formula <-
      as.formula(
        paste0(
          "~ . + ",
          as.character(parse(text = deparse(fcall)))
        ),
        env = .envir
      )

    component$fcall <- fcall
  }

  class(component) <- c("component", "list")
  component
}



#' Methods for inlabru component lists
#'
#' Constructor methods for inlabru component lists. Syntax details are given in
#' [component()].
#'
#' @param \dots Parameters passed on to other methods. Also see Details.
#' @family component constructors
#' @aliases bru_component_list
#' @param object The object to operate on
#' @param lhoods A `bru_like_list` object
#' @param .envir An evaluation environment for non-formula input
#' @export
#' @rdname component_list
component_list <- function(object,
                           lhoods = NULL,
                           .envir = parent.frame(),
                           ...) {
  UseMethod("component_list")
}

#' @details
#' * `component_list.formula`: Convert a component formula
#' into a `component_list` object
#'
#' @export
#' @family component constructors
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname component_list
#'
#' @examples
#' # As an example, let us create a linear component. Here, the component is
#' # called "myLinearEffectOfX" while the covariate the component acts on is
#' # called "x". Note that a list of components is returned because the
#' # formula may define multiple components
#'
#' eff <- component_list(~ myLinearEffectOfX(main = x, model = "linear"))
#' summary(eff[[1]])
#' # Equivalent shortcuts:
#' eff <- component_list(~ myLinearEffectOfX(x, model = "linear"))
#' eff <- component_list(~ myLinearEffectOfX(x))
#' # Individual component
#' eff <- component("myLinearEffectOfX", main = x, model = "linear")
component_list.formula <- function(object,
                                   lhoods = NULL,
                                   .envir = parent.frame(), ...) {
  if (!is.null(environment(object))) {
    .envir <- environment(object)
  }
  # Automatically add Intercept and remove -1 unless -1
  # or -Intercept is in components formula
  object <- auto_intercept(object)

  code <- code.components(object)
  parsed <- lapply(code, function(x) parse(text = x))
  components <- lapply(
    parsed,
    function(component.expression) {
      eval(component.expression,
        envir = .envir
      )
    }
  )
  environment(object) <- .envir
  component_list(components, lhoods = lhoods, .envir = .envir)
}





#' @details
#' * `component_list.list`: Combine a list of components and/or component formulas
#' into a `component_list` object
#' @export
#' @rdname component_list
component_list.list <- function(object,
                                lhoods = NULL,
                                .envir = parent.frame(),
                                ...) {
  # Maybe the list has been given an environment?
  if (!is.null(environment(object))) {
    .envir <- environment(object)
  } else if (is.null(.envir)) {
    # Later code needs an actual environment
    .envir <- new.env()
  }
  if (any(vapply(object, function(x) inherits(x, "formula"), TRUE))) {
    object <-
      do.call(
        c,
        lapply(
          object,
          function(x) {
            if (inherits(x, "formula")) {
              component_list(x, lhoods = lhoods, .envir = .envir)
            } else {
              list(x)
            }
          }
        )
      )
  }
  stopifnot(all(vapply(object, function(x) inherits(x, "component"), TRUE)))
  names(object) <- lapply(object, function(x) x$label)
  if (anyDuplicated(names(object))) {
    stop(paste0(
      "Duplicated component labels detected: ",
      paste0(
        "'",
        sort(unique(names(object)[duplicated(names(object))])),
        "'",
        collapse = ", "
      )
    ))
  }
  class(object) <- c("component_list", "list")
  environment(object) <- .envir
  if (!is.null(lhoods)) {
    lhoods <- bru_used_update(lhoods, names(object))
    object <- add_mappers(object, lhoods = lhoods)
  }
  object
}






#' @export
#' @details * `c.component_list`: The `...` arguments should be `component_list`
#' objects. The environment from the first argument will be applied to the
#' resulting `component_list`.
#' @rdname component_list
`c.component_list` <- function(...) {
  stopifnot(all(vapply(
    list(...),
    function(x) inherits(x, "component_list"),
    TRUE
  )))
  env <- environment(list(...)[[1]])
  object <- NextMethod()
  class(object) <- c("component_list", "list")
  environment(object) <- env
  object
}

#' @export
#' @details * `c.component`: The `...` arguments should be `component`
#' objects. The environment from the first argument will be applied to the
#' resulting ``component_list`.
#' @rdname component_list
`c.component` <- function(...) {
  stopifnot(all(vapply(
    list(...),
    function(x) inherits(x, "component"),
    TRUE
  )))
  env <- environment(list(...)[[1]])
  object <- list(...)
  class(object) <- c("component_list", "list")
  environment(object) <- env
  object
}

#' @export
#' @param x `component_list` object from which to extract a sub-list
#' @param i indices specifying elements to extract
#' @rdname component_list
`[.component_list` <- function(x, i) {
  env <- environment(x)
  object <- NextMethod()
  class(object) <- c("component_list", "list")
  environment(object) <- env
  object
}





#' @title Equip components with mappers
#' @description Equip component(s) with mappers for subcomponents that do not
#' have predefined mappers. When needed, the data in `lhoods` is used to determine
#' the appropriate mapper(s).
#' @param component A `component` object
#' @param lhoods A `bru_like_list` object
#' @return A `component` object with completed mapper information
#' @examples
#' \dontrun{
#' if (interactive()) {
#' }
#' }
#' @rdname add_mappers
#' @keywords internal
#' @export

add_mappers.component <- function(component, lhoods, ...) {
  # Filter out lhoods that don't use/support the component
  keep_lh <-
    vapply(lhoods,
      function(lh, label) {
        label %in% bru_used(lh)[["effect"]]
      },
      TRUE,
      label = component$label
    )
  lh <- lhoods[keep_lh]

  component$main <- add_mapper(
    component$main,
    label = component$label,
    lhoods = lh,
    env = component$env,
    require_indexed = FALSE
  )
  component$group <- add_mapper(
    component$group,
    label = component$label,
    lhoods = lh,
    env = component$env,
    require_indexed = TRUE
  )
  component$replicate <- add_mapper(
    component$replicate,
    label = component$label,
    lhoods = lh,
    env = component$env,
    require_indexed = TRUE
  )
  # Add scalable multi-mapper
  component[["mapper"]] <-
    bru_mapper_pipe(
      list(
        mapper = bru_mapper_multi(list(
          main = component$main$mapper,
          group = component$group$mapper,
          replicate = component$replicate$mapper
        )),
        scale = bru_mapper_scale()
      )
    )

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

  if (!(component[["main"]][["type"]] %in% c("offset", "const"))) {
    # Update the formula that will be presented to INLA
    component$inla.formula <-
      as.formula(
        paste0(
          "~ . + ",
          paste0(deparse(fcall),
            collapse = "\n"
          )
        ),
        env = component$env
      )
  }

  component$fcall <- fcall

  component
}

#' @param components A `component_list` object
#' @export
#' @rdname add_mappers
add_mappers.component_list <- function(components, lhoods, ...) {
  is_copy <- vapply(components, function(x) !is.null(x[["copy"]]), TRUE)
  for (k in which(!is_copy)) {
    components[[k]] <- add_mappers(components[[k]], lhoods)
  }
  for (k in which(is_copy)) {
    if (is.null(components[[k]][["copy"]])) {
      stop("Internal error: copy model detected, but no copy information available.")
    }
    if (is.null(components[[components[[k]][["copy"]]]])) {
      stop(paste0(
        "Could not find component '",
        components[[k]][["copy"]],
        "' to use as copy for component '",
        components[[k]][["label"]],
        "'."
      ))
    }
    components[[k]]$mapper <- components[[components[[k]][["copy"]]]][["mapper"]]
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
                             season.length = NULL) {
  type <- model
  factor_mapping <- NULL
  if (inherits(model, "inla.spde")) {
    type <- "spde"
  } else if (inherits(model, "inla.rgeneric")) {
    type <- "rgeneric"
  } else if (inherits(model, "inla.cgeneric")) {
    type <- "cgeneric"
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
      warning(
        paste0(
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
    } else if (identical(model, "fixed")) {
      model <- "iid"
      type <- "fixed"
    } else if (model %in% c("offset", "const")) {
      model <- "const"
      type <- "const"
    } else {
      type <- model
    }
  } else {
    type <- "unknown"
  }
  if (!is.null(mapper)) {
    if (!inherits(mapper, "bru_mapper")) {
      stop(
        "Unknown mapper class '",
        paste0(class(mapper), collapse = ", "),
        "'"
      )
    }
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
      factor_mapping = factor_mapping
    )
  class(subcomponent) <- c("bru_subcomponent", "list")

  subcomponent
}




make_unique_inputs <- function(inp, allow_list = FALSE) {
  is_spatial <- vapply(inp, function(x) inherits(x, "Spatial"), TRUE)
  is_sfc <- vapply(inp, function(x) inherits(x, "sfc"), TRUE)
  is_matrix <- vapply(inp, function(x) is.matrix(x), TRUE)
  is_Matrix <- vapply(inp, function(x) inherits(x, "Matrix"), TRUE)
  is_factor <- vapply(inp, function(x) is.factor(x), TRUE)
  is_data_frame <- vapply(inp, function(x) is.data.frame(x), TRUE)
  is_list <- vapply(inp, function(x) is.list(x), TRUE) &
    !is_data_frame
  if (any(is_spatial)) {
    if (!all(is_spatial)) {
      stop("Inconsistent spatial/non-spatial input. Unable to infer mapper information.")
    }
    inconsistent_crs <- FALSE
    inp_crs <- lapply(inp, fm_CRS)
    crs_info <- lapply(inp_crs, fm_wkt)
    null_crs <- vapply(crs_info, is.null, logical(1))
    inconsistent_crs <-
      (length(unique(unlist(crs_info))) > 1) ||
        (any(null_crs) && !all(null_crs))
    if (inconsistent_crs) {
      stop("Inconsistent spatial CRS information. Unable to infer mapper information.")
    }
    inp_values <- unique(do.call(
      rbind,
      lapply(
        inp,
        function(x) coordinates(x)
      )
    ))
    n_values <- nrow(inp_values)
  } else if (any(is_sfc)) {
    if (!all(is_sfc)) {
      stop("Inconsistent spatial/non-spatial input. Unable to infer mapper information.")
    }
    inconsistent_crs <- FALSE
    inp_crs <- lapply(inp, fm_crs)
    null_crs <- vapply(inp_crs, is.na, logical(1))
    inconsistent_crs <-
      (length(unique(inp_crs)) > 1) ||
        (any(null_crs) && !all(null_crs))
    if (inconsistent_crs) {
      stop("Inconsistent spatial crs information. Unable to infer mapper information.")
    }
    inp_values <- sf::st_sfc(unique(do.call(c, inp)),
      crs = inp_crs[[1]]
    )
    n_values <- NROW(inp_values)
  } else if (any(is_matrix | is_Matrix)) {
    if (!all(is_matrix | is_Matrix)) {
      stop("Inconsistent input types; matrix and non-matrix")
    }
    # Add extra column to work around bug in unique.matrix for single-column
    # Matrix and ModelMatrix matrices:
    inp_values <- unique.matrix(cbind(1, do.call(rbind, inp)))[, -1, drop = FALSE]
    n_values <- nrow(inp_values)
  } else if (any(is_data_frame)) {
    if (!all(is_data_frame)) {
      stop("Inconsistent input types; data.frame and non-data.frame")
    }
    inp_values <- unique(do.call(rbind, inp))
    n_values <- nrow(inp_values)
  } else if (any(is_factor)) {
    if (!all(is_factor)) {
      stop("Inconsistent input types; factor and non-factor")
    }
    inp_values <- sort(unique(do.call(c, lapply(inp, as.character))),
      na.last = NA
    )
    n_values <- length(inp_values)
  } else if (any(is_list)) {
    if (!all(is_list)) {
      stop("Inconsistent input types; list and non-list")
    }
    if (!allow_list) {
      stop("Unable to create automatic mapper. List data at this level requires an explicit mapper definition.")
    }
    inp_values <- list()
    n_values <- list()
    for (i in seq_along(inp[[1]])) {
      inp_ <- list()
      for (j in seq_along(inp)) {
        inp_[[j]] <- inp[[j]][[i]]
      }
      result <- make_unique_inputs(inp_, allow_list = FALSE)
      inp_values[[i]] <- result$inp_values
      n_values[[i]] <- result$n_values
    }
  } else {
    inp_values <- sort(unique(unlist(inp)), na.last = NA)
    n_values <- length(inp_values)
  }

  return(list(
    inp_values = inp_values,
    n_values = n_values,
    is_list
  ))
}


add_mapper <- function(subcomp, label, lhoods = NULL, env = NULL,
                       require_indexed = FALSE) {
  if (is.null(subcomp[["mapper"]])) {
    if (!inherits(subcomp[["model"]], "character")) {
      subcomp[["mapper"]] <- bru_get_mapper_safely(subcomp[["model"]])
    }
  }
  if (!is.null(subcomp[["mapper"]])) {
    if (!inherits(subcomp[["mapper"]], "bru_mapper")) {
      stop(paste0(
        "Unknown mapper of type '",
        paste0(class(subcomp[["mapper"]]), collapse = ", "),
        "' for ", label
      ))
    }
  } else if (is.null(subcomp[["input"]][["input"]])) {
    # Used for automatic group and replicate mappers
    subcomp <- make_mapper(subcomp,
      label,
      input_values = 1,
      strict = TRUE,
      require_indexed = require_indexed
    )
  } else {
    if (!is.null(lhoods)) {
      if (length(lhoods) > 0) {
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
      } else {
        # Component not directly used in any likelihood.
        # Attempt to evaluate with no data;
        # useful for intercept-like components only used via the
        # *_latent technique.
        inp <- list(
          input_eval(subcomp$input,
            data = NULL,
            env = env,
            label = subcomp$input$label,
            null.on.fail = TRUE
          )
        )
      }
      # Check for
      # 1) All NULL; Deprecated unless input is NULL. Since version 2.1.14,
      #              intercepts should be notated explicitly with label(1)
      # 2) Some NULL; exclude NULL results
      # TODO: Check for vector/matrix/coordinate inconsistency
      null.results <- vapply(inp, function(x) is.null(x), TRUE)
      if (all(null.results)) {
        warning(
          paste0(
            "All covariate evaluations for '", label,
            "' are NULL; an intercept component was likely intended.\n",
            "  Implicit latent intercept component specification is deprecated since version 2.1.14.\n",
            "  Use explicit notation '+ ", label, "(1)' instead",
            if (identical(label, "Intercept")) {
              " (or '+1' for '+ Intercept(1)')"
            },
            "."
          ),
          immediate. = TRUE
        )
        unique_inputs <- list(
          inp_values = 1,
          n_values = 1
        )
      } else {
        if (any(null.results)) {
          inp_ <- inp[!null.results]
        } else {
          inp_ <- inp
        }

        unique_inputs <- make_unique_inputs(inp_, allow_list = TRUE)
      }
      if (sum(unlist(unique_inputs$n_values)) < 1) {
        subcomp$n <- 1
        subcomp$values <- NULL
        inp_values <- NULL
      }
      subcomp <- make_mapper(subcomp,
        label,
        input_values = unique_inputs$inp_values,
        strict = TRUE,
        require_indexed = require_indexed
      )
    }
  }
  if (!is.null(subcomp[["mapper"]])) {
    # Check internal consistency of user specified n and the mapper:
    mapper_n <- ibm_n(subcomp[["mapper"]], inla_f = TRUE)
    if (!is.null(subcomp[["n"]]) &&
      subcomp[["n"]] != mapper_n) {
      stop(paste0(
        "Size mismatch, n=", subcomp[["n"]], " != ibm_n(inla_f = TRUE)=",
        mapper_n, " mapper for label ", label
      ))
    }
    subcomp[["n"]] <- mapper_n
    subcomp[["values"]] <- ibm_values(subcomp[["mapper"]], inla_f = TRUE)
  }
  subcomp
}



make_submapper <- function(subcomp_n,
                           subcomp_values,
                           input_values,
                           label,
                           subcomp_type,
                           subcomp_factor_mapping,
                           require_indexed,
                           allow_interpolation = TRUE) {
  if (!is.null(subcomp_n)) {
    if (!is.null(subcomp_values)) {
      warning(
        "Both 'n' and 'values' provided. Ignoring 'values'",
        immediate. = TRUE
      )
    }
    values <- seq_len(subcomp_n)
  } else if (!is.null(subcomp_values)) {
    values <- subcomp_values
  } else if (!is.null(input_values)) {
    values <- input_values
  } else {
    stop(paste0("No mapper, no n, and no values given for ", label))
  }

  if (is.factor(values) ||
    is.character(values) ||
    (!is.null(subcomp_type) && (subcomp_type %in% "factor"))) {
    mapper <- bru_mapper_factor(
      values,
      factor_mapping = subcomp_factor_mapping,
      indexed = require_indexed
    )
  } else {
    values <- sort(unique(values), na.last = NA)
    if (length(values) > 1) {
      if (allow_interpolation) {
        mapper <-
          bru_mapper(
            fm_mesh_1d(values),
            indexed = require_indexed
          )
      } else {
        mapper <- bru_mapper_index(n = ceiling(max(values)))
      }
    } else if (all(values == 1)) {
      if (require_indexed) {
        mapper <- bru_mapper_index(n = 1)
      } else {
        mapper <- bru_mapper_linear()
      }
    } else {
      if (require_indexed) {
        mapper <- bru_mapper_factor(values,
          factor_mapping = "full",
          indexed = TRUE
        )
      } else {
        mapper <- bru_mapper_linear()
      }
    }
  }

  mapper
}


#' Extract mapper information from INLA model component objects
#'
#' The component definitions will automatically attempt to extract mapper
#' information from any model object by calling the generic `bru_get_mapper`.
#' Any class method implementation should return a [bru_mapper] object suitable
#' for the given latent model.
#'
#' @param model A model component object
#' @param \dots Arguments passed on to other methods
#' @return A [bru_mapper] object defined by the model component
#' @seealso [bru_mapper] for mapper constructor methods, and
#' [bru_mapper_methods] for method generics and specific implementations.
#' @export
#' @examples
#' if (bru_safe_inla(quietly = TRUE)) {
#'   library(INLA)
#'   mesh <- fmesher::fm_rcdt_2d_inla(globe = 2)
#'   spde <- inla.spde2.pcmatern(mesh,
#'     prior.range = c(1, 0.5),
#'     prior.sigma = c(1, 0.5)
#'   )
#'   mapper <- bru_get_mapper(spde)
#'   ibm_n(mapper)
#' }
bru_get_mapper <- function(model, ...) {
  UseMethod("bru_get_mapper", model)
}


#' @rdname bru_get_mapper
#' @details * `bru_get_mapper.inla.spde` extract an indexed mapper for
#' the `model$mesh` object contained in the model object.
#' It returns `NULL` gives a warning
#' if no known mesh type is found in the model object.
#' @export
bru_get_mapper.inla.spde <- function(model, ...) {
  if (inherits(model$mesh, c("fm_mesh_2d", "inla.mesh"))) {
    mapper <- bru_mapper(model$mesh)
  } else if (inherits(model$mesh, c("fm_mesh_1d", "inla.mesh.1d"))) {
    mapper <- bru_mapper(model$mesh, indexed = TRUE)
  } else {
    mapper <- NULL
    warning(
      paste0(
        "Unknown SPDE mesh class '",
        paste0(class(model$mesh), collapse = ", "),
        "' for bru_get_mapper.inla.spde. Please specify a mapper manually instead."
      ),
      immediate. = TRUE
    )
  }
  mapper
}

#' @rdname bru_get_mapper
#' @details * `bru_get_mapper.inla.rgeneric` returns the mapper given by a call to
#' `model$f$rgeneric$definition("mapper")`. To support this for your own
#' `inla.rgeneric` models, add a `"mapper"` option to the `cmd` argument
#' of your rgeneric definition function. You will need to store the mapper
#' in your object as well.  Alternative, define your model using a subclass
#' and define a corresponding `bru_get_mapper.subclass` method that should return
#' the corresponding `bru_mapper` object.
#' @export
bru_get_mapper.inla.rgeneric <- function(model, ...) {
  if (is.null(model[["f"]][["rgeneric"]][["definition"]])) {
    NULL
  } else {
    model[["f"]][["rgeneric"]][["definition"]]("mapper")
  }
}

#' @rdname bru_get_mapper
#' @details * `bru_get_mapper_safely` tries to call the `bru_get_mapper`,
#' and returns `NULL` if it fails (e.g. due to no available class method).
#' If the call succeeds and returns non-`NULL`, it checks that the object
#' inherits from the `bru_mapper` class, and gives an error if it does not.
#' @export
bru_get_mapper_safely <- function(model, ...) {
  m <- tryCatch(
    bru_get_mapper(model, ...),
    error = function(e) {
    }
  )
  if (!is.null(m) && !inherits(m, "bru_mapper")) {
    stop(paste0(
      "The bru_get_mapper method for model class '",
      paste0(class(model), collapse = ", "),
      "' did not return a bru_mapper object"
    ))
  }
  m
}


# Defines a default mapper given the type of model and parameters provided
# Checks subcomp$mapper, subcomp$model (for "bym2" and other multicomponent
# models), subcomp$model$mesh (for spde models), subcomp$n,
# subcomp$values, or input_values, in that order.
make_mapper <- function(subcomp,
                        label,
                        input_values = NULL,
                        strict = TRUE,
                        require_indexed = FALSE) {
  if (is.null(subcomp[["mapper"]])) {
    if (!inherits(subcomp[["model"]], "character")) {
      subcomp[["mapper"]] <- bru_get_mapper_safely(subcomp[["model"]])
    }
  }
  if (!is.null(subcomp[["mapper"]])) {
    if (!inherits(subcomp[["mapper"]], "bru_mapper")) {
      stop(paste0(
        "Unknown mapper of type '",
        paste0(class(subcomp[["mapper"]]), collapse = ", "),
        "' for ", label
      ))
    }
  } else if (subcomp[["type"]] %in% c("linear", "clinear")) {
    subcomp[["mapper"]] <- bru_mapper_linear()
  } else if (subcomp[["type"]] %in% c("offset", "const")) {
    subcomp[["mapper"]] <- bru_mapper_const()
  } else if (subcomp[["type"]] %in% c("fixed")) {
    if (!is.null(subcomp[["values"]])) {
      labels <- subcomp[["values"]]
    } else if (!is.null(input_values)) {
      tmp <- as(input_values, "Matrix")
      labels <- colnames(tmp)
      if (is.null(labels)) {
        labels <- as.character(seq_len(ncol(tmp)))
      }
    } else if (!is.null(subcomp[["n"]])) {
      labels <- as.character(seq_len(subcomp[["n"]]))
    } else {
      stop(paste0(
        "Need to specify at least one of values (labels), input, or n for '",
        subcomp[["label"]],
        "' of component '", label, "'."
      ))
    }
    subcomp[["mapper"]] <- bru_mapper_matrix(labels)
  } else if (subcomp[["model"]] %in% c("bym", "bym2")) {
    # No mapper; construct based on input values
    # Default mapper will be integer or factor indexed, not
    # interpolated
    mappers <- list(
      make_submapper(
        subcomp_n = subcomp[["n"]],
        subcomp_values = subcomp[["values"]],
        input_values = input_values,
        label = paste0(subcomp[["label"]], " part 1"),
        subcomp_type = subcomp[["type"]],
        subcomp_factor_mapping = subcomp[["factor_mapping"]],
        require_indexed = require_indexed,
        allow_interpolation = FALSE
      )
    )
    mappers[[2]] <- mappers[[1]]
    names(mappers) <- c("u", "v")
    subcomp[["mapper"]] <-
      bru_mapper_collect(mappers, hidden = TRUE)
  } else {
    # No mapper; construct based on input values
    # Interpolation on by default
    subcomp[["mapper"]] <-
      make_submapper(
        subcomp_n = subcomp[["n"]],
        subcomp_values = subcomp[["values"]],
        input_values = input_values,
        label = subcomp[["label"]],
        subcomp_type = subcomp[["type"]],
        subcomp_factor_mapping = subcomp[["factor_mapping"]],
        require_indexed = require_indexed,
        allow_interpolation = TRUE
      )
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
  # If rhs is "~1", make sure there's at least one component to parse
  # and that offsets show up in "factors"
  tms <- terms(update.formula(components, . ~ . + BRU_DUMMY_COMPONENT + 1))
  codes <- attr(tms, "term.labels")

  # Check for offset()
  offset_idx <- attr(tms, "offset")
  if (length(offset_idx) > 0) {
    isoff <- as.vector(unlist(lapply(
      rownames(attr(tms, "factors")),
      function(s) substr(s, 1, 6) == "offset"
    )))
    if (!any(isoff)) {
      stop("Internal error: multiple offsets indicated but none extracted")
    }
    codes <- c(codes, rownames(attr(tms, "factors"))[isoff])
  }
  codes <- codes[codes != "BRU_DUMMY_COMPONENT"]

  for (k in seq_along(codes)) {
    code <- codes[[k]]

    # Function syntax or linear component?
    ix <- regexpr("(", text = code, fixed = TRUE)
    is.offset <- FALSE
    if (ix > 0) {
      label <- substr(code, 1, ix - 1)
      is.fixed <- FALSE
      if (label %in% c("offset", "const")) {
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
    } else {
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
          sub(
            paste0(label, "("),
            paste0(
              fname, '("', label, '"',
              ', model = "const", main = '
            ),
            code,
            fixed = TRUE
          )
      } else {
        codes[[k]] <- sub(paste0(label, "("),
          paste0(fname, "(\"", label, "\", "),
          code,
          fixed = TRUE
        )
      }
    }
  }

  codes
}




# OPERATORS ----


#' Summarise components
#'
#' @export
#' @method summary component
#' @keywords internal
#' @param object A `component` or `component_list`.
#' @param ... ignored.
#' @param depth The depth of which to expand the component mapper.
#' Default `Inf`, to traverse the entire mapper tree.
#' @param verbose logical; If `TRUE`, includes more details of the
#' component definitions. When `FALSE`, only show basic component
#' definition information.  Default `TRUE`.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#'

summary.component <- function(object, ..., depth = Inf, verbose = TRUE) {
  result <- list(
    "Label" = object[["label"]],
    "Copy_of" = object[["copy"]],
    "Type" = paste0(
      unlist(lapply(
        c("main", "group", "replicate"),
        function(x) {
          if (is.null(object[[x]][["input"]][["input"]])) {
            NULL
          } else {
            paste0(x, " = ", object[[x]][["type"]])
          }
        }
      )),
      collapse = ", "
    ),
    "Input" = paste0(
      unlist(lapply(
        c("main", "group", "replicate", "weights"),
        function(x) {
          if (is.null(object[[x]][["input"]][["input"]])) {
            NULL
          } else {
            paste0(
              x, " = ",
              deparse(object[[x]][["input"]][["input"]])
            )
          }
        }
      )),
      collapse = ", "
    ),
    "Mapper" =
      if (is.null(object[["mapper"]])) {
        "Not yet initialised"
      } else {
        strwrap(
          summary(object[["mapper"]],
            prefix = "    ",
            initial = "",
            depth = depth
          ),
          prefix = "    ",
          initial = ""
        )
      },
    "INLA formula" =
      paste0(
        "\n",
        strwrap(
          paste(as.character(object$inla.formula), collapse = " "),
          width = 0.9 * getOption("width"),
          indent = 4,
          exdent = 6
        )
      )
  )
  if (!verbose) {
    Summary <- paste0(
      object[["label"]],
      if (is.null(object[["copy"]])) NULL else paste0("(=", object[["copy"]], ")"),
      ": ",
      paste0(
        unlist(lapply(
          c("main", "group", "replicate", "weights"),
          function(x) {
            if (is.null(object[[x]][["input"]][["input"]])) {
              NULL
            } else {
              paste0(
                x, " = ", object[[x]][["type"]],
                "(",
                deparse(object[[x]][["input"]][["input"]]),
                ")"
              )
            }
          }
        )),
        collapse = ", "
      )
    )
    result$Summary <- Summary
  }
  class(result) <- c("summary_component", "list")
  result
}

#' @export
#' @method summary component_list
#' @param ... passed on to other summary methods
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname summary.component

summary.component_list <- function(object, verbose = TRUE, ...) {
  result <- lapply(
    object,
    function(x) {
      summary(x, verbose = verbose, ...)
    }
  )
  for (cp in which(vapply(result, function(x) !is.null(x$Copy_of), FALSE))) {
    result[[cp]]$Type <- result[[result[[cp]][["Copy_of"]]]][["Type"]]
  }
  class(result) <- c("summary_component_list", "list")
  result
}

#' @export
#' @param x A 'summary_component' object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.component

print.summary_component <- function(x, ...) {
  if (is.null(x[["Summary"]])) {
    for (name in names(x)) {
      if (!is.null(x[[name]])) {
        if (name %in% "Label") {
          cat("Label:", "\t", x[[name]], "\n", sep = "")
        } else if (name %in% "Mapper") {
          cat("  ", "Map: ", x[[name]], "\n", sep = "")
        } else {
          cat("  ", name, ":", "\t", x[[name]], "\n", sep = "")
        }
      }
    }
  } else {
    cat(x[["Summary"]], "\n", sep = "")
  }
  invisible(x)
}

#' @export
#' @param x A summary object to be printed.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.component

print.summary_component_list <- function(x, ...) {
  lapply(x, print)
  invisible(x)
}




#' @export
#' @keywords internal
#' @param component A component.
#' @param input Component inputs, from `input_eval()`
#' @param state linearisation evaluation state
#' @param ... Optional parameters passed on to `ibm_eval`
#' and `ibm_jacobian.
#' @return A `bru_mapper_taylor` or `comp_simple_list` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname comp_lin_eval

comp_lin_eval.component <- function(component,
                                    input = NULL,
                                    state = NULL,
                                    ...) {
  if (is.null(state)) {
    state <- rep(0, ibm_n(component[["mapper"]]))
  }
  ibm_linear(component[["mapper"]], input = input, state = state, ...)
}

#' @export
#' @rdname comp_lin_eval

comp_lin_eval.component_list <- function(components, input, state, ...) {
  # Note: Make sure the list element names carry over!
  mappers <-
    lapply(
      components,
      function(x) {
        label <- x[["label"]]
        comp_lin_eval(x,
          input = input[[label]],
          state = state[[label]],
          ...
        )
      }
    )

  class(mappers) <- c("comp_simple_list", class(mappers))
  mappers
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
#' In inlabru this can be achieved in several ways of using the `main` parameter
#' (`map` in version 2.1.13 and earlier), which does not need to be named.
#'
#' \itemize{
#' \item{`components = y ~ psi(main = x^2, model = "linear")`}
#' \item{`components = y ~ psi(x^2, model = "linear")`}
#' \item{`components = y ~ psi(mySquareFun(x), model = "linear")`,}
#' \item{`components = y ~ psi(myOtherSquareFun, model = "linear")`,}
#'
#' }
#'
#' In the first example inlabru will interpret the map parameter as an expression to be evaluated within
#' the data provided. Since \eqn{x} is a known covariate it will know how to calculate it. The second
#' example is an expression as well but it uses a function called `mySquareFun`. This function is
#' defined by user but has to be accessible within the work space when setting up the components.
#' The third example provides the function `myOtherSquareFun`. In this case,
#'  inlabru will call the function as `myOtherSquareFun(.data.)`, where `.data.`
#'  is the data provided via the [like()] `data` parameter.
#' The function needs to know what parts of the data to use to construct the
#' needed output. For example,
#' ```
#' myOtherSquareFun <- function(data) {
#'   data[ ,"x"]^2
#' }
#' ```
#'
#' @section Spatial Covariates:
#'
#' When fitting spatial models it is common to work with covariates that depend on space, e.g. sea
#' surface temperature or elevation. Although it is straightforward to add this data to the input
#' data frame or write a covariate function like in the previous section there is an even more
#' convenient way in inlabru. Spatial covariates are often stored as `SpatialPixelsDataFrame`,
#' `SpatialPixelsDataFrame` or `RasterLayer` objects. These can be provided directly via
#' the input expressions if they are supported by [eval_spatial()], and
#' the [like()] data is an `sf` or `SpatialPointsDataFrame` object.
#' `inlabru` will then automatically
#' evaluate and/or interpolate the covariate at your data locations when using code like
#' ```
#' components = y ~ psi(mySpatialPixels, model = "linear")
#' ```
#' For more precise control, use the the `layer` and `selector` arguments (see [component()]),
#' or call `eval_spatial()` directly, e.g.:
#' ```
#' components = y ~ psi(eval_spatial(mySpatialPixels, where = .data.), model = "linear")
#' ```
#'
#' @section Coordinates:
#'
#' A common spatial modelling component when using inla are SPDE models. An important feature of
#' inlabru is that it will automatically calculate the so called A-matrix (a component model matrix)
#' which maps SPDE
#' values at the mesh vertices to values at the data locations. For this purpose, the input
#' can be set to `coordinates`, which is the `sp` package function that extracts point
#' coordinates from the `SpatialPointsDataFrame` that was provided as input to [like()]. The code for
#' this would look as follows:
#' ```
#' components = y ~ field(coordinates, model = inla.spde2.matern(...))
#' ```
#' Since `coordinates` is a function from the `sp` package, this results in
#' evaluation of `sp::coordinates(.data.)`, which loses any CRS information
#' from the data object.
#'
#' For `sf` data with a geometry column (by default named `geometry`), use
#' ```
#' components = y ~ field(geometry, model = inla.spde2.matern(...))
#' ```
#' Since the CRS information is part of the geometry column of the `sf` object,
#' this retains CRS information, so this is more robust, and allows the model
#' to be built on a different CRS than the observation data.
#'
#'
#' @export
#' @keywords internal
#' @param component A component.
#' @param data A `data.frame`, `tibble`, `sf`, `list`, or `Spatial*` object of
#' covariates and/or point locations.
#' If `NULL`, return the component's map.
#' @param ... Unused.
#' @return An list of mapper input values, formatted for the full component mapper
#' (of type `bru_mapper_pipe`)
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}, Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname input_eval

input_eval.component <- function(component,
                                 data,
                                 ...) {
  stopifnot(inherits(component[["mapper"]], "bru_mapper_pipe"))

  # The names should be a subset of main, group, replicate
  part_names <- ibm_names(component[["mapper"]][["mappers"]][[1]])
  mapper_val <- list()
  for (part in part_names) {
    mapper_val[[part]] <-
      input_eval(
        component[[part]]$input,
        data,
        env = component$env,
        label = part,
        ...
      )
  }
  if (is.null(component[["weights"]])) {
    scale_val <- NULL
  } else {
    scale_val <-
      input_eval(
        component[["weights"]],
        data,
        env = component$env,
        label = paste(component$label, "(weights)"),
        ...
      )
  }
  # Any potential length mismatches must be handled by bru_mapper_multi and
  # bru_mapper_collect, since e.g. 'main' might take a list() input object.
  list(mapper = mapper_val, scale = scale_val)
}

#' @export
#' @rdname input_eval

input_eval.component_list <-
  function(components,
           data,
           ...) {
    lapply(components, function(x) input_eval(x, data = data, ...))
  }



input_eval_layer <- function(layer, selector = NULL, envir, enclos,
                             label,
                             e_input) {
  input_layer <- tryCatch(
    eval(layer, envir = envir, enclos = enclos),
    error = function(e) {
      e
    }
  )
  if (inherits(input_layer, "error")) {
    stop(paste0(
      "Failed to evaluate 'layer' input '",
      deparse(layer),
      "' for '",
      paste0(label, ":layer"),
      "'."
    ))
  }
  if (is.null(input_layer) && is.null(selector)) {
    if (label %in% names(e_input)) {
      input_layer <- label
    }
  }
  input_layer
}


#' @describeIn input_eval Attempts to evaluate a component input (e.g. `main`,
#' `group`, `replicate`, or `weight`), and process the results:
#' 1. If `eval()` failed, return NULL or map everything to 1
#'    (see the `null.on.fail` argument). This should normally not
#'    happen, unless the component use logic is incorrect,
#'    (e.g. via `include`/`exclude`)
#'    leading to missing columns for a certain likelihood in a
#'    multi-`like()` model.
#' 2. If we obtain a function, apply the function to the data object
#' 3. If we obtain an object supported by [eval_spatial()], extract the values
#'    of that data frame at the point locations
#' 4. Else we obtain a vector and return as-is. This happens when input
#'    references a column of the data points, or some other complete expression
#'
#' @seealso [component()]
#' @export
input_eval.bru_input <- function(input, data, env = NULL,
                                 null.on.fail = FALSE, ...) {
  # Evaluate the map with the data in an environment
  enclos <-
    if (is.null(env)) {
      parent.frame()
    } else {
      env
    }
  envir <- new.env(parent = enclos)
  if (is.list(data)) {
    for (nm in names(data)) {
      assign(nm, data[[nm]], envir = envir)
    }
  } else {
    data_df <- as.data.frame(data)
    for (nm in names(data_df)) {
      assign(nm, data_df[[nm]], envir = envir)
    }
  }
  assign(".data.", data, envir = envir)

  e_input <- tryCatch(
    eval(input$input, envir = envir, enclos = enclos),
    error = function(e) {
      e
    }
  )

  # ## Need to handle varying input lengths.
  # ## auto-expansion of scalars needs to happen elsewhere, where it's needed,
  # ## e.g. when using component A-matrices to construct a default linear model
  # ## A matrix.
  #  n <- nrow(as.data.frame(data))

  handle_problems <- function(val) {
    if (is.null(val)) {
      return(NULL)
    }
    if (inherits(e_input, "simpleError")) {
      if (null.on.fail) {
        return(NULL)
      }

      val <- 1
      input_string <- deparse(input$input)
      if (identical(input_string, "coordinates")) {
        warning(
          paste0(
            "The input evaluation '",
            input_string,
            "' for '", input$label,
            "' failed. Perhaps you need to load the 'sp' package with 'library(sp)'?",
            " Attempting 'sp::coordinates'."
          ),
          immediate. = TRUE
        )
        input$input <- expression(sp::coordinates)
        return(input_eval(
          input,
          data = data,
          env = env,
          null.on.fail = null.on.fail,
          ...
        ))
      } else {
        warning(
          paste0(
            "The input evaluation '",
            input_string,
            "' for '", input$label,
            "' failed. Perhaps the data object doesn't contain the needed variables?",
            " Falling back to '1'."
          ),
          immediate. = TRUE
        )
      }
      return(val)
    }
    return(val)
  }

  e_input <- handle_problems(e_input)
  if (is.null(e_input)) {
    return(NULL)
  }

  if (is.function(e_input)) {
    # Allow but detect failures:
    val <- tryCatch(
      e_input(data),
      error = function(e) {
        e
      }
    )

    val <- handle_problems(val)
    if (is.null(val)) {
      return(NULL)
    }

    if (identical(as.character(input$input), "coordinates")) {
      tryCatch(
        expr = {
          # Return SpatialPoints instead of a matrix
          val <- as.data.frame(val)
          coordinates(val) <- seq_len(ncol(val))
          # Allow proj4string failures:
          data_crs <- tryCatch(fm_CRS(data),
            error = function(e) {
            }
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
  } else if (inherits(e_input, "formula")) {
    # Allow but detect failures:
    val <- tryCatch(
      MatrixModels::model.Matrix(e_input, data = data, sparse = TRUE),
      error = function(e) {
        e
      }
    )
    val <- handle_problems(val)
    if (is.null(val)) {
      return(NULL)
    }
    val <- as(val, "Matrix")
  } else if (inherits(
    e_input,
    gsub(
      pattern = "^eval_spatial\\.([^*]*)\\*?",
      replacement = "\\1",
      x = format(utils::.S3methods("eval_spatial"))
    )
  )) {
    input_layer <-
      input_eval_layer(
        layer = input[["layer"]],
        selector = input[["selector"]],
        envir = envir,
        enclos = enclos,
        label = input[["label"]],
        e_input = e_input
      )
    layer <- extract_layer(
      data,
      input_layer,
      input[["selector"]]
    )
    check_layer(e_input, data, layer)
    val <- eval_spatial(
      e_input,
      data,
      layer = layer,
      selector = NULL
    )
    #  } else if ((input$label %in% c("offset", "const")) &&
    #    is.numeric(e_input) &&
    #    (length(e_input) == 1)) {
    #    val <- rep(e_input, n)
  } else {
    val <- e_input
  }

  # ## Need to allow different sizes; K %*% effect needs to have same length as response, but the input and effect effect itself doesn't.
  #  # Always return as many rows as data has
  #  if (is.vector(val) && (length(val) == 1) && (n > 1)) {
  #    val <- rep(val, n)
  #  }

  # Check if any of the locations are NA. If we are dealing with SpatialGridDataFrame try
  # to fix that by filling in nearest neighbour values.
  # # TODO: Check how to deal with this fully in the case of multilikelihood models
  # Answer: should respect the lhood "include/exclude" info for the component list
  if ((inherits(e_input, c(
    "SpatialGridDataFrame",
    "SpatialPixelsDataFrame",
    "SpatRaster"
  ))) &&
    any(is.na(as.data.frame(val)))) {
    warning(
      paste0(
        "Model input '",
        deparse(input$input),
        "' for '", input$label,
        "' returned some NA values.\n",
        "Attempting to fill in spatially by nearest available value.\n",
        "To avoid this basic covariate imputation, supply complete data."
      ),
      immediate. = TRUE
    )

    val <- bru_fill_missing(
      data = e_input, where = data, values = val,
      layer = layer, selector = NULL
    )
  }

  # Check for NA values.
  # TODO: update 2022-11-04:
  # iid models now by default are given _index mappers, which treat NAs
  # as zero.  An option to turn on a warning could be useful for checking
  # models that aren't supposed to have NA inputs.
  #
  #  if (any(is.na(unlist(as.vector(val))))) {
  #    # TODO: remove this check and make sure NAs are handled properly elsewhere,
  #    # if possible. Problem: treating NA as "no effect" can have negative side
  #    # effects. For spatial covariates, can be handled by infill, but a general
  #    # solution doesn't exist, so not sure how to deal with this.
  #    stop(sprintf(
  #      "Input '%s' of component '%s' has returned NA values. Please design your
  #                 argument as to return non-NA for all points in your model domain/mesh.",
  #      as.character(input$input)[[1]], label
  #    ))
  #  }

  val
}




#' @export
#' @keywords internal
#' @param component A component.
#' @param inla_f logical; when `TRUE`, must result in
#' values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' @param ... Unused.
#' @return a list of indices into the latent variables compatible with the
#' component mapper.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com},
#'   Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname index_eval

index_eval.component <- function(component, inla_f, ...) {
  idx <- ibm_values(component[["mapper"]], inla_f = inla_f, multi = TRUE)
  names(idx) <- paste0(component[["label"]], c("", ".group", ".repl"))
  idx
}


#' @export
#' @rdname index_eval

index_eval.component_list <- function(components, inla_f, ...) {
  lapply(components, function(x) index_eval(x, inla_f = inla_f, ...))
}

#' @export
#' @keywords internal
#' @param components A component list.
#' @param ... Unused.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname inla_subset_eval

inla_subset_eval.component_list <- function(components, ...) {
  lapply(components, function(x) ibm_inla_subset(x[["mapper"]], multi = TRUE))
}
