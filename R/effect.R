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

#' Methods for `bru_mapper` objects
#'
#' @details
#' * `bru_mapper` Generic mapper S3 constructor. See below for details of the
#' default constructor that can be used to define new mappers in user code.
#' @param \dots Arguments passed on to other methods
#' @param mapper A mapper S3 object, normally inheriting from `bru_mapper`.
#' For the default `bru_mapper` method, a list that will be converted to a
#' `bru_mapper` object by adding class information and (optional) methods.
#' @export
#' @seealso [bru_mapper_methods] for specific method implementations.
#' @rdname bru_mapper
#' @examples
#' mapper <- bru_mapper_index(5)
#' ibm_amatrix(mapper, c(1, 3, 4, 5, 2))
bru_mapper <- function(...) {
  UseMethod("bru_mapper")
}
#' @details
#' * `ibm_n` Generic. Implementations must return the size of the latent vector
#' being mapped to.
#' @param inla_f logical; when `TRUE` in `ibm_n`, `ibm_values`, and
#' `ibm_amatrix` methods, these must result in values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#'
#' @export
#' @rdname bru_mapper
ibm_n <- function(mapper, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_n"]])) {
    mapper[[".envir"]][["ibm_n"]](mapper, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_n")
  }
}
#' @details
#' * `ibm_values` Generic. Implementations must return a vector that
#' would be interpretable by an `INLA::f(..., values = ...)` specification.
#' The exception is the method for `bru_mapper_multi`, that returns a
#' multi-column data frame
#' @export
#' @rdname bru_mapper
ibm_values <- function(mapper, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_values"]])) {
    mapper[[".envir"]][["ibm_values"]](mapper, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_values")
  }
}
#' @details
#' * `ibm_amatrix` Generic.
#' Implementations must return a (sparse) matrix of size `NROW(input)`
#' (except for the `bru_mapper_multi` and `bru_mapper_collect` methods,
#' that require `list()` inputs, and the input size is determined by the
#' combined inputs)
#' by `ibm_n(mapper, inla_f = FALSE)`. The `inla_f=TRUE` argument should only affect
#' the allowed type of input format.
#' @param input The values for which to produce a mapping matrix
#' @export
#' @rdname bru_mapper
ibm_amatrix <- function(mapper, input, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_amatrix"]])) {
    mapper[[".envir"]][["ibm_amatrix"]](
      mapper, input, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_amatrix")
  }
}
#' @details
#' * `ibm_inla_subset` Generic.
#' Implementations must return a logical vector of `TRUE/FALSE` for
#' the subset such that, given the full A matrix and values output,
#' `A[, subset, drop = FALSE]` and `values[subset]`
#' (or `values[subset, , drop = FALSE]` for data.frame values) are equal
#' to the `inla_f = TRUE` version of A and values. The default method uses
#' the `ibm_values` output to construct the subset indexing.
#' @export
#' @rdname bru_mapper
ibm_inla_subset <- function(mapper, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_inla_subset"]])) {
    mapper[[".envir"]][["ibm_inla_subset"]](
      mapper, ...)
  } else {
    UseMethod("ibm_inla_subset")
  }
}
#' @details
#' * `ibm_valid_input` Generic.
#' Implementations must return a logical vector of length `NROW(input)` (
#' or for `bru_mapper_multi` and `bru_mapper_collect` a list of such
#' vectors)
#' @param input The values for which to produce validity information
#' @export
#' @rdname bru_mapper
ibm_valid_input <- function(mapper, input, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_valid_input"]])) {
    mapper[[".envir"]][["ibm_valid_input"]](
      mapper, input, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_valid_input")
  }
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
#' \item{\eqn{\beta} }{is a *random variable* associated with \eqn{x} and}
#' \item{\eqn{\psi = \beta * x }}{ is called the *random effect* of \eqn{x}}
#' }
#'
#' A problem that arises when using this kind of R formula is that it does not
#' clearly reflect the mathematical
#' formula. For instance, when providing the formula to inla, the resulting
#' object will refer to the random
#' effect \eqn{\psi = \beta * x } as `x`. Hence, it is not clear if `x` refers to the covariate
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
#' @param model Either one of "offset", "factor_full", "factor_contrast", "linear",
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
# Deprecated parameters
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
#' if (require("INLA", quietly = TRUE)) {
#'
#'   # As an example, let us create a linear component. Here, the component is
#'   # called "myEffectOfX" while the covariate the component acts on is called "x":
#'
#'   cmp <- component("myEffectOfX", main = x, model = "linear")
#'   summary(cmp)
#'
#'   # A more complicated component:
#'   cmp <- component("myEffectOfX",
#'     main = x,
#'     model = inla.spde2.matern(inla.mesh.1d(1:10))
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
                                group = NULL,
                                group_mapper = NULL,
                                group_layer = NULL,
                                group_selector = NULL,
                                ngroup = NULL,
                                control.group = NULL,
                                # Replicate model parameters
                                replicate = NULL,
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
  # btypes = c("offset", "factor_full", "factor_contrast", "linear", "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")

  # The label
  label <- object

  if (is.null(model) && is.null(copy)) {
    # TODO: may need a special marker to autodetect factor and matrix variables,
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
        layer = main_layer,
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
        layer = group_layer,
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
        layer = replicate_layer,
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
          layer = weights_layer,
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
  if (component$main$type %in% c("offset")) {
    # The offset is included either automatically for ~ . linear models,
    # or explicitly by name in the predictor expression, so no INLA formula
    # component is needed.
    component$inla.formula <- as.formula(paste0("~ ."),
      env = .envir
    )
    component$main$mapper <- bru_mapper_offset()
    component$group$mapper <- bru_mapper_index(1L)
    component$replicate$mapper <- bru_mapper_index(1L)
    # Add multi-mapper
    component[["mapper"]] <-
      bru_mapper_multi(list(
        main = component$main$mapper,
        group = component$group$mapper,
        replicate = component$replicate$mapper
      ))
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
      as.formula(paste0(
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
    warning(paste0(
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
    object <- add_mappers(object, lhoods = lhoods)
  }
  object
}






#' @export
#' @details * `c.component_list`: The `...` arguments should be `component_list`
#' objects. The environment from the first argument will be applied to the resulting list.
#' @rdname component_list
`c.component_list` <- function(...) {
  env <- environment(..1)
  object <- NextMethod()
  class(object) <- c("component_list", "list")
  environment(object) <- env
  object
}

#' @export
#' @param x `component_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @rdname component_list
`[.component_list` <- function(x, i) {
  env <- environment(x)
  object <- NextMethod()
  class(object) <- c("component_list", "list")
  environment(object) <- env
  object
}





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
#' @rdname add_mappers
#' @keywords internal
#' @export

add_mappers.component <- function(component, lhoods, ...) {
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

  component$main <- add_mapper(component$main,
    label = component$label,
    lhoods = lh,
    env = component$env,
    require_indexed = FALSE
  )
  component$group <- add_mapper(component$group,
    label = component$label,
    lhoods = lh,
    env = component$env,
    require_indexed = TRUE
  )
  component$replicate <- add_mapper(component$replicate,
    label = component$label,
    lhoods = lh,
    env = component$env,
    require_indexed = TRUE
  )
  # Add multi-mapper
  component[["mapper"]] <-
    bru_mapper_multi(list(
      main = component$main$mapper,
      group = component$group$mapper,
      replicate = component$replicate$mapper
    ))

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

  if (!identical(component[["main"]][["type"]], "offset")) {
    # Update the formula that will be presented to INLA
    component$inla.formula <-
      as.formula(paste0(
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
    } else if (identical(model, "fixed")) {
      model <- "iid"
      type <- "fixed"
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
    inp_crs <- lapply(inp, fm_sp_get_crs)
    if (fm_has_PROJ6()) {
      crs_info <- lapply(inp_crs, fm_crs_get_wkt)
      null_crs <- vapply(crs_info, is.null, logical(1))
      inconsistent_crs <-
        (length(unique(unlist(crs_info))) > 1) ||
          (any(null_crs) && !all(null_crs))
    } else {
      crs_info <- vapply(inp_crs, fm_CRSargs, "")
      inconsistent_crs <- length(unique(crs_info)) > 1
    }
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
  } else if (any(is_matrix | is_Matrix)) {
    if (!all(is_matrix | is_Matrix)) {
      stop("Inconsistent input types; matrix and non-matrix")
    }
    inp_values <- unique.matrix(do.call(rbind, inp))
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
      # 1) All NULL; Deprecated unless input is NULL. Since version 2.1.14,
      #              intercepts should be notated explicitly with label(1)
      # 2) Some NULL; exclude NULL results
      # TODO: Check for vector/matrix/coordinate inconsistency
      null.results <- vapply(inp, function(x) is.null(x), TRUE)
      if (all(null.results)) {
        warning(paste0(
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
      if (unique_inputs$n_values < 1) {
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
                           require_indexed) {
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
      factor_mapping = subcomp_factor_mapping
    )
  } else {
    values <- sort(unique(values), na.last = NA)
    if (length(values) > 1) {
      mapper <-
        bru_mapper(
          INLA::inla.mesh.1d(values),
          indexed = require_indexed
        )
    } else if (all(values == 1)) {
      mapper <- bru_mapper_linear()
    } else {
      mapper <- bru_mapper_linear()
    }
  }

  mapper
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
    if (subcomp[["type"]] %in% c("spde")) {
      if (inherits(subcomp[["model"]]$mesh, "inla.mesh")) {
        subcomp[["mapper"]] <- bru_mapper(subcomp[["model"]]$mesh)
      } else if (inherits(subcomp[["model"]]$mesh, "inla.mesh.1d")) {
        subcomp[["mapper"]] <-
          bru_mapper(subcomp[["model"]]$mesh, indexed = TRUE)
      } else {
        stop(paste0(
          "Unknown SPDE mesh class '",
          paste0(class(subcomp[["model"]]$mesh), collapse = ", "),
          "' for ", label, ". Please specify a mapper manually instead."
        ))
      }
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
    subcomp[["mapper"]] <- bru_mapper(subcomp[["mapper"]])
  } else if (subcomp[["type"]] %in% c("linear", "clinear")) {
    subcomp[["mapper"]] <- bru_mapper_linear()
  } else if (subcomp[["type"]] %in% c("offset")) {
    subcomp[["mapper"]] <- bru_mapper_offset()
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
    mappers <- list(
      make_submapper(
        subcomp_n = subcomp[["n"]],
        subcomp_values = subcomp[["values"]],
        input_values = input_values,
        label = paste0(subcomp[["label"]], " part 1"),
        subcomp_type = subcomp[["type"]],
        subcomp_factor_mapping = subcomp[["factor_mapping"]],
        require_indexed = require_indexed
      )
    )
    mappers[[2]] <- mappers[[1]]
    names(mappers) <- c("u", "v")
    subcomp[["mapper"]] <-
      bru_mapper_collect(mappers, hidden = TRUE)
  } else {
    # No mapper; construct based on input values
    subcomp[["mapper"]] <-
      make_submapper(
        subcomp_n = subcomp[["n"]],
        subcomp_values = subcomp[["values"]],
        input_values = input_values,
        label = subcomp[["label"]],
        subcomp_type = subcomp[["type"]],
        subcomp_factor_mapping = subcomp[["factor_mapping"]],
        require_indexed = require_indexed
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
              ', model = "offset", main = '
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

# MAPPERS ----

#' @details * `bru_mapper.default` adds the "bru_mapper" class and `new_class`
#' to an object. If provided, mapper method functions are added to an environment
#' `.envir` in the object.  The generic methods `ibm_n`, `ibm_n_inla`,
#' `ibm_values`, `ibm_values_inla`,
#' `ibm_amatrix`, `ibm_amatrix_inla`,
#' `ibm_valid_input`, and `ibm_valid_input_inla` look for these
#' functions first,
#' and otherwise call `UseMethod()`.  This is an alternative to using `.S3method()`
#' to register the methods, e.g.
#' `.S3method("ibm_amatrix", "my_mapper_class", ibm_amatrix.my_mapper_class)`.
#' @param new_class If non-`NULL`, this is added at the front of the class definition
#' @param methods, optional `list` of named method definitions; See Details.
#'
#' @export
#' @rdname bru_mapper
bru_mapper.default <- function(mapper,
                               new_class = NULL,
                               methods = NULL,
                               ...) {
  if (any(c(
    "ibm_n", "ibm_n_inla",
    "ibm_values", "ibm_values_inla",
    "ibm_amatrix", "ibm_amatrix_inla",
    "ibm_valid_input", "ibm_valid_input_inla"
  ) %in%
    ...names())) {
    warning(
      paste0(
        "Deprecated use of named method arguments for 'bru_mapper'.\n",
        "Use methods = list(methodname = ..., ...) instead."
      )
    )
    method_names <- intersect(
      c(
        "ibm_n",
        "ibm_values",
        "ibm_amatrix",
        "ibm_inla_subset",
        "ibm_valid_input"
      ),
      ...names()
    )
    methods <- list(...)[method_names]
  }
  if (!inherits(mapper, "bru_mapper")) {
    class(mapper) <- c("bru_mapper", class(mapper))
  }
  if (!is.null(new_class) && !inherits(mapper, new_class)) {
    class(mapper) <- c(new_class, class(mapper))
  }
  if (is.null(mapper[[".envir"]])) {
    mapper$.envir <- new.env()
  }
  if (!is.null(methods)) {
    for (method in names(methods)) {
      if (!is.null(methods[[method]])) {
        assign(method, methods[[method]], envir = mapper[[".envir"]])
      }
    }
  }
  mapper
}


#' Implementation methods for mapper objects
#'
#' A `bru_mapper` sub-class implementation must provide an
#' `ibm_matrix()` method. If the model size 'n' and definition
#' values 'values' are stored in the object itself, default methods are
#' available (see Details). Otherwise the
#' `ibm_n()` and `ibm_values()` methods also need to be provided.
#'
#' @param \dots Arguments passed on to other methods
#' @param mapper A mapper S3 object, normally inheriting from `bru_mapper`
#' @param inla_f logical; when `TRUE` in `ibm_n` and `ibm_values`, these must result in values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' For `ibm_amatrix` methods, it may influence how the input data is interpreted.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#' @seealso [bru_mapper()]
#' @name bru_mapper_methods
#'
#' @details
#' * The default `ibm_n()` method returns a non-null element 'n' from the
#' mapper object, and gives an error if it doesn't exist. If `inla_f=TRUE`,
#' first checks for a 'n_inla' element.
#' @export
#' @rdname bru_mapper_methods
ibm_n.default <- function(mapper, inla_f = FALSE, ...) {
  if (inla_f && !is.null(mapper[["n_inla"]])) {
    mapper[["n_inla"]]
  } else if (!is.null(mapper[["n"]])) {
    mapper[["n"]]
  } else {
    stop("Default 'ibm_n()' method called but mapper doesn't have an 'n' element.")
  }
}


#' @details
#' * The default `ibm_values()` method returns a non-null element
#' 'values' from the mapper object, and `seq_len(ibm_n(mapper))` if
#' it doesn't exist.
#' @export
#' @rdname bru_mapper_methods
ibm_values.default <- function(mapper, inla_f = FALSE, ...) {
  if (inla_f && !is.null(mapper[["values_inla"]])) {
    mapper[["values_inla"]]
  } else if (!inla_f && !is.null(mapper[["values"]])) {
    mapper[["values"]]
  } else {
    seq_len(ibm_n(mapper, inla_f = inla_f))
  }
}

#' @details
#' * The default `ibm_amatrix()` gives an error message.
#' Mapper classes must implement their own `ibm_amatrix` method.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.default <- function(mapper, input, inla_f = FALSE, ...) {
  stop(paste0(
    "Missing implementation of 'ibm_amatrix()' for mapper of class '",
    class(mapper)[1], "'."
  ))
}

#' @details
#' * The default `ibm_inla_subset` method uses
#' the `ibm_values` output to construct the inla subset indexing, passing
#' extra arguments such as `multi` on to the methods (this means it supports
#' both regular vector values and `multi=1` data.frame values).
#' @export
#' @rdname bru_mapper
ibm_inla_subset.default <- function(mapper, ...) {
  values_full <- ibm_values(mapper, inla_f = FALSE, ...)
  values_inla <- ibm_values(mapper, inla_f = TRUE, ...)
  if (is.data.frame(values_full)) {
    subset <- logical(NROW(values_full))
    if (length(subset) > 0) {
      subset[
        plyr::match_df(
          cbind(
            .inla_subset = seq_len(NROW(values_full)),
            values_full
          ),
          values_inla,
          on = names(values_full)
        )[[".inla_subset"]]
      ] <- TRUE
    }
  } else {
    subset <- base::match(values_full, values_inla, nomatch = 0) > 0
  }
  subset
}


#' @details
#' * The default `ibm_valid_input()` method returns an all-TRUE logical vector.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.default <- function(mapper, input, ...) {
  rep(TRUE, NROW(input))
}





## inla.mesh ####

#' @param mesh An `inla.mesh.1d` or `inla.mesh.2d` object to use as a mapper
#' @param indexed logical; If `TRUE`, the `ibm_values()` output will be the
#' integer indexing sequence for the latent variables. If `FALSE`, the knot
#' locations are returned (useful as an interpolator for `rw2` models
#' and similar).
#' Default: `TRUE`
#' @export
#' @rdname bru_mapper
bru_mapper.inla.mesh <- function(mesh, ...) {
  mapper <- list(mesh = mesh)
  class(mapper) <- c("bru_mapper_inla_mesh_2d", "list")
  bru_mapper.default(mapper)
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_inla_mesh_2d <- function(mapper, ...) {
  mapper[["mesh"]]$n
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_inla_mesh_2d <- function(mapper, ...) {
  seq_len(mapper[["mesh"]]$n)
}
#' @param input The values for which to produce a mapping matrix
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_inla_mesh_2d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  if (!is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

## inla.mesh.1d ####

#' @param indexed logical; If `TRUE`, the `ibm_values()` output will be the
#' integer indexing sequence for the latent variables (needed for `spde` models).
#' If `FALSE`, the knot
#' locations are returned (useful as an interpolator for `rw2` models
#' and similar).
#' Default: `NULL`, to force user specification of this parameter
#' @export
#' @rdname bru_mapper
bru_mapper.inla.mesh.1d <- function(mesh, indexed = NULL, ...) {
  if (is.null(indexed)) {
    stop("indexed=TRUE/FALSE needs to be specified to convert inla.mesh.1d to a bru_mapper")
  }
  mapper <- list(
    mesh = mesh,
    indexed = indexed
  )
  class(mapper) <- c("bru_mapper_inla_mesh_1d", "list")
  bru_mapper.default(mapper)
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_inla_mesh_1d <- function(mapper, ...) {
  mapper[["mesh"]]$m
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_inla_mesh_1d <- function(mapper, ...) {
  if (mapper[["indexed"]]) {
    seq_len(mapper[["mesh"]]$m)
  } else {
    mapper[["mesh"]]$loc
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_inla_mesh_1d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

## _index ####

#' @param n Size of a model for `bru_mapper_index`
#' @export
#' @rdname bru_mapper
bru_mapper_index <- function(n = 1L, ...) {
  mapper <- list(
    n = n
  )
  class(mapper) <- c("bru_mapper_index", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_index <- function(mapper, input, ...) {
  ok <- !is.na(input)
  ok[ok] <- (input[ok] >= 1) & (input[ok] <= ibm_n(mapper))
  ok
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_index <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  ok <- which(ibm_valid_input(mapper, input, ...))
  Matrix::sparseMatrix(
    i = ok,
    j = input[ok],
    x = rep(1, length(ok)),
    dims = c(NROW(input), ibm_n(mapper))
  )
}

## _linear ####

#' @export
#' @rdname bru_mapper
bru_mapper_linear <- function(...) {
  mapper <- list()
  class(mapper) <- c("bru_mapper_linear", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_linear <- function(mapper, ...) {
  1L
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_linear <- function(mapper, ...) {
  1.0
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_linear <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  ok <- !is.na(input)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = rep(1, sum(ok)),
    x = input[ok],
    dims = c(NROW(input), ibm_n(mapper))
  )
  A
}

## _matrix ####

#' @param labels Column labels for matrix mappings
#' @export
#' @rdname bru_mapper
bru_mapper_matrix <- function(labels, ...) {
  if (is.factor(labels)) {
    mapper <- list(
      labels = levels(labels)
    )
  } else {
    mapper <- list(
      labels = as.character(labels)
    )
  }
  class(mapper) <- c("bru_mapper_matrix", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_matrix <- function(mapper, ...) {
  length(mapper$labels)
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_matrix <- function(mapper, ...) {
  factor(x = mapper$labels, levels = mapper$labels)
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_matrix <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  } else if (is.matrix(input)) {
    A <- as(input, "Matrix")
  } else if (inherits(input, "Matrix")) {
    A <- input
  } else if (inherits(input, "Spatial")) {
    A <- sp::coordinates(input)
    A <- as(A, "Matrix")
  } else {
    A <- as(input, "Matrix")
  }
  if (ncol(A) != ibm_n(mapper)) {
    stop(paste0(
      "Input to matrix mapper has ", ncol(A),
      " columns but should have ", ibm_n(mapper),
      " columns."
    ))
  }
  A
}



## _factor ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param factor_mapping character; selects the type of factor mapping.
#' * `'contrast'` for leaving out the first factor level.
#' * `'full'` for keeping all levels.
#' @export
#' @rdname bru_mapper
bru_mapper_factor <- function(values, factor_mapping, ...) {
  factor_mapping <- match.arg(factor_mapping, c("full", "contrast"))
  if (is.factor(values)) {
    mapper <- list(
      levels = levels(values),
      factor_mapping = factor_mapping
    )
  } else if (is.character(values)) {
    mapper <- list(
      levels = unique(values),
      factor_mapping = factor_mapping
    )
  } else {
    mapper <- list(
      levels = as.character(sort(unique(values))),
      factor_mapping = factor_mapping
    )
  }
  class(mapper) <- c("bru_mapper_factor", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_factor <- function(mapper, ...) {
  length(ibm_values(mapper))
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_factor <- function(mapper, ...) {
  if (identical(mapper$factor_mapping, "contrast")) {
    mapper$levels[-1L]
  } else {
    mapper$levels
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_factor <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
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



## _offset ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @export
#' @rdname bru_mapper
bru_mapper_offset <- function(...) {
  mapper <- list()
  class(mapper) <- c("bru_mapper_offset", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_offset <- function(mapper, ...) {
  0L
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_offset <- function(mapper, ...) {
  NULL
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_offset <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, 1L))
  }
  ok <- !is.na(input)
  ok[which(ok)] <- (input[ok] > 0L)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = rep(1L, sum(ok)),
    x = input[ok],
    dims = c(NROW(input), 1L)
  )
  A
}


## _multi ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param mappers A list of `bru_mapper` objects
#' @details * `bru_mapper_multi` constructs a kronecker product mapping
#' @export
#' @rdname bru_mapper
bru_mapper_multi <- function(mappers, ...) {
  mapper <- list(
    mappers = mappers,
    n_multi = lapply(mappers, function(x) {
      ibm_n(x)
    }),
    n_inla_multi = lapply(mappers, function(x) {
      ibm_n(x, inla_f = TRUE)
    }),
    values_multi = lapply(mappers, function(x) {
      ibm_values(x)
    }),
    values_inla_multi = lapply(mappers, function(x) {
      ibm_values(x, inla_f = TRUE)
    })
  )
  mapper[["n"]] <- prod(unlist(mapper[["n_multi"]]))
  mapper[["n_inla"]] <- prod(unlist(mapper[["n_inla_multi"]]))
  class(mapper) <- c("bru_mapper_multi", "bru_mapper", "list")
  mapper
}

#' @param multi integer or logical;
#' If positive, the number of levels to recurse in a `bru_multi_mapper`.
#' If `TRUE`, equivalent to `1L`. If `FALSE`, equivalent to `0L`.
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_multi <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (multi > 1) {
    lapply(mapper[["mappers"]], function(x) {
      ibm_n(x, multi = multi - 1)
    })
  } else if (multi == 1) {
    if (inla_f) {
      mapper[["n_inla_multi"]]
    } else {
      mapper[["n_multi"]]
    }
  } else {
    if (inla_f) {
      mapper[["n_inla"]]
    } else {
      mapper[["n"]]
    }
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_multi <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (multi > 1) {
    lapply(mapper[["mappers"]], function(x) {
      ibm_values(x, multi = multi - 1)
    })
  } else if (multi == 1) {
    # Expand indices/values. First sub-mapper varies fastest
    as.data.frame(
      do.call(
        expand.grid,
        c(
          if (inla_f) {
            mapper[["values_inla_multi"]]
          } else {
            mapper[["values_multi"]]
          },
          list(
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
          )
        )
      )
    )
  } else {
    if (inla_f) {
      seq_len(mapper[["n_inla"]])
    } else {
      seq_len(mapper[["n"]])
    }
  }
}


#' @details
#' * `ibm_amatrix` for `bru_mapper_multi` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_multi()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_multi <- function(mapper, input,
                                         inla_f = FALSE, multi = 0L, ...) {
  by_number <- FALSE
  if (is.matrix(input)) {
    if (is.null(colnames(input))) {
      by_number <- TRUE
    }
    input <- as.data.frame(input)
  }
  if (by_number || is.null(names(input))) {
    indexing <- seq_along(mapper[["mappers"]])
  } else {
    indexing <- names(mapper[["mappers"]])
  }
  A <-
    lapply(
      indexing,
      function(x) {
        ibm_amatrix(mapper[["mappers"]][[x]],
          input = input[[x]],
          inla_f = inla_f,
          multi = multi - 1
        )
      }
    )
  if (multi < 1) {
    # Combine the matrices (A1, A2, A3) -> rowkron(A3, rowkron(A2, A1))
    A_ <- A[[1]]
    for (k in seq_len(length(mapper[["mappers"]]) - 1)) {
      A_ <- row_kron(A[[k + 1]], A_)
    }
    return(A_)
  }
  A
}

#' @details
#' * `ibm_valid_input` for `bru_mapper_multi` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_multi()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_multi <- function(mapper, input,
                                             inla_f = FALSE, multi = 0L, ...) {
  by_number <- FALSE
  if (is.matrix(input)) {
    if (is.null(colnames(input))) {
      by_number <- TRUE
    }
    input <- as.data.frame(input)
  }
  if (!is.list(input)) {
    validity <- as.list(rep(FALSE, length(mapper[["mappers"]])))
  } else {
    if (by_number || is.null(names(input))) {
      indexing <- seq_along(mapper[["mappers"]])
    } else {
      indexing <- names(mapper[["mappers"]])
    }
    validity <-
      lapply(
        indexing,
        function(x) {
          ibm_valid_input(mapper[["mappers"]][[x]],
            input = input[[x]],
            inla_f = inla_f,
            multi = 0
          )
        }
      )
  }
  if (multi < 1) {
    # Combine the vectors (v1, v2, v3) -> v1 & v2 & v3
    validity_ <- validity[[1]]
    for (k in seq_len(length(mapper[["mappers"]]) - 1)) {
      validity_ <- validity_ & validity[[k + 1]]
    }
    return(validity_)
  }
  validity
}

#' @return
#' * `[`-indexing a `bru_mapper_multi` extracts a subset
#'   `bru_mapper_multi` object (for drop `FALSE`) or an individual sub-mapper
#'   (for drop `TRUE`, and `i` identifies a single element)
#' @export
#' @param x object from which to extract element(s)
#' @param i indices specifying element(s) to extract
#' @param drop logical;
#' For `[.bru_mapper_multi`, whether to extract an individual mapper when
#' `i` identifies a single element. If `FALSE`, a list of sub-mappers is
#' returned (suitable e.g. for creating a new `bru_mapper_multi` object).
#' Default: `TRUE`
#' @rdname bru_mapper_methods
`[.bru_mapper_multi` <- function(x, i, drop = TRUE) {
  if (is.logical(i)) {
    i <- which(i)
  }
  mapper <- x[["mappers"]][i]
  if (drop) {
    if (length(mapper) == 1) {
      mapper <- mapper[[1]]
    } else if (length(mapper) == 0) {
      mapper <- NULL
    }
  }
  mapper
}

#' @return
#' * The `names()` method for `bru_mapper_multi` returns the names from the
#' sub-mappers list
#' @export
#' @rdname bru_mapper_methods
`names.bru_mapper_multi` <- function(x) {
  names(x[["mappers"]])
}

#' @param value a character vector of up to the same length as the number
#' of mappers in the multi-mapper x
#' @export
#' @rdname bru_mapper_methods
`names<-.bru_mapper_multi` <- function(x, value) {
  names(x[["mappers"]]) <- value
  names(x[["n_multi"]]) <- value
  names(x[["n_inla_multi"]]) <- value
  names(x[["values_multi"]]) <- value
  names(x[["values_inla_multi"]]) <- value
  x
}



## _collect ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param mappers A list of `bru_mapper` objects
#' @param hidden `logical`, set to `TRUE` to flag that the mapper is to be used
#' as a first level input mapper for `INLA::f()` in a model that requires making
#' only the first mapper visible to `INLA::f()` and `INLA::inla.stack()`, such
#' as for "bym2" models, as activated by the `inla_f` argument to `ibm_n`,
#' `ibm_values`, and `ibm_amatrix`. Set to `FALSE` to always access the full
#' mapper, e.g. for `rgeneric` models
#' @details * `bru_mapper_collect` constructs concatenated collection mapping
#' @export
#' @rdname bru_mapper
bru_mapper_collect <- function(mappers, hidden = FALSE, ...) {
  mapper <- list(
    mappers = mappers,
    n_multi = lapply(mappers, function(x) {
      ibm_n(x)
    }),
    values_multi = lapply(mappers, function(x) {
      ibm_values(x)
    }),
    hidden = hidden
  )
  mapper[["n"]] <- sum(unlist(mapper[["n_multi"]]))
  mapper[["values"]] <- seq_len(mapper[["n"]])
  class(mapper) <- c("bru_mapper_collect", "bru_mapper", "list")
  mapper
}

#' @param multi integer or logical;
#' If positive, the number of levels to recurse in a `bru_collect_mapper`.
#' If `TRUE`, equivalent to `1L`. If `FALSE`, equivalent to `0L`.
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_collect <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    if (multi > 0) {
      ibm_n(mapper, multi = multi, ...)
    } else {
      ibm_n(mapper[["mappers"]][[1]])
    }
  } else {
    if (multi > 1) {
      lapply(mapper[["mappers"]], function(x) {
        ibm_n(x, multi = multi - 1)
      })
    } else if (multi == 1) {
      mapper[["n_multi"]]
    } else {
      mapper[["n"]]
    }
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_collect <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    if (multi > 0) {
      ibm_values(mapper, multi = multi, ...)
    } else {
      ibm_values(mapper[["mappers"]][[1]])
    }
  } else {
    if (multi > 1) {
      lapply(mapper[["mappers"]], function(x) {
        ibm_values(x, multi = multi - 1)
      })
    } else if (multi == 1) {
      mapper[["values_multi"]]
    } else {
      mapper[["values"]]
    }
  }
}
#' @details
#' * `ibm_amatrix` for `bru_mapper_collect` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns. When `inla_f=TRUE` and `hidden=TRUE` in
#' the mapper definition, the input format should instead match that of
#' the first, non-hidden, sub-mapper.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_collect <- function(mapper, input, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    input <- list(input)
  }
  if (is.null(names(input))) {
    indexing <- seq_len(length(mapper[["mappers"]]))
    A <-
      lapply(
        indexing,
        function(x) {
          ibm_amatrix(
            mapper[["mappers"]][[x]],
            input = if (x <= length(input)) {
              input[[x]]
            } else {
              NULL
            },
            multi = multi - 1
          )
        }
      )
  } else {
    indexing <- names(mapper[["mappers"]])
    A <-
      lapply(
        indexing,
        function(x) {
          ibm_amatrix(
            mapper[["mappers"]][[x]],
            input = input[[x]],
            multi = multi - 1
          )
        }
      )
  }
  if (multi < 1) {
    # Combine the matrices (A1, A2, A3, ...) -> bdiag(A1, A2, A3, ...)
    A_ <- Matrix::.bdiag(A)
    return(A_)
  }
  A
}

#' @details
#' * `ibm_valid_input` for `bru_mapper_collect` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_collect <- function(mapper, input, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    validity <- ibm_valid_input(mapper[["mappers"]][[1]],
      input = input,
      multi = multi - 1
    )
  } else {
    if (!is.list(input)) {
      validity <- as.list(rep(FALSE, mapper[["mappers"]]))
    } else {
      if (is.null(names(input))) {
        indexing <- seq_along(mapper[["mappers"]])
      } else {
        indexing <- names(mapper[["mappers"]])
      }
      validity <-
        lapply(
          indexing,
          function(x) {
            ibm_valid_input(mapper[["mappers"]][[x]],
              input = input[[x]],
              multi = 0L
            )
          }
        )
    }
    if (multi < 1) {
      # Combine the vectors (v1, v2, v3) -> c(v1, v2, v3)
      validity_ <- do.call(c, validity)
      return(validity_)
    }
  }
  validity
}

#' @return
#' * `[`-indexing a `bru_mapper_collect` extracts a subset
#'   `bru_mapper_collect` object (for drop `FALSE`) or an individual sub-mapper
#'   (for drop `TRUE`, and `i` identifies a single element)
#' @export
#' @param x object from which to extract element(s)
#' @param i indices specifying element(s) to extract
#' @param drop logical;
#' For `[.bru_mapper_collect`, whether to extract an individual mapper when
#' `i` identifies a single element. If `FALSE`, a list of sub-mappers is
#' returned (suitable e.g. for creating a new `bru_mapper_collect` object).
#' Default: `TRUE`
#' @rdname bru_mapper_methods
`[.bru_mapper_collect` <- function(x, i, drop = TRUE) {
  if (is.logical(i)) {
    i <- which(i)
  }
  mapper <- x[["mappers"]][i]
  if (drop) {
    if (length(mapper) == 1) {
      mapper <- mapper[[1]]
    } else if (length(mapper) == 0) {
      mapper <- NULL
    }
  }
  mapper
}

#' @return
#' * The `names()` method for `bru_mapper_collect` returns the names from the
#' sub-mappers list
#' @export
#' @rdname bru_mapper_methods
`names.bru_mapper_collect` <- function(x) {
  names(x[["mappers"]])
}

#' @param value a character vector of up to the same length as x
#' @export
#' @rdname bru_mapper_methods
`names<-.bru_mapper_collect` <- function(x, value) {
  names(x[["mappers"]]) <- value
  names(x[["n_multi"]]) <- value
  names(x[["values_multi"]]) <- value
  x
}



# OPERATORS ----


#' Summarise components
#'
#' @export
#' @method summary component
#' @keywords internal
#' @param object A `component` or `component_list`.
#' @param ... ignored.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'

summary.component <- function(object, ...) {
  result <- list(
    "Label" = object$label,
    "  Main" = sprintf(
      "input '%s',\ttype '%s'",
      deparse(object[["main"]][["input"]][["input"]]),
      object[["main"]][["type"]]
    ),
    "  Group" =
      if (!is.null(object[["group"]][["input"]][["input"]])) {
        sprintf(
          "input '%s',\ttype '%s'",
          deparse(object[["group"]][["input"]][["input"]]),
          object[["group"]][["type"]]
        )
      } else {
        NULL
      },
    "  Replicate" =
      if (!is.null(object[["replicate"]][["input"]][["input"]])) {
        sprintf(
          "input '%s',\ttype '%s'",
          deparse(object[["replicate"]][["input"]][["input"]]),
          object[["replicate"]][["type"]]
        )
      } else {
        NULL
      },
    "  INLA formula" = as.character(object$inla.formula)
  )
  if (!is.null(object[["copy"]])) {
    result[["Copy"]] <- sprintf("Copy of component '%s'", object[["copy"]])
  }
  class(result) <- c("summary_component", "list")
  result
}

#' @export
#' @method summary component_list
#' @param ... ignored.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname summary.component

summary.component_list <- function(object, ...) {
  result <- lapply(
    object,
    summary
  )
  class(result) <- c("summary_component_list", "list")
  result
}

#' @export
#' @param x A 'summary_component' object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.component

print.summary_component <- function(x, ...) {
  for (name in names(x)) {
    if (!is.null(x[[name]])) {
      cat(name, ":", "\t", x[[name]], "\n", sep = "")
    }
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
#' @param data A `data.frame` or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname amatrix_eval

amatrix_eval.component <- function(component, data, ...) {
  val <- input_eval(component, data)
  A <- ibm_amatrix(component$mapper, input = val, ...)

  if (!is.null(val[["weights"]])) {
    A <- val[["weights"]] * A
  }

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
#' the data provided. Since \eqn{x} is a known covariate it will know how to calculate it. The second
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
#' convenient way in inlabru. Spatial covariates are often stored as `SpatialPixelsDataFrame`,
#' `SpatialPixelsDataFrame` or `RasterLayer` objects. These can be provided directly via
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
#' @param component A component.
#' @param data A `data.frame` or Spatial* object of covariates and/or point locations. If null, return the component's map.
#' @param ... Unused.
#' @return An vector or a coordinate matrix
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname input_eval

input_eval.component <- function(component,
                                 data,
                                 ...) {
  val <- list()
  # The names should be a subset of main, group, replicate
  part_names <- names(component[["mapper"]])
  for (part in part_names) {
    val[[part]] <-
      input_eval(
        component[[part]]$input,
        data,
        env = component$env,
        label = part,
        ...
      )
  }
  if (!is.null(component[["weights"]])) {
    val[["weights"]] <-
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
  val
}

#' @export
#' @rdname input_eval

input_eval.component_list <-
  function(components,
           data,
           ...) {
    part <- match.arg(part)
    result <- lapply(components, function(x) input_eval(x, data = data, ...))
    result
  }



#' @export
#' @rdname input_eval

input_eval.bru_input <- function(input, data, env = NULL, label = NULL,
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

  emap <- tryCatch(eval(input$input, envir = envir, enclos = enclos),
    error = function(e) {
    }
  )

  # 0) Eval failed. map everything to 1. This happens for automatically
  #    added Intercept, and for some components that cannot be evaluated
  #    based on the data, e.g. missing columns for a certain likelihood in a
  #    multilikelihood model.
  # 1) If we obtain a function, apply the function to the data
  # 2) If we obtain a SpatialGridDataFrame extract the values of that data
  #    frame at the point locations using the over() function
  # 3) Else we obtain a vector and return as-is. This happens when input
  #    references a column of the data points, or some other complete expression

  # ## Need to handle varying input lengths.
  # ## auto-expansion of scalars needs to happen elsewhere, where it's needed,
  # ## e.g. when using component A-matrices to construct a default linear model
  # ## A matrix.
  #  n <- nrow(as.data.frame(data))

  if (is.null(emap)) {
    if (null.on.fail) {
      return(NULL)
    }
    #    val <- rep(1, n)
    val <- 1
  } else if (is.function(emap)) {
    # Allow but detect failures:
    val <- tryCatch(emap(data),
      error = function(e) {
      }
    )
    if (is.null(val)) {
      # TODO: figure out if we need to do something else in this case.
      # Returning NULL seems like a good way to keep track of these failures
      # that are likely to happen for multilikelihood models; A component only
      # needs to be evaluable for at least one of the likelihoods.
    } else if (identical(as.character(input$input), "coordinates")) {
      tryCatch(
        expr = {
          # Return SpatialPoints instead of a matrix
          val <- as.data.frame(val)
          coordinates(val) <- seq_len(ncol(val))
          # Allow proj4string failures:
          data_crs <- tryCatch(fm_sp_get_crs(data),
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
  } else if (inherits(emap, "formula")) {
    # Allow but detect failures:
    val <- tryCatch(MatrixModels::model.Matrix(emap, data = data, sparse = TRUE),
      error = function(e) {
      }
    )
    if (is.null(val)) {
      # TODO: figure out if we need to do something else in this case.
      # Returning NULL seems like a good way to keep track of these failures
      # that are likely to happen for multilikelihood models; A component only
      # needs to be evaluable for at least one of the likelihoods.
    }
  } else if (inherits(emap, "SpatialGridDataFrame") |
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
    #  } else if ((input$label == "offset") &&
    #    is.numeric(emap) &&
    #    (length(emap) == 1)) {
    #    val <- rep(emap, n)
  } else {
    val <- emap
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
  if ((inherits(emap, "SpatialGridDataFrame") ||
    inherits(emap, "SpatialPixelsDataFrame")) &&
    any(is.na(as.data.frame(val)))) {
    warning(
      paste0(
        "Model input '",
        deparse(input$input),
        "' for '", label, "' returned some NA values.\n",
        "Attempting to fill in spatially by nearest available value.\n",
        "To avoid this basic covariate imputation, supply complete data."
      ),
      immediate. = TRUE
    )

    val <- bru_fill_missing(
      data = emap, where = data, values = val,
      layer = layer, selector = input[["selector"]]
    )
  }

  # Check for NA values.
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
  idx <- ibm_values(component[["mapper"]], inla_f = inla_f, multi = 1)
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
  lapply(components, function(x) ibm_inla_subset(x[["mapper"]], multi = 1))
}
