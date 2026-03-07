# GENERICS ----

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

#' @title Latent model component construction
#'
#' @description
#' Similar to `glm()`, `gam()` and `inla()`, [bru()] models can be constructed
#' via a formula-like syntax, where each latent effect is specified. However, in
#' addition to the parts of the syntax compatible with `INLA::inla`, `bru`
#' components offer additional functionality which facilitates modelling, and
#' the predictor expression can be specified separately, allowing more complex
#' and non-linear predictors to be defined. The formula syntax is just a way to
#' allow all model components to be defined in a single line of code, but the
#' definitions can optionally be split up into separate component definitions.
#' See Details for more information.
#'
#' The `bru_comp` methods all rely on the [bru_comp.character()]
#' method, that defines a model component with a given label/name. The user
#' usually doesn't need to call these methods directly, but can instead supply a
#' formula expression that can be interpreted by the
#' [bru_comp_list.formula()] method, called inside [bru()].
#'
#' @details
#' As shorthand, [bru()] will understand basic additive formulae describing
#' fixed effect models. For instance, the components specification `y ~ x` will
#' define the linear combination of an effect named `x` and an intercept to the
#' response `y` with respect to the likelihood family stated when calling
#' [bru()]. Mathematically, the linear predictor \eqn{\eta}{eta} would be
#' written as
#'
#' \deqn{\eta = \beta x + c,}{eta = beta * x + c,}
#'
#' where:
#'
#' \describe{
#' \item{\eqn{c}{c}}{is the *intercept*}
#' \item{\eqn{x}{x}}{is a *covariate*}
#' \item{\eqn{\beta}{beta}}{is a *latent variable* associated with \eqn{x}{x}
#' and} \item{\eqn{\psi = \beta x}{psi = beta * x}}{ is called the *effect* of
#' \eqn{x}{x}}
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
#' where `f()` is the inla specific function to set up random effects of all
#' kinds. The underlying predictor would again be \eqn{\eta = \beta * x + c} but
#' the result of fitting the model would state `x` as the random effect's name.
#' bru allows rewriting this formula in order to explicitly state the name of
#' the random effect and the name of the associated covariate. This is achieved
#' by replacing `f` with an arbitrary name that we wish to assign to the effect,
#' e.g.
#'
#' \itemize{\item{`components = y ~ psi(x, model = "linear")`.}}
#'
#' Being able to discriminate between \eqn{x} and \eqn{\psi} is relevant because
#' of two functionalities bru offers. The formula parameters of both [bru()] and
#' the prediction method [predict.bru] are interpreted in the mathematical
#' sense. For instance, `predict` may be used to analyze the analytical
#' combination of the covariate \eqn{x} and the intercept using
#'
#' \itemize{\item{`predict(fit, data.frame(x=2)), ~ exp(psi + Intercept)`.}}
#'
#' which corresponds to the mathematical expression \ifelse{html}{\out{e
#' <sup>&#946; + c</sup>}}{\eqn{e^{x \beta + c}}}.
#'
#' On the other hand, predict may be used to only look at a transformation of
#' the latent variable \eqn{\beta_\psi}
#'
#' \itemize{\item{`predict(fit, NULL, ~ exp(psi_latent))`.}}
#'
#' which corresponds to the mathematical expression \ifelse{html}{\out{e
#' <sup>&#946;</sup>}}{\eqn{e^{\beta}}}.
#'
#' @param \dots Parameters passed on to other methods
#'
#' @rdname bru_comp
#' @aliases component
#' @seealso [bru_input()], [summary.bru_comp()]
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
#' cmp <- bru_comp_list(~ myLinearEffectOfX(main = x, model = "linear"))
#' summary(cmp)
#' # Equivalent shortcuts:
#' cmp <- bru_comp_list(~ myLinearEffectOfX(x, model = "linear"))
#' cmp <- bru_comp_list(~ myLinearEffectOfX(x))
#' # Individual component
#' cmp <- bru_comp("myLinearEffectOfX", main = x, model = "linear")
#' summary(cmp)
bru_comp <- function(...) {
  UseMethod("bru_comp")
}

#' @describeIn bru_comp
#' Backwards compatibility alias for `bru_comp()`
bru_component <- function(...) {
  bru_comp(...)
}

#' @export
#' @param object A character label for the component
#' @param main
#' `main` takes an R expression that evaluates to where the latent variables
#' should be evaluated (coordinates, indices, continuous scalar (for rw2 etc)).
#' Arguments starting with weights, group, replicate behave similarly to main,
#' but for the corresponding features of `INLA::f()`. `main` supports `tidy`
#' evaluation expressions, with `.data` and `.env` pronouns.
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
#' for most or all models. Default: NULL, for auto-detection. Models with
#' matrix input for `Cmatrix` or `graph` will use the matrix size.
#' An error is given if it realises it can't figure it out by itself.
#' @param values Specifies for what covariate/index values INLA should build
#' the latent model. Normally generated internally based on the mapping details.
#' (Default: NULL, for auto-determination)
#' @param season.length Passed on to `INLA::f()` for model `"seasonal"`
#' (TODO: check if this parameter is still fully handled)
#' @param nrow,ncol Number of rows and columns for `model` types
#'  "rw2d", "rw2diid", and "matern2d". Default is `NULL`.
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
#' @param marginal May specify a `bm_marginal()` mapper,
#' that is applied before scaling by `weights`.
#' @param A.msk `r lifecycle::badge("deprecated")` and has no effect.
#' @param .envir Evaluation environment. Can later be accessed via
#'   [bru_comp_env()]. Default: `parent.frame()`
#' @param envir_extra Environment for storing `INLA::f()` argument values.
#'   If NULL, an new environment is created. Can be accessed via
#'  [bru_comp_env_extra()]
#'
#' @details The `bru_comp.character` method is inlabru's equivalent to
#'   `INLA`'s `f()` function but adds functionality that is unique to inlabru.
#'
#'   The `main`, `weights`, `group`, `replicate`, and `*_layer` arguments
#'   support `tidy` evaluation expressions, with `.data` and `.env` pronouns.
#'
#' @rdname bru_comp
#' @aliases bru_comp
#' @aliases bru_component
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla()) {
#'   # As an example, let us create a linear component. Here, the component is
#'   # called "myEffectOfX" while the covariate the component acts on is called
#'   # "x":
#'
#'   cmp <- bru_comp("myEffectOfX", main = x, model = "linear")
#'   summary(cmp)
#'
#'   # A more complicated component:
#'   cmp <- bru_comp("myEffectOfX",
#'     main = x,
#'     model = INLA::inla.spde2.matern(fm_mesh_1d(1:10))
#'   )
#'
#'   # Compound fixed effect component, where x and z are in the input data.
#'   # The formula will be passed on to MatrixModels::model.Matrix:
#'   cmp <- bru_comp("eff", ~ -1 + x:z, model = "fixed")
#'   summary(cmp)
#' }
#' }
#'
bru_comp.character <- function(object,
                               # Main model parameters
                               main = NULL, # This must be kept as 1st arg.
                               weights = NULL, # This must be 2nd arg.
                               ..., # Prevent partial matching
                               model = NULL,
                               mapper = NULL,
                               main_layer = NULL,
                               main_selector = NULL,
                               n = NULL,
                               values = NULL,
                               season.length = NULL,
                               nrow = NULL,
                               ncol = NULL,
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
                               # Marginal transformation
                               marginal = NULL,
                               A.msk = deprecated(),
                               .envir = parent.frame(),
                               envir_extra = NULL) {
  # INLA models:
  # itypes = c(linear, iid, mec, meb, rgeneric, rw1, rw2, crw2, seasonal, besag,
  # besag2, bym, bym2, besagproper, besagproper2, fgn, fgn2, ar1, ar1c, ar, ou,
  # generic, generic0, generic1, generic2, generic3, spde, spde2, spde3, iid1d,
  # iid2d, iid3d, iid4d, iid5d, 2diid, z, rw2d, rw2diid, slm, matern2d, copy,
  # clinear, sigm, revsigm, log1exp, logdist)
  #
  # Supported:
  # btypes = c("const", "offset", "factor_full", "factor_contrast", "linear",
  # "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")

  # The label
  label <- object

  if (is.null(model) && is.null(copy)) {
    # TODO: may need a special marker to autodetect factor and matrix variables,
    #       which can only be autodetected in a data-aware pass.
    model <- "linear"
  }

  # Force evaluation of explicit inputs
  force(values)

  include_weights <- !missing(weights) &&
    (!is.null(substitute(weights)))
  include_group <- !missing(group) ||
    !missing(ngroup) || !missing(group_mapper)
  include_repl <- !missing(replicate) ||
    !missing(nrep) || !missing(replicate_mapper)
  if (include_group) {
    if (!is.null(substitute(group)) &&
      !identical(deparse(substitute(group)), "1L")) {
      if (is.null(control.group)) {
        control.group <- INLA::inla.set.control.group.default()
      }
      group_model <- control.group$model
    } else {
      group_model <- "exchangeable"
    }
  }

  if (is.null(envir_extra)) {
    envir_extra <- new.env(parent = .envir)
  }

  if (is.null(n)) {
    arg_names <- names(list(...))
    if ("Z" %in% arg_names) {
      Zmat <- list(...)[["Z"]]
      if (is.matrix(Zmat) || inherits(Zmat, "Matrix")) {
        n <- nrow(Zmat) + ncol(Zmat)
      }
    } else if ("Cmatrix" %in% arg_names) {
      Cmat <- list(...)[["Cmatrix"]]
      if (is.matrix(Cmat) || inherits(Cmat, "Matrix")) {
        n <- nrow(Cmat)
      }
    } else if ("graph" %in% arg_names) {
      n <- INLA::inla.read.graph(list(...)[["graph"]], size.only = TRUE)
    }
  }

  # Convert weights, ngroup, and nrep to bru_mapper info
  if (include_group) {
    if (!is.null(ngroup)) {
      if (!is.null(group_mapper)) {
        stop("At most one of 'ngroup' and 'group_mapper' should be supplied.")
      }
      group_mapper <- bm_index(ngroup)
      ngroup <- NULL
    }
  }
  if (include_repl) {
    if (!is.null(nrep)) {
      if (!is.null(replicate_mapper)) {
        stop("At most one of 'nrep' and 'replicate_mapper' should be supplied.")
      }
      replicate_mapper <- bm_index(nrep)
      nrep <- NULL
    }
  }

  inp <- new_bru_input(
    {{ main }},
    label = label,
    layer = {{ main_layer }},
    selector = main_selector
  )
  subcomp <- list(
    main = bru_subcomp(
      input = inp,
      mapper = mapper,
      n = n,
      model = model,
      values = values,
      season.length = season.length,
      nrow = nrow,
      ncol = ncol
    )
  )
  if (include_group) {
    subcomp$group <-
      bru_subcomp(
        input = new_bru_input(
          {{ group }},
          label = glue("{label}.group"),
          layer = {{ group_layer }},
          selector = group_selector
        ),
        mapper = group_mapper,
        n = NULL,
        model = group_model
      )
  }
  if (include_repl) {
    subcomp$replicate <-
      bru_subcomp(
        input = new_bru_input(
          {{ replicate }},
          label = glue("{label}.repl"),
          layer = {{ replicate_layer }},
          selector = replicate_selector
        ),
        mapper = replicate_mapper,
        n = NULL,
        model = "iid"
      )
  }
  if (include_weights) {
    weights_input <- new_bru_input(
      {{ weights }},
      label = glue("{label}.weights"),
      layer = {{ weights_layer }},
      selector = weights_selector
    )
    subcomp$weights <- ibm_input_set(bm_scale(), weights_input)
  }

  # Default component (to be filled)
  component <- c(
    list(
      label = label,
      inla.formula = NULL
    ),
    subcomp,
    list(
      copy = copy,
      marginal =
        if (is.null(marginal)) {
          NULL
        } else if (ibm_input_available(marginal)) {
          marginal
        } else {
          ibm_input_new(
            marginal,
            NULL,
            label = glue("{label}.marginal"),
            layer = NULL,
            selector = NULL
          )
        },
      env = .envir,
      env_extra = envir_extra
    )
  )

  # Main bit
  # Construct a call to the f function from the parameters provided
  # Ultimately, this call will be converted to the model formula presented to
  # INLA
  fcall <- sys.call()
  fcall[[1]] <- "f"
  fcall[[2]] <- as.symbol(label)
  names(fcall)[2] <- ""
  # 'main' and 'weights' are the only regular parameter allowed to be nameless,
  # and only if they are the first parameters (position 3 and 4 in fcall)
  if ((length(fcall) >= 3L) &&
    (is.null(names(fcall)) || (names(fcall)[3] %in% c("", NA_character_)))) {
    names(fcall)[3] <- "main"
  }
  if ((length(fcall) >= 4L) &&
    (is.null(names(fcall)) || (names(fcall)[4] %in% c("", NA_character_)))) {
    names(fcall)[4] <- "weights"
  }
  unnamed_arguments <- which(names(fcall[-c(1, 2)]) %in% c("", NA_character_))
  if (length(unnamed_arguments) > 0) {
    # Without this check, R gives the error
    #   'In str2lang(s) : parsing result not of length one, but 0'
    # in the INLA call instead, which isn't very informative.
    stop(glue(
      "Unnamed arguments detected in component '{label}'.\n",
      "  Only 'main' and 'weights' parameters may be unnamed.\n",
      "  Unnamed arguments at position(s) ",
      glue_collapse(unnamed_arguments, sep = ", ")
    ))
  }

  # Special and general cases:
  if (component$main$type %in% c("offset", "const")) {
    # The offset is included either automatically for ~ . linear models,
    # or explicitly by name in the predictor expression, so no INLA formula
    # component is needed.
    component$inla.formula <- as.formula("~ .", env = .envir)
    component$main$mapper <-
      ibm_input_set(bm_const(), component$main$mapper)
    component$group <- NULL
    component$replicate <- NULL
    # Add scalable multi-mapper
    component[["mapper"]] <-
      bm_pipe(
        c(
          list(
            core = bm_multi(list(
              main = component[["main"]][["mapper"]]
            ))
          ),
          if (include_weights) {
            list(scale = component[["weights"]])
          } else {
            NULL
          }
        )
      )
  } else {
    if (!is.null(copy)) {
      # Store copy-model name or object in the environment
      #    model_name <- paste0("BRU_", label, "_copy_model")
      #    fcall[["copy"]] <- component$copy # as.symbol(model_name)
      #    bru_comp_env_extra(component, model_name) <- component$copy
    } else {
      # Store model name or object in the environment
      model_name <- glue("BRU_{label}_main_model")
      fcall[["model"]] <- as.symbol(model_name)
      bru_comp_env_extra(component, model_name) <- component$main$model
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
      "weights_selector",
      "marginal"
    ))]

    # Replace arguments that will be evaluated by a mapper
    # TODO: make a more general system
    suffixes <- list(
      "group" = "group",
      "replicate" = "repl"
    )
    for (arg in names(suffixes)) {
      if (arg %in% names(fcall)) {
        fcall[[arg]] <- as.symbol(glue("{label}.{suffixes[[arg]]}"))
      }
    }

    if (is.null(copy)) {
      # These values were previously initialised by an INLA::f call, but ::f
      # should only be called by INLA, and never by inlabru, since it sets up
      # temporary files that will not be removed, and also requires data not
      # available at this point!
      if (!is.null(component$main[["n"]])) {
        fcall[["n"]] <- component$main$n
      }
      if (!is.null(season.length)) {
        fcall[["season.length"]] <- season.length
      }
      if (!is.null(component$main[["nrow"]])) {
        fcall[["nrow"]] <- component$main$nrow
      }
      if (!is.null(component$main[["ncol"]])) {
        fcall[["ncol"]] <- component$main$ncol
      }

      # Make sure 'values' is setup properly.
      if (is.null(component$main$values)) {
        fcall <- fcall[!("values" %in% names(fcall))]
      } else {
        values_name <- glue("BRU_{label}_values")
        fcall[["values"]] <- as.symbol(values_name)
        bru_comp_env_extra(component, values_name) <- component$main$values
      }

      # Setup factor precision parameter
      if (component$main$type %in% c("factor", "fixed")) {
        if (is.null(fcall[["hyper"]])) {
          # TODO: allow configuration of the precision via prec.linear
          fixed_hyper_name <- glue("BRU_{label}_main_fixed_hyper")
          fcall[["hyper"]] <- as.symbol(fixed_hyper_name)
          bru_comp_env_extra(component, fixed_hyper_name) <-
            list(prec = list(
              initial = log(INLA::inla.set.control.fixed.default()$prec),
              fixed = TRUE
            ))
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

  class(component) <- "bru_comp"

  component <- bru_comp_update_mapper(component)

  component <- bru_compat_pre_2_14_bru_comp(component)

  component
}

bru_compat_pre_2_14_bru_comp <- function(comp) {
  if (!isTRUE(bru_options_get("bru_compat_pre_2_14_enable"))) {
    return(comp)
  }
  if (!inherits(comp, "bru_comp")) {
    class(comp) <- c("bru_comp", class(comp))
  }
  comp$main$input <- ibm_input_get(comp$main$mapper)
  if (!is.null(comp$group)) {
    comp$group$input <- ibm_input_get(comp$group$mapper)
  }
  if (!is.null(comp$replicate)) {
    comp$replicate$input <- ibm_input_get(comp$replicate$mapper)
  }
  comp
}


#' Methods for inlabru component lists
#'
#' Constructor methods for inlabru component lists. Syntax details are given in
#' [bru_comp()].
#'
#' @param \dots Parameters passed on to other methods. Also see Details.
#' @family component constructors
#' @param object The object to operate on
#' @param lhoods A [bru_obs_list] object
#' @param .envir An evaluation environment for non-formula input
#' @export
#' @rdname bru_comp_list
#' @aliases bru_component_list
#' @aliases component_list
bru_comp_list <- function(object,
                          ...,
                          .envir = parent.frame()) {
  UseMethod("bru_comp_list")
}

#' @describeIn bru_comp_list Convert a component formula
#' into a `bru_comp_list` object
#'
#' @export
#' @family component constructors
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' # As an example, let us create a linear component. Here, the component is
#' # called "myLinearEffectOfX" while the covariate the component acts on is
#' # called "x". Note that a list of components is returned because the
#' # formula may define multiple components
#'
#' eff <- bru_comp_list(~ myLinearEffectOfX(main = x, model = "linear"))
#' summary(eff[[1]])
#' # Equivalent shortcuts:
#' eff <- bru_comp_list(~ myLinearEffectOfX(x, model = "linear"))
#' eff <- bru_comp_list(~ myLinearEffectOfX(x))
#' # Individual component
#' eff <- bru_comp("myLinearEffectOfX", main = x, model = "linear")
bru_comp_list.formula <- function(object,
                                  ...,
                                  lhoods = NULL,
                                  .envir = parent.frame()) {
  if (!is.null(environment(object))) {
    .envir <- environment(object)
  }
  # Automatically add Intercept and remove -1 unless -1
  # or -Intercept is in components formula
  object <- auto_intercept(object)

  code <- bru_formula_to_bru_obs_code(object)
  parsed <- lapply(code, function(x) rlang::parse_expr(x))
  components <- lapply(
    parsed,
    function(component.expression) {
      rlang::eval_bare(component.expression,
        env = .envir
      )
    }
  )
  environment(object) <- .envir
  bru_comp_list(components, lhoods = lhoods, .envir = .envir)
}


#' @describeIn bru_comp_list Combine a list of components, component lists,
#'   and/or component formulas into a single `bru_comp_list` object
#' @param inputs A tree-like list of component input evaluations,
#' from [bru_input.bru_obs_list()].
#' @export
bru_comp_list.list <- function(object,
                               ...,
                               lhoods = NULL,
                               .envir = parent.frame(),
                               inputs = NULL) {
  # Maybe the list has been given an environment?
  if (!is.null(environment(object))) {
    .envir <- environment(object)
  } else if (is.null(.envir)) {
    # Later code needs an actual environment
    .envir <- new.env()
  }
  is_comp <- vapply(object, function(x) inherits(x, "bru_comp"), TRUE)
  if (!all(is_comp)) {
    is_comp_list <- vapply(object, function(x) {
      inherits(x, "bru_comp_list")
    }, TRUE)
    if (!all(is_comp_list)) {
      object <- lapply(seq_along(object), function(k) {
        bru_comp_list(object[[k]], .envir = .envir)
      })
    }
    object <- unlist(object, recursive = FALSE)
  }
  object <- structure(
    object,
    class = c("bru_comp_list", "list")
  )
  environment(object) <- .envir

  object <- set_list_names(object, tag = "label", priority = "name")
  if (anyDuplicated(names(object))) {
    stop(paste0(
      "Duplicated component labels detected: ",
      glue_collapse(
        glue(
          "'{sort(unique(names(object)[duplicated(names(object))]))}'"
        ),
        sep = ", "
      )
    ))
  }

  if (!is.null(lhoods)) {
    lhoods <- bru_used_update(lhoods, names(object))
    object <- add_mappers(object, lhoods = lhoods, inputs = inputs)
  }
  object
}

#' @export
#' @describeIn bru_comp_list Place a single `bru_comp` object into a
#'   `bru_comp_list` object.
bru_comp_list.bru_comp <- function(object,
                                   ...,
                                   .envir = parent.frame()) {
  bru_comp_list(list(object), .envir = .envir, ...)
}

#' @export
#' @describeIn bru_comp_list Make sure a `bru_comp_list` object is fully
#'   configured.
bru_comp_list.bru_comp_list <- function(object,
                                        ...,
                                        .envir = parent.frame()) {
  bru_comp_list(list(object), .envir = .envir, ...)
}


#' @export
#' @describeIn bru_comp_list The `...` arguments should be `bru_comp_list`
#' objects. The environment from the first argument will be applied to the
#' resulting `bru_comp_list`.
`c.bru_comp_list` <- function(...) {
  bru_comp_list(list(...), .envir = environment(list(...)[[1]]))
}

#' @export
#' @describeIn bru_comp_list The `...` arguments should be `component`
#'   objects from [bru_comp()]. The environment from the first argument
#'   will be applied to the resulting `bru_comp_list`.
`c.bru_comp` <- function(...) {
  bru_comp_list(list(...), .envir = environment(list(...)[[1]]))
}

#' @export
#' @param x `bru_comp_list` object from which to extract a sub-list
#' @param i indices specifying elements to extract
#' @rdname bru_comp_list
`[.bru_comp_list` <- function(x, i) {
  env <- environment(x)
  object <- NextMethod()
  class(object) <- c("bru_comp_list", "list")
  environment(object) <- env
  object
}


#' Get/set component environment data
#'
#' @description Get or set data in the component's `env_extra` and `env`
#'   environments. If `name` is `NULL`, the entire environment is returned.
#' @param x A `bru_comp` object
#' @param name The name of the variable to get or set. Names
#'   prefixed by `"BRU_"` are reserved for internal use and should not be set
#'   by users or external package authors.
#' @param value The value to set for `name` in `env_extra` (for the setter
#'   function).
#' @return For the getter, either the entire environment or the
#'   value of `name` in the environment. For the setters, the modified
#'   `bru_comp` object.
#' @export
#' @rdname bru_comp_env_extra
#' @seealso bru_comp
#' @examples
#' if (bru_safe_inla()) {
#'   cmp <- bru_comp("x", x)
#'   as.list(bru_comp_env_extra(cmp))
#'
#'   cmp <- bru_comp("x", x, model = "fixed")
#'   as.list(bru_comp_env_extra(cmp))
#'
#'   bru_comp_env(cmp)
#' }
#'
bru_comp_env_extra <- function(x, name = NULL) {
  if (is.null(name)) {
    return(x[["env_extra"]])
  }
  get(name, envir = x[["env_extra"]], inherits = FALSE)
}

#' @describeIn bru_comp_env_extra Set data in the component's `env_extra`
#' environment
#' @export
`bru_comp_env_extra<-` <- function(x, name, value) {
  assign(name, value, envir = x[["env_extra"]])
  x
}

#' @describeIn bru_comp_env_extra Get the component's `env` environment, or an
#'   element from that environment.
#' Note that in most cases this environment is the global R environment.
#' @export
bru_comp_env <- function(x, name = NULL) {
  if (is.null(name)) {
    return(x[["env"]])
  }
  get(name, envir = x[["env"]], inherits = FALSE)
}

#' @describeIn bru_comp_env_extra Set data in the component's `env` environment.
#' Note that in most cases this environment is the global R environment, so
#' modifying it is normally not recommended.
#' @export
`bru_comp_env<-` <- function(x, name, value) {
  assign(name, value, envir = x[["env"]])
  x
}


bru_comp_update_mapper <- function(component) {
  # Core component multi-mapper
  core_mapper <- list(
    core = bm_multi(
      lapply(
        component[intersect(
          c("main", "group", "replicate"),
          names(component)
        )],
        function(subcomp) {
          subcomp[["mapper"]]
        }
      )
    )
  )
  # Add marginal transformation mapper
  if (!is.null(component[["marginal"]])) {
    marginal_mapper <- list(marginal = component[["marginal"]])
  } else {
    marginal_mapper <- NULL
  }
  # Add weights-mapper
  if (!is.null(component[["weights"]])) {
    weights_mapper <- list(scale = component[["weights"]])
  } else {
    weights_mapper <- NULL
  }
  # Combine the mappers
  component[["mapper"]] <-
    bm_pipe(
      c(
        core_mapper,
        marginal_mapper,
        weights_mapper
      )
    )
  component
}


#' @title Equip components with mappers
#' @description Equip component(s) with mappers for subcomponents that do not
#'   have predefined mappers. When needed, the data in `lhoods` is used to
#'   determine the appropriate mapper(s).
#' @param component A [component] object
#' @param lhoods A [bru_obs_list] object
#' @param inputs A precomputed list of inputs, as returned by
#'   [bru_input.bru_comp_list()]
#' @param lh_data A list of data object, one for each likelihood in `lhoods`
#'   that is used to determine the mapper(s). If `NULL`, the data is
#'   or inputs are used instead.
#' @return A `component` object with completed mapper information
#' @examples
#' \dontrun{
#' if (interactive() && bru_safe_inla()) {}
#' }
#' @rdname add_mappers
#' @keywords internal
#' @export

add_mappers.bru_comp <- function(component,
                                 lhoods,
                                 inputs = NULL,
                                 lh_data = NULL,
                                 ...) {
  if (is.null(lh_data)) {
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
    inputs <- inputs[keep_lh]

    lh_data <- lapply(lh, function(x) x[["data"]])
  }

  component$main <- add_mapper(
    component$main,
    label = component$label,
    data = lh_data,
    inputs = lapply(
      inputs,
      function(x) x[[component$label]][["core"]][["main"]]
    ),
    env = component$env,
    require_indexed = FALSE
  )
  if (!is.null(component[["group"]])) {
    component$group <- add_mapper(
      component$group,
      label = component$label,
      data = lh_data,
      inputs = lapply(
        inputs,
        function(x) x[[component$label]][["core"]][["group"]]
      ),
      env = component$env,
      require_indexed = TRUE
    )
  }
  if (!is.null(component[["replicate"]])) {
    component$replicate <- add_mapper(
      component$replicate,
      label = component$label,
      data = lh_data,
      inputs = lapply(
        inputs,
        function(x) x[[component$label]][["core"]][["replicate"]]
      ),
      env = component$env,
      require_indexed = TRUE
    )
  }
  component <- bru_comp_update_mapper(component)

  fcall <- component$fcall

  # Set ngroup and nrep defaults
  if (!is.null(component[["group"]])) {
    if (is.null(component$group[["n"]])) {
      fcall[["ngroup"]] <- 1L
    } else {
      fcall[["ngroup"]] <- component$group$n
    }
  }
  if (!is.null(component[["replicate"]])) {
    if (is.null(component$replicate[["n"]])) {
      fcall[["nrep"]] <- 1L
    } else {
      fcall[["nrep"]] <- component$replicate$n
    }
  }

  if (is.null(component$main$values)) {
    fcall <- fcall[!("values" %in% names(fcall))]
  } else {
    values_name <- glue("BRU_{component$label}_values")
    fcall[["values"]] <- as.symbol(values_name)
    bru_comp_env_extra(component, values_name) <- component$main$values
  }

  if (!(component[["main"]][["type"]] %in% c("offset", "const"))) {
    # Update the formula that will be presented to INLA
    component$inla.formula <-
      as.formula(
        glue(
          "~ . + ",
          glue_collapse(deparse(fcall), sep = "\n")
        ),
        env = component$env
      )
  }

  component$fcall <- fcall

  component
}

#' @param components A `bru_comp_list` object
#' @export
#' @rdname add_mappers
add_mappers.bru_comp_list <- function(components,
                                      lhoods,
                                      inputs = NULL,
                                      lh_data = NULL,
                                      ...) {
  if (is.null(inputs)) {
    inputs <- bru_input(lhoods, components = components)
  }
  is_copy <- vapply(components, function(x) !is.null(x[["copy"]]), TRUE)
  for (k in which(!is_copy)) {
    components[[k]] <- add_mappers(
      components[[k]], lhoods,
      inputs = inputs, lh_data = lh_data
    )
  }
  for (k in which(is_copy)) {
    if (is.null(components[[k]][["copy"]])) {
      stop(paste0(
        "Internal error: copy model detected, ",
        "but no copy information available."
      ))
    }
    if (is.null(components[[components[[k]][["copy"]]]])) {
      stop(glue(
        "Could not find component '{components[[k]][['copy']]}'",
        " to use as copy for component '{components[[k]][['label']]}'."
      ))
    }
    components[[k]]$mapper <-
      components[[components[[k]][["copy"]]]][["mapper"]]
  }
  components
}

bru_subcomp <- function(input = NULL,
                        mapper = NULL,
                        model = NULL,
                        n = NULL,
                        values = NULL,
                        season.length = NULL,
                        nrow = NULL,
                        ncol = NULL) {
  type <- model
  factor_mapping <- NULL
  if (inherits(model, "inla.spde")) {
    type <- "spde"
  } else if (inherits(model, "inla.rgeneric")) {
    type <- "rgeneric"
  } else if (inherits(model, "inla.cgeneric")) {
    type <- "cgeneric"
  } else if (is.character(model)) {
    if (identical(model, "factor")) {
      model <- "factor_contrast"
      bru_log_warn(
        paste0(
          "Deprecated model 'factor'. Please use 'factor_full' or ",
          "'factor_contrast' instead.\n",
          "Defaulting to 'factor_contrast' that matches the old 'factor' model."
        )
      )
    }
    if (model %in%
      c(
        "clinear", "sigm", "revsigm",
        "log1exp", "logdist"
      )
    ) {
      type <- "specialnonlinear"
    } else if (model %in%
      c(
        "rw2d", "rw2diid", "matern2d"
      )
    ) {
      type <- "inlalattice"
      if (is.null(n)) {
        n <- nrow * ncol
      } else {
        if (n != nrow * ncol) {
          stop(
            "n = ",
            n,
            " must be equal to nrow * ncol = ",
            nrow,
            " * ",
            ncol,
            " = ",
            nrow * ncol,
            " for inla lattice model '",
            model,
            "'."
          )
        }
      }
    } else if (identical(model, "factor_full")) {
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
  if (is.null(mapper)) {
    mapper <- bm_autodetect()
  }
  if (!inherits(mapper, "bru_mapper")) {
    stop(
      "Unknown mapper class '",
      glue_collapse(class(mapper), sep = ", "),
      "'"
    )
  }
  mapper <- ibm_input_set(mapper, input)
  subcomponent <-
    structure(
      list(
        mapper = mapper,
        model = model,
        type = type,
        n = n,
        nrow = nrow,
        ncol = ncol,
        values = values,
        season.length = season.length,
        factor_mapping = factor_mapping
      ),
      class = "bru_subcomp"
    )

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
    bru_safe_sp(force = TRUE)
    if (!all(is_spatial)) {
      stop(paste0(
        "Inconsistent spatial/non-spatial input. ",
        "Unable to infer mapper information."
      ))
    }
    inconsistent_crs <- FALSE
    inp_crs <- lapply(inp, fm_CRS)
    crs_info <- lapply(inp_crs, fm_wkt)
    null_crs <- vapply(crs_info, is.null, logical(1))
    inconsistent_crs <-
      (length(unique(unlist(crs_info))) > 1) ||
        (any(null_crs) && !all(null_crs))
    if (inconsistent_crs) {
      stop(paste0(
        "Inconsistent spatial CRS information. ",
        "Unable to infer mapper information."
      ))
    }
    inp_values <- unique(do.call(
      rbind,
      lapply(
        inp,
        function(x) sp::coordinates(x)
      )
    ))
    n_values <- nrow(inp_values)
  } else if (any(is_sfc)) {
    if (!all(is_sfc)) {
      stop(paste0(
        "Inconsistent spatial/non-spatial input. ",
        "Unable to infer mapper information."
      ))
    }
    inconsistent_crs <- FALSE
    inp_crs <- lapply(inp, fm_crs)
    null_crs <- vapply(inp_crs, is.na, logical(1))
    inconsistent_crs <-
      (length(unique(inp_crs)) > 1) ||
        (any(null_crs) && !all(null_crs))
    if (inconsistent_crs) {
      stop(paste0(
        "Inconsistent spatial crs information. ",
        "Unable to infer mapper information."
      ))
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
    inp_values <-
      unique.matrix(cbind(1, do.call(rbind, inp)))[, -1, drop = FALSE]
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
      stop(paste0(
        "Unable to create automatic mapper. List data at this level requires ",
        "an explicit mapper definition."
      ))
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

  list(
    inp_values = inp_values,
    n_values = n_values,
    is_list
  )
}


# @param inputs list with precomputed inputs, for each lhoods element;
# @param data list of data objects, one for each lhoods element
# lapply(full_inputs,
#        function(x) x[[component$label]][["core"]][[subcomponent_label]])
add_mapper <- function(subcomp, label, data = NULL, env = NULL,
                       inputs = NULL,
                       require_indexed = FALSE) {
  if (inherits(subcomp[["mapper"]], "bm_autodetect") &&
    ibm_input_available(subcomp[["mapper"]])) {
    mapper_input <- ibm_input_get(subcomp[["mapper"]])
  } else {
    mapper_input <- NULL
  }
  if (is.null(subcomp[["mapper"]]) ||
    inherits(subcomp[["mapper"]], "bm_autodetect")) {
    if (!inherits(subcomp[["model"]], "character")) {
      subcomp[["mapper"]] <- bru_get_mapper_safely(subcomp[["model"]])
      if (inherits(subcomp[["mapper"]], "bru_mapper") &&
        ibm_input_available(subcomp[["mapper"]])) {
        mapper_input <- ibm_input_get(subcomp[["mapper"]])
      }
    }
  }
  if (!is.null(subcomp[["mapper"]]) &&
    !inherits(subcomp[["mapper"]], "bm_autodetect")) {
    if (!inherits(subcomp[["mapper"]], "bru_mapper")) {
      stop(glue(
        "Unknown mapper of type '",
        glue_collapse(class(subcomp[["mapper"]]), sep = ", "),
        "' for {label}"
      ))
    }
  } else if (is.null(mapper_input)) {
    stop("This code should be unused?")
    # Used for automatic group and replicate mappers
    subcomp[["mapper"]] <- make_mapper(
      subcomp,
      label,
      input_values = 1,
      strict = TRUE,
      require_indexed = require_indexed
    )
  } else {
    if (!is.null(data) || !is.null(inputs)) {
      if (is.null(inputs) || (length(inputs) == 0)) {
        if (!is.null(data) && (length(data) > 0)) {
          inputs <- lapply(
            data,
            function(dat) {
              bru_input(
                mapper_input,
                data = dat,
                .envir = env,
                label = mapper_input$label,
                null.on.fail = TRUE
              )
            }
          )
        } else {
          # Component not directly used in any likelihood.
          # Attempt to evaluate with no data;
          # useful for intercept-like components only used via the
          # *_latent technique.
          inputs <- list(
            bru_input(
              mapper_input,
              data = NULL,
              .envir = env,
              label = mapper_input$label,
              null.on.fail = TRUE
            )
          )
        }
      }
      # Check for
      # 1) All NULL; Deprecated unless input is NULL. Since version 2.1.14,
      #              intercepts should be notated explicitly with label(1)
      # 2) Some NULL; exclude NULL results
      # TODO: Check for vector/matrix/coordinate inconsistency
      null.results <- vapply(inputs, function(x) is.null(x), TRUE)
      if (all(null.results)) {
        msg <- glue(
          "All covariate evaluations for '{label}'",
          " are NULL; an intercept component was likely intended.\n",
          "  Implicit latent intercept component specification is ",
          "deprecated since version 2.1.14.\n",
          "  Use explicit notation '+ {label}(1)' instead",
          if (identical(label, "Intercept")) {
            " (or '+1' for '+ Intercept(1)')"
          },
          "."
        )
        bru_log_message(msg)
        unique_inputs <- list(
          inp_values = 1,
          n_values = 1
        )
      } else {
        if (any(null.results)) {
          inp_ <- inputs[!null.results]
        } else {
          inp_ <- inputs
        }

        unique_inputs <- make_unique_inputs(inp_, allow_list = TRUE)
      }
      if (sum(unlist(unique_inputs$n_values)) < 1) {
        subcomp$n <- 1
        subcomp$values <- NULL
        inp_values <- NULL
      }
      subcomp[["mapper"]] <- make_mapper(
        subcomp,
        label,
        input_values = unique_inputs$inp_values,
        strict = TRUE,
        require_indexed = require_indexed
      )
    }
  }
  if (!is.null(subcomp[["mapper"]]) &&
    !inherits(subcomp[["mapper"]], "bm_autodetect")) {
    if (!is.null(mapper_input) && !ibm_input_available(subcomp[["mapper"]])) {
      subcomp[["mapper"]] <- ibm_input_set(subcomp[["mapper"]], mapper_input)
    }
    # Check internal consistency of user specified n and the mapper:
    mapper_n <- ibm_n(subcomp[["mapper"]], inla_f = TRUE)
    if (!is.null(subcomp[["n"]]) &&
      subcomp[["n"]] != mapper_n) {
      stop(glue(
        "Size mismatch, n={subcomp[['n']]}  != ",
        "ibm_n(inla_f = TRUE)={mapper_n}, ",
        " in mapper for label '{label}'."
      ))
    }
    subcomp[["n"]] <- mapper_n
    subcomp[["values"]] <- ibm_values(subcomp[["mapper"]], inla_f = TRUE)
  }
  subcomp
}


make_values <- function(subcomp_n, subcomp_values, input_values, label) {
  if (!is.null(subcomp_n)) {
    if (!is.null(subcomp_values)) {
      warning(
        "Both 'n' and 'values' provided. Ignoring 'values'",
        immediate. = TRUE
      )
    }
    return(seq_len(subcomp_n))
  }
  if (!is.null(subcomp_values)) {
    return(subcomp_values)
  }
  if (!is.null(input_values)) {
    return(input_values)
  }

  stop(glue("No mapper, no n, and no values given for '{label}'"))
}

make_submapper <- function(subcomp_n,
                           subcomp_nrow_ncol,
                           subcomp_values,
                           input_values,
                           label,
                           subcomp_type,
                           subcomp_factor_mapping,
                           require_indexed,
                           allow_interpolation = TRUE) {
  if (identical(subcomp_type, "inlalattice")) {
    if (is.null(subcomp_nrow_ncol)) {
      stop(
        glue("nrow and ncol must be provided for inla lattice model '{label}'.")
      )
    }
    n <- prod(subcomp_nrow_ncol)
    if (!is.null(subcomp_n)) {
      stopifnot(n == subcomp_n)
    }
    subcomp_n <- n
    return(
      bm_index(n = subcomp_n)
    )
  }

  values <- make_values(subcomp_n, subcomp_values, input_values, label)

  if (is.factor(values) ||
    is.character(values) ||
    (!is.null(subcomp_type) && (subcomp_type %in% "factor"))) {
    return(
      bm_factor(
        values,
        factor_mapping = subcomp_factor_mapping,
        indexed = require_indexed
      )
    )
  }

  values <- sort(unique(values), na.last = NA)
  if (length(values) > 1) {
    if (allow_interpolation) {
      return(
        bru_mapper(
          fm_mesh_1d(values),
          indexed = require_indexed
        )
      )
    }
    mapper <- bm_index(n = ceiling(max(values)))
    return(mapper)
  }
  if (all(values == 1)) {
    if (require_indexed) {
      return(bm_index(n = 1L))
    }
    return(bm_linear())
  }
  if (require_indexed) {
    return(
      bm_factor(values,
        factor_mapping = "full",
        indexed = TRUE
      )
    )
  }

  bm_linear()
}


#' @title
#' Extract mapper information from INLA model component objects
#'
#' @description
#' The component definitions will automatically attempt to extract mapper
#' information from any model object by calling the generic `bru_get_mapper`.
#' Any class method implementation should return a [bru_mapper] object suitable
#' for the given latent model.
#'
#' @details
#' Before implementing your own `bru_get_mapper` method, check if there
#' is already a general method available that handles your model class, such as
#' [bru_get_mapper.inla.spde()], [bru_get_mapper.inla.rgeneric()], and
#' [bru_get_mapper.inla.cgeneric()].
#'
#' @param model A model component object
#' @param \dots Arguments passed on to other methods
#' @return A [bru_mapper] object defined by the model component
#' @seealso [bru_mapper] for mapper constructor methods, and
#' the individual mappers for specific implementation details.
#' @family mappers
#' @export
#' @examples
#' if (bru_safe_inla()) {
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


#' @describeIn bru_get_mapper Extract an indexed mapper for
#' the `model$mesh` object contained in the model object,
#' which is assumed to be of a class supporting relevant `fmesher` methods.
#' @export
bru_get_mapper.inla.spde <- function(model, ...) {
  bm_fmesher(model[["mesh"]])
}

#' @describeIn bru_get_mapper Returns the mapper given by a pre-computed mapper,
#'   or an index mapper mapping the size of the model graph.
#'
#'   The easiest method to define a mapper for an `inla.rgeneric` model is to
#'   store the mapper in the object.
#'   Alternatively, define your model using a subclass and define a
#'   `bru_get_mapper.<subclass>` method that should return the
#'   corresponding `bru_mapper` object.
#'
#' The order of precedence for the mapper construction when calling
#' `bru_get_mapper(model)` has the following precedence:
#'
#' 1. `bru_get_mapper.<subclass>`, if `model` has a subclass, otherwise
#' 2. `model[["mapper"]]` if that is `NULL`, and otherwise
#' 3. [bm_index()] using the size of the graph returned by `model[["f"]][["n"]]`
#' @export
bru_get_mapper.inla.rgeneric <- function(model, ...) {
  if (!is.null(model[["mapper"]])) {
    return(model[["mapper"]])
  }
  if (!is.null(model[["f"]][["n"]])) {
    return(bm_index(n = model[["f"]][["n"]]))
  }
  NULL
}

#' @describeIn bru_get_mapper Works the same as the method of `inla.rgeneric`
#' @export
bru_get_mapper.inla.cgeneric <- function(model, ...) {
  if (!is.null(model[["mapper"]])) {
    return(model[["mapper"]])
  }
  if (!is.null(model[["f"]][["n"]])) {
    return(bm_index(n = model[["f"]][["n"]]))
  }
  NULL
}

#' @describeIn bru_get_mapper Tries to call the `bru_get_mapper`,
#' and returns `NULL` if it fails (e.g. due to no available class method).
#' If the call succeeds and returns non-`NULL`, it checks that the object
#' inherits from the `bru_mapper` class, and gives an error if it does not.
#' @export
bru_get_mapper_safely <- function(model, ...) {
  m <- tryCatch(
    bru_get_mapper(model, ...),
    error = function(e) {
      NULL
    }
  )
  if (!is.null(m) && !inherits(m, "bru_mapper")) {
    stop(paste0(
      "The bru_get_mapper method for model class '",
      glue_collapse(class(model), sep = ", "),
      "' did not return a bru_mapper object."
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
  mapper <- subcomp[["mapper"]]
  if ((is.null(mapper) ||
    inherits(mapper, "bm_autodetect")) &&
    !inherits(subcomp[["model"]], "character")) {
    mapper <- bru_get_mapper_safely(subcomp[["model"]])
  }
  if (!is.null(mapper) && !inherits(mapper, "bm_autodetect")) {
    if (inherits(mapper, "bru_mapper")) {
      return(mapper)
    }
    stop(paste0(
      "Unknown mapper of type '",
      glue_collapse(class(mapper), sep = ", "),
      "' for ", label
    ))
  }
  if (subcomp[["type"]] %in% c("linear", "specialnonlinear")) {
    return(bm_linear())
  }
  if (subcomp[["type"]] %in% c("offset", "const")) {
    return(bm_const())
  }
  if (subcomp[["type"]] %in% c("fixed")) {
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
      stop(glue(
        "Need to specify at least one of values (labels), input, or n for ",
        "'{subcomp[['label']]}'",
        " of component '{label}'."
      ))
    }
    return(bm_matrix(labels))
  }
  allow_interpolation <- TRUE
  stopifnot(bru_safe_inla(multicore = TRUE))
  inla_model <- INLA::inla.models()[["latent"]][[subcomp[["model"]]]]
  if (is.null(inla_model) ||
    !is.integer(inla_model[["aug.factor"]]) ||
    (inla_model[["aug.factor"]] == 1L)) {
    the_mapper <- subcomp[["mapper"]]
    the_input <- ibm_input_get(the_mapper)
    mapper_names <- the_input[["label"]]
    labels <- the_input[["label"]]
  } else {
    if (subcomp[["model"]] %in% c("bym", "bym2")) {
      mapper_names <- c("u", "v")
      allow_interpolation <- FALSE
    } else {
      mapper_names <- glue("u{seq_len(inla_model[['aug.factor']])}")
    }
    the_mapper <- subcomp[["mapper"]]
    the_input <- ibm_input_get(the_mapper)
    labels <- glue("{the_input[['label']]}_{mapper_names}")
  }

  mappers <- lapply(
    labels,
    function(lab) {
      make_submapper(
        subcomp_n = subcomp[["n"]],
        subcomp_nrow_ncol = c(subcomp[["nrow"]], subcomp[["ncol"]]),
        subcomp_values = subcomp[["values"]],
        input_values = input_values,
        label = lab,
        subcomp_type = subcomp[["type"]],
        subcomp_factor_mapping = subcomp[["factor_mapping"]],
        require_indexed = require_indexed,
        allow_interpolation = allow_interpolation
      )
    }
  )

  if (length(mapper_names) > 1L) {
    names(mappers) <- mapper_names
    return(bm_collect(mappers, hidden = TRUE))
  }
  mappers[[1]]
}

#' Convert components to R code
#'
#' @aliases bru_formula_to_bru_obs_code
#' @keywords internal
#' @param components A [formula] describing latent model components.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'

bru_formula_to_bru_obs_code <- function(components, add = "") {
  fname <- "inlabru:::bru_comp.character"
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
      codes[[k]] <- glue(
        '{fname}("{label}", main = {label}, model = "linear"{add})'
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
            glue("{label}("),
            glue('{fname}("{label}", model = "const", main = '),
            code,
            fixed = TRUE
          )
      } else {
        codes[[k]] <- sub(
          glue("{label}("),
          glue('{fname}("{label}", '),
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
#' @keywords internal
#' @param object Object to be summarised.
#' @param \dots Passed on to other summary methods.
#' @param depth The depth of which to expand the component mapper.
#' Default `Inf`, to traverse the entire mapper tree.
#' @param verbose logical; If `TRUE`, includes more details of the
#' component definitions. When `FALSE`, only show basic component
#' definition information.  Default `TRUE`.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [bru_comp()], [summary.bru_input()]
#' @export
#' @method summary bru_comp

summary.bru_comp <- function(object, ..., depth = Inf, verbose = TRUE) {
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
            glue("{x} = {object[[x]][['type']]}")
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
      if (is.null(object[["copy"]])) {
        NULL
      } else {
        glue("(={object[['copy']]})")
      },
      ": ",
      paste0(
        unlist(lapply(
          c("main", "group", "replicate", "marginal", "weights"),
          function(x) {
            if (!is.null(object[[x]])) {
              format(object[[x]],
                verbose = verbose,
                ...,
                label.override = x
              )
            }
          }
        )),
        collapse = ", "
      )
    )
    result$Summary <- Summary
  }
  class(result) <- "summary_component"
  result
}

#' @export
#' @method summary bru_comp_list
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @rdname summary.bru_comp

summary.bru_comp_list <- function(object, verbose = TRUE, ...) {
  result <- lapply(
    object,
    function(x) {
      summary(x, verbose = verbose, ...)
    }
  )
  for (cp in which(vapply(result, function(x) !is.null(x$Copy_of), FALSE))) {
    result[[cp]]$Type <- result[[result[[cp]][["Copy_of"]]]][["Type"]]
  }
  class(result) <- c("summary_bru_comp_list", "list")
  result
}

#' @export
#' @method print bru_comp
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_comp
print.bru_comp <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

#' @export
#' @method print bru_comp_list
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_comp
print.bru_comp_list <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}


#' @export
#' @method format bru_subcomp
#' @param label.override character; If not `NULL`, use this label instead of
#' the object's label.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_comp
format.bru_subcomp <- function(x,
                               verbose = TRUE,
                               ...,
                               label.override = NULL) {
  text <- format(
    ibm_input_get(x[["mapper"]]),
    verbose = verbose,
    ...,
    label.override = label.override,
    type = if (verbose) NULL else x[["type"]]
  )
  text
}

#' @export
#' @method summary bru_subcomp
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_comp
summary.bru_subcomp <- function(object,
                                verbose = TRUE,
                                ...,
                                label.override = NULL) {
  text <- format(
    object,
    verbose = verbose,
    ...,
    label.override = label.override
  )
  res <- structure(
    list(
      text = text
    ),
    class = "summary_bru_subcomp"
  )
  res
}

#' @export
#' @method print bru_subcomp
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_comp
print.bru_subcomp <- function(x,
                              verbose = TRUE,
                              ...,
                              label.override = NULL) {
  text <- format(x, verbose = verbose, ..., label.override = label.override)
  cat(text, "\n", sep = "")
  invisible(x)
}


#' @export
#' @param x Object to be printed.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_comp

print.summary_component <- function(x, ...) {
  if (is.null(x[["Summary"]])) {
    for (name in names(x)) {
      if (!is.null(x[[name]])) {
        if (name %in% "Label") {
          cat(glue("Label:\t{x[[name]]}"), "\n")
        } else if (name %in% "Mapper") {
          cat(glue("  Map:\t{glue::glue_collapse(x[[name]])}"), "\n")
        } else {
          cat(glue("  {name}:\t{glue::glue_collapse(x[[name]])}"), "\n")
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
#' @rdname summary.bru_comp

print.summary_bru_comp_list <- function(x, ...) {
  lapply(x, print)
  invisible(x)
}

#' @export
#' @rdname summary.bru_comp

print.summary_bru_subcomp <- function(x, ...) {
  if (!is.null(x[["text"]])) {
    cat(x[["text"]], "\n", sep = "")
  }
  invisible(x)
}
