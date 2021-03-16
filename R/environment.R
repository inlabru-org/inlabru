#' @include 0_inlabru_envir.R

#' @title Get access to the internal environment
#' @details The environment is defined in aaaaa.R which is loaded first.
#' @keywords internal
bru_env_get <- function() {
  pkg_envir <- parent.env(environment())
  envir <- get0(".inlabru_envir", envir = pkg_envir)
  if (!is.environment(envir)) {
    stop("Something went wrong: cannot find internal .inlabru_envir environment.")
  }
  envir
}


# Documentation would clash with base .onLoad documentation
# @title Initialise log storage and global options
# @param libname a character string giving the library directory where the
#   package defining the namespace was found.
# @param pkgname a character string giving the name of the package.
# @aliases namespace_hooks
# @keywords internal
# @rdname namespace_hooks
.onLoad <- function(libname, pkgname) {
  bru_log_reset()
  bru_log_message("inlabru loaded", allow_verbose = FALSE)
  bru_log_message("Clear override options", allow_verbose = FALSE)
  bru_options_reset()
}




# inlabru log methods ----

#' @title inlabru log message methods
#' @description Resets the inlabru log object
#' @details `bru_log_reset()` clears the log contents.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname bru_log

bru_log_reset <- function() {
  envir <- bru_env_get()
  envir$log <- character(0)
  invisible(NULL)
}


#' @param pretty logical; If `TRUE`, return a single string with the log
#' messages separated and terminated by line feeds, suitable for `cat(...)`.
#' If `FALSE`, return the raw log as a vector of strings, suitable for
#' `cat(..., sep = "\n")`. Default: `FALSE`
#' @return `bru_log_get` RETURN_VALUE
#' @export
#' @rdname bru_log

bru_log_get <- function(pretty = FALSE) {
  if (pretty) {
    paste0(paste0(bru_env_get()[["log"]], collapse = "\n"), "\n")
  } else {
    bru_env_get()[["log"]]
  }
}


#' @param ... Zero or more objects passed on to [`base::.makeMessage()`]
#' @param domain Domain for translations, passed on to [`base::.makeMessage()`]
#' @param appendLF logical; whether to add a newline to the message. Only
#'   used for verbose output.
#' @param verbosity numeric value describing the verbosity level of the message
#' @param allow_verbose Whether to allow verbose output. Must be set to FALSE
#' until the options object has been initialised.
#' @param verbose logical, numeric, or `NULL`; local override for verbose
#' output. If `NULL`, the global option `bru_verbose` or default value is used.
#' If `FALSE`, no messages are printed. If `TRUE`, messages with `verbosity`
#' \eqn{\le 1}{<=1}
#' are printed. If numeric, messages with `verbosity` \eqn{\le}{<=} `verbose` are
#' printed.
#' @param verbose_store Same as `verbose`, but controlling what messages are
#' stored in the global log object. Can be controlled via the `bru_verbose_store`
#' with [bru_options_set()].
#' @return
#' * `bru_log_message` OUTPUT_DESCRIPTION
#' @details
#' * `bru_log_message` DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname bru_log

bru_log_message <- function(..., domain = NULL, appendLF = TRUE,
                            verbosity = 1,
                            allow_verbose = TRUE, verbose = NULL,
                            verbose_store = NULL) {
  if (allow_verbose) {
    if ((!is.null(verbose) && (verbose >= verbosity)) ||
      (is.null(verbose) &&
        bru_options_get("bru_verbose", include_default = TRUE) >= verbosity)) {
      message(..., domain = domain, appendLF = appendLF)
    }
  }
  if ((!is.null(verbose_store) && (verbose_store >= verbosity)) ||
    !allow_verbose ||
    (is.null(verbose_store) &&
      bru_options_get("bru_verbose_store", include_default = TRUE) >= verbosity)) {
    envir <- bru_env_get()
    envir$log <- c(
      envir$log,
      .makeMessage(Sys.time(), ": ", ...,
        domain = domain,
        appendLF = FALSE
      )
    )
  }
  invisible()
}



# Options methods ----

#' @title Create or update an options objects
#' @description Create a new options object, or merge information from several
#'   objects.
#'
#'   The `_get`, `_set`, and `_reset` functions operate on a global
#'   package options override object. In many cases, setting options in
#'   specific calls to [bru()] is recommended instead.
#' @param ... A collection of named options, optionally including one or more
#' [`bru_options`] objects. Options specified later override the previous options.
#' @return `bru_options()` returns an `bru_options` object.
#' @section Valid options:
#' For `bru_options` and `bru_options_set`, recognised options are:
#' \describe{
#' \item{bru_verbose}{logical or numeric; if `TRUE`, log messages of `verbosity`
#' \eqn{\le 1} are printed by [bru_log_message()]. If numeric, log messages
#' of
#' verbosity \eqn{\le}`bru_verbose` are printed.
#' For line search details, set `bru_verbose=2` or `3`.
#' Default: 0, to not print any messages}
#' \item{bru_verbose_store}{logical or numeric; if `TRUE`, log messages of
#' `verbosity` \eqn{\le 1} are stored by [bru_log_message()]. If numeric,
#' log messages of verbosity \eqn{\le} are stored. Default: Inf, to store all messages.}
#' \item{bru_run}{If TRUE, run inference. Otherwise only return configuration needed
#'   to run inference.}
#' \item{bru_max_iter}{maximum number of inla iterations}
#' \item{bru_initial}{An `inla` object returned from previous calls of
#'   `INLA::inla`, [bru()] or [lgcp()], or a list of named vectors of starting
#'   values for the latent variables. This will be used as a
#'   starting point for further improvement of the approximate posterior.}
#' \item{bru_int_args}{List of arguments passed all the way to the
#' integration method `ipoints` and `int.polygon` for 'cp' family models;
#' \describe{
#' \item{method}{"stable" or "direct". For "stable" (default) integration points
#' are aggregated to mesh vertices.}
#' \item{nsub1}{Number of integration points per knot interval in 1D. Default 30.}
#' \item{nsub2}{Number of integration points along a triangle edge for 2D. Default 9.}
#' \item{nsub}{Deprecated parameter that overrides `nsub1` and `nsub2` if set. Default `NULL`.}
#' }
#' }
#' \item{bru_method}{List of arguments controlling the iterative inlabru method:
#' \describe{
#' \item{taylor}{Either 'legacy' (for the pre-2.1.15 method) or 'pandemic'
#' (default, from version 2.1.15).}
#' \item{search}{Either 'all' (default), to use all available line search
#' methods, or one or more of
#' 'finite' (reduce step size until predictor is finite),
#' 'contract' (decrease step size until trust hypersphere reached)
#' 'expand' (increase step size until no improvement),
#' 'optimise' (fast approximate error norm minimisation).
#' To disable line search, set to an empty vector. Line search is not
#' available for `taylor="legacy"`.}
#' \item{factor}{Numeric, \eqn{> 1} determining the line search step scaling
#' multiplier. Default \eqn{(1 + \sqrt{5})/2}{(1+sqrt(5))/2}.}
#' }
#' }
#' \item{`inla()` options}{
#' All options not starting with `bru_` are passed on to `inla()`, sometimes
#' after altering according to the needs of the inlabru method.

#' Warning:
#'   Due to how inlabru currently constructs the `inla()` call, the `mean`,
#'   `prec`, `mean.intercept`, and `prec.intercept` settings in
#'   `control.fixed` will have no
#'   effect. Until a more elegant alternative has been implemented, use explicit
#'   `mean.linear` and `prec.linear` specifications in each
#'   `model="linear"` component instead.
#' }
#' }
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Combine global and user options:
#'   options1 <- bru_options(bru_options_get(), bru_verbose = TRUE)
#'   # Create a proto-options object in two equivalent ways:
#'   options2 <- as.bru_options(bru_verbose = TRUE)
#'   options2 <- as.bru_options(list(bru_verbose = TRUE))
#'   # Combine options objects:
#'   options3 <- bru_options(options1, options2)
#' }
#' }
#' @export
#' @rdname bru_options

bru_options <- function(...) {
  new_bru_options <- function() {
    options <- list()
    class(options) <- c("bru_options", "list")
    options
  }

  input_options <- list(...)

  options <- new_bru_options()
  for (k in seq_along(input_options)) {
    if (inherits(input_options[[k]], "bru_options")) {
      options <- modifyList(options, input_options[[k]], keep.null = TRUE)
    } else {
      options <- modifyList(options, input_options[k], keep.null = TRUE)
    }
  }

  options
}

#' @param x An object to be converted to an `bru_options` object.
#' @return For `as.bru_options()`, `NULL` or no input returns an empty
#' `bru_options` object, a `list` is converted via `bru_options(...)`,
#' and `bru_options` input is passed through. Other types of input generates
#' an error.
#'
#' @export
#' @rdname bru_options

as.bru_options <- function(x = NULL) {
  if (inherits(x, "bru_options")) {
    x
  } else if (is.null(x)) {
    bru_options()
  } else if (is.list(x)) {
    do.call(bru_options, x)
  } else {
    stop("Unable coerce object to 'bru_options'")
  }
}

#' @return `bru_options_default()` returns an `bru_options` object containing
#'   default options.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname bru_options

bru_options_default <- function() {
  bru_options(
    # inlabru options
    bru_verbose = 0,
    bru_verbose_store = Inf,
    bru_max_iter = 10,
    bru_run = TRUE,
    bru_int_args = list(method = "stable", nsub1 = 30, nsub2 = 9),
    bru_method = list(
      taylor = "pandemic",
      search = "all",
      factor = (1 + sqrt(5)) / 2,
      stop_at_max_rel_deviation = 0.01
    ),
    # bru_initial: NULL
    # inla options
    E = 1,
    Ntrials = 1,
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
    control.inla = list(int.strategy = "auto"),
    control.fixed = list(expand.factor.strategy = "inla")
  )
}


bru_options_deprecated <- function(args) {
  stopifnot(inherits(args, "bru_options"))
  cl <- class(args)
  deprecated_args <- c(
    mesh = "",
    run = "bru_run",
    max.iter = "bru_max_iter",
    result = "bru_initial",
    bru_result = "bru_initial",
    int.args = "int.args"
  )
  names_args <- names(args)
  deprecated_args <-
    deprecated_args[intersect(names_args, names(deprecated_args))]
  if (any(nzchar(names(deprecated_args)) == 0)) {
    warning(paste0(
      "Ignoring deprecated global options '",
      paste0(names(deprecated_args)[nzchar(names(deprecated_args)) == 0], collapse = "', '"),
      "'."
    ))
    names_args <- setdiff(
      names_args,
      names(deprecated_args)[nzchar(names(deprecated_args)) == 0]
    )
    deprecated_args <- deprecated_args[nzchar(names(deprecated_args)) > 0]
    args <- args[names_args]
  }
  if (length(deprecated_args) > 0) {
    warning(paste0(
      "Converting deprecated global option(s) '",
      paste0(names(deprecated_args), collapse = "', '"),
      "' to new option(s) '",
      paste0(deprecated_args, collapse = "', '"),
      "'."
    ))
    names_args[names_args %in% names(deprecated_args)] <- deprecated_args
    names(args) <- names_args
  }
  class(args) <- cl
  args
}


#' Additional bru options
#'
#' Construct a `bru_options` object including the default and global options,
#' and converting deprecated option names.
#'
#' @param \dots Options passed on to [as.bru_options()]
#'
#' @aliases bru_call_options
#' @export
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#'
#' opts <- bru_call_options()
#'
#' # Print them:
#' opts
#' }
#'
bru_call_options <- function(...) {
  opt <- as.bru_options(...)
  opt <- bru_options_deprecated(opt)
  opt <- bru_options(bru_options_get(), opt)
  bru_options_check(opt)
  opt
}


# Extract bru options
bru_options_bru <- function(options) {
  stopifnot(inherits(options, "bru_options"))
  cl <- class(options)
  options <- options[grepl("^bru_", names(options))]
  class(options) <- cl
  options
}

# Extract non-bru options
bru_options_inla <- function(options) {
  stopifnot(inherits(options, "bru_options"))
  cl <- class(options)
  options <- options[!grepl("^bru_", names(options))]
  class(options) <- cl
  options
}



#' @details `bru_options_check` checks for valid contents of an `bru_options`
#' object
#' @param options An `bru_options` object to be checked
#' @param ignore_null Ignore missing or NULL options.
#' @return `bru_options_check()` returns a `logical`; `TRUE` if the object
#'   contains valid options for use by other functions
#' @details `bru_options_check()` produces warnings for invalid options.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   bru_options_check(bru_options(bru_max_iter = "text"))
#' }
#' }
#' @export
#' @rdname bru_options

bru_options_check <- function(options, ignore_null = TRUE) {
  options <- as.bru_options(options)
  ok <- TRUE
  are_null <- vapply(options, is.null, TRUE)
  if (!ignore_null) {
    disallowed_null <-
      intersect(
        names(options),
        c("bru_max_iter")
      )
    disallowed_null <- disallowed_null[are_null[disallowed_null]]
    if (length(disallowed_null) > 0) {
      warning(paste0(
        paste0("'", disallowed_null, "'", collapse = ", "),
        " should not be set to NULL."
      ))
    }
  }
  for (name in names(options)[!are_null]) {
    # Check valid max_iter
    opt <- options[[name]]
    if (name == "bru_max_iter") {
      if (!is.numeric(opt) || !(opt > 0)) {
        ok <- FALSE
        warning("'bru_max_iter' should be a positive integer.")
      }
    }
  }

  ok
}



#' @param name Either `NULL`, or single option name string, or character vector
#'  or list with option names,
#'   Default: NULL
#' @param include_default logical; If `TRUE`, the default options are included
#'   together with the global override options. Default: `TRUE`
#' @return `bru_options_get` returns either an [`bru_options`] object, for
#'   `name == NULL`, the contents of single option, if `name` is a options name
#'   string, or a named list of option contents, if `name` is a list of option
#'   name strings.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname bru_options

bru_options_get <- function(name = NULL, include_default = TRUE) {
  if (include_default) {
    default <- bru_options_default()
  } else {
    default <- bru_options()
  }
  global <- bru_env_get()$options
  options <- bru_options(default, global)
  if (is.null(name)) {
    return(options)
  }
  if (is.list(name)) {
    mget(unlist(name), as.environment(options))
  } else {
    options[[name]]
  }
}


#' @details `bru_options_set()` is used to set global package options.
#' @return `bru_options_set()` returns a copy of the global override options,
#' invisibly (as `bru_options_get(include_default = FALSE)`).
#' @seealso [bru_options()], [bru_options_default()], [bru_options_get()]
#' @param .reset For `bru_options_set`, logical indicating if the global override
#' options list should be emptied before setting the new option(s).
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   bru_options_set(
#'     bru_verbose = TRUE,
#'     verbose = TRUE
#'   )
#' }
#' }
#' @export
#' @rdname bru_options

bru_options_set <- function(..., .reset = FALSE) {
  envir <- bru_env_get()
  if (.reset) {
    envir$options <- bru_options(...)
  } else {
    envir$options <- bru_options(envir$options, ...)
  }
  invisible(bru_options_get(include_default = FALSE))
}

#' @details `bru_options_reset()` clears the global option overrides.
#' @export
#' @rdname bru_options

bru_options_reset <- function() {
  envir <- bru_env_get()
  envir$options <- bru_options()
  invisible(bru_options_get(include_default = FALSE))
}



#' @title Print inlabru options
#' @param object A [bru_options] object to be summarised
#' @param x A `summary_bru_options` object to be printed
#' @param legend logical; If `TRUE`, include explanatory text, Default: `TRUE`
#' @param include_global logical; If `TRUE`, include global override options
#' @param include_default logical; If `TRUE`, include default options
#' @param ... Further parameters, currently ignored
#'
#' @examples
#' if (interactive()) {
#'   options <- bru_options(verbose = TRUE)
#'
#'   # Don't print options only set in default:
#'   print(options, include_default = FALSE)
#'
#'   # Only include options set in the object:
#'   print(options, include_default = FALSE, include_global = FALSE)
#' }
#' @method summary bru_options
#' @export
#' @rdname summary.bru_options
summary.bru_options <- function(object,
                                legend = TRUE,
                                include_global = TRUE,
                                include_default = TRUE,
                                ...) {
  traverse <- function(combined, default, global, object) {
    result <- list()
    for (name in sort(names(combined))) {
      if (is.list(combined[[name]])) {
        result[[name]] <- list(
          is_list = TRUE,
          value = traverse(
            combined[[name]], default[[name]],
            global[[name]], object[[name]]
          )
        )
      } else {
        result[[name]] <- list(
          is_list = FALSE,
          value = if (is.null(combined[[name]])) {
            "NULL"
          } else {
            combined[[name]]
          },
          origin =
            if (
              !is.null(default[[name]]) &&
                identical(default[[name]], combined[[name]])
            ) {
              "default"
            } else if (
              !is.null(global[[name]]) &&
                identical(global[[name]], combined[[name]])
            ) {
              "global"
            } else if (
              !is.null(object[[name]]) &&
                identical(object[[name]], combined[[name]])
            ) {
              "user"
            } else {
              "unknown"
            }
        )
      }
    }
    result
  }
  if (include_default) {
    default <- bru_options_default()
  } else {
    default <- bru_options()
  }
  if (include_global) {
    global <- bru_options_get(include_default = FALSE)
  } else {
    global <- bru_options()
  }
  combined <- bru_options(default, global, object)

  if (legend) {
    legend <- c(
      "default = value from the default options",
      "global  = value from the global override object",
      "user    = value from the user override object"
    )
  } else {
    legend <- NULL
  }
  result <- list(
    legend = legend,
    value = traverse(combined, default, global, object)
  )
  class(result) <- c("summary_bru_options", "list")
  result
}

#' @export
#' @rdname summary.bru_options
print.summary_bru_options <- function(x, ...) {
  traverse <- function(tree, prefix = "") {
    for (name in sort(names(tree))) {
      if (tree[[name]]$is_list) {
        cat(paste0(prefix, name, " =\n"))
        traverse(
          tree[[name]]$value,
          prefix = paste0(prefix, "\t")
        )
      } else {
        cat(paste0(
          prefix,
          name, " =\t",
          tree[[name]]$value,
          "\t(",
          tree[[name]]$origin,
          ")\n"
        ))
      }
    }
  }

  if (!is.null(x[["legend"]])) {
    cat("Legend:\n")
    cat(paste0("  ", x[["legend"]], collapse = "\n"))
  }
  cat("Options for inlabru:\n")
  traverse(x[["value"]], prefix = "  ")
  invisible(x)
}







# Old log methods ----

#' @title inlabru log messages
#'
#' @description Retrieve, add, and/or print log messages
#'
#' @param txt character; log message.
#' @param verbose logical; if `TRUE`, print the log message on screen with
#' `message(txt)`. Default: `bru_options_get("bru_verbose")`
#' @details The log message is stored if the log is active, see
#' [bru_log_active()]
#' @return `bru_log` invisibly returns the added log message.
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and Finn Lindgren
#' \email{finn.lindgren@@gmail.com}

bru_log <- function(txt, verbose = NULL) {
  if (is.null(verbose)) {
    verbose <- bru_options_get("bru_verbose")
  }
  bru_log_message(txt, verbose = verbose)
}


#' @param activation logical; whether to activate (`TRUE`) or deactivate
#' (`FALSE`) the inlabru logging system. Default: NULL, to keep the current
#' activation state
#' @return `bru_log_active` returns the previous activation state
#' @export
#' @examples
#' code_runner <- function() {
#'   oa <- bru_log_active(TRUE)
#'   on.exit(bru_log_active(oa))
#'   bru_log("Test message")
#' }
#' bru_log_active()
#' code_runner()
#' cat(bru_log_get())
#' bru_log_active()
#' @rdname bru_log

bru_log_active <- function(activation = NULL) {
  current <- bru_options_get("bru_verbose_store")
  if (!is.null(activation)) {
    bru_options_set("bru_verbose_store" = activation)
  }
  current
}


# Old options methods ----

#' @describeIn inlabru-deprecated Use `bru_option_get` instead.
#' @param name character; an option name
#' @export
iinla.getOption <- function(name = NULL) {
  if (identical(name, "iinla.verbose")) {
    new_option <- "bru_verbose"
  } else {
    new_option <- name
  }
  .Deprecated(
    paste0('bru_options_get("', new_option, '")'),
    old = deparse(sys.call(sys.parent()))
  )
  res <- bru_options_get(new_option)
  res
}


#' @describeIn inlabru-deprecated Use `bru_option_set` instead.
#' @export
iinla.setOption <- function(...) {
  iinla.setOption.core <- function(option = c(
                                     "control.compute",
                                     "control.inla",
                                     "iinla.verbose",
                                     "control.fixed"
                                   ),
                                   value,
                                   option_value_alternating = TRUE) {
    if (identical(option, "iinla.verbose")) {
      new_option <- "bru_verbose"
    } else {
      new_option <- option
    }
    if (option_value_alternating) {
      .Deprecated(
        paste0("bru_options_get(", new_option, " = value)"),
        old = paste0('iinla.setOption("', option, '", value)')
      )
    } else {
      .Deprecated(
        paste0("bru_options_get(", new_option, " = value)"),
        old = paste0("iinla.setOption(", option, " = value)")
      )
    }
    setting <- list(value)
    names(setting) <- new_option
    do.call(bru_options_set, setting)
  }
  called <- list(...)
  len <- length(names(called))
  if (len > 0L) {
    for (i in 1L:len) {
      do.call(
        iinla.setOption.core,
        args = list(
          names(called)[i],
          called[[i]],
          option_value_alternating = FALSE
        )
      )
    }
  } else {
    iinla.setOption.core(..., option_value_alternating = TRUE)
  }
  return(invisible())
}



# Utils ----



#' @describeIn inlabru-deprecated Global setting for tutorial sessions.
#'
#' Use [bru_options_set()] to set specific
#' options instead instead.  In versions <= 2.1.15, this function set the INLA
#' integration strategy to "eb" to speed up calculations. This is normally not
#' needed since version 2.2.0, since the only the final iteration will use
#' other than "eb".
#'
#' @aliases init.tutorial
#' @export
#'
#' @return NULL
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # Note: Only run this if you want to change the inlabru options for this session
#'
#' # Determine current bru defaults:
#' bo <- bru_options_get()
#'
#' init.tutorial()
#'
#' # Check if it worked:
#' bru_options_get("control.inla")
#' }
#'
init.tutorial <- function() {
  .Deprecated("bru_options_set(bru_verbose = TRUE, control.compute = list(dic = TRUE, waic = TRUE))")
  message("Setting defaults for tutorial session.")
  local_bru_options_set(
    bru_verbose = TRUE,
    control.compute = list(dic = TRUE, waic = TRUE)
  )
}
