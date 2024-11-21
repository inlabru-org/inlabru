#' @include 0_inlabru_envir.R

#' @title Get access to the internal environment
#' @details The environment is defined in `0_inlabru_envir.R` which is loaded
#'   first.
#' @keywords internal
bru_env_get <- function() {
  pkg_envir <- parent.env(environment())
  envir <- get0(".inlabru_envir", envir = pkg_envir)
  if (!is.environment(envir)) {
    stop(
      "Something went wrong: cannot find internal .inlabru_envir environment."
    )
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
  # For Matrix coercion deprecation testing: 1=warn, 2=stop, NA=something else
  #  options(Matrix.warnDeprecatedCoerce = 2)
}





# inlabru log methods ----

#' @title Clear log contents
#' @description
#' Clears the log contents up to
#' a given `offset` or `bookmark`. Default: clear the entire log.
#' When `x` is NULL, the global `inlabru` log is updated, and `invisible(NULL)`
#' is returned. Otherwise the updated object is returned (invisibly).
#' @param x A `bru_log` object, or in some cases, and object that can be
#' converted/extracted to a `bru_log` object. `NULL` denotes the global
#' `inlabru` log object.
#' @inheritParams bru_log_offset
#' @returns Returns (invisibly) the modified `bru_log` object, or `NULL` (when
#'   `x` is `NULL`)
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   bru_log_reset()
#' }
#' }
#' @export
#' @family inlabru log methods

bru_log_reset <- function(x = NULL, bookmark = NULL, offset = NULL) {
  offset <- bru_log_offset(x = x, bookmark = bookmark, offset = offset)
  if (is.null(x)) {
    # Clear log up to the given offset
    envir <- bru_env_get()
    log_length <- length(envir[["log"]])
    if (offset >= log_length) {
      envir[["log"]] <- bru_log_new()
      return(invisible(NULL))
    }
    if (offset == 0L) {
      return(invisible(NULL))
    }
    index <- offset + seq_len(log_length - offset)
    envir[["log"]] <- envir[["log"]][index]
    return(invisible(NULL))
  }
  # Clear log up to the given offset
  log_length <- length(x)
  if (offset >= log_length) {
    x <- bru_log_new()
    return(invisible(x))
  }
  if (offset == 0L) {
    return(invisible(x))
  }
  index <- offset + seq_len(log_length - offset)
  x <- x[index]
  return(invisible(x))
}


#' @export
#' @title Create a `bru_log` object
#' @description
#' Create a `bru_log` object, by default empty.
#' @param x An optional character vector of log messages, or `data.frame`
#' with columns `message`, `timestamp`, and `verbosity`.
#' @param bookmarks An optional `integer` vector of named bookmarks
#' message in `x`.
#' @family inlabru log methods
#' @examples
#' x <- bru_log_new()
#' x <- bru_log_message("Test message", x = x)
#' print(x)
bru_log_new <- function(x = NULL, bookmarks = NULL) {
  if (is.null(x)) {
    x <- character(0)
  }
  if (!is.data.frame(x)) {
    x <- data.frame(
      message = as.character(x)
    )
  }
  if (is.null(x[["timestamp"]])) {
    ts <- Sys.time()
    is.na(ts) <- TRUE
    ts <- rep(ts, nrow(x))
    x[["timestamp"]] <- ts
  }
  if (is.null(x[["verbosity"]])) {
    x[["verbosity"]] <- integer(nrow(x))
  }
  bookmarks <- if (is.null(bookmarks)) {
    integer(0)
  } else {
    storage.mode(bookmarks) <- "integer"
    bookmarks
  }
  x <- structure(
    list(
      log = x,
      bookmarks = bookmarks
    ),
    class = "bru_log"
  )
  x
}


#' @title Methods for `bru_log` bookmarks
#' @description
#' Methods for `bru_log` bookmarks.
#'
#' @param bookmark character; The label for a bookmark with a stored offset.
#' @param offset integer; a position offset in the log, with `0L` pointing at
#' the start of the log. If negative, denotes the point `abs(offset)` elements
#' from tail of the log. When `bookmark` is non-NULL, the `offset` applies a
#' shift (forwards or backwards) to the bookmark list.
#' @param x A `bru_log` object. If `NULL`, the global `inlabru` log is used.
#' @returns `bru_log_bookmark()`: Returns the modified `bru_log` object if `x`
#'   is non-NULL.
#' @export
#' @describeIn bru_log_bookmark
#' Set a log bookmark. If `offset` is `NULL` (the default),
#' the bookmark will point to the current end of the log.
#' @family inlabru log methods
bru_log_bookmark <- function(bookmark = "", offset = NULL, x = NULL) {
  if (nchar(bookmark) == 0) {
    stop("Bookmark labels must have at least one character.")
  }
  offset <- bru_log_offset(x = x, bookmark = NULL, offset = offset)
  if (is.null(x)) {
    envir <- bru_env_get()
    envir[["log"]][["bookmarks"]] <- c(
      envir[["log"]][["bookmarks"]],
      structure(offset, names = bookmark)
    )
    return(invisible(NULL))
  }
  x[["bookmarks"]] <- c(
    x[["bookmarks"]],
    structure(offset, names = bookmark)
  )
  return(invisible(NULL))
  invisible(x)
}

#' @export
#' @describeIn bru_log_bookmark Return a integer vector with named elements
#'   being bookmarks into the global `inlabru` log with associated log position
#'   offsets.
#' @returns `bru_log_bookmarks()`: Returns the bookmark vector associated with
#'   `x`
bru_log_bookmarks <- function(x = NULL) {
  if (is.null(x)) {
    bru_env_get()[["log"]][["bookmarks"]]
  } else {
    x[["bookmarks"]]
  }
}

#' @title Position methods for `bru_log` objects
#' @export
#' @description
#' Position methods for `bru_log` objects.
#'
#' @describeIn bru_log_offset
#' Utility function for computing log position offsets.
#' @family inlabru log methods
#' @inheritParams bru_log_bookmark
bru_log_offset <- function(x = NULL,
                           bookmark = NULL,
                           offset = NULL) {
  if (is.null(x)) {
    log_length <- length(bru_env_get()[["log"]])
  } else {
    log_length <- length(x)
  }
  if (is.null(bookmark)) {
    if (is.null(offset) || (offset > log_length)) {
      return(log_length)
    }
    if (offset < 0L) {
      offset <- max(0L, log_length + offset)
    }
    return(offset)
  }
  marks <- c(
    "STARTOFLOG" = 0L,
    bru_log_bookmarks(x = x),
    "ENDOFLOG" = log_length
  )
  if (bookmark == "") {
    bookmark <- names(marks)[length(marks) - 1L]
  }
  found <- names(marks) %in% bookmark
  if (any(found)) {
    which_found <- max(which(found))
    found <- marks[[which_found]]
  } else {
    warning(paste0(
      "Log bookmark '",
      bookmark,
      "' not found; assuming start of log."
    ))
    marks <- c(0L, log_length)
    which_found <- 1L
  }
  which_found <- max(1L, min(length(marks), which_found + offset))
  return(marks[[which_found]])
}

#' @export
#' @describeIn bru_log_offset Utility function for computing index vectors
#' for `bru_log` objects.
#' @param i indices specifying elements to extract. If `character`, denotes
#' the sequence between bookmark `i` and the next bookmark (or the end of the
#' log if `i` is the last bookmark)
#' @param verbosity integer value for limiting the highest verbosity level being
#'   returned.
bru_log_index <- function(x = NULL, i, verbosity = NULL) {
  if (is.null(x)) {
    x <- bru_env_get()[["log"]]
  }
  log_length <- length(x)
  if (is.null(i)) {
    i <- seq_len(log_length)
  } else if (is.character(i)) {
    offset0 <- bru_log_offset(x, bookmark = i, offset = 0L)
    offset1 <- bru_log_offset(x, bookmark = i, offset = 1L)
    i <- offset0 + seq_len(offset1 - offset0)
  } else if (is.logical(i)) {
    stopifnot(length(i) == log_length)
    i <- which(i)
  }
  if (any(i < 0L)) {
    stopifnot(all(i < 0L))
    i <- setdiff(seq_len(log_length), -i)
  }
  if ((length(i) > 0) && !is.null(verbosity)) {
    i <- i[x[["log"]][["verbosity"]][i] <= verbosity]
  }
  i
}


#' @title Access methods for `bru_log` objects
#' @description Access method for `bru_log` objects.
#' Note: Up to version `2.8.0`, `bru_log()` was a deprecated alias for
#' `bru_log_message()`. When running on `2.8.0` or earlier, use `bru_log_get()`
#' to access the global log, and `cat(fit$bru_iinla$log, sep = "\n")` to print a
#' stored estimation object log.
#' After version `2.8.0`, use `bru_log()` to access the global log, and
#' `bru_log(fit)` to access a stored estimation log.
#' @param x An object that is, contains, or can be converted to,
#' a `bru_log` object. If `NULL`, refers to the global `inlabru` log.
#' @param verbosity integer value for limiting the highest verbosity level being
#'   returned.
#' @return `bru_log` A `bru_log` object, containing a
#' character vector of log messages, and potentially a vector of bookmarks.
#' @export
#' @family inlabru log methods
#' @describeIn bru_log Extract stored log messages. If non-`NULL`, the
#'   `verbosity` argument determines the maximum verbosity level of the messages
#'   to extract.
bru_log <- function(x = NULL, verbosity = NULL) {
  if (is.null(x)) {
    x <- bru_env_get()[["log"]]
  }
  if (inherits(x, "bru_log")) {
    if (!is.null(verbosity)) {
      i <- bru_log_index(x, i = NULL, verbosity = verbosity)
      x <- x[i]
    }
    return(x)
  }
  UseMethod("bru_log")
}

#' @rdname bru_log
#' @export
bru_log.character <- function(x, verbosity = NULL) {
  y <- bru_log_new(x = x)
  if (!is.null(verbosity)) {
    i <- bru_log_index(y, i = NULL, verbosity = verbosity)
    y <- y[i]
  }
  y
}

#' @rdname bru_log
#' @export
bru_log.bru_log <- function(x, verbosity = NULL) {
  if (!is.null(verbosity)) {
    i <- bru_log_index(x, i = NULL, verbosity = verbosity)
    x <- x[i]
  }
  x
}

#' @rdname bru_log
#' @export
bru_log.iinla <- function(x, verbosity = NULL) {
  if (is.null(x[["log"]])) {
    return(bru_log_new(character(0)))
  }
  bru_log(x[["log"]], verbosity = verbosity)
}

#' @rdname bru_log
#' @export
bru_log.bru <- function(x, verbosity = NULL) {
  x <- bru_check_object_bru(x)
  if (is.null(x[["bru_iinla"]][["log"]])) {
    return(bru_log_new(character(0)))
  }
  bru_log(x[["bru_iinla"]][["log"]], verbosity = verbosity)
}

#' @describeIn bru_log Print a `bru_log` object with `cat(x, sep = "\n")`.
#' If `verbosity` is `TRUE`, include the verbosity level of each message.
#' @param ... further arguments passed to or from other methods.
#' @param timestamp If `TRUE`, include the timestamp of each message. Default
#'   `TRUE`.
#' @export
#' @examples
#' bru_log(verbosity = 2L)
#' print(bru_log(), timestamp = TRUE, verbosity = TRUE)
#'
print.bru_log <- function(x, ..., timestamp = TRUE, verbosity = FALSE) {
  msg <- x[["log"]][["message"]]
  if (timestamp) {
    msg <- paste0(
      x[["log"]][["timestamp"]], ": ", msg
    )
  }
  if (verbosity) {
    msg <- paste0(
      msg, " (level ", x[["log"]][["verbosity"]], ")"
    )
  }
  cat(msg, sep = "\n")
  invisible(x)
}


#' @export
#' @describeIn bru_log Convert `bru_log` object to a plain `character` vector
as.character.bru_log <- function(x, ...) {
  x[["log"]][["message"]]
}


#' @export
#' @param i indices specifying elements to extract. If `character`, denotes
#' the sequence between bookmark `i` and the next bookmark (or the end of the
#' log if `i` is the last bookmark)
#' @describeIn bru_log Extract a subset of a `bru_log` object
`[.bru_log` <- function(x, i) {
  i <- bru_log_index(x, i)
  if (length(i) == 0L) {
    marks <- integer(0)
  } else {
    marks <- x[["bookmarks"]]
    marks <- marks[(marks >= min(i) - 1L) & (marks <= max(i))]
    for (k in seq_along(marks)) {
      marks[k] <- sum(i <= marks[k])
    }
  }
  bru_log_new(
    x = x[["log"]][i, , drop = FALSE],
    bookmarks = marks
  )
}

#' @export
#' @describeIn bru_log Concatenate several `bru_log` or `character` objects
#' into a `bru_log` object.
`c.bru_log` <- function(...) {
  obj <- list(...)
  is_ch <- vapply(
    obj,
    function(x) inherits(x, "character"),
    TRUE
  )
  if (any(is_ch)) {
    for (k in which(is_ch)) {
      obj[[k]] <- bru_log_new(x = obj[[k]])
    }
  }
  stopifnot(all(vapply(
    obj,
    function(x) inherits(x, "bru_log"),
    TRUE
  )))
  offset <- 0L
  for (k in seq_along(obj)) {
    obj[[k]][["bookmarks"]] <- obj[[k]][["bookmarks"]] + offset
    offset <- offset + length(obj[[k]])
  }
  bru_log_new(
    x = do.call(rbind, lapply(obj, function(x) x[["log"]])),
    bookmarks = do.call(c, lapply(obj, function(x) x[["bookmarks"]]))
  )
}

#' @export
#' @describeIn bru_log Obtain the number of log entries
#' into a `bru_log` object.
`length.bru_log` <- function(x) {
  NROW(x[["log"]])
}


#' @title Add a log message
#' @description Adds a log message.
#' @param ... For `bru_log_message()`, zero or more objects passed on to
#' [`base::.makeMessage()`]
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
#' are printed. If numeric, messages with `verbosity` \eqn{\le}{<=} `verbose`
#' are printed.
#' @param verbose_store Same as `verbose`, but controlling what messages are
#'   stored in the global log object. Can be controlled via the
#'   `bru_verbose_store` with [bru_options_set()].
#' @param x A `bru_log` object. If `NULL`, refers to the global `inlabru` log.
#' @return
#' `bru_log_message` returns `invisible(x)`, where `x` is the updated `bru_log`
#' object, or `NULL`.
#' @family inlabru log methods
#' @examples
#' if (interactive()) {
#'   code_runner <- function() {
#'     local_bru_options_set(
#'       # Show messages up to and including level 2 (default 0)
#'       bru_verbose = 2,
#'       # Store messages to an including level 3 (default Inf, storing all)
#'       bru_verbose_store = 3
#'     )
#'
#'     bru_log_bookmark("bookmark 1")
#'     bru_log_message("Test message 1", verbosity = 1)
#'     bru_log_message("Test message 2", verbosity = 2)
#'     bru_log_bookmark("bookmark 2")
#'     bru_log_message("Test message 3", verbosity = 3)
#'     bru_log_message("Test message 4", verbosity = 4)
#'
#'     invisible()
#'   }
#'   message("Run code")
#'   code_runner()
#'   message("Check log from bookmark 1")
#'   print(bru_log()["bookmark 1"])
#'   message("Check log from bookmark 2")
#'   print(bru_log()["bookmark 2"])
#' }
#' @export

bru_log_message <- function(..., domain = NULL, appendLF = TRUE,
                            verbosity = 1L,
                            allow_verbose = TRUE, verbose = NULL,
                            verbose_store = NULL,
                            x = NULL) {
  new_x <- data.frame(
    message = .makeMessage(...,
      domain = domain,
      appendLF = FALSE
    ),
    timestamp = Sys.time(),
    verbosity = as.integer(verbosity)
  )
  if (allow_verbose) {
    if ((!is.null(verbose) && (verbose >= verbosity)) ||
      (is.null(verbose) &&
        bru_options_get("bru_verbose", include_default = TRUE) >= verbosity)) {
      message(new_x[["message"]], domain = NA, appendLF = appendLF)
    }
  }
  if ((!is.null(verbose_store) && (verbose_store >= verbosity)) ||
    !allow_verbose ||
    (is.null(verbose_store) &&
      (bru_options_get("bru_verbose_store",
        include_default = TRUE
      ) >= verbosity))) {
    if (is.null(x)) {
      envir <- bru_env_get()
      envir[["log"]][["log"]] <- rbind(envir[["log"]][["log"]], new_x)
    } else {
      x[["log"]] <- rbind(x[["log"]], new_x)
    }
  }
  if (is.null(x)) {
    invisible(NULL)
  } else {
    invisible(x)
  }
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
#'   [`bru_options`] objects. Options specified later override the previous
#'   options.
#' @return `bru_options()` returns a `bru_options` object.
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
#' log messages of verbosity \eqn{\le} are stored. Default: Inf, to store all
#' messages.}
#' \item{bru_run}{If TRUE, run inference. Otherwise only return configuration
#' needed to run inference.}
#' \item{bru_max_iter}{maximum number of inla iterations, default 10.
#'  Also see the `bru_method$rel_tol` and related options below.}
#' \item{bru_initial}{An `inla` object returned from previous calls of
#'   `INLA::inla`, [bru()] or [lgcp()], or a list of named vectors of starting
#'   values for the latent variables. This will be used as a
#'   starting point for further improvement of the approximate posterior.}
#' \item{bru_int_args}{List of arguments passed all the way to the
#' integration method `ipoints` and `int.polygon` for 'cp' family models;
#' \describe{
#' \item{method}{"stable" or "direct". For "stable" (default) integration points
#' are aggregated to mesh vertices.}
#' \item{nsub1}{Number of integration points per knot interval in 1D.
#'   Default 30.}
#' \item{nsub2}{Number of integration points along a triangle edge for 2D.
#'   Default 9.}
#' \item{nsub}{Deprecated parameter that overrides `nsub1` and `nsub2` if set.
#'   Default `NULL`.}
#' }
#' }
#' \item{bru_method}{List of arguments controlling the iterative inlabru method:
#' \describe{
#' \item{taylor}{'pandemic'
#' (default, from version 2.1.15).}
#' \item{search}{Either 'all' (default), to use all available line search
#' methods, or one or more of
#' \describe{
#' \item{'finite'}{(reduce step size until predictor is finite)}
#' \item{'contract'}{(decrease step size until trust hypersphere reached)}
#' \item{'expand'}{(increase step size until no improvement)}
#' \item{'optimise'}{(fast approximate error norm minimisation)}
#' }
#' To disable line search, set to an empty vector. Line search is not
#' available for `taylor="legacy"`.}
#' \item{factor}{Numeric, \eqn{> 1} determining the line search step scaling
#' multiplier. Default \eqn{(1 + \sqrt{5})/2}{(1+sqrt(5))/2}.}
#' \item{rel_tol}{Stop the iterations when the largest change in linearisation
#' point
#' (the conditional latent state mode) in relation to the estimated posterior
#' standard deviation is less than `rel_tol`. Default 0.1 (ten percent).}
#' \item{max_step}{The largest allowed line search step factor. Factor 1 is the
#' full INLA step. Default is 2.}
#' \item{line_opt_method}{Which method to use for the line search optimisation
#' step.
#' Default "onestep", using a quadratic approximation based on the value and
#' gradient at zero, and the value at the current best step length guess.
#' The method "full" does line optimisation on the full nonlinear predictor;
#' this is slow and intended for debugging purposes only.}
#' }
#' }
#' \item{bru_compress_cp}{logical; when `TRUE`, compress the
#' \eqn{\sum_{i=1}^n \eta_i}{sum_i=1^n eta_i}
#' part of the Poisson process likelihood (`family="cp"`) into a single term,
#' with \eqn{y=n}{y=n}, and predictor `mean(eta)`. Default: `TRUE`}
#' \item{bru_debug}{logical; when `TRUE`, activate temporary debug features for
#' package development. Default: `FALSE`}
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
#' @describeIn bru_options Coerces inputs to a `bru_options` object.

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

#' @describeIn bru_options Returns the default options.
#' @return `bru_options_default()` returns an `bru_options` object containing
#'   default options.
#' @export

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
      rel_tol = 0.1,
      max_step = 2,
      line_opt_method = "onestep"
    ),
    bru_compress_cp = TRUE,
    bru_debug = FALSE,
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
  handle_args <- function(args, replacements) {
    names_args <- names(args)
    deprecated_args <- replacements[names(replacements) %in% names_args]
    depr_list <- vapply(deprecated_args, is.list, TRUE)
    if (any(depr_list)) {
      for (k in which(depr_list)) {
        args[[names(deprecated_args[k])]] <-
          handle_args(args[[names(deprecated_args[k])]], deprecated_args[[k]])
      }
    }
    deprecated_args <- deprecated_args[!depr_list]
    if (any(nzchar(names(deprecated_args)) == 0)) {
      warning(paste0(
        "Ignoring deprecated global options '",
        paste0(names(deprecated_args)[nzchar(names(deprecated_args)) == 0],
          collapse = "', '"
        ),
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
      names_args <- names(args)
      names_args[names_args %in% names(deprecated_args)] <-
        deprecated_args[names_args[names_args %in% names(deprecated_args)]]
      names(args) <- names_args
    }
    args
  }

  stopifnot(inherits(args, "bru_options"))
  cl <- class(args)
  deprecated_args <- list(
    mesh = "",
    run = "bru_run",
    max.iter = "bru_max_iter",
    result = "bru_initial",
    bru_result = "bru_initial",
    int.args = "bru_int_args",
    bru_method = list(stop_at_max_rel_deviation = "rel_tol")
  )
  args <- handle_args(args, deprecated_args)
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
#' @export
#' @returns A `bru_options` object
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' \donttest{
#'
#' opts <- bru_call_options()
#'
#' # Print them:
#' opts
#' }
#'
#' @keywords internal
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



#' @describeIn bru_options Checks for valid contents of a `bru_options`
#' object, and produces warnings for invalid options.
#' @param options An `bru_options` object to be checked
#' @param ignore_null Ignore missing or NULL options.
#' @return `bru_options_check()` returns a `logical`; `TRUE` if the object
#'   contains valid options for use by other functions
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   bru_options_check(bru_options(bru_max_iter = "text"))
#' }
#' }
#' @export

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
#' bru_options_get("bru_verbose")
#' @export
#' @describeIn bru_options Used to access global package options.

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


#' @describeIn bru_options Used to set global package options.
#' @return `bru_options_set()` returns a copy of the global override options,
#' invisibly (as `bru_options_get(include_default = FALSE)`).
#' @seealso [bru_options()], [bru_options_default()], [bru_options_get()]
#' @param .reset For `bru_options_set`, logical indicating if the global
#'   override options list should be emptied before setting the new option(s).
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

bru_options_set <- function(..., .reset = FALSE) {
  envir <- bru_env_get()
  if (.reset) {
    envir$options <- bru_options(...)
  } else {
    envir$options <- bru_options(envir$options, ...)
  }
  invisible(bru_options_get(include_default = FALSE))
}

#' @describeIn bru_options Clears the global option overrides.
#' @export

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
  result <- structure(
    list(
      legend = legend,
      value = traverse(combined, default, global, object)
    ),
    class = "summary_bru_options"
  )
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
