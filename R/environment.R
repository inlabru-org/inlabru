iinla.env <- new.env()
iinla.env$log <- paste0(Sys.time(), ": inlabru start")

.onAttach <- function(libname, pkgname) {
  iinla.env$log_active <- FALSE
  if (length(find.package("INLA", quiet = TRUE)) == 0) {
    iinla.env$inla.installed <- FALSE
  } else {
    iinla.env$inla.installed <- TRUE
  }
}

requireINLA <- function() {
  tc <- tryCatch(attachNamespace("INLA"), error = function(x) {})
  # if ( iinla.env$inla.installed ) {
  # if (!isNamespaceLoaded("INLA")) {
  #   attachNamespace("INLA")
  # }

  # } else {
  #  stop("This function requires INLA to be installed. Please consult www.r-inla.org on how to install INLA. After installing INLA please reload inlabru.")
  # }
}



#' @title Global setting for tutorial sessions
#'
#' @description Increases verbosity and sets the inference strategy to empirical Bayes.
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
#' # Determine current bru default:
#' bo <- bru.options()
#'
#' # By default, INLA's integration strategy is set to the INLA default 'auto':
#' bo$inla.options$control.inla
#'
#' # Now, let's run init.tutorial() to make empirical Bayes the default
#' # integration method when `bru` calls `inla`
#'
#' init.tutorial()
#'
#' # Check if it worked:
#' bru.options()$inla.options$control.inla
#' }
#'
init.tutorial <- function() {
  cat("Setting defaults for tutorial session. \n")
  iinla.setOption("iinla.verbose", list(TRUE))
  iinla.setOption("control.inla", list(list(int.strategy = "eb")))
  iinla.setOption("control.compute", list(list(config = TRUE, dic = TRUE, waic = TRUE)))
}

#' @title inlabru log messages
#'
#' @description Retrieve, add, and/or print log messages
#'
#' @param txt character; log message.
#' @param verbose logical; if `TRUE`, print the log message on screen with
#' `message(txt)`. Default: `iinla.getOptions("iinla.verbose")`
#' @details The log message is stored if the log is active, see
#' [bru_log_active()]
#' @return `bru_log` invisibly returns the added log message.
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and Finn Lindgren
#' \email{finn.lindgren@@gmail.com}

bru_log <- function(txt, verbose = NULL) {
  if (is.null(verbose)) {
    verbose <- iinla.getOption("iinla.verbose")
  }
  if (isTRUE(verbose)) {
    message(txt)
  }
  msg <- paste0(Sys.time(), ": ", txt)
  if (iinla.env$log_active) {
    iinla.env$log <- c(iinla.env$log, msg)
  }
  invisible(msg)
}


#' @param pretty logical; If `TRUE`, return a single string with the log
#' messages separated and terminated by line feeds, suitable for `cat()`.
#' If `FALSE`, return the raw log as a vector of strings. Default: `TRUE`
#' @export
#' @rdname bru_log

bru_log_get <- function(pretty = TRUE) {
  if (pretty) {
    paste0(paste0(iinla.env$log, collapse = "\n"), "\n")
  } else {
    iinla.env
  }
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
  current <- iinla.env$log_active
  if (!is.null(activation)) {
    stopifnot(is.logical(activation))
    iinla.env$log_active <- activation
  }
  current
}
#' @details `bru_log_reset` clears the log
#' @export
#' @rdname bru_log

bru_log_reset <- function() {
  iinla.env$log <- paste0(Sys.time(), ": inlabru log reset")
  invisible(iinla.env$log)
}

#' Merge defaults with overriding options
#'
#' Helper function for setting option variables and lists.
#'
#' @details
#' \itemize{
#'   \item Atomic values override defaults.
#'   \item `NULL` values are replaced by defaults.
#'   \item Missing elements of an option list are set to the default values.
#' }
#'
#' @param options An atomic element or a list with named entries.
#' @param defaults An atomic element or a list with named entries.
#' @return The merged option(s).
#'
#' @keywords internal
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' def <- list(iterations = 10, method = "newton")
#' inlabru:::override_config_defaults(list(iterations = 20, verbose = TRUE), def)
override_config_defaults <- function(options, defaults) {
  if (is.null(options)) {
    return(defaults)
  }
  if (!is.list(options)) {
    return(options)
  }
  for (name in names(options)) {
    defaults[[name]] <- options[[name]]
  }
  defaults
}

iinla.getOption <- function(
                            option = c(
                              "control.compute",
                              "control.inla",
                              "iinla.verbose",
                              "control.fixed"
                            )) {
  if (missing(option)) {
    stop("argument is required.")
  }
  envir <- iinla.env
  option <- match.arg(option, several.ok = TRUE)
  if (exists("iinla.options", envir = envir)) {
    opt <- get("iinla.options", envir = envir)
  } else {
    opt <- list()
  }

  default.opt <-
    list(
      iinla.verbose = FALSE,
      control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
      control.inla = list(int.strategy = "auto"),
      control.fixed = list(expand.factor.strategy = "inla")
    )
  res <- c()
  for (i in 1:length(option)) {
    if (option[i] %in% names(opt)) {
      res[[option[i]]] <- override_config_defaults(
        opt[[option[i]]],
        default.opt[[option[i]]]
      )
    } else {
      res[[option[i]]] <- default.opt[[option[i]]]
    }
  }
  if (length(res) == 1) {
    res <- res[[1]]
  }
  return(res)
}


iinla.setOption <- function(...) {
  iinla.setOption.core <- function(option = c(
                                     "control.compute",
                                     "control.inla",
                                     "iinla.verbose",
                                     "control.fixed"
                                   ),
                                   value) {
    envir <- iinla.env
    option <- match.arg(option, several.ok = FALSE)
    if (!exists("iinla.options", envir = envir)) {
      assign("iinla.options", list(), envir = envir)
    }
    if (is.character(value)) {
      eval(parse(text = paste("iinla.options$", option,
        "=", shQuote(value),
        sep = ""
      )), envir = envir)
    } else {
      eval(parse(text = paste("iinla.options$", option,
        "=", ifelse(is.null(value), "NULL", value),
        sep = ""
      )), envir = envir)
    }
    return(invisible())
  }
  called <- list(...)
  len <- length(names(called))
  if (len > 0L) {
    for (i in 1L:len) {
      do.call(iinla.setOption.core, args = list(
        names(called)[i],
        called[[i]]
      ))
    }
  } else {
    iinla.setOption.core(...)
  }
  return(invisible())
}
