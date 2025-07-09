#' @title Unit test helpers
#' @name local_testthat
#' @description Local helper functions for package unit tests
#' @param envir environment for exit handlers
#' @rdname local_testthat
#' @keywords internal
NULL

#' @param x character; Name of variable to assign to
#' @param values the object to assign to `x`
#' @export
#' @describeIn local_testthat Assign local variable. Useful for easy cleanup
#' of global workspace with `withr::deferred_run()` when running tests
#' interactively.
local_bru_testthat_assign <- function(x, values, envir = parent.frame()) {
  exist <- exists(x, envir = envir)
  if (exist) {
    old_value <- envir[[x]]
    withr::defer(assign(x, old_value, envir = envir), envir = envir)
  } else {
    withr::defer(rm(list = x, envir = envir), envir = envir)
  }
  assign(x, values, envir = envir)
}

#' @param tolerances numeric vector of length 3; `[lowtol, midtol, hitol]`
#' @export
#' @describeIn local_testthat Assign test tolerances
#' Assign local tolerance variables. Useful for easy cleanup
#' of global workspace with `withr::deferred_run()` when running tests
#' interactively.
local_bru_testthat_tolerances <- function(tolerances = c(1e-4, 1e-2, 1e-1),
                                          envir = parent.frame()) {
  local_bru_testthat_assign("lowtol", tolerances[1], envir = envir)
  local_bru_testthat_assign("midtol", tolerances[2], envir = envir)
  local_bru_testthat_assign("hitol", tolerances[3], envir = envir)
}



#' @describeIn local_testthat Wrapper for [bru_options_set_local()],
#' to locally override the global package options.
#' @return `local_bru_options_set()` returns a copy of the global override
#' options (not including the defaults), invisibly.
#' @seealso [bru_options_set_local()], [bru_options_default()],
#'   [bru_options_get()]
#' @param .reset For `local_bru_options_set`, logical indicating if the global
#' override options list should be emptied before setting the new option(s).
#'
#' @export

local_bru_options_set <- function(...,
                                  .reset = FALSE,
                                  envir = parent.frame()) {
  bru_options_set_local(..., .reset = .reset, .envir = envir)
}



#' @export
#' @rdname local_testthat
local_basic_intercept_testdata <- function() {
  withr::local_seed(123)
  data.frame(
    Intercept = 1,
    y = rnorm(100)
  )
}

#' @export
#' @rdname local_testthat
local_basic_fixed_effect_testdata <- function() {
  withr::local_seed(123)
  cbind(
    local_basic_intercept_testdata(),
    data.frame(x1 = rnorm(100))
  )
}


# @returns logical; If the option setting was successful, `TRUE` is returned,
# otherwise `FALSE`.
local_inla_options_set <- function(...,
                                   envir = parent.frame(),
                                   .save_only = FALSE) {
  # Set INLA options
  inla_options <- list(...)
  old_inla_options <- list()
  if (length(inla_options) > 0) {
    for (name in names(inla_options)) {
      old_inla_options[[name]] <- tryCatch(
        INLA::inla.getOption(name),
        error = function(e) {
          e
        }
      )

      if (inherits(old_inla_options[[name]], "simpleError")) {
        return(FALSE)
      }

      if (!.save_only) {
        e <- tryCatch(
          INLA::inla.setOption(name, inla_options[[name]]),
          error = function(e) {
            e
          }
        )
        if (inherits(e, "simpleError")) {
          return(FALSE)
        }
      }

      withr::defer(
        INLA::inla.setOption(name, old_inla_options[[name]]),
        envir
      )
    }
  }
  TRUE
}


#' @describeIn local_testthat Tests should set num.threads = "1:1" to ensure
#'   within-system repeatability by calling `local_bru_safe_inla()`; see also
#'   [bru_safe_inla()]
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the INLA
#'   `num.threads` option is not checked or altered. Default: `FALSE`, multicore
#'   not allowed (used for examples and unit tests).
#' @param quietly logical; if `TRUE`, prints diagnostic messages. A message is
#'   always printed if the INLA `num.threads` option is altered, regardless of
#'   the `quietly` argument. Default: TRUE.
#' @export
local_bru_safe_inla <- function(multicore = FALSE,
                                quietly = TRUE,
                                envir = parent.frame()) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    inla.call <- tryCatch(
      INLA::inla.getOption("inla.call"),
      error = function(e) {
        e
      }
    )
    if (inherits(inla.call, "simpleError")) {
      return(testthat::skip(
        "inla.getOption('inla.call') failed, skip INLA tests."
      ))
    }

    # Save the num.threads option so it can be restored
    local_inla_options_set(
      num.threads = NULL,
      envir = envir,
      .save_only = TRUE
    )

    local_inla_options_set(
      inla.timeout = 60,
      fmesher.timeout = 30,
      fmesher.evolution = 2L,
      fmesher.evolution.warn = TRUE,
      fmesher.evolution.verbosity = "stop",
      envir = envir,
      .save_only = FALSE
    )

    # withr::local_options(lifecycle_verbosity = "quiet", .local_envir = envir)
  }

  if (!multicore) {
    local_bru_options_set(num.threads = "1:1", envir = envir)
  }
  testthat::skip_if_not(bru_safe_inla(multicore = multicore, quietly = quietly))
  # Ignore spurious warnings for INLA 23.12.17. Is fixed in later INLA versions:
  assign(
    "processed.status.for.model.scopy.in.section.latent",
    TRUE,
    INLA::inla.get.inlaEnv()
  )
}


#' @describeIn local_testthat Initialise environment for tests.
#' Assigns tolerance variables.
#' To be called either at the top of a testfile, or inside tests.
#' Does *not* call [local_bru_safe_inla()], since that may invoke a skip and
#' should be called inside each test that relies on INLA.
#' @export
local_bru_testthat_setup <- function(envir = parent.frame()) {
  local_bru_testthat_tolerances(envir = envir)
  local_bru_options_set(
    # Need to specify specific smtp to ensure consistent tests.
    # To specifically test pardiso, need to override locally
    control.compute = list(smtp = "taucs"),
    inla.mode = "compact",
    envir = envir
  )
  sp_version <- tryCatch(
    getNamespaceVersion("sp"),
    error = function(e) {
      NULL
    }
  )
  if (!is.null(sp_version) &&
    (utils::compareVersion(sp_version, "1.6-0") >= 0)) {
    old_sp_evolution_status <- tryCatch(
      sp::get_evolution_status(),
      error = function(e) {
        2L
      },
      warning = function(e) {
        2L
      }
    )
    withr::defer(
      tryCatch(sp::set_evolution_status(old_sp_evolution_status),
        warning = function(e) invisible(NULL),
        error = function(e) invisible(NULL)
      ),
      envir = envir
    )
    bru_safe_sp(quietly = TRUE, force = TRUE)
  }

  invisible()
}
