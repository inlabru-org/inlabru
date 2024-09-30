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



#' @details `local_bru_options_set()` is used to set global package options.
#' @return `local_bru_options_set()` returns a copy of the global override
#' options (not including the defaults), invisibly.
#' @seealso [bru_options_set()], [bru_options_default()], [bru_options_get()]
#' @param .reset For `local_bru_options_set`, logical indicating if the global
#' override options list should be emptied before setting the new option(s).
#'
#' @examples
#' my_fun <- function(val) {
#'   local_bru_options_set(bru_verbose = val)
#'   bru_options_get("bru_verbose")
#' }
#' # Inside the function, the bru_verbose option is changed.
#' # Outside the function, the bru_verbose option is unchanged.
#' print(my_fun(TRUE))
#' print(bru_options_get("bru_verbose"))
#' print(my_fun(FALSE))
#' print(bru_options_get("bru_verbose"))
#' @export
#' @describeIn local_testthat Calls [bru_options_set()] in a reversible way

local_bru_options_set <- function(...,
                                  .reset = FALSE,
                                  envir = parent.frame()) {
  old_opt <- bru_options_get(include_default = FALSE)
  withr::defer(bru_options_set(old_opt, .reset = TRUE), envir = envir)
  bru_options_set(..., .reset = .reset)
  invisible(old_opt)
}



#' @export
#' @rdname local_testthat
local_basic_intercept_testdata <- function() {
  set.seed(123)
  data.frame(
    Intercept = 1,
    y = rnorm(100)
  )
}

#' @export
#' @rdname local_testthat
local_basic_fixed_effect_testdata <- function() {
  set.seed(123)
  cbind(
    local_basic_intercept_testdata(),
    data.frame(x1 = rnorm(100))
  )
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
    old_threads <- tryCatch(
      INLA::inla.getOption("num.threads"),
      error = function(e) {
        e
      }
    )
    if (inherits(old_threads, "simpleError")) {
      return(testthat::skip("inla.getOption() failed, skip INLA tests."))
    }
    withr::defer(
      INLA::inla.setOption(num.threads = old_threads),
      envir
    )

    # Save the fmesher.timeout option so it can be restored
    old_fmesher_timeout <- INLA::inla.getOption("fmesher.timeout")
    withr::defer(
      INLA::inla.setOption(fmesher.timeout = old_fmesher_timeout),
      envir
    )
    INLA::inla.setOption(fmesher.timeout = 30)

    if ("fmesher.evolution" %in% names(INLA::inla.getOption())) {
      # Save the fmesher.evolution option so it can be restored
      old_fmesher_evolution <- INLA::inla.getOption("fmesher.evolution")
      withr::defer(
        INLA::inla.setOption(fmesher.evolution = old_fmesher_evolution),
        envir
      )
      INLA::inla.setOption(fmesher.evolution = max(old_fmesher_evolution, 2L))
    }

    if ("fmesher.evolution.warn" %in% names(INLA::inla.getOption())) {
      # Save the fmesher.evolution.warn option so it can be restored
      old_fmesher_evolution_warn <-
        INLA::inla.getOption("fmesher.evolution.warn")
      withr::defer(
        INLA::inla.setOption(
          fmesher.evolution.warn = old_fmesher_evolution_warn
        ),
        envir
      )
      INLA::inla.setOption(fmesher.evolution.warn = TRUE)
    }

    if ("fmesher.evolution.verbosity" %in% names(INLA::inla.getOption())) {
      # Save the fmesher.evolution.verbosity option so it can be restored
      old_fmesher_evolution_verbosity <-
        INLA::inla.getOption("fmesher.evolution.verbosity")
      withr::defer(
        INLA::inla.setOption(
          fmesher.evolution.verbosity = old_fmesher_evolution_verbosity
        ),
        envir
      )
      INLA::inla.setOption(fmesher.evolution.verbosity = "stop")
    }
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
