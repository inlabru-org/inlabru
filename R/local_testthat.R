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
local_testthat_assign <- function(x, values, envir = parent.frame()) {
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
local_testthat_tolerances <- function(tolerances = c(1e-4, 1e-2, 1e-1),
                                      envir = parent.frame()) {
  local_testthat_assign("lowtol", tolerances[1], envir = envir)
  local_testthat_assign("midtol", tolerances[2], envir = envir)
  local_testthat_assign("hitol", tolerances[3], envir = envir)
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



#' @describeIn local_testthat Disable PROJ4/6 warnings.
#' To be used within package tests. Restores state on exit.
#'
#' @param proj4 logical; whether to show PROJ4 conversion warnings. Default `FALSE`
#' @param thin logical; whether to show only a thinned version of rgdal PROJ6
#' warnings. Default `TRUE`
#' @export
local_set_PROJ6_warnings <- function(proj4 = FALSE,
                                     thin = TRUE,
                                     envir = parent.frame()) {
  withr::local_options(
    list(
      "rgdal_show_exportToProj4_warnings" =
        if (!proj4) {
          "none"
        } else if (thin) {
          "thin"
        } else {
          "all"
        }
    ),
    .local_envir = envir
  )
  requireNamespace("rgdal", quietly = TRUE)
  if (fm_has_PROJ6()) {
    old1 <- rgdal::get_rgdal_show_exportToProj4_warnings()
    withr::defer(
      rgdal::set_rgdal_show_exportToProj4_warnings(old1),
      envir = envir
    )
    rgdal::set_rgdal_show_exportToProj4_warnings(proj4)

    old2 <- rgdal::get_thin_PROJ6_warnings()
    withr::defer(
      rgdal::set_thin_PROJ6_warnings(old2),
      envir = envir
    )
    rgdal::set_thin_PROJ6_warnings(thin)
  }
}


#' @export
#' @describeIn local_testthat Return a list of the current rgdal warning options
local_get_rgdal_options <- function() {
  requireNamespace("rgdal", quietly = TRUE)
  list(
    option_rgdal_show_exportToProj4_warnings =
      getOption("rgdal_show_exportToProj4_warnings"),
    rgdal_show_exportToProj4_warnings = rgdal::get_rgdal_show_exportToProj4_warnings(),
    thin_PROJ6_warnings = rgdal::get_thin_PROJ6_warnings()
  )
}

#' @export
#' @describeIn local_testthat Disable rgdal PROJ4 conversion warnings and thin
#' PROJ6 warnings.
local_disable_PROJ6_warnings <- function(envir = parent.frame()) {
  local_set_PROJ6_warnings(proj4 = FALSE, thin = TRUE, envir = envir)
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




#' @export
#' @rdname local_testthat
local_mrsea_convert <- function(x, use_km = FALSE) {
  # The estimation is numerically unreliable when the spatial
  # domain is represented in metres, and has been seen to produce
  # different results on different systems (e.g. Travis CI).

  # The data is stored in km scale
  if (!use_km) {
    # Transform km to m:
    crs_m <- fm_crs_set_lengthunit(x$mesh$crs, "m")
    x$mesh <- fm_spTransform(x$mesh, crs_m)
    x$samplers <- sp::spTransform(x$samplers, crs_m)
    x$samplers$weight <- x$samplers$weight * 1000
    x$points <- sp::spTransform(x$points, crs_m)
    x$boundary <- sp::spTransform(x$boundary, crs_m)
    x$covar <- sp::spTransform(x$covar, crs_m)
    x$points$Effort <- x$points$Effort * 1000
    x$points$mid.x <- x$points$mid.x * 1000
    x$points$mid.y <- x$points$mid.y * 1000
    x$points$start.x <- x$points$start.x * 1000
    x$points$start.y <- x$points$start.y * 1000
    x$points$end.x <- x$points$end.x * 1000
    x$points$end.y <- x$points$end.y * 1000
    x$points$distance <- x$points$distance * 1000
    x$samplers$Effort <- x$samplers$Effort * 1000
    x$samplers$mid.x <- x$samplers$mid.x * 1000
    x$samplers$mid.y <- x$samplers$mid.y * 1000
  }
  x
}




#' @describeIn local_testthat Tests should set num.threads = "1:1" to ensure
#' within-system repeatability by calling `local_bru_safe_inla()`;
#' see also [bru_safe_inla()]
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered. Default: `FALSE`, multicore
#' not allowed (used for examples and unit tests).
#' @param quietly logical; if `TRUE`, prints diagnostic messages. A message is
#' always printed if the INLA `num.threads` option is altered, regardless of the
#' `quietly` argument. Default: TRUE.
#' @export
local_bru_safe_inla <- function(multicore = FALSE,
                                quietly = TRUE,
                                envir = parent.frame()) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    # Save the num.threads option so it can be restored
    old_threads <- tryCatch(
      INLA::inla.getOption("num.threads"),
      error = function(e) { e }
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
  }
  if (!multicore) {
    local_bru_options_set(num.threads = "1:1", envir = envir)
  }
  testthat::skip_if_not(bru_safe_inla(multicore = multicore, quietly = quietly))
}


#' @describeIn local_testthat Initialise environment for tests.
#' Disables PROJ4/PROJ6 warnings, and assigns tolerance variables.
#' To be called either at the top of a testfile, or inside tests.
#' Does *not* call [local_bru_safe_inla()], since that may invoke a skip and
#' should be called inside each test that relies on INLA.
#' @export
local_bru_testthat_setup <- function(envir = parent.frame()) {
  local_disable_PROJ6_warnings(envir = envir)
  local_testthat_tolerances(envir = envir)
  local_bru_options_set(
    control.compute = list(dic = FALSE, waic = FALSE),
    inla.mode = "experimental",
    envir = envir
  )
}
