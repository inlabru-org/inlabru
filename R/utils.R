#' Load INLA safely for examples and tests
#'
#' Loads the INLA package with `requireNamespace("INLA", quietly = TRUE)`, and
#' optionally checks and sets the multicore `num.threads` INLA option.
#'
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered. Default: FALSE, multicore
#' not allowed (used for examples and unit tests).
#' @param quietly logical; if `TRUE`, prints diagnostic messages. A message is
#' always printed if the INLA `num.threads` option is altered, regardless of the
#' `quietly` argument. Default: TRUE.
#' @export
#' @return logical; `TRUE` if INLA was loaded safely, otherwise FALSE
#'
#' @examples
#' \dontrun{
#' if (bru_safe_inla()) {
#'   # Run inla dependent calculations
#' }
#' }
#'
bru_safe_inla <- function(multicore = FALSE, quietly = TRUE) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    if (!multicore) {
      n.t <- INLA::inla.getOption("num.threads")
      if (!quietly) {
        message(paste0("Current num.threads is '", n.t, "'."))
      }
      if (!identical(n.t, "1:1")) {
        message(paste0(
          "Setting INLA option num.threads to '1:1'.",
          " Previous value '", n.t, "'."
        ))
        INLA::inla.setOption(num.threads = "1:1")
      } else {
        if (!quietly) {
          message("No num.threads change needed.")
        }
      }
    }
    TRUE
  } else {
    if (!quietly) {
      print("INLA not loaded safely.")
    }
    FALSE
  }
}


#' Expand labels
#' 
#' @param labels character vector; original labels
#' @param expand character vector; subset of labels to expand
#' @param suffix character; the suffix to add to the labels selected by `expand`
#' @return a vector of labels with suffix appended to the selected labels
expand_labels <- function(labels, expand, suffix) {
  labels[labels %in% expand] <- paste0(labels[labels %in% expand], suffix)
}

extract_matrixlist_column <- function(thelist, col) {
  vapply(
    names(thelist),
    function(x) {
      list(as.vector(thelist[[x]][, col]))
    },
    list(1)
  )
}

extract_vectorlist_column <- function(thelist) {
  vapply(
    names(thelist),
    function(x) {
      list(as.vector(thelist[[x]]))
    },
    list(1)
  )
}
