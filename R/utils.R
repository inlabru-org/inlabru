#' Load INLA safely for examples and tests
#'
#' Loads the INLA package with `requireNamespace("INLA", quietly = TRUE)`, and
#' optionally checks and sets the multicore `num.threads` INLA option.
#'
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered.
#' If `FALSE`, forces `num.threads="1:1"`. Default: NULL, checks
#' if running in testthat or non-interactively, in which case sets
#' `multicore=FALSE`, otherwise `TRUE`.
#' @param quietly logical; if `TRUE`, prints diagnostic messages. A message is
#' always printed if the INLA `num.threads` option is altered, regardless of the
#' `quietly` argument. Default: FALSE.
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
bru_safe_inla <- function(multicore = NULL,
                          quietly = FALSE) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    if (is.null(multicore)) {
      multicore <-
        !identical(Sys.getenv("TESTTHAT"), "true") ||
          interactive()
    }
    if (!multicore) {
      n.t <- INLA::inla.getOption("num.threads")
      if (!quietly) {
        message(paste0("Current num.threads is '", n.t, "'."))
      }
      if (!identical(n.t, "1:1")) {
        if (!quietly) {
          message(paste0(
            "Setting INLA option num.threads to '1:1'.",
            " Previous value '", n.t, "'."
          ))
        }
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
      message("INLA not loaded safely.")
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
  labels
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







check_selector <- function(data, where, layer, selector) {
  if (is.null(selector)) {
    if (is.character(layer)) {
      if (!(layer %in% names(data))) {
        return(
          paste0(
            "Input layer name '",
            layer,
            "' doesn't match available variable names.\n",
            "Available names are '",
            paste0(names(data), collapse = "', '"),
            "'.\n",
            "Use *_layer for the input component to specify a valid name."
          )
        )
      }
    } else if (is.numeric(layer)) {
      if ((layer < 1) || (layer > NCOL(data))) {
        return(
          paste0(
            "Input layer nr ", layer,
            " is not in the valid range, [",
            1, ", ", NCOL(data)
          )
        )
      }
    } else {
      return(
        paste0(
          "Unable to identify the spatial data frame layer to evaluate.\n",
          "Available names are '",
          paste0(names(data), collapse = "', '"),
          "'.\n",
          "Use *_layer for the input component to specify a valid name."
        )
      )
    }
  } else {
    if (is.null(where[[selector]])) {
      return(
        paste0("'selector' is non-null, but not such label found in the 'where' object")
      )
    }
  }
  TRUE
}

eval_SpatialDF <- function(data, where, layer = NULL, selector = NULL) {
  stopifnot(inherits(
    data,
    c(
      "SpatialPixelsDataFrame",
      "SpatialGridDataFrame"
    )
  ))
  if (is.null(layer) && is.null(selector)) {
    layer <- 1
  }
  if (!isTRUE({
    msg <- check_selector(data, where, layer, selector)
  })) {
    stop(msg)
  }
  if (is.null(selector)) {
    val <- sp::over(
      where,
      data
    )[, layer, drop = TRUE]
  } else {
    layer <- where[[selector]]
    val <- numeric(NROW(where))
    for (l in unique(layer)) {
      val[layer == l] <- over(
        where[layer == l, , drop = FALSE],
        data
      )[, l, drop = TRUE]
    }
  }
  val
}


# Fill in missing values in Spatial grids
# @param data A SpatialPointsDataFrame, SpatialPixelsDataFrame, or a
# SpatialGridDataFrame containg data to use for filling
# @param data A, matrix, data.frame, or SpatialPoints or
# SpatialPointsDataFrame, containing the locations of the evaluated values
# @param values A vector of values to be filled in where `is.na(values)` is
# `FALSE`
# @param layer,selector Specifies what data column or columns from which to
# extract data, see [component()] for details.
# @param batch_size
# @export
bru_fill_missing <- function(data, where, values,
                             layer = NULL, selector = NULL,
                             batch_size = 1000) {
  stopifnot(inherits(
    data,
    c(
      "SpatialPixelsDataFrame",
      "SpatialGridDataFrame"
    )
  ))
  if (is.null(layer) && is.null(selector)) {
    layer <- 1
  }
  if (!isTRUE({
    msg <- check_selector(data, where, layer, selector)
  })) {
    stop(msg)
  }
  if (inherits(where, "Spatial")) {
    data_crs <- INLA::inla.sp_get_crs(data)
    where_crs <- INLA::inla.sp_get_crs(where)
    if (!fm_identical_CRS(data_crs, where_crs)) {
      warning("'data' and 'where' for spatial infilling have different CRS")
    }
  }

  if (!is.null(selector)) {
    selection <- where[[selector]]
    selector_notok <- is.na(selection)
    if (any(selector_notok)) {
      # Only works if the selector is also in the data object.
      selection <-
        bru_fill_missing(
          data = data,
          where = where,
          values = selection,
          layer = selection,
          batch_size = batch_size
        )
    }
    layers <- unique(selection)
    for (l in layers) {
      values[selection == l] <-
        bru_fill_missing(
          data = data,
          where = where[selection == l, , drop = FALSE],
          values = values[selection == l],
          layer = l,
          batch_size = batch_size
        )
    }
    return(values)
  }

  notok <- is.na(values)
  ok <- which(!notok)
  notok <- which(notok)
  data_notok <- is.na(data[[layer]])
  data_ok <- which(!data_notok)
  data_notok <- which(data_notok)

  batch_size <- 200
  for (batch in seq_len(ceiling(length(notok) / batch_size))) {
    subset <- notok[seq((batch - 1) * batch_size,
      min(length(notok), batch * batch_size),
      by = 1
    )]
    dst <- rgeos::gDistance(
      sp::SpatialPoints(data[data_ok, , drop = FALSE]),
      sp::SpatialPoints(where[subset, , drop = FALSE]),
      byid = TRUE
    )

    nn <- apply(dst, MARGIN = 1, function(row) which.min(row)[[1]])
    values[subset] <- data[[layer]][data_ok[nn]]
  }
  values
}
