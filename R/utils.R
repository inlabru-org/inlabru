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
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default: FALSE.
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
      n.t <- tryCatch(
        INLA::inla.getOption("num.threads"),
        error = function(e) { e }
      )
      if (inherits(n.t, "simpleError")) {
        if (!quietly) {
          message("inla.getOption() failed. INLA not installed correctly.")
        }
        return(FALSE)
      }
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


#' Fill in missing values in Spatial grids
#'
#' Computes nearest-available-value imputation for missing values in space
#'
#' @param data A SpatialPointsDataFrame, SpatialPixelsDataFrame, or a
#' SpatialGridDataFrame containg data to use for filling
#' @param where A, matrix, data.frame, or SpatialPoints or
#' SpatialPointsDataFrame, containing the locations of the evaluated values
#' @param values A vector of values to be filled in where `is.na(values)` is
#' `TRUE`
#' @param layer,selector Specifies what data column or columns from which to
#' extract data, see [component()] for details.
#' @param batch_size Size of nearest-neighbour calculation blocks, to limit the
#' memory and computational complexity.
#' @return An infilled vector of values
#' @export
#' @examples
#' \dontrun{
#' if (bru_safe_inla()) {
#'   points <-
#'     sp::SpatialPointsDataFrame(
#'       matrix(1:6, 3, 2),
#'       data = data.frame(val = c(NA, NA, NA))
#'     )
#'   input_coord <- expand.grid(x = 0:7, y = 0:7)
#'   input <-
#'     sp::SpatialPixelsDataFrame(
#'       input_coord,
#'       data = data.frame(val = as.vector(input_coord$y))
#'     )
#'   points$val <- bru_fill_missing(input, points, points$val)
#'   print(points)
#'
#'   # To fill in missing values in a grid:
#'   print(input$val[c(3, 30)])
#'   input$val[c(3, 30)] <- NA # Introduce missing values
#'   input$val <- bru_fill_missing(input, input, input$val)
#'   print(input$val[c(3, 30)])
#' }
#' }
bru_fill_missing <- function(data, where, values,
                             layer = NULL, selector = NULL,
                             batch_size = 500) {
  stopifnot(inherits(
    data,
    c(
      "SpatialPointsDataFrame",
      "SpatialPixelsDataFrame",
      "SpatialGridDataFrame"
    )
  ))
  if (inherits(data, "SpatialGridDataFrame")) {
    data <- as(data, "SpatialPixelsDataFrame")
  }
  if (is.null(layer) && is.null(selector)) {
    layer <- 1
  }
  if (!isTRUE({
    msg <- check_selector(data, where, layer, selector)
  })) {
    stop(msg)
  }
  if (inherits(where, "Spatial")) {
    data_crs <- fm_sp_get_crs(data)
    where_crs <- fm_sp_get_crs(where)
    if (!fm_identical_CRS(data_crs, where_crs)) {
      warning("'data' and 'where' for spatial infilling have different CRS")
    }
    where_coord <- sp::coordinates(where)
  } else {
    data_crs <- sp::CRS(NA_character_)
    where_crs <- sp::CRS(NA_character_)
    where_coord <- where
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
          where = sp::SpatialPoints(
            where_coord[selection == l, , drop = FALSE],
            proj4string = where_crs
          ),
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

  for (batch in seq_len(ceiling(length(notok) / batch_size))) {
    subset <- notok[seq((batch - 1) * batch_size,
      min(length(notok), batch * batch_size),
      by = 1
    )]
    dst <- rgeos::gDistance(
      sp::SpatialPoints(data[data_ok, , drop = FALSE], proj4string = data_crs),
      sp::SpatialPoints(where_coord[subset, , drop = FALSE], proj4string = where_crs),
      byid = TRUE
    )

    nn <- apply(dst, MARGIN = 1, function(row) which.min(row)[[1]])
    values[subset] <- data[[layer]][data_ok[nn]]
  }
  values
}



# Resave data
resave_package_data <- function() {
  name_list <- c(
    "gorillas", "mexdolphin", "mrsea",
    "Poisson1_1D", "Poisson2_1D", "Poisson3_1D",
    "seals", "shrimp", "toygroups"
  )
  for (name in name_list) {
    message(paste0("Data: ", name))
    env <- new.env()
    data(list = name, package = "inlabru", envir = env)

    # Find paths
    new_path <- file.path("data", paste0(name, ".rda"))
    old_path <- file.path("data", paste0(name, ".RData"))
    if (!file.exists(old_path)) {
      old_path <- new_path
    }

    old_info <- file.info(old_path)
    if (length(names(env)) == 1) {
      eval(
        parse(text = paste0(
          "usethis::use_data(",
          paste0(names(env), collapse = ", "),
          ", compress = 'xz', overwrite = TRUE)"
        )),
        envir = env
      )
    } else {
      eval(
        parse(text = paste0(
          "save(",
          paste0(names(env), collapse = ", "),
          ", file = '",
          new_path,
          "', compress = 'xz')"
        )),
        envir = env
      )
    }
    new_info <- file.info(new_path)
    browser()
  }
}


#' Row-wise Kronecker products
#'
#' Takes two Matrices and computes the row-wise Kronecker product.  Optionally
#' applies row-wise weights and/or applies an additional 0/1 row-wise Kronecker
#' matrix product.
#'
#' @param M1 A matrix that can be transformed into a sparse Matrix.
#' @param M2 A matrix that can be transformed into a sparse Matrix.
#' @param repl An optional index vector.  For each entry, specifies which
#' replicate the row belongs to, in the sense used in
#' `INLA::inla.spde.make.A`
#' @param n.repl The maximum replicate index, in the sense used in
#' `INLA::inla.spde.make.A()`.
#' @param weights Optional scaling weights to be applied row-wise to the
#' resulting matrix.
#' @return A `Matrix::sparseMatrix` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export row_kron
row_kron <- function(M1, M2, repl = NULL, n.repl = NULL, weights = NULL) {
  M1 <- as(as(as(M1, "CsparseMatrix"), "dgCMatrix"), "dgTMatrix")
  M2 <- as(as(as(M2, "CsparseMatrix"), "dgCMatrix"), "dgTMatrix")
  n1 <- nrow(M1)
  n2 <- nrow(M2)
  if ((n1 == 1) && (n2 > 1)) {
    M1 <- Matrix::kronecker(rep(1, n2), M1)
    n <- n2
  } else if ((n1 > 1) && (n2 == 1)) {
    M2 <- Matrix::kronecker(rep(1, n1), M2)
    n <- n1
  } else if (n1 != n2) {
    stop(paste0("Size mismatch for row.kron, (n1, n2) = (", n1, ", ", n2, ")"))
  } else {
    n <- n1
  }
  if (is.null(repl)) {
    repl <- rep(1L, n)
  }
  if (is.null(n.repl)) {
    n.repl <- max(repl)
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else if (length(weights) == 1L) {
    weights <- rep(weights[1], n)
  }

  ## TODO: Check robustness for all-zero rows.
  ## TODO: Maybe move big sparseMatrix call outside the loop.
  ## TODO: Automatically choose M1 or M2 for looping.

  n1 <- (as.vector(Matrix::sparseMatrix(
    i = 1L + M1@i, j = rep(1L, length(M1@i)),
    x = 1L, dims = c(n, 1)
  )))
  n2 <- (as.vector(Matrix::sparseMatrix(
    i = 1L + M2@i, j = rep(1L, length(M2@i)),
    x = 1L, dims = c(n, 1)
  )))

  M <- (Matrix::sparseMatrix(
    i = integer(0), j = integer(0), x = numeric(0),
    dims = c(n, ncol(M2) * ncol(M1) * n.repl)
  ))
  n1 <- n1[1L + M1@i]
  for (k in unique(n1)) {
    sub <- which(n1 == k)
    n.sub <- length(sub)

    i.sub <- 1L + M1@i[sub]
    j.sub <- 1L + M1@j[sub]
    o1 <- order(i.sub, j.sub)
    jj <- rep(seq_len(k), times = n.sub / k)
    i.sub <- i.sub[o1]
    j.sub <- (Matrix::sparseMatrix(
      i = i.sub,
      j = jj,
      x = j.sub[o1],
      dims = c(n, k)
    ))
    x.sub <- (Matrix::sparseMatrix(
      i = i.sub,
      j = jj,
      x = weights[i.sub] * M1@x[sub][o1],
      dims = c(n, k)
    ))
    sub2 <- which(is.element(1L + M2@i, i.sub))

    if (length(sub2) > 0) {
      i <- 1L + M2@i[sub2]
      ii <- rep(i, times = k)
      repl.i <- repl[ii]

      M <- (M +
        Matrix::sparseMatrix(
          i = ii,
          j = (1L + rep(M2@j[sub2], times = k) +
            ncol(M2) * (as.vector(j.sub[i, ]) - 1L) +
            ncol(M2) * ncol(M1) * (repl.i - 1L)),
          x = (rep(M2@x[sub2], times = k) *
            as.vector(x.sub[i, ])),
          dims = c(n, ncol(M2) * ncol(M1) * n.repl)
        ))
    }
  }

  return(M)
}
