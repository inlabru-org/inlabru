#' @include deprecated.R

#' @title Load INLA safely for examples and tests
#'
#' @description
#' Loads the INLA package with `requireNamespace("INLA", quietly = TRUE)`, and
#' optionally checks and sets the multicore `num.threads` INLA option.
#'
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered.
#' If `FALSE`, forces `num.threads="1:1:1"`. Default: NULL, checks
#' if running in testthat or non-interactively, in which case sets
#' `multicore=FALSE`, otherwise `TRUE`.
#' @param quietly logical; if `FALSE` and `multicore` is `FALSE`,
#' prints a message if the `num.threads` option
#' isn't already "1.1:1" to alert the user to the change.
#' Default: FALSE.
#' @param minimum_version character; the minimum required INLA version.
#' Default 23.1.31 (should always match the requirement in the package
#' DESCRIPTION)
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
#' @keywords internal
bru_safe_inla <- function(multicore = NULL,
                          quietly = FALSE,
                          minimum_version = "23.1.31") {
  inla_version <-
    check_package_version_and_load(
      pkg = "INLA",
      minimum_version = minimum_version,
      quietly = quietly
    )
  if (is.na(inla_version)) {
    return(FALSE)
  }

  inla.call <- tryCatch(
    INLA::inla.getOption("inla.call"),
    error = function(e) {
      e
    }
  )
  if (inherits(inla.call, "error")) {
    if (!quietly) {
      message(
        "inla.getOption('inla.call') failed. INLA not installed correctly."
      )
    }
    return(FALSE)
  }
  if (all(grepl("\\.run$", inla.call))) {
    inla.binary <- gsub("\\.run$", "", inla.call)
    if (!file.exists(inla.binary)) {
      if (!quietly) {
        message(
          paste0(
            "INLA binary '", inla.binary, "' not found. ",
            "INLA not installed correctly, or with platform mismatch."
          )
        )
      }
      return(FALSE)
    }
  }

  if (is.null(multicore)) {
    multicore <-
      interactive() && !identical(Sys.getenv("TESTTHAT"), "true")
  }
  if (!multicore) {
    n.t <- tryCatch(
      INLA::inla.getOption("num.threads"),
      error = function(e) {
        e
      }
    )
    if (inherits(n.t, "error")) {
      if (!quietly) {
        message("inla.getOption() failed. INLA not installed correctly.")
      }
      return(FALSE)
    }
    if (!identical(n.t, "1:1:1")) {
      if (!quietly) {
        message(paste0(
          "Changing INLA option num.threads from '", n.t, "' to '1:1:1'."
        ))
      }
      INLA::inla.setOption(num.threads = "1:1:1")
    }
  }
  TRUE
}





check_package_version_and_load <-
  function(pkg, minimum_version, quietly = FALSE) {
    version <- tryCatch(utils::packageVersion(pkg),
      error = function(e) NA_character_
    )
    if (is.na(version)) {
      if (!quietly) {
        message(paste0("Package '", pkg, "' is not installed."))
      }
      return(NA_character_)
    }
    if (version < minimum_version) {
      if (!quietly) {
        message(paste0(
          "Installed '", pkg, "' version is ", version, " but ",
          "version >= ", minimum_version, " is required."
        ))
      }
      return(NA_character_)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!quietly) {
        message("Package '", pkg, "' not loaded safely.")
      }
      return(NA_character_)
    }
    version
  }


#' Check for potential `sp` version compatibility issues
#'
#' Loads the sp package with `requireNamespace("sp", quietly = TRUE)`, and
#' checks and optionally sets the `sp` evolution status flag if `rgdal` is
#' unavailable.
#'
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default
#'   `FALSE`
#' @param force logical; If `rgdal` is unavailable
#' and evolution status is less that `2L`, return `FALSE` if `force` is `FALSE`.
#' If `force` is `TRUE`, return `TRUE` if the package configuration is safe,
#' potentially after forcing the evolution status to `2L`.
#' Default `FALSE`
#' @param minimum_version character; the minimum required version.
#' Default 2.1 (should always match the requirement in the package
#' DESCRIPTION)
#' @return Returns (invisibly) `FALSE` if a potential issue is detected, and
#'   give a message if `quietly` is `FALSE`. Otherwise returns `TRUE`
#' @export
#' @examples
#' \dontrun{
#' if (bru_safe_sp()) {
#'   # Run sp dependent calculations
#' }
#' }
#'
#' @keywords internal
bru_safe_sp <- function(quietly = FALSE,
                        force = FALSE,
                        minimum_version = "2.1") {
  sp_version <-
    check_package_version_and_load(
      pkg = "sp",
      minimum_version = minimum_version,
      quietly = quietly
    )
  if (is.na(sp_version)) {
    return(invisible(FALSE))
  }

  if (sp_version >= "2.1.4") {
    return(invisible(TRUE))
  }

  if (sp_version >= "1.6-0") {
    # Default to 2L to allow future sp to stop supporting
    # get_evolution_status; assume everything is fine if it fails.
    evolution_status <- tryCatch(
      sp::get_evolution_status(),
      error = function(e) 2L,
      warning = function(e) 2L
    )
    rgdal_version <- tryCatch(utils::packageVersion("rgdal"),
      error = function(e) NA_character_
    )
    if ((evolution_status < 2L) && is.na(rgdal_version)) {
      if (!quietly) {
        message(
          "'sp' version >= 1.6-0 detected, rgdal isn't installed, and ",
          "evolution status is < 2L."
        )
      }
      if (!force) {
        if (!quietly) {
          message(
            "This may cause issues with some CRS handling code. ",
            "To avoid this, use 'sp::set_evolution_status(2L)'"
          )
        }
        return(invisible(FALSE))
      }

      sp::set_evolution_status(2L)
      if (!quietly) {
        message(
          "Ran 'sp::set_evolution_status(2L)' to avoid ",
          "issues with some CRS handling code."
        )
      }
    }
  }
  invisible(TRUE)
}


#' @title Check for potential `terra` installation issues
#'
#' @description
#' Loads the terra package with `requireNamespace("terra", quietly = TRUE)`,
#' and checks whether `eval_spatial()` successfully extracts non-NA values.
#' NA values may appear when there is a GDAL/PROJ installation issue that causes
#' it to not find the `proj.db` database file.
#'
#' For use in package tests and examples that depend on `terra`.
#'
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default
#'   `FALSE`
#' @param minimum_version character; the minimum required version.
#' Default 1.7-66 (should always match the requirement in the package
#' DESCRIPTION)
#' @return Returns (invisibly) `FALSE` if a potential issue is detected, and
#'   give a message if `quietly` is `FALSE`. Otherwise returns `TRUE`
#' @export
#' @examples
#' \dontrun{
#' if (bru_safe_terra()) {
#'   # Run terra dependent calculations
#' }
#' }
#'
#' @keywords internal
bru_safe_terra <- function(quietly = FALSE,
                           minimum_version = "1.7-66") {
  terra_version <-
    check_package_version_and_load(
      pkg = "terra",
      minimum_version = minimum_version,
      quietly = quietly
    )
  if (is.na(terra_version)) {
    return(invisible(FALSE))
  }

  gcov <- gorillas_sf_gcov()
  where <- inlabru::gorillas_sf$nests[seq_len(5), , drop = FALSE]
  values <- eval_spatial(gcov, where = where)
  if (anyNA(values)) {
    if (!quietly) {
      message(
        paste0(
          "Potential issue detected with 'terra' package spatial evaluation,\n",
          "likely due to GDAL/PROJ installation issues."
        )
      )
    }
    return(invisible(FALSE))
  }

  invisible(TRUE)
}


#' Expand labels
#'
#' @param labels character vector; original labels
#' @param expand character vector; subset of labels to expand
#' @param suffix character; the suffix to add to the labels selected by `expand`
#' @return a vector of labels with suffix appended to the selected labels
#' @keywords internal
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


check_layer <- function(data, where, layer) {
  # This works for both SpatialGrid/PixelDataFrame and SpatRaster
  names_data <- names(data)
  if (is.character(layer)) {
    unique_layer <- unique(layer)
    if (!all((unique_layer %in% names_data))) {
      stop(
        paste0(
          "Input layer name(s) '",
          paste0(unique_layer, collapse = "', '"),
          "' mismatch with available variable names.\n",
          "Available names are '",
          paste0(names_data, collapse = "', '"),
          "'.\n",
          "Use *_layer for the input component to specify a valid name."
        )
      )
    }
  } else if (is.numeric(layer)) {
    ok_layer <- (layer >= 1) & (layer <= length(names_data))
    if (!all(ok_layer)) {
      stop(
        paste0(
          "Input layer(s) nr ",
          paste0(unique(layer[!ok_layer]), collapse = ", "),
          " not in the valid range, [",
          1, ", ", length(names_data), "]"
        )
      )
    }
  } else {
    stop(
      paste0(
        "Unable to identify the spatial data frame layer to evaluate.\n",
        "Available names are '",
        paste0(names(data), collapse = "', '"),
        "'.\n",
        "Use *_layer for the input component to specify a valid name."
      )
    )
  }
  TRUE
}


extract_selector <- function(where, selector) {
  if (is.null(selector)) {
    return(NULL)
  }
  if (inherits(where, "SpatVector")) {
    requireNamespace("terra")
    layer <- terra::values(where)[[selector]]
  } else {
    layer <- where[[selector]]
  }
  if (is.null(layer)) {
    stop(
      "'selector' is non-null, but no such label found in the 'where' object"
    )
  }
  layer
}

extract_layer <- function(where, layer, selector) {
  if (!is.null(layer)) {
    if (!is.null(selector)) {
      warning("Both layer and selector specified. Ignoring selector",
        immediate. = TRUE
      )
    }
    selector <- NULL
  } else if (!is.null(selector)) {
    layer <- extract_selector(where, selector)
  } else if (is.null(layer)) {
    layer <- 1
  }
  if (length(layer) == 1) {
    layer <- rep(layer, NROW(where))
  }
  layer
}


#' Evaluate spatial covariates
#'
#' @param data Spatial data
#' @param where Where to evaluate the data
#' @param layer Which `data` layer to extract (as integer or character).
#' May be a vector, specifying a separate layer for each `where` item.
#' @param selector The name of a variable in `where` specifying the `layer`
#' information.
#'
#' @export
eval_spatial <- function(data, where, layer = NULL, selector = NULL) {
  UseMethod("eval_spatial")
}


#' @export
#' @describeIn eval_spatial Compatibility wrapper for `eval_spatial.sf`
eval_spatial.SpatialPolygonsDataFrame <- function(data,
                                                  where,
                                                  layer = NULL,
                                                  selector = NULL) {
  eval_spatial(
    sf::st_as_sf(data),
    where = where,
    layer = layer,
    selector = selector
  )
}

#' @export
#' @rdname eval_spatial
eval_spatial.SpatialPixelsDataFrame <- function(data,
                                                where,
                                                layer = NULL,
                                                selector = NULL) {
  eval_spatial_Spatial(
    data = data,
    where = where,
    layer = layer,
    selector = selector
  )
}

#' @export
#' @rdname eval_spatial
eval_spatial.SpatialGridDataFrame <- function(data,
                                              where,
                                              layer = NULL,
                                              selector = NULL) {
  eval_spatial_Spatial(
    data = data,
    where = where,
    layer = layer,
    selector = selector
  )
}

eval_spatial_Spatial <- function(data, where, layer = NULL, selector = NULL) {
  stopifnot(inherits(
    data,
    c(
      "SpatialPixelsDataFrame",
      "SpatialGridDataFrame"
    )
  ))
  if (inherits(where, "sf")) {
    where <- sf::as_Spatial(where)
  }

  if (inherits(where, "SpatialPoints")) {
    where <- fm_transform(where, crs = fm_CRS(data), passthrough = TRUE)
    if (ncol(sp::coordinates(where)) >= 3) {
      where <- sp::SpatialPoints(
        coords = sp::coordinates(where)[, 1:2, drop = FALSE],
        proj4string = fm_CRS(where)
      )
    }
  }
  layer <- extract_layer(where, layer, selector)
  check_layer(data, where, layer)
  unique_layer <- unique(layer)
  if (length(unique_layer) == 1) {
    val <- sp::over(
      where,
      data
    )[, unique_layer, drop = TRUE]
  } else {
    val <- numeric(NROW(where))
    for (l in unique(layer)) {
      val[layer == l] <- sp::over(
        where[layer == l, , drop = FALSE],
        data
      )[, l, drop = TRUE]
    }
  }
  val
}


#' @export
#' @describeIn eval_spatial Supports point-in-polygon information lookup.
#' Other combinations are untested.
eval_spatial.sf <- function(data, where, layer = NULL, selector = NULL) {
  if (!inherits(where, c("sf", "sfc", "sfg"))) {
    where <- sf::st_as_sf(where)
  }
  where <- fm_transform(where, crs = sf::st_crs(data), passthrough = TRUE)

  layer <- extract_layer(where, layer, selector)
  check_layer(data, where, layer)
  unique_layer <- unique(layer)
  data_example <- data[1, unique_layer[1], drop = TRUE]
  data_NA <- data_example
  is.na(data_NA) <- TRUE
  if (length(unique_layer) == 1) {
    idx <- sf::st_intersects(where, data, sparse = TRUE)
    val <- vapply(
      idx,
      function(i) {
        if (is.null(i) || (length(i) == 0)) {
          data_NA
        } else {
          data[min(i), unique_layer, drop = TRUE]
        }
      },
      data_example
    )
  } else {
    val <- numeric(NROW(where))
    idx <- sf::st_intersects(where, data, sparse = TRUE)
    for (l in unique(layer)) {
      val[layer == l] <- vapply(
        idx[layer == l],
        function(i) {
          if (is.null(i) || (length(i) == 0)) {
            data_NA
          } else {
            data[min(i), unique_layer, drop = TRUE]
          }
        },
        data_example
      )
    }
  }
  val
}



terra_factor_levels <- function(data) {
  layers <- names(data)
  is_factor <- terra::is.factor(data)
  if (!any(is_factor)) {
    return(NULL)
  }
  if (!all(is_factor)) {
    stop(
      glue::glue(
        "Some layers are factors ",
        "({paste(sort(layers[is_factor]), collapse = ', ')})",
        " and some are not ",
        "({paste(sort(layers[!is_factor]), collapse = ', ')})."
      )
    )
  }
  factor_levels <- lapply(
    terra::levels(data),
    function(l) {
      if (is.data.frame(l)) {
        l[, 2]
      } else {
        NULL
      }
    }
  )
  unique_levels_idx <- which(!duplicated(factor_levels))
  factor_levels <- factor_levels[unique_levels_idx]
  layers <- layers[unique_levels_idx]
  if (length(unique_levels_idx) > 1) {
    stop(
      glue::glue(
        "Factor levels differ between layers. ",
        "Found factor levels:\n",
        glue::glue_collapse(
          glue::glue("{layers}: {factor_levels}", sep = ", "),
          sep = "\n"
        )
      )
    )
  }

  factor_levels[[1]]
}


#' @export
#' @rdname eval_spatial
eval_spatial.SpatRaster <- function(data,
                                    where,
                                    layer = NULL,
                                    selector = NULL) {
  requireNamespace("terra")
  layer <- extract_layer(where, layer, selector)
  check_layer(data, where, layer)
  if (!inherits(where, "SpatVector")) {
    handle_crs <- (!fm_crs_is_null(fm_crs(where)) &&
      !fm_crs_is_null(fm_crs(data)))
    if (handle_crs) {
      where <- fmesher::fm_transform(where, fm_crs(data))
    }
    where <- terra::vect(where)
    if (handle_crs) {
      # Terra's data crs might not be detected as the same as fm_crs(data)
      where <- terra::project(where, data)
    }
  } else {
    if (!fm_crs_is_null(fm_crs(where)) &&
      !fm_crs_is_null(fm_crs(data))) {
      where <- terra::project(where, data)
    }
  }
  val <- terra::extract(
    data,
    where,
    ID = FALSE,
    layer = layer
  )
  val <- val[["value"]]

  layer_set <- sort(unique(layer))
  factor_levels <- terra_factor_levels(data[[layer_set]])
  if (!is.null(factor_levels)) {
    if (is.factor(val)) {
      val <- as.character(val)
      levs <- factor_levels
    } else {
      levs <- seq_along(factor_levels)
    }
    val <- factor(val, levels = levs, labels = factor_levels)
  }
  val
}


#' @export
#' @rdname eval_spatial
eval_spatial.stars <- function(data, where, layer = NULL, selector = NULL) {
  layer <- extract_layer(where, layer, selector)
  check_layer(data, where, layer)
  if (!inherits(where, c("sf", "sfc"))) {
    where <- sf::st_as_sf(where)
  }
  val <- stars::st_extract(
    data,
    sf::st_geometry(where)
  )
  val <- sf::st_as_sf(val)

  unique_layer <- unique(layer)
  if (length(unique_layer) == 1) {
    val <- val[, unique_layer, drop = TRUE]
  } else {
    val_ <- numeric(NROW(where))
    for (l in unique(layer)) {
      val_[layer == l] <- val[layer == l, l, drop = TRUE]
    }
    val <- val_
  }

  val
}


#' Fill in missing values in Spatial grids
#'
#' Computes nearest-available-value imputation for missing values in space
#'
#' @param data A SpatialPointsDataFrame, SpatialPixelsDataFrame,
#' SpatialGridDataFrame, SpatRaster, Raster, or sf object
#' containing data to use for filling
#' @param where A, matrix, data.frame, or SpatialPoints or
#'   SpatialPointsDataFrame, or sf object, containing the locations of the
#'   evaluated values
#' @param values A vector of values to be filled in where `is.na(values)` is
#' `TRUE`
#' @param layer,selector Specifies what data column or columns from which to
#' extract data, see [bru_comp()] for details.
#' @param batch_size `r lifecycle::badge("deprecated")` due to improved
#'   algorithm.
#' Size of nearest-neighbour calculation blocks, to limit the
#' memory and computational complexity.
#' @return An infilled vector of values
#' @export
#' @examples
#' \dontrun{
#' if (require("sf", quietly = TRUE)) {
#'   points <-
#'     sf::st_as_sf(
#'       data.frame(
#'         x = 1:3,
#'         y = 4:6,
#'         val = c(NA, NA, NA)
#'       ),
#'       coords = c("x", "y")
#'     )
#'   input_coord <- expand.grid(x = 0:7, y = 0:7)
#'   input <-
#'     sf::st_as_sf(
#'       cbind(input_coord, val = as.vector(input_coord$y)),
#'       coords = c("x", "y")
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
                             batch_size = deprecated()) {
  stopifnot(inherits(
    data,
    c(
      "SpatialPointsDataFrame",
      "SpatialPixelsDataFrame",
      "SpatialGridDataFrame",
      "Raster",
      "SpatRaster",
      "sf"
    )
  ))
  # Convert to sf and terra
  if (inherits(data, c("SpatialGrid", "SpatialPixelsDataFrame", "Raster"))) {
    requireNamespace("terra")
    data <- terra::rast(data)
  } else if (inherits(data, "SpatialPointsDataFrame")) {
    data <- sf::st_as_sf(data)
  }
  # Convert to sf
  if (inherits(where, c("Spatial", "SpatVector"))) {
    where <- sf::st_as_sf(where)
  } else if (!inherits(where, "sf")) {
    where <- sf::st_as_sf(where, coords = seq_len(min(3, NCOL(where))))
  }

  layer <- extract_layer(where, layer, selector)
  check_layer(data, where, layer)

  if (anyNA(layer)) {
    stop("NAs detected in the 'layer' information.")
  }

  layers <- unique(layer)
  if (length(layers) > 1) {
    for (l in layers) {
      values[layer == l] <-
        bru_fill_missing(
          data = data,
          where = where[layer == l, , drop = FALSE],
          values = values[layer == l],
          layer = l,
          batch_size = batch_size
        )
    }
    return(values)
  }

  # Only one layer from here on.
  layer <- layers

  data_crs <- fm_crs(data)
  if (inherits(data, "SpatRaster")) {
    requireNamespace("terra")
    data_values <- terra::values(
      data[[layer]],
      dataframe = TRUE,
      na.rm = TRUE
    )[[layer]]
    data_coord <- as.data.frame(terra::crds(data[[layer]], na.rm = TRUE))
    data_coord <- sf::st_as_sf(data_coord,
      coords = seq_len(ncol(data_coord)),
      crs = data_crs
    )
  } else {
    data_values <- data[[layer]]
    data_coord <- data
  }
  data_notok <- is.na(data_values)
  if (all(data_notok)) {
    stop("All data values are NA; unable to fill missing values.")
  }
  data_values <- data_values[!data_notok]
  if (inherits(data_coord, "sf")) {
    data_coord <- data_coord[!data_notok, , drop = FALSE]
  } else {
    data_coord <- data_coord[!data_notok, , drop = FALSE]
  }

  where <- fm_transform(where, crs = data_crs, passthrough = TRUE)

  values_notok <- is.na(values)

  nn <- sf::st_nearest_feature(
    where[values_notok, , drop = FALSE],
    data_coord
  )

  values[values_notok] <- data_values[nn]

  values
}


# Resave data
resave_package_data <- function() {
  name_list <- c(
    "gorillas_sf", "mexdolphin_sf",
    "mrsea",
    "Poisson1_1D", "Poisson2_1D", "Poisson3_1D",
    "shrimp", "toygroups",
    "robins_subset",
    "toypoints"
  )
  compress_values <- c("gzip", "bzip2", "xz")

  do_compress <- function(env, compress, the_path) {
    if (length(names(env)) == 1) {
      thetext <- paste0(
        "usethis::use_data(",
        paste0(names(env), collapse = ", "),
        ", compress = '", compress, "', overwrite = TRUE)"
      )
    } else {
      thetext <- paste0(
        "save(",
        paste0(names(env), collapse = ", "),
        ", file = '",
        the_path,
        "', compress = '", compress, "')"
      )
    }
    eval(
      parse(text = thetext),
      envir = env
    )
  }

  new_info <- NULL
  old_info <- NULL
  for (name in name_list) {
    message(paste0("Data: ", name))
    env <- new.env()
    utils::data(list = name, package = "inlabru", envir = env)

    # Find paths
    the_path <- file.path("data", paste0(name, ".rda"))

    old_info <- rbind(old_info, file.info(the_path))

    smallest_size <- Inf
    smallest_compress <- NULL
    for (compress in compress_values) {
      do_compress(env, compress, the_path)
      new_info_ <- file.info(the_path)

      if (new_info_$size < smallest_size) {
        smallest_size <- new_info_$size
        smallest_compress <- compress
      }
    }

    do_compress(env, smallest_compress, the_path)

    new_info <- rbind(
      new_info,
      cbind(
        file.info(the_path),
        method = smallest_compress
      )
    )
  }
  df <- data.frame(
    path = rownames(old_info),
    size_old = old_info$size,
    size_new = new_info$size,
    method = new_info$method
  )
  df$"new/old" <- round(df$size_new / df$size_old * 1000) / 10
  df
}


#' @describeIn inlabru-deprecated Row-wise Kronecker products
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_row_kron()].
#'
#' Takes two Matrices and computes the row-wise Kronecker product.  Optionally
#' applies row-wise weights and/or applies an additional 0/1 row-wise Kronecker
#' matrix product.
#'
#' Returns a `Matrix::sparseMatrix` object.
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
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export row_kron
#' @keywords internal
row_kron <- function(M1, M2, repl = NULL, n.repl = NULL, weights = NULL) {
  lifecycle::deprecate_warn(
    "2.10.0",
    "row_kron()",
    "fmesher::fm_row_kron()"
  )
  fm_row_kron(M1 = M1, M2 = M2, repl = repl, n.repl = n.repl, weights = weights)
}
