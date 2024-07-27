# More on special bru likelihoods - aggregated poisson count (apc)
if (family == "apc") {
  if (is.null(data)) {
    stop("You called like() with family='apc' but no 'data' argument was supplied.")
  }

  if (is.null(ips)) {
    ips <- apmaker(
      samplers = samplers,
      domain = domain,
      dnames = response,
      int.args = options[["bru_int_args"]]
    )
  }

  if (length(E) > 1) {
    warning("Exposure/effort parameter E should be a scalar for likelihood 'apc'.")
  }

  ips_is_Spatial <- inherits(ips, "Spatial")
  if (ips_is_Spatial) {
    ips_coordnames <- sp::coordnames(ips)
    ips_crs <- fm_sp_get_crs(ips)
  }
  data_is_Spatial <- inherits(data, "Spatial")
  if (data_is_Spatial) {
    data_coordnames <- sp::coordnames(data)
    data_crs <- fm_sp_get_crs(data)
    if (ips_is_Spatial) {
      new_coordnames <- data_coordnames[seq_len(min(
        length(ips_coordnames),
        length(data_coordnames)
      ))]
      new_coordnames[new_coordnames %in% ""] <-
        paste0(
          "BRU_dummy_coordinate_",
          seq_along(new_coordnames)
        )[new_coordnames %in% ""]
      ips_coordnames <- paste0(
        "BRU_dummy_coordinate_",
        seq_along(ips_coordnames)
      )
      data_coordnames <- paste(
        "BRU_dummy_coordinate_",
        seq_along(data_coordnames)
      )
      ips_coordnames[seq_along(new_coordnames)] <- new_coordnames
      data_coordnames[seq_along(new_coordnames)] <- new_coordnames
      sp::coordnames(ips) <- ips_coordnames
      sp::coordnames(data) <- data_coordnames

      # TODO: check that the crs info is the same
    }
  }
  data <- as.data.frame(data)
  if (!is.null(response_data)) {
    warning("Ignoring non-null response_data input for 'apc' likelihood")
  }
  ips <- as.data.frame(ips)
  dim_names <- intersect(names(data), names(ips))
  if (identical(options[["bru_compress_apc"]], TRUE)) {
    allow_combine <- TRUE
    response_data <- data.frame(
      BRU_E = c(
        0,
        E * ips[["weight"]]
      ),
      BRU_response_apc = c(
        NROW(data),
        rep(0, NROW(ips))
      )
    )
    if (!linear) {
      expr_text <- as.character(formula)[length(as.character(formula))]
      expr_text <- paste0(
        "{BRU_eta <- ", expr_text, "\n",
        " c(mean(BRU_eta[BRU_aggregate]), BRU_eta[!BRU_aggregate])}"
      )
    } else {
      expr_text <- paste0(
        "{BRU_eta <- BRU_EXPRESSION\n",
        " c(mean(BRU_eta[BRU_aggregate]), BRU_eta[!BRU_aggregate])}"
      )
    }
    expr <- parse(text = expr_text)
    data <- rbind(
      cbind(data[dim_names], BRU_aggregate = TRUE),
      cbind(ips[dim_names], BRU_aggregate = FALSE)
    )
    formula
  } else {
    response_data <- data.frame(
      BRU_E = c(
        rep(0, NROW(data)),
        E * ips[["weight"]]
      ),
      BRU_response_apc = c(
        rep(1, NROW(data)),
        rep(0, NROW(ips))
      )
    )
    data <- rbind(
      data[dim_names],
      ips[dim_names]
    )
  }
  if (ips_is_Spatial) {
    non_coordnames <- setdiff(names(data), data_coordnames)
    data <- sp::SpatialPointsDataFrame(
      coords = data[new_coordnames],
      data = data[non_coordnames],
      proj4string = data_crs,
      match.ID = FALSE
    )
  }

  response <- "BRU_response_apc"
  inla.family <- "poisson"
  E <- response_data[["BRU_E"]]
}
