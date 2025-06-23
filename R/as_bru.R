as_bru <- function(x, ...) {
  UseMethod("as_bru", x)
}

as_bru.bru <- function(x, ...) {
  x
}

as_bru.inla <- function(x, ...) {
  if (!inherits(x, "inla")) {
    stop("Input must be an 'inla' object.")
  }

  # Convert component definitions
  form <- x[[".args"]][["formula"]]
  form_char <- as.character(form)
  response_char <- form_char[2]
  response_expr <- parse(text = response_char)
  response_value <- eval(response_expr, envir = x[[".args"]][["data"]])

  family <- x[[".args"]][["family"]]
  link <- x[[".args"]][["control.predictor"]][["link"]]
  if (is.null(link)) {
    # Detect family data split
    if (is.matrix(response_value) || is.data.frame(response_value)) {
      family_idx <- vapply(
        seq_len(NROW(response_value)),
        function(i) {
          fam <- which(!is.na(response_value[i, ]))
          if (length(fam) == 0) {
            return(NA_integer_)
          } else if (length(fam) > 1) {
            return(-fam[1])
          } else {
            return(fam)
          }
        },
        integer(1)
      )
      if (any(is.na(family_idx))) {
        warning("Some rows in the response value matrix/data frame have no family defined.")
        family_idx[is.na(family_idx)] <- 1L
      }
      if (any(family_idx < 0)) {
        warning("Some rows in the response value matrix/data frame have ambiguous family.")
        family_idx[family_idx < 0] <- -family_idx[family_idx < 0]
      }
    } else {
      stop("Unhandled response value class: ", paste0(class(response_value), collapse = ", "))
    }
  }

  # Split response into each family.
  # Easy for tibbles/data.frame, need special handling for inla.mdata, inla.surv, etc
  if (is.matrix(response_value) || is.data.frame(response_value)) {
    response_value <- lapply(sort(unique(family_idx)),
                             function(i) {
                               response_value[family_idx == i, i, drop = TRUE]
                             })
  } else {
    stop("Unhandled response value class: ", paste0(class(response_value), collapse = ", "))
  }

  # Split data into each family.
  # Easy for tibbles/data.frame, difficult for list() data
  if (is.data.frame(x[[".args"]][["data"]])) {
    data <- lapply(sort(unique(family_idx)),
                             function(i) {
                               x[[".args"]][["data"]][family_idx == i, , drop = FALSE]
                             })
  } else {
    stop("Unhandled data class: ", paste0(class(x[[".args"]][["data"]]), collapse = ", "))
  }

  x
}
