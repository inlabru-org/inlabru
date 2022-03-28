#' Standardise inla hyperparameter names
#'
#' The inla hyperparameter output uses parameter names that can include
#' whitespace and special characters. This function replaces those characters
#' with underscores.
#'
#' @param x character vector; names to be standardised
#' @returns A character vector with standardised names
#' @examples
#' bru_standardise_names("Precision for the Gaussian observations")
#' @export
bru_standardise_names <- function(x) {
  new_names <- vapply(
    x,
    function(x) {
      gsub("[-() ]", "_", x = x, fixed = FALSE)
    },
    "name"
  )
  not_ok <- grepl("__", x = new_names)
  while (any(not_ok)) {
    new_names[not_ok] <- vapply(
      new_names[not_ok],
      function(x) {
        gsub("__", "_", x = x, fixed = FALSE)
      },
      "name"
    )
    not_ok <- grepl("__", x = new_names)
  }
  new_names
}



#' Prediction from fitted inla model
#'
#' This method is not supported for plain inla objects. Please see [predict.bru()] instead.
#'
#' @aliases predict.inla
#' @export
#' @param object A `bru` object obtained by calling [bru()] or [lgcp()].
#' @param ... Arguments passed on to [predict.bru()].
#' @return A `prediction` object.
#' @seealso [predict.inla()]
#' @keywords internal
predict.inla <- function(object, ...) {
  stop("predict() is not supported for plain inla objects. See https://github.com/fbachl/inlabru/issues/78")
}



#' Sampling based on inla posteriors
#'
#' @description
#' This method is not supported for plain inla objects. Please see [generate.bru()] instead.
#' See https://github.com/fbachl/inlabru/issues/78 for more information.
#' @aliases generate.inla
#' @rdname predict.inla
#' @export
#' @seealso [generate.inla()]
#' @param object An `inla` object obtained by calling `INLA::inla()`.
#' @param ... additional arguments passed on to[generate.bru].
#'
#' @return List of generated samples
#' @seealso [predict.bru()]
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
generate.inla <- function(object,
                          ...) {
  stop("generate() is no longer supported for plain inla objects. See https://github.com/fbachl/inlabru/issues/78")
}





inla_result_latent_idx <- function(result) {
  do.call(
    c,
    lapply(
      result$misc$configs$contents$tag,
      function(x) {
        idx <-
          list(
            result$misc$configs$contents$start[
              result$misc$configs$contents$tag == x
            ] - 1 + seq_len(
              result$misc$configs$contents$length[
                result$misc$configs$contents$tag == x
              ]
            )
          )
        names(idx) <- x
        idx
      }
    )
  )
}

#' Extract a summary property from all results of an inla result
#'
#' @param result an `inla` result object
#' @param property character; "mean", "sd", "mode", or some other column
#' identifier for inla result `$summary.fixed`, `$summary.random$label`, and
#' `$summary.hyperpar`, or "joint_mode". For "joint_mode", the joint latent mode
#' is extracted, and the joint hyperparameter mode, in the internal scale.
#' For "predictor_sd" the posterior standard deviations of the linear predictor
#' are returned.
#' @param internal_hyperpar logical; if `TRUE`, use internal scale for
#' hyperparamter properties. Default is `FALSE`, except when `property` is
#' "joint_mode" which forces `internal_hyperpar=TRUE`.
#' @return named list for each estimated fixed effect coefficient,
#' random effect vector, and hyperparameter. The hyperparameter names are
#' standardised with [bru_standardise_names()]
#' @keywords internal
extract_property <- function(result, property,
                             internal_hyperpar = FALSE) {
  stopifnot(inherits(result, "inla"))
  ret <- list()

  if (property == "joint_mode") {
    mode_idx <- inla_result_latent_idx(result)
    for (label in c(
      rownames(result$summary.fixed),
      names(result$summary.random)
    )) {
      ret[[label]] <- result$mode$x[mode_idx[[label]]]
    }
    theta_names <- bru_standardise_names(result$mode$theta.tags)
    for (idx in seq_along(theta_names)) {
      ret[[theta_names[idx]]] <- result$mode$theta[idx]
    }
    return(ret)
  } else if (property == "predictor_sd") {
    idx <- grepl(
      pattern = "^APredictor",
      x = rownames(result$summary.linear.predict)
    )
    ret <- result$summary.linear.predictor$sd[idx]
    return(ret)
  }

  for (label in rownames(result$summary.fixed)) {
    ret[[label]] <- result$summary.fixed[label, property]
  }

  for (label in names(result$summary.random)) {
    ret[[label]] <- result$summary.random[[label]][, property]
  }

  if (internal_hyperpar) {
    for (label in rownames(result$internal.summary.hyperpar)) {
      new.label <- bru_standardise_names(label)
      ret[[new.label]] <- result$internal.summary.hyperpar[label, property]
    }
  } else {
    for (label in rownames(result$summary.hyperpar)) {
      new.label <- bru_standardise_names(label)
      ret[[new.label]] <- result$summary.hyperpar[label, property]
    }
  }

  fac.names <- names(result$model$effects)[
    vapply(
      result$model$effects,
      function(e) {
        identical(e$type, "factor")
      },
      TRUE
    )
  ]
  # TODO: Consider whether to extract/convert effect$type == "factor" models
  # from random effects into fixed effects
  #
  # Old code:
  # For factors we add a data.frame with column names equivalent to the factor levels
  #    for (name in fac.names) {
  #      tmp <- unlist(ret[startsWith(names(ret), name)])
  #      names(tmp) <- lapply(names(tmp), function(nm) {
  #        substring(nm, nchar(name) + 1)
  #      })
  #      ret[[name]] <- tmp
  #    }

  ret
}



##
## A wrapper for inla.posterior.sample()
##
## Converts each sample into a list of sub-samples representing the
## latent variables
##
## Example: samples[[1]]$somefield is the value of the field "somefield"
##          in the first sample
##

inla.posterior.sample.structured <- function(result, n, seed = NULL,
                                             num.threads = NULL, ...) {
  if (!is.null(seed) && (seed != 0L)) {
    num.threads <- "1:1"
  }
  # Workaround for older versions of INLA
  if ("hyper.user.scale" %in% formalArgs(INLA::inla.posterior.sample)) {
    samples <- INLA::inla.posterior.sample(
      n = n,
      result = result,
      seed = seed,
      num.threads = num.threads
    )
  } else if ("parallel.configs" %in% formalArgs(INLA::inla.posterior.sample)) {
    samples <- INLA::inla.posterior.sample(
      n = n,
      result = result,
      seed = seed,
      intern = FALSE,
      num.threads = num.threads,
      parallel.configs = FALSE
    )
  } else {
    samples <- INLA::inla.posterior.sample(
      n = n,
      result = result,
      seed = seed,
      intern = FALSE,
      num.threads = num.threads
    )
  }

  ssmpl <- list()
  for (i in seq_along(samples)) {
    smpl.latent <- samples[[i]]$latent
    smpl.hyperpar <- samples[[i]]$hyperpar
    vals <- list()

    # Extract simulated predictor and fixed effects
    for (name in unique(c("Predictor", result$names.fixed))) {
      vals[[name]] <- extract.entries(name, smpl.latent)
    }

    # For fixed effects that were modeled via factors we attach an extra vector holding the samples
    fac.names <- names(result$model$effects)[
      vapply(
        result$model$effects,
        function(e) {
          identical(e$model, "factor")
        },
        TRUE
      )
    ]
    for (name in fac.names) {
      vals[[name]] <- smpl.latent[startsWith(rownames(smpl.latent), name), ]
      names(vals[[name]]) <- lapply(names(vals[[name]]), function(nm) {
        substring(nm, nchar(name) + 1)
      })
    }


    # Extract simulated latent variables.
    # If the model is "clinear", however, we might extract the realisations
    # from the hyperpar field. TODO: check if all the special models now have
    # their results available as latent random effects, and avoid special code,
    # since the hyperpar name definition has changed
    if (length(result$summary.random) > 0) {
      for (k in seq_along(result$summary.random)) {
        name <- unlist(names(result$summary.random[k]))
        model <- result$model.random[k]
        #        if (!(model == "Constrained linear")) {
        vals[[name]] <- extract.entries(name, smpl.latent)
        #        }
        #        else {
        #         vals[[name]] <- smpl.hyperpar[paste0("Beta for ", name)]
        #        }
      }
    }
    if (length(smpl.hyperpar) > 0) {
      ## Sanitize the variable names; replace problems with "_".
      ## Needs to handle whatever INLA uses to describe the hyperparameters.
      ## Known to include " " and "-" and potentially "(" and ")".
      names(smpl.hyperpar) <- bru_standardise_names(names(smpl.hyperpar))
    }
    ssmpl[[i]] <- c(vals, smpl.hyperpar)
  }

  #
  # Return
  #
  return(ssmpl)
}

extract.entries <- function(name, smpl) {
  ename <- gsub("\\.", "\\\\.", name)
  ename <- gsub("\\(", "\\\\(", ename)
  ename <- gsub("\\)", "\\\\)", ename)
  ptn <- paste("^", ename, "[\\:]*[\\.]*[0-9]*[\\.]*[0-9]*$", sep = "")
  return(smpl[grep(ptn, rownames(smpl))])
}

# Expand observation vectors/matrices in stacks into to a multicolumn matrix for multiple likelihoods
#
# @aliases inla.stack.mexpand
# @name inla.stack.mexpand
# @export
# @param ... List of stacks that contain vector observations
#            (existing multilikelihood observation matrices are also permitted)
# @param old.names A vector of strings with the names of the observation vector/matrix for each stack.
#        If a single string, this is assumed for all the stacks. (default "BRU.response")
# @param new.name The name to be used for the expanded observation matrix,
#        possibly the same as an old name. (default "BRU.response")
# @return a list of modified stacks with multicolumn observations
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk} and Finn Lindgren \email{finn.lindgren@@gmail.com}
#

inla.stack.mexpand <- function(...,
                               old.names = "BRU.response",
                               new.name = "BRU.response") {
  stacks <- list(...)
  if (length(old.names) == 1) {
    old.names <- rep(old.names, length(stacks))
  }
  y.cols <- unlist(lapply(seq_along(stacks),
    function(x, stacks, old.names) {
      LHS <- INLA::inla.stack.LHS(stacks[[x]])[[old.names[x]]]
      ifelse(is.vector(LHS), 1, ncol(LHS))
    },
    stacks = stacks, old.names = old.names
  ))
  y.offset <- c(0, cumsum(y.cols))
  y.cols.total <- sum(y.cols)
  for (j in seq_along(stacks)) {
    LHS <- INLA::inla.stack.LHS(stacks[[j]])
    RHS <- INLA::inla.stack.RHS(stacks[[j]])
    A <- INLA::inla.stack.A(stacks[[j]])
    # Access the raw tag indexing information
    tags <- list(
      data = stacks[[j]]$data$index,
      effects = stacks[[j]]$effects$index
    )

    # Expand the observation vector/matrix into a multilikelihood observation matrix:
    y.rows <- ifelse(is.vector(LHS[[old.names[j]]]),
      length(LHS[[old.names[j]]]),
      nrow(LHS[[old.names[j]]])
    )
    LHS[[new.name]] <-
      cbind(
        matrix(NA, nrow = y.rows, ncol = y.offset[j]),
        LHS[[old.names[j]]],
        matrix(NA, nrow = y.rows, ncol = y.cols.total - y.offset[j + 1])
      )

    # Create the modified stack, with model compression disabled to prevent modifications:
    stacks[[j]] <-
      INLA::inla.stack.sum(data = LHS, A = A, effects = RHS, compress = FALSE, remove.unused = FALSE)
    # Since the row indexing is unchanged, copy the tag index information:
    stacks[[j]]$data$index <- tags$data
    stacks[[j]]$effects$index <- tags$effects
  }
  stacks
}

# Stack multiple exposures
# Obsolete. Do not use. Exposure vector stacking is handled automatically by inla.stack
#
# @aliases inla.stack.e
# @param ... observation vectors
# @return e observation vector
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}

inla.stack.e <- function(...) {
  stop("This function is obsolete and should not be used.")
}


#' Join stacks intended to be run with different likelihoods
#'
#' @param ... List of stacks that contain vector observations
#'            (existing multi-likelihood observation matrices are also permitted)
#' @param compress If `TRUE`, compress the model by removing duplicated rows of
#' effects, replacing the corresponding A-matrix columns with a single column
#' containing the sum.
#' @param remove.unused	If `TRUE`, compress the model by removing rows of
#' effects corresponding to all-zero columns in the A matrix (and removing those columns).
#' @param old.names A vector of strings with the names of the observation
#'   vector/matrix for each stack.
#'   If a single string, this is assumed for all the stacks. (default "BRU.response")
#' @param new.name The name to be used for the expanded observation matrix,
#'        possibly the same as an old name. (default "BRU.response")
#' @export
#' @rdname inla.stack.mjoin
#'

inla.stack.mjoin <- function(..., compress = TRUE, remove.unused = TRUE,
                             old.names = "BRU.response", new.name = "BRU.response") {
  stacks <- inla.stack.mexpand(..., old.names = old.names, new.name = new.name)
  do.call(
    INLA::inla.stack.join,
    c(stacks, list(compress = compress, remove.unused = remove.unused))
  )
}





plotmarginal.inla <- function(result,
                              varname = NULL,
                              index = NULL,
                              link = function(x) {
                                x
                              },
                              add = FALSE,
                              ggp = TRUE,
                              lwd = 3,
                              ...) {
  requireNamespace("ggplot2")
  vars <- variables.inla(result)
  ovarname <- varname

  if (varname %in% c(result$names.fixed, rownames(result$summary.hyperpar)) ||
    (!is.null(index) && (varname %in% names(result$summary.random)))) {
    if (varname %in% rownames(vars) &&
      vars[varname, "type"] == "fixed") {
      marg <- INLA::inla.tmarginal(link, result$marginals.fixed[[varname]])
    } else if (varname %in% names(result$summary.random)) {
      marg <- INLA::inla.tmarginal(link, result$marginals.random[[varname]][[index]])
      ovarname <- paste0(ovarname, " ", index)
      if (ovarname %in% rownames(vars)) {
        if (!identical(vars[ovarname, "ID"], as.character(index - 1))) {
          # Use factor level name
          ovarname <- vars[ovarname, "ID"]
        }
      }
    } else if (varname %in% rownames(vars) &&
      vars[varname, "type"] == "hyperpar") {
      marg <- INLA::inla.tmarginal(link, result$marginals.hyperpar[[varname]])
    }
    uq <- INLA::inla.qmarginal(0.975, marg)
    uqy <- INLA::inla.dmarginal(uq, marg)
    lq <- INLA::inla.qmarginal(0.025, marg)
    lqy <- INLA::inla.dmarginal(lq, marg)
    inner.x <- seq(lq, uq, length.out = 100)
    inner.marg <- data.frame(x = inner.x, y = INLA::inla.dmarginal(inner.x, marg))

    df <- data.frame(marg)
    ggplot2::ggplot(data = df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
      ggplot2::geom_path() +
      ggplot2::geom_ribbon(
        ymin = 0,
        ggplot2::aes(ymax = .data[["y"]]),
        alpha = 0.1
      ) +
      ggplot2::geom_segment(x = lq, y = 0, xend = lq, yend = lqy) +
      ggplot2::geom_segment(x = uq, y = 0, xend = uq, yend = uqy) +
      ggplot2::geom_ribbon(
        data = inner.marg, ymin = 0,
        ggplot2::aes(ymax = .data[["y"]]), alpha = 0.1
      ) +
      ggplot2::xlab(ovarname) +
      ggplot2::ylab("pdf")
  } else {
    df <- result$summary.random[[varname]]
    colnames(df) <- c("ID", "mean", "sd", "lower", "mid", "upper", "mode", "kld")
    df$mean <- link(df$mean)
    h <- 1e-6
    df$sd <- link(df$sd) * abs((link(df$mean + h) - link(df$mean - h)) / (2 * h))
    df$lower <- link(df$lower)
    df$mid <- link(df$mid)
    df$upper <- link(df$upper)
    df$mode <- link(df$mode)
    p <- ggplot2::ggplot(df, ggplot2::aes(.data[["ID"]], .data[["mode"]]))
    p +
      ggplot2::geom_ribbon(ggplot2::aes(
        ymin = .data[["lower"]],
        ymax = .data[["upper"]]
      ), alpha = 0.1) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(y = .data[["mean"]]), col = 2) +
      ggplot2::geom_line(ggplot2::aes(y = .data[["mid"]]), col = 2, lty = 2) +
      ggplot2::ylab("mode and quantiles") +
      ggplot2::xlab(paste0(varname, " ID"))
  }
}

variables.inla <- function(result, include.random = TRUE) {
  handle.missing <- function(col.names) {
    cbind(data.frame(
      type = character(0),
      model = character(0),
      as.data.frame(
        matrix(NA, 0, length(col.names),
          dimnames = list(c(), col.names)
        )
      )
    ))
  }

  handle.missing.columns <- function(data, col.names) {
    missing.names <- setdiff(col.names, colnames(data))
    if (length(missing.names) > 0) {
      df <- as.data.frame(matrix(NA, nrow(data), length(missing.names),
        dimnames = list(NULL, missing.names)
      ))
      return(cbind(data, df))
    } else {
      return(data)
    }
  }

  handle.data.frame <- function(data, type, model, col.names) {
    ##    rownames(data) ???
    cbind(
      data.frame(
        type = rep(type, nrow(data)),
        model = rep(model, nrow(data))
      ),
      handle.missing.columns(data, col.names)
    )
  }

  ## Get column names, handling possibly empty output:
  fixed.missing <- (is.null(result$summary.fixed) ||
    (ncol(result$summary.fixed) == 0))
  hyperpar.missing <- (is.null(result$summary.hyperpar) ||
    (ncol(result$summary.hyperpar) == 0))
  random.missing <- (is.null(result$summary.random) ||
    (length(result$summary.random) == 0) ||
    (ncol(result$summary.random[[1]]) == 0))
  col.names <- c()
  if (!fixed.missing) {
    col.names <- union(col.names, colnames(result$summary.fixed))
  }
  if (!hyperpar.missing) {
    col.names <- union(col.names, colnames(result$summary.hyperpar))
  }
  if (!random.missing) {
    col.names <- union(col.names, colnames(result$summary.random[[1]]))
  }
  if (fixed.missing) {
    fixed <- handle.missing(col.names)
  } else {
    fixed <- handle.data.frame(
      result$summary.fixed,
      "fixed", "fixed", col.names
    )
  }
  if (hyperpar.missing) {
    hyperpar <- handle.missing(col.names)
  } else {
    hyperpar <- handle.data.frame(
      result$summary.hyperpar,
      "hyperpar", NA, col.names
    )
  }
  if (random.missing) {
    random <- list(handle.missing(col.names))
  } else {
    if (!include.random) {
      random <-
        lapply(
          seq_along(result$summary.random),
          function(x) {
            dat <- data.frame(
              mean = mean(result$summary.random[[x]]$mean),
              sd = mean(result$summary.random[[x]]$sd)
            )
            output <- handle.data.frame(
              dat,
              "random",
              result$model.random[x],
              col.names
            )
            rownames(output) <-
              names(result$summary.random)[x]
            output
          }
        )
    } else {
      random <-
        lapply(
          seq_along(result$summary.random),
          function(x) {
            output <- handle.data.frame(
              result$summary.random[[x]],
              "random",
              result$model.random[x],
              col.names
            )
            rownames(output) <-
              paste(
                names(result$summary.random)[x],
                seq_len(nrow(result$summary.random[[x]]))
              )
            output
          }
        )
    }
  }

  variables <- do.call(rbind, c(list(fixed, hyperpar), random))
  return(variables)
}
