bru_check_object_bru <- function(object,
                                 new_version = getNamespaceVersion("inlabru")) {
  object <-
    bru_info_upgrade(
      object,
      new_version = new_version
    )
  object
}

# Used for upgrading from versions <= 2.7.0.9017 to >= 2.7.0.2021
bru_used_upgrade_2.7.0.9017_to_2.7.0.9021 <- function(lhoods, labels) {
  for (k in seq_along(lhoods)) {
    if (is.null(lhoods[[k]][["used"]]) &&
      is.null(lhoods[[k]][["used_components"]])) {
      used <- bru_used(
        NULL,
        effect = lhoods[[k]][["include_components"]],
        latent = lhoods[[k]][["include_latent"]],
        effect_exclude = lhoods[[k]][["exclude_components"]],
      )

      lhoods[[k]][["include_components"]] <- NULL
      lhoods[[k]][["exclude_components"]] <- NULL
      lhoods[[k]][["include_latent"]] <- NULL
    } else if (!is.null(lhoods[[k]][["used_components"]])) {
      used <- bru_used(
        NULL,
        effect = lhoods[[k]][["effect"]],
        latent = lhoods[[k]][["latent"]]
      )
      lhoods[[k]][["used_components"]] <- NULL
    } else {
      used <- bru_used(lhoods[[k]])
    }

    used <- bru_used_update(used, labels = labels)
    lhoods[[k]][["used"]] <- used
  }
  lhoods
}


bru_info_upgrade_functions <- function() {
  list(
    "2.1.14.901" = function(object) {
      # Ensure INLA_version stored
      if (is.null(object[["INLA_version"]])) {
        object[["INLA_version"]] <- "0.0.0"
      }
      # Check for component$mapper as a bm_multi
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        cmp[["mapper"]] <-
          bm_multi(
            list(
              main = cmp$main$mapper,
              group = cmp$group$mapper,
              replicate = cmp$replicate$mapper
            ),
            simplify = TRUE
          )
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object
    },
    "2.5.3.9003" = function(object) {
      # Check that likelihoods store 'weights'
      for (k in seq_along(object[["lhoods"]])) {
        lhood <- object[["lhoods"]][[k]]
        if (is.null(lhood[["weights"]])) {
          lhood[["weights"]] <- 1
          object[["lhoods"]][[k]] <- lhood
        }
      }
      object
    },
    "2.5.3.9005" = function(object) {
      # Make sure component$mapper is a bm_pipe
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        # Convert offset mappers to const mappers
        if (inherits(cmp$main$mapper, c("bm_offset", "bru_mapper_offset"))) {
          cmp$main$mapper <- bm_const()
        }
        cmp[["mapper"]] <-
          bm_pipe(
            list(
              mapper = bm_multi(
                list(
                  main = cmp$main$mapper,
                  group = cmp$group$mapper,
                  replicate = cmp$replicate$mapper
                ),
                simplify = TRUE
              ),
              scale = bm_scale()
            )
          )
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object
    },
    "2.6.0.9000" = function(object) {
      # Make sure component$mapper is a bm_pipe
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        cmp[["mapper"]] <-
          bm_pipe(
            list(
              mapper = bm_multi(
                list(
                  main = cmp$main$mapper,
                  group = cmp$group$mapper,
                  replicate = cmp$replicate$mapper
                ),
                simplify = TRUE
              ),
              scale = bm_scale()
            )
          )
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object
    },
    "2.7.0.9010" = function(object) {
      # Make sure component$group/replicate$input isn't NULL
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        if (identical(deparse(cmp[["group"]][["input"]][["input"]]), "NULL")) {
          cmp[["group"]][["input"]][["input"]] <- expression(1L)
        }
        if (
          identical(
            deparse(cmp[["replicate"]][["input"]][["input"]]),
            "NULL"
          )
        ) {
          cmp[["replicate"]][["input"]][["input"]] <- expression(1L)
        }
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object
    },
    "2.7.0.9016" = function(object) {
      # Convert old style include_component/allow_latent to intermediate
      # include_component/include_latent format

      # Update include/exclude information to limit it to existing components
      for (k in seq_along(object[["lhoods"]])) {
        if (isTRUE(object[["lhoods"]][[k]][["allow_latent"]])) {
          # Force inclusion of all components
          object[["lhoods"]][[k]][["include_latent"]] <- NULL
        } else {
          if (is.null(object[["lhoods"]][[k]][["include_latent"]])) {
            object[["lhoods"]][[k]][["include_latent"]] <- character(0)
          }
        }
      }
      object[["lhoods"]] <-
        bru_used_upgrade_2.7.0.9017_to_2.7.0.9021(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object
    },
    "2.7.0.9017" = function(object) {
      # Convert old style include_component/allow_latent to new
      # bru_used format

      # Update include/exclude information to limit it to existing components
      for (k in seq_along(object[["lhoods"]])) {
        if (is.null(object[["lhoods"]][[k]][["used_components"]])) {
          if (isTRUE(object[["lhoods"]][[k]][["allow_latent"]])) {
            # Force inclusion of all components
            object[["lhoods"]][[k]][["include_latent"]] <- NULL
          } else {
            if (is.null(object[["lhoods"]][[k]][["include_latent"]])) {
              object[["lhoods"]][[k]][["include_latent"]] <- character(0)
            }
          }
        }
      }
      object[["lhoods"]] <-
        bru_used_upgrade_2.7.0.9017_to_2.7.0.9021(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object
    },
    "2.7.0.9021" = function(object) {
      # Make sure 'used' components format is properly stored

      if (is.null(object[["lhoods"]][["is_additive"]])) {
        object[["lhoods"]][["is_additive"]] <-
          object[["lhoods"]][["linear"]]
      }

      object[["lhoods"]] <-
        bru_used_upgrade_2.7.0.9017_to_2.7.0.9021(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object
    },
    "2.10.1.9007" = function(object_full) {
      # Update timings info to difftime format

      warning(
        glue(
          "
           From 2.10.1.9007, elapsed time is in Elapsed, Time is CPU time.
             Copying old elapsed time to both Elapsed and Time,
             setting System to zero;
             Do not over-interpret."
        ),
        immediate. = TRUE
      )
      conversion <- as.difftime(
        object_full[["bru_timings"]][["Time"]],
        units = "secs"
      )
      object_full[["bru_timings"]][["Time"]] <- conversion
      object_full[["bru_timings"]][["System"]] <-
        as.difftime(rep(0.0, length(conversion)), units = "secs")
      object_full[["bru_timings"]][["Elapsed"]] <- conversion

      object_full
    },
    "2.10.1.9012" = function(object_full) {
      # Update log format to new format

      object_full[["bru_iinla"]][["log"]] <-
        bru_log_new(
          object_full[["bru_iinla"]][["log"]][["log"]],
          bookmarks = object_full[["bru_iinla"]][["log"]][["bookmarks"]]
        )

      object_full
    },
    "2.12.0.9014" = function(object) {
      # Add bru_obs class name, to support old bru_like objects.
      # Technically for 2.12.0.9017, but the upgrade code needs it for 2014
      if (!is.null(object[["lhoods"]])) {
        object[["lhoods"]] <-
          lapply(object[["lhoods"]], function(x) {
            class(x) <- c("bru_obs", class(x))
            x
          })
        class(object[["lhoods"]]) <- c("bru_obs_list", "list")
      }

      # Update is_additive/linear storage
      object[["lhoods"]][["is_additive"]] <-
        object[["lhoods"]][["linear"]]

      # Update predictor expression
      # object[["lhoods"]] <-
      #   bru_used_update(
      #     object[["lhoods"]],
      #     labels = names(object[["model"]][["effects"]])
      #   )

      object
    },
    "2.12.0.9017" = function(object) {
      # Update bru_like class names to bru_obs
      if (!is.null(object[["lhoods"]])) {
        object[["lhoods"]] <-
          lapply(object[["lhoods"]], function(x) {
            class(x) <- "bru_obs"
            x
          })
        class(object[["lhoods"]]) <- c("bru_obs_list", "list")
      }

      eff <- object[["model"]][["effects"]]
      for (k in seq_along(object[["model"]][["effects"]])) {
        class(eff[[k]][["main"]]) <- "bru_subcomp"
        class(eff[[k]][["group"]]) <- "bru_subcomp"
        class(eff[[k]][["replicate"]]) <- "bru_subcomp"
        class(eff[[k]]) <- "bru_comp"
      }
      class(eff) <- c("bru_comp_list", "list")
      object[["model"]][["effects"]] <- eff

      object
    },
    "2.13.0.9012" = function(object) {
      # Update bru_input input&layer to quosures
      quosure_update <- function(x, env) {
        if (is.null(x)) {
          return(NULL)
        }
        rlang::as_quosure(x, env = env)
      }
      eff <- object[["model"]][["effects"]]
      for (k in seq_along(object[["model"]][["effects"]])) {
        for (part in c("main", "group", "replicate")) {
          if (!is.null(eff[[k]][[part]][["input"]])) {
            for (input_part in c("input", "layer")) {
              eff[[k]][[part]][["input"]][[input_part]] <-
                quosure_update(
                  eff[[k]][[part]][["input"]][[input_part]],
                  eff[[k]][["env_extra"]]
                )
            }
          }
        }
        part <- "weights"
        if (!is.null(eff[[k]][[part]])) {
          for (input_part in c("input", "layer")) {
            eff[[k]][[part]][[input_part]] <-
              quosure_update(
                eff[[k]][[part]][[input_part]],
                eff[[k]][["env_extra"]]
              )
          }
        }
      }
      object[["model"]][["effects"]] <- eff

      object
    },
    "2.13.0.9015" = function(object) {
      # Update lhoods to streamlined storage format
      if (!is.null(object[["lhoods"]])) {
        object[["lhoods"]] <-
          lapply(object[["lhoods"]], function(x) {
            for (what in c("E", "Ntrials", "weights", "scale")) {
              BRU_what <- paste0("BRU_", what)
              if (is.null(x[["response_data"]][[BRU_what]])) {
                x[["response_data"]][[BRU_what]] <- x[[what]]
              }
              x[[what]] <- NULL
            }
            x
          })
        class(object[["lhoods"]]) <- c("bru_obs_list", "list")
      }

      object
    },
    "2.13.0.9017" = function(object) {
      # Update lhoods to new bru_pred_expr object storage format
      if (!is.null(object[["lhoods"]])) {
        lhood_upgrader <- function(x) {
          pred_text <- rlang::expr_deparse(x[["expr"]])
          pred_text <- gsub(
            pattern = "^expression\\(",
            replacement = "(",
            x = pred_text,
            fixed = FALSE
          )
          pred_expr <- x[["expr"]]
          x[["pred_expr"]] <- structure(
            list(
              pred_text = pred_text,
              pred_expr = rlang::parse_expr(pred_text),
              is_additive = x[["is_additive"]],
              is_rowwise = !x[["allow_combine"]],
              is_linear = x[["linear"]],
              used = x[["used"]],
              .envir = x[[".envir"]],
              resp_text = as.character(x[["formula"]])[2]
            ),
            class = "bru_pred_expr"
          )
          x
        }
        object[["lhoods"]] <-
          lapply(object[["lhoods"]], lhood_upgrader)
        class(object[["lhoods"]]) <- c("bru_obs_list", "list")
      }

      object
    },
    "2.13.0.9018" = function(object) {
      # Update component input storage to mapper-attached storage
      if (!is.null(object[["model"]][["effects"]])) {
        object[["model"]][["effects"]] <-
          lapply(object[["model"]][["effects"]], function(x) {
            for (subcomp in c("main", "group", "replicate")) {
              if (
                !is.null(x[[subcomp]][["input"]]) &&
                  inherits(x[[subcomp]][["input"]], "bru_input") &&
                  is.null(x[[subcomp]][["mapper"]][[".input"]])
              ) {
                x[[subcomp]][["mapper"]][[".input"]] <- x[[subcomp]][["input"]]
                x[[subcomp]][["input"]] <- NULL
              }
            }
            if (
              !is.null(x[["marginal"]]) &&
                is.null(x[["marginal"]][[".input"]])
            ) {
              x[["marginal"]][[".input"]] <- new_bru_input(NULL)
            }
            if (
              !is.null(x[["weights"]]) &&
                inherits(x[["weights"]], "bru_input")
            ) {
              weights_input <- x[["weights"]]
              x[["weights"]] <- bm_scale()
              x[["weights"]][[".input"]] <- weights_input
            }
            x <- bru_comp_update_mapper(x)
            x <- bru_compat_pre_2_14_bru_comp(x)
            x
          })
        class(object[["model"]][["effects"]]) <- c("bru_comp_list", "list")
      }

      object
    },
    "2.13.0.9036" = function(object) {
      # Update component input storage to mapper-attached storage
      if (!is.null(object[["model"]][["effects"]])) {
        object[["model"]][["effects"]] <-
          lapply(object[["model"]][["effects"]], function(x) {
            x <- bru_comp_update_mapper(x)
            x <- bru_compat_pre_2_14_bru_comp(x)
            x
          })
        class(object[["model"]][["effects"]]) <- c("bru_comp_list", "list")
      }

      object
    },
    "2.14.0.9002" = function(object) {
      # Update bm_factor mappers to proper levels for cases where previously
      # incorrect data reading resulted in an incorrect (but matching) mapper
      the_inputs <- bru_input(
        object[["lhoods"]],
        object[["model"]][["effects"]]
      )
      if (!is.null(object[["model"]][["effects"]])) {
        object[["model"]][["effects"]] <-
          lapply(object[["model"]][["effects"]], function(x) {
            lab <- x[["label"]]
            mapper <- x[["main"]][["mapper"]]
            if (inherits(mapper, "bm_factor")) {
              inp_levels <- lapply(the_inputs, function(xx) {
                if (is.factor(xx[[lab]][["core"]][["main"]])) {
                  levels(xx[[lab]][["core"]][["main"]])
                } else {
                  NULL
                }
              })
              inp_levels_non_null <-
                vapply(inp_levels, function(x) !is.null(x), logical(1))
              if (any(inp_levels_non_null)) {
                inp_levels <- unique(inp_levels[inp_levels_non_null])
                if (length(unique(inp_levels)) > 1) {
                  warning(glue::glue(
                    "Inconsistent input factor levels for {lab}.
                   Levels found: {paste(unique(inp_levels), collapse = '; ')}.
                   Using levels from first input with non-null levels."
                  ))
                }
                inp_levels <- inp_levels[[1]]
                if (length(inp_levels) != length(mapper$levels)) {
                  warning(glue::glue(
                    "Number of levels for {lab} in mapper ",
                    "({length(mapper$levels)}) does ",
                    "not match number of levels ",
                    "in input ({length(inp_levels)}). ",
                    "Using levels from input."
                  ))
                }
                mapper[["levels"]] <- inp_levels
                mapper[["n"]] <- length(inp_levels)
                if (identical(mapper[["factor_mapping"]], "contrast")) {
                  mapper[["n"]] <- mapper[["n"]] - 1L
                }
                x[["main"]][["mapper"]] <- mapper
              }
            }
            x <- bru_comp_update_mapper(x)
            x
          })
        class(object[["model"]][["effects"]]) <- c("bru_comp_list", "list")
      }

      object
    }
  )
}
bru_info_upgrade <- function(object,
                             new_version = getNamespaceVersion("inlabru")) {
  object_full <- object
  object <- object[["bru_info"]]
  msg <- NULL
  if (!is.list(object)) {
    msg <- "Not a list"
  } else if (!inherits(object, "bru_info")) {
    msg <- "Not a bru_info object"
  } else if (is.null(object[["inlabru_version"]])) {
    msg <- "`inlabru_version` is missing"
  }
  if (!is.null(msg)) {
    stop(glue(
      "bru_info part of the object can't be converted to `bru_info`; {msg}"
    ))
  }

  old_ver <- object[["inlabru_version"]]
  if (utils::compareVersion(new_version, old_ver) > 0) {
    warning(
      glue("Old bru_info object version {old_ver} detected.
            Attempting upgrade to version {new_version}.")
    )

    upgrade_functions <- bru_info_upgrade_functions()

    message("Detected bru_info version ", old_ver)
    for (ver in names(upgrade_functions)) {
      if (utils::compareVersion(ver, old_ver) > 0) {
        message("Upgrading bru_info to ", ver)
        if ("object_full" %in% names(formals(upgrade_functions[[ver]]))) {
          object_full <- upgrade_functions[[ver]](object_full = object_full)
        } else {
          object <- upgrade_functions[[ver]](object)
        }
        object[["inlabru_version"]] <- ver
        object_full[["bru_info"]] <- object
      }
    }

    object[["inlabru_version"]] <- new_version
    object_full[["bru_info"]] <- object
    message(glue("Upgraded bru_info to {new_version}"))
  }
  object_full
}
