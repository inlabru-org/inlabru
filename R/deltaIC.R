#' @title Summarise DIC from `lgcp` estimates.
#'
#' @description
#' Calculates DIC differences and produces an ordered summary.
#'
#' @param \dots Comma-separated objects inheriting from class `inla` and
#'   obtained from a run of `INLA::inla()`, [bru()] or [lgcp()]
#' @param criterion character vector.
#' If it includes 'DIC', computes DIC differences;
#' 'WAIC' is also allowed, but note that plain WAIC values from inla are not
#' well-defined for point process models. Default 'WAIC'
#'
#' @return A data frame with each row containing the Model name,
#' DIC and Delta.DIC.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla()) {
#'   # Generate some data
#'   input.df <- data.frame(idx = 1:10) |>
#'     dplyr::mutate(
#'       x = cos(idx),
#'       y = rpois(length(x), 5 + 2 * x + rnorm(length(x), mean = 0, sd = 0.1))
#'     )
#'
#'   # Fit two models
#'   fit1 <- bru(
#'     y ~ x,
#'     family = "poisson",
#'     data = input.df,
#'     options = list(control.compute = list(dic = TRUE))
#'   )
#'   fit2 <- bru(
#'     y ~ x + rand(idx, model = "iid"),
#'     family = "poisson",
#'     data = input.df,
#'     options = list(control.compute = list(dic = TRUE))
#'   )
#'
#'   # Compare DIC
#'
#'   deltaIC(fit1, fit2)
#' }
#' }
deltaIC <- function(..., criterion = "DIC") {
  criterion <- match.arg(criterion, c("DIC", "WAIC"), several.ok = TRUE)

  if ("WAIC" %in% criterion) {
    warning(
      paste0(
        "WAIC values from INLA are not well-defined for point process models.",
        " Use with caution, or better, not at all."
      ),
      .immediately = TRUE
    )
  }

  names <- as.character(substitute(list(...)))[-1L]
  model <- eval(list(...))
  nmod <- length(model)
  dic <- waic <- rep(NA, nmod)
  for (i in 1:nmod) {
    mod <- model[[i]]
    if ("DIC" %in% criterion) {
      if (is.null(mod$dic$dic)) {
        stop("Object ", i, " does not have DIC information.")
      }
      dic[i] <- mod$dic$dic
    }
    if ("WAIC" %in% criterion) {
      if (is.null(mod$waic$waic)) {
        stop("Object ", i, " does not have WAIC information.")
      }
      waic[i] <- mod$waic$waic
    }
  }
  if ("DIC" %in% criterion) {
    ord <- order(dic)
  } else if ("WAIC" %in% criterion) {
    ord <- order(waic)
  }
  dic <- dic[ord]
  waic <- waic[ord]
  names <- names[ord]
  ddic <- dic - min(dic)
  dwaic <- waic - min(waic)
  result <- data.frame(Model = names)
  if ("DIC" %in% criterion) {
    result <- cbind(result, data.frame(DIC = dic, Delta.DIC = ddic))
  }
  if ("WAIC" %in% criterion) {
    result <- cbind(result, data.frame(WAIC = waic, Delta.WAIC = dwaic))
  }
  result
}


#' @title Information criteria for point processes
#' @description `r lifecycle::badge("experimental")`
#'   Calculate information criteria for point process models, using a fitted
#'   model and optionally new data. The interface is still under development
#'   and may change without warning in future versions.
#' @param fit A fitted `bru` model object, from [bru()] or [lgcp()].
#' @param data A data frame containing the observed points, with columns for
#'  coordinates and any covariates used in the model. The data frame may also
#'  include a column named `.block` that indicates the block membership of each
#'  point, if blockwise criteria are to be calculated.
#' @param predictor An expression that evaluates to the linear predictor of the
#'  model, given the data. This should be an expression that can be evaluated
#'  in the context of the data frame and model components, just as in the
#'  `formula` argument of [bru_obs()] and [lgcp()].
#' @param domain,samplers Arguments defining the integration domain(s) and
#'   sampling regions, as for [bru_obs()] and [lgcp()].
#' @param n.samples Integer. The number of posterior samples to draw for
#'   estimating the criteria. Default is 2000.
#' @returns A `tibble` with columns for the `Criterion` name, `Method`,
#' log pointwise predictive density (`lppd`), effective number of parameters
#' (`p_eff`), and the information criterion value (`IC`)
#' @export
#' @keywords internal
# lgcp_IC(fit_veg,nests,vegetation,list(geometry=mesh_int),boundary)
# lgcp_IC(fit_veg,nests,vegetation,list(geometry=mesh_int),cvpart)
lgcp_IC <- function(fit, data, predictor, domain, samplers, n.samples = 2000L) {
  int_points <- fm_int(domain = domain, samplers = samplers)
  agg <- bm_aggregate(type = "sum", n_block = 1L)
  agg_block <- bm_aggregate(type = "sum", n_block = max(int_points$.block))
  if (is.null(data[[".block"]]) || (max(int_points$.block) == 1L)) {
    # TODO: implement automated linking of points to integration blocks.
    data$.block <- 1L
  }
  pred_quo <- rlang::enquo(predictor)
  pred_text <- as.character(as.formula(deparse1(pred_quo)))
  pred_text <- pred_text[length(pred_text)]
  form <- as.formula(
    glue::glue(
      '~ {{
      eta <- {{
               {pred_text}
             }}
      log_lambda_y <- eta[part == "Y"]
      lambda_y <- exp(log_lambda_y)
      # Full domain
      lambda_int <- sum(weight[part == "int"] * exp(eta[part == "int"]))
      log_lambda_int <- log(lambda_int)
      area <- sum(weight[part == "int"])
      log_like <- sum(log_lambda_y) + area - lambda_int
      # Blockwise
      lambda_block_int <-
        ibm_eval(agg_block,
                 input = list(
                   block = .block[part == "int"],
                   weights = weight[part == "int"]
                 ),
                 state = exp(eta[part == "int"])
                )
      block_area <-
        ibm_eval(agg_block,
                 input = list(
                   block = .block[part == "int"],
                   weights = weight[part == "int"]
                 ),
                 state = rep(1, sum(part == "int"))
                )
      block_log_lambda_y <-
        ibm_eval(agg_block,
                 input = list(block = .block[part == "Y"]),
                 state = log(lambda_block_int / block_area)[.block[part == "Y"]]
                )
      block_n <-
        ibm_eval(agg_block,
                 input = list(block = .block[part == "Y"]),
                 state = rep(1, sum(part == "Y"))
                )

      block_log_like <- block_area - lambda_block_int + block_log_lambda_y
      block_like <-
        exp(block_log_like - block_area +
              block_n * log(block_area) - lgamma(block_n + 1L))
      list(
        lambda_y = lambda_y,
        log_lambda_y = log_lambda_y,
        lambda_int = lambda_int,
        log_lambda_int = log_lambda_int,
        area = area,
        log_like = log_like,
        block_log_like = block_log_like,
        block_like = block_like
      )
    }}'
    )
  )
  newdata <- cbind(
    dplyr::bind_rows(data, int_points),
    part = c(
      rep("Y", nrow(data)),
      rep("int", nrow(int_points))
    )
  )
  pred <- predict(fit, newdata = newdata, formula = form, n.samples = n.samples)
  # Limit:
  lppd_limit <- sum(log(pred$lambda_y$mean)) +
    pred$area$mean -
    pred$lambda_int$mean
  p_WAIC_limit <- sum(pred$log_lambda_y$sd^2)
  p_DIC <- 2 * pred$log_like$sd^2
  WAIC_limit <- -2 * (lppd_limit - p_WAIC_limit)
  DIC <- -2 * (lppd_limit - p_DIC)

  lppd_limit.se <-
    sqrt(
      sum(pred$lambda_y$mean.mc_std_err^2 / pred$lambda_y$mean^2) +
        pred$lambda_int$mean.mc_std_err^2
    )
  p_DIC.se <- sqrt(4 * (2 * pred$log_like$sd)^2 * pred$log_like$sd.mc_std_err^2)
  p_WAIC_limit.se <- sqrt(sum(
    (2 * pred$log_lambda_y$sd)^2 *
      pred$log_lambda_y$sd.mc_std_err^2
  ))

  # True Blockwise counts:
  block_area <- ibm_eval(
    agg_block,
    input = list(
      block = newdata$.block[newdata$part == "int"],
      weights = newdata$weight[newdata$part == "int"]
    ),
    state = rep(1, sum(newdata$part == "int"))
  )
  block_n <-
    ibm_eval(
      agg_block,
      input = list(block = newdata$.block[newdata$part == "Y"]),
      state = rep(1, sum(newdata$part == "Y"))
    )
  adj <- newdata$.block[newdata$part == "Y"]

  ok <- block_area > 0
  lppd_blockwise <- sum((log(pred$block_like$mean) +
    block_area -
    block_n * log(block_area) +
    lgamma(block_n + 1L))[ok])
  p_WAIC_blockwise <- sum(pred$block_log_like$sd[ok]^2)
  WAIC_blockwise <- -2 * (lppd_blockwise - p_WAIC_blockwise)

  lppd_blockwise.se <-
    sqrt(
      sum(pred$block_like$mean.mc_std_err[ok]^2 / pred$block_like$mean[ok]^2)
    )
  p_WAIC_blockwise.se <- sqrt(sum(
    (2 * pred$block_log_like$sd[ok])^2 *
      pred$block_log_like$sd.mc_std_err[ok]^2
  ))

  # INLA ICs:
  lppd_INLA_DIC <- fit$dic$dic / (-2) + fit$dic$p.eff
  lppd_INLA_WAIC <- fit$waic$waic / (-2) + fit$waic$p.eff
  lppd_INLA_DIC_adjusted <- lppd_INLA_DIC + sum(int_points$weight)
  lppd_INLA_WAIC_adjusted <- lppd_INLA_WAIC + sum(int_points$weight)

  # fmt: skip
  df <- tibble::tribble(
    ~Criterion, ~Method, ~lppd, ~p_eff, ~IC,
    ~lppd.std_err, ~p_eff.std_err, ~IC.std_err,
    "DIC", "Limit", lppd_limit, p_DIC, DIC, lppd_limit.se, p_DIC.se, NA,
    "DIC", "INLA", lppd_INLA_DIC, fit$dic$p.eff, fit$dic$dic, NA, NA, NA,
    "DIC_adjusted", "INLA", lppd_INLA_DIC_adjusted, fit$dic$p.eff,
    -2 * (lppd_INLA_DIC_adjusted - fit$dic$p.eff), NA, NA, NA,
    "WAIC", "Limit", lppd_limit, p_WAIC_limit, WAIC_limit,
    lppd_limit.se, p_WAIC_limit.se, NA,
    "WAIC", "Blockwise", lppd_blockwise, p_WAIC_blockwise, WAIC_blockwise,
    lppd_blockwise.se, p_WAIC_blockwise.se, NA,
    "WAIC", "INLA", lppd_INLA_WAIC, fit$waic$p.eff, fit$waic$waic, NA, NA, NA,
    "WAIC_adjusted", "INLA", lppd_INLA_WAIC_adjusted, fit$waic$p.eff,
    -2 * (lppd_INLA_WAIC_adjusted - fit$waic$p.eff), NA, NA, NA
  )
  df$IC.std_err <- 2 * (df$lppd.std_err + df$p_eff.std_err)
  df
}

#' @export
#' @describeIn lgcp_IC Calculate IC differences based on joined output from
#' several `lgcp_IC` calls;
#' ```
#' delta_IC(rbind(cbind(lgcp_IC(...), Model = "A"),
#'                cbind(lgcp_IC(...), Model = "B")))
#' ```
delta_IC <- function(data) {
  data <- data |>
    dplyr::group_by(.data$Criterion, .data$Method) |>
    dplyr::mutate(Delta_IC = .data$IC - min(.data$IC)) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$Criterion, .data$Method, .data$IC)
  data
}
