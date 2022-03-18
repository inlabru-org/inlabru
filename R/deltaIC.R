#' @title Summarise DIC and WAIC from `lgcp` objects.
#'
#' @description
#' Calculates DIC and/or WAIC differences and produces an ordered summary.
#'
#' @param ... Comma-separated objects inheriting from class `inla` and obtained
#' from a run of `INLA::inla()`, [bru()] or [lgcp()]
#' @param criterion character vector. If it includes 'DIC', computes DIC differences;
#' If it contains 'WAIC', computes WAIC differences. Default: 'DIC'
#'
#' @return A data frame with each row containing the Model name, DIC and Delta.DIC,
#' and/or WAIC and Delta.WAIC.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (bru_safe_inla(multicore = FALSE)) {
#'   # Generate some data
#'   input.df <- data.frame(idx = 1:10, x = cos(1:10))
#'   input.df <- within(
#'     input.df,
#'     y <- rpois(10, 5 + 2 * cos(1:10) + rnorm(10, mean = 0, sd = 0.1))
#'   )
#'
#'   # Fit two models
#'   fit1 <- bru(y ~ x, family = "poisson", data = input.df)
#'   fit2 <- bru(y ~ x + rand(idx, model = "iid"), family = "poisson", data = input.df)
#'
#'   # Compare DIC
#'
#'   deltaIC(fit1, fit2)
#' }
#' }
deltaIC <- function(..., criterion = "DIC") {
  criterion <- match.arg(criterion, c("DIC", "WAIC"), several.ok = TRUE)

  names <- as.character(substitute(list(...)))[-1L]
  model <- eval(list(...))
  nmod <- length(model)
  dic <- waic <- rep(NA, nmod)
  for (i in 1:nmod) {
    mod <- model[[i]]
    if (("DIC" %in% criterion) && is.null(mod$dic$dic)) {
      stop("Object ", i, " does not have DIC information.")
    }
    if (("WAIC" %in% criterion) && is.null(mod$waic$waic)) {
      stop("Object ", i, " does not have WAIC information.")
    }
    dic[i] <- mod$dic$dic
    waic[i] <- mod$waic$waic
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
  return(result)
}
