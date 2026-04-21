#' @include deprecated.R
#' @include mappers.R

## _expr ####

#' @title Mapper for general expressions
#' @description
#' Constructs an expression mapping
#' @export
#' @param expr An expression object, typically created using `rlang::quo()`,
#'   that can be evaluated in a data mask containing the root variables, derived
#'   variables, and data variables.
#' @param root_label,derived_label Character strings specifying the pronoun
#'   labels for the data mask for root and derived variables. Default "root" and
#'   "derived", respectively. For example, if `root_label = "latent"`, then the
#'   root variables will be accessible in the expression as `.latent$varname`.
#' @param root_suffix If non-NULL, aharacter string specifying a suffix to add
#'   to the root variable names when making them directly available to the
#'   expression. Default is `""`. For example, if `root_suffix = "_latent"`,
#'   then a root variable named `"x"` will be available as `"x_latent"` in the
#'   expression, in addition to the data mask pronoun version.
#' @inheritParams bru_mapper_generics
#' @details
#' The `input` should be a list with elements
#' \describe{
#' \item{data}{Defaults to `NULL`. If not `NULL`, should be a list of
#' data.frame or similar objects, with `data$data` being the main data
#' container.
#' If `is_rowwise == TRUE`, the number of rows in the `data$data` data.frame
#' determines the number of rows in the output, and the columns can be used as
#' constants in the expression, accessed via `.data$colname`,
#' `.data.$colname`, or `.data.[["colname"]]`.}
#' \item{roots}{If `NULL` or missing, defaults to `names(state)`, and must
#' otherwise be a subset of `names(state)`.
#' If not `NULL`, should be a character vector of variable names.
#' By definition, the pre-Jacobian for each root variable is an identity
#' matrix.}
#' \item{derived}{The state vectors of variables derived from the root
#'   variables. If `NULL` or missing, defaults to an empty list.}
#' \item{jacobians}{If `NULL` or missing, defaults to an empty list.
#' If not `NULL`, should be a list with named entries, one for variable derived
#' from the root variables. Each list element should be a named list of Jacobian
#' matrices, with names matching the root variables.
#' Missing entries are treated as all-zero matrices.
#' }
#' }
#'
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' # Basic expression with only root variables ("x"). Implies
#' # input$roots = names(state) = "x" and identity Jacobian for "x".
#' (m <- bm_expr(rlang::quo(cos(x))))
#' ibm_eval(m, list(), list(x = 1:5))
#' ibm_eval2(m, list(), list(x = 1:5))
#'
#' # Expression with data
#' (m <- bm_expr(rlang::quo(cos(x) * .data$z)))
#' ibm_eval(m, list(data = list(data = data.frame(z = 11:15))), list(x = 1:5))
#' ibm_eval2(m, list(data = data.frame(z = 11:15)), list(x = 1:5))
#'
#' # Expression with data, root variables, and derived variables.
#' (m <- bm_expr(
#'   rlang::quo(sin(x_latent) + cos(.effect$y) * .data$z),
#'   root_label = "latent",
#'   derived_label = "effect",
#'   root_suffix = "_latent"
#' )
#' )
#' ibm_eval(
#'   m,
#'   list(
#'     data = list(data = data.frame(z = 11:15)),
#'     derived = list(y = 2:6), # y = x + 1
#'     jacobians = list(y = list(x = Matrix::Diagonal(1.0, 5)))
#'   ),
#'   state = list(x = 1:5)
#' )
#'
#' @rdname bm_expr
#'
bm_expr <- function(
  expr,
  ...,
  root_label = "root",
  derived_label = "derived",
  root_suffix = NULL,
  .envir = parent.frame()
) {
  mapper <- list(
    expr = rlang::as_quosure(expr, env = .envir),
    is_linear = FALSE,
    is_additive = FALSE,
    is_rowwise = FALSE,
    root_label = root_label,
    derived_label = derived_label,
    root_suffix = root_suffix
  )
  bru_mapper_define(mapper, new_class = "bm_expr")
}

#' @export
#' @rdname ibm_n
#'
ibm_n.bm_expr <- function(
  mapper,
  ...,
  input = NULL,
  state = NULL,
  multi = FALSE
) {
  if (is.null(state)) {
    return(NA_integer_)
  }

  n <- lengths(state)
  if (!multi) {
    n <- sum(n)
  }

  n
}


#' @export
#' @rdname ibm_n_output
#'
ibm_n_output.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  ...
) {
  if (
    !mapper[["is_rowwise"]] ||
      is.null(input[["data"]]) ||
      is.null(input[["data"]][["data"]])
  ) {
    return(NA_integer_)
  }
  NROW(input[["data"]][["data"]][[1]])
}

#' @export
#' @rdname ibm_values
#'
ibm_values.bm_expr <- function(mapper, inla_f = FALSE, ...) {
  NULL
}

#' @export
#' @rdname ibm_is_linear
#'
ibm_is_linear.bm_expr <- function(mapper, ...) {
  mapper[["is_linear"]]
}


#' @describeIn ibm_jacobian
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` format should match the description given for [bm_expr()].
#' @export
#'
ibm_jacobian.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  env = rlang::caller_env()
) {



  if (multi) {
    return(A)
  }

  # Combine the matrices (A1, A2, A3, ...) -> rbind(A1, A2, A3, ...)
  A <- do.call(rbind, A)
  A
}

bm_expr_data_mask <- function(
  mapper,
  input,
  state = NULL,
  env = rlang::caller_env()
) {
  if (
    is.null(mapper[["root_suffix"]]) ||
      (is.character(mapper[["root_suffix"]]) &&
        identical(mapper[["root_suffix"]], ""))
  ) {
    state_with_suffix <- NULL
  } else {
    state_names_with_suffix <-
      expand_labels(names(state), names(state), mapper[["root_suffix"]])
    state_with_suffix <- stats::setNames(state, state_names_with_suffix)
  }
  data_mask <- bru_data_mask(
    stats::setNames(
      c(
        list(
          derived = input[["derived"]],
          state_with_suffix,
          root = state
        ),
        input[["data"]]
      ),
      c(
        mapper[["derived_label"]],
        "",
        mapper[["root_label"]],
        names(input[["data"]])
      )
    )
  )
}

#' @export
#' @describeIn ibm_eval
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` format should match the description given for [bm_expr()].
#'
ibm_eval.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  ...,
  data_mask = NULL,
  env = rlang::caller_env()
) {
  if (is.null(data_mask)) {
    data_mask <- bm_expr_data_mask(mapper, input, state = state, env = env)
  }
  val <- rlang::eval_tidy(
    mapper[["expr"]],
    data = data_mask,
    env = env
  )
  val
}

#' @export
#' @rdname ibm_linear
#'
ibm_linear.bm_expr <- function(mapper, input, state, inla_f = FALSE, ...) {
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    inla_f = FALSE,
    multi = TRUE,
    ...,
  )
  bm_taylor(
    offset = eval2$offset,
    jacobian = eval2$jacobian,
    state0 = state,
    values_mapper = mapper
  )
}
