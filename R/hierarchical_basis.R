#' @include deprecated.R
#' @include effect.R
#' @include mappers.R

# pcmatern_B ####

#' @title Make hierarchical mesh basis functions
#' @description `r lifecycle::badge("experimental")` Construct hierarchical
#'   basis functions. This is highly experimental and may change or be removed
#'   in future versions.
#' @export
#' @keywords internal
#' @rdname pcmatern_B
#' @examples
#' m <- fmesher::fm_subdivide(fmesher::fmexample$mesh, 1))
#' B <- list()
#' for (method in c("distance", "laplace", "graph")) {
#'   B[[method]] <- make_hierarchical_mesh_basis(m, method = method)
#' }
#' \donttest{
#' if (require("ggplot2")) {
#'   print(
#'     ggplot(
#'       data =
#'         data.frame(
#'           idx = seq_len(ncol(B$distance)),
#'           nnz = as.vector(Matrix::colSums(0 != B$distance)),
#'           distance = as.vector(Matrix::colSums(B$distance)),
#'           laplace = as.vector(Matrix::colSums(B$laplace)),
#'           graph = as.vector(Matrix::colSums(B$graph))
#'         )
#'     ) +
#'       geom_point(aes(x = idx, y = nnz, color = "nnz")) +
#'       geom_point(aes(x = idx, y = distance, color = "distance")) +
#'       geom_point(aes(x = idx, y = laplace, color = "laplace")) +
#'       geom_point(aes(x = idx, y = graph, color = "graph")) +
#'       scale_y_log10() +
#'       scale_x_log10()
#'   )
#'
#'   idx <- seq_len(200)
#'   theta <- qr.solve(B$laplace[, idx, drop = FALSE], m$loc[, 1])
#'   print(
#'     ggplot() +
#'       gg(m, col = B$laplace[, idx, drop = FALSE] %*% theta, nx = 60, ny = 60) +
#'       geom_fm(data = m, color = ggplot2::alpha("black", 0.1), alpha = 0) +
#'       scale_fill_distiller(palette = "RdBu")
#'   )
#' }
#' }
#'
make_hierarchical_mesh_basis <- function(mesh, forward = TRUE, method = NULL) {
  # Construct neighbour matrix in a way that doesn't involve the mesh specifics;
  # only the computational neighbourhood structure:
  fem <- fm_fem(mesh, order = 2)
  G1 <- fem$g1
  G2 <- fem$g2
  G <- G1
  # Adjacency matrix:
  Gadj <- (G1 != 0) * 1.0
  Gadj <- Gadj - Matrix::Diagonal(nrow(G), diag(Gadj))

  # First point for each disconnected mesh component
  # Calculate graph distances
  ii <- list()
  jj <- list()
  xx <- list()
  D <- rep(Inf, nrow(G))
  while (!all(is.finite(D))) {
    set <- rep(FALSE, nrow(G))
    front <- rep(FALSE, nrow(G))
    start <- min(which(!is.finite(D)))
    front[start] <- TRUE
    max_dist <- -1
    while (any(front)) {
      max_dist <- max_dist + 1
      D[front] <- max_dist
      set <- set | front
      front <- (as.vector(Gadj %*% front) > 0.5) & !set
    }
    set[set] <- (D[set] < max_dist)
    ijx <- make_basis_fcn(
      start = start,
      inner = set,
      D = D,
      max_dist = max_dist,
      G = G,
      Gadj = Gadj,
      method = method
    )
    ii[[length(ii) + 1]] <- ijx[["ii"]]
    jj[[length(jj) + 1]] <- length(jj) + ijx[["jj"]]
    xx[[length(xx) + 1]] <- ijx[["xx"]]
  }

  # Iteratively add basis functions for the point furthest away from the the
  # previous core points, i.e. where D is maximal.  The radius of each is equal
  # to the initial D-value for the new point.
  while (any(D > 0)) {
    D_local <- rep(Inf, nrow(G))
    set <- rep(FALSE, nrow(G))
    front <- rep(FALSE, nrow(G))
    start <- which.max(D) # The first maximal distance point
    front[start] <- TRUE
    max_dist <- D[start]
    for (the_dist in c(0, seq_len(max_dist))) {
      D[front] <- pmin(D[front], the_dist)
      D_local[front] <- the_dist
      set <- set | front
      front <- (as.vector(Gadj %*% front) > 0.5) & !set
    }
    set[set] <- (D_local[set] < max_dist)
    ijx <- make_basis_fcn(
      start = start,
      inner = set,
      D = D_local,
      max_dist = max_dist,
      G = G,
      Gadj = Gadj,
      method = method
    )
    ii[[length(ii) + 1]] <- ijx[["ii"]]
    jj[[length(jj) + 1]] <- length(jj) + ijx[["jj"]]
    xx[[length(xx) + 1]] <- ijx[["xx"]]
  }

  if (forward) {
    B <- Matrix::sparseMatrix(
      i = unlist(ii),
      j = unlist(jj),
      x = unlist(xx),
      dims = c(nrow(G), length(ii))
    )
  } else {
    B <- Matrix::sparseMatrix(
      i = unlist(ii),
      j = length(ii) + 1 - unlist(jj),
      x = unlist(xx),
      dims = c(nrow(G), length(ii))
    )
  }
  B
}

make_basis_fcn <- function(start, inner, D, max_dist, G, Gadj, method = NULL) {
  method <- match.arg(method, c("distance", "laplace", "graphlaplace"))
  front <- D == max_dist
  if (is.logical(inner)) {
    inner <- which(inner)
  }
  # Make sure the start and end are not part of inner
  inner <- setdiff(inner, start)
  if (length(inner) == 0L) {
    return(list(ii = start, jj = 1L, xx = 1.0))
  }
  ii <- c(start, inner)
  jj <- rep(1, length(ii))
  xx <- switch(method,
    "distance" = {
      c(1, (max_dist - D[inner]) / max_dist)
    },
    "laplace" = {
      b <- G[inner, start, drop = FALSE]
      L <- G[inner, inner, drop = FALSE]
      c(1, -as.vector(Matrix::solve(L, b)))
    },
    "graphlaplace" = {
      Glaplace <- Matrix::Diagonal(nrow(G), Matrix::rowSums(Gadj)) - Gadj
      b <- Glaplace[inner, start, drop = FALSE]
      L <- Glaplace[inner, inner, drop = FALSE]
      c(1, -as.vector(Matrix::solve(L, b)))
    }
  )
  list(ii = ii, jj = jj, xx = xx)
}


#' @export
#' @keywords internal
#' @describeIn pcmatern_B Construct a pcmatern model with basis change
#' @seealso [bru_get_mapper.inla_model_reparam()]
#' `r lifecycle::badge("experimental")`
inla.spde2.pcmatern_B <- function(mesh, ..., B) {
  model <- INLA::inla.spde2.pcmatern(mesh, ...)
  model$n.spde <- ncol(B)
  model$f$n <- ncol(B)
  if (nrow(B) != ncol(B)) {
    stop("Rectangular B not supported")
  }
  # TODO: check that it's a stationary model, since non-stationary would need a
  # different precision structure (should use rgeneric or cgeneric) and
  # different B0, B1, B2 matrices
  model$param.inla$M0 <- Matrix::t(B) %*% model$param.inla$M0 %*% B
  model$param.inla$M1 <- Matrix::t(B) %*% model$param.inla$M1 %*% B
  model$param.inla$M2 <- Matrix::t(B) %*% model$param.inla$M2 %*% B
  model$B <- B
  class(model) <- c("inla_model_reparam", class(model))
  model
}

#' @describeIn bru_get_mapper Reparameterised inla model mapper, see
#'  [inla.spde2.pcmatern_B()]
#' `r lifecycle::badge("experimental")`
#' @export
bru_get_mapper.inla_model_reparam <- function(model, ...) {
  map <- NextMethod()
  bm_reparam(map, model[["B"]])
}
