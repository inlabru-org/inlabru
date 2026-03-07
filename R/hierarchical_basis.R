#' @include deprecated.R
#' @include effect.R
#' @include mappers.R

# pcmatern_B ####

#' @title Make hierarchical mesh basis functions
#' @description `r lifecycle::badge("experimental")` Construct hierarchical
#'   basis functions. This is highly experimental and may change or be removed
#'   in future versions.
#' @param mesh A `fmesher` object.
#' @param forward logical; if `TRUE` (default), the basis functions are
#'   constructed in forward order; otherwise in reverse order.
#' @param method Character; one of `"laplace"` (default), `"graphdistance"` or
#'   `"graphlaplace"`. The method used to construct the basis functions. See
#'   Details.
#' @param alpha Numeric; only used if `method = "laplace"`. The approximate
#'   Laplacian power used to construct the basis functions. Default is
#'   `(1 + fm_manifold_dim(mesh)) / 2`, which gives an approximately linear
#'   function along the shortest path. Must be in the range `[1, 2]`. See
#'   Details.
#' @details The hierarchical basis functions are constructed by iteratively
#'   adding basis functions centered at the point furthest away from the
#'   previously added basis functions. The radius of each basis function is
#'   equal to the distance to the previously added basis functions at the time
#'   of addition. Three methods are available to construct the basis functions:
#'   \describe{
#'     \item{`laplace`}{The basis functions are constructed using a combination
#'       of
#'  two operators, the Laplacian and the squared Laplacian (the biharmonic
#'  operator). This approximates the Laplacian raised to the power `alpha`. The
#'  default value for `alpha` gives an approximately linear function along the
#'  shortest path.one based on the standard FEM stiffness matrix, and one based
#'  on the graph Laplacian. The `alpha` parameter controls the weight between
#'  the two operators, with `alpha = 1` using only the Laplacian and `alpha = 2`
#'  using only the squared Laplacian.}
#'    \item{`graphdistance`}{The basis functions are constructed using a simple
#'  neighbour distance-based approach, where the basis function value
#'  decreases by equal amounts for each traversed edge.}
#'    \item{`graphlaplace`}{The basis functions are constructed using the
#'  graph Laplacian only.}
#'  }
#' @returns A sparse matrix of class `dgCMatrix`, where each column
#'   corresponds to a basis function.
#' @export
#' @keywords internal
#' @rdname pcmatern_B
#' @examples
#' m <- fmesher::fmexample$mesh
#' B <- make_hierarchical_mesh_basis(m, method = "laplace2", alpha = 1.5)
#'
#' \donttest{
#' m <- fmesher::fm_subdivide(fmesher::fmexample$mesh, 1)
#' B <- list()
#' for (method in c(
#'   "laplace", "laplace2",
#'   "graphdistance",
#'   "graphlaplace", "graphlaplace2"
#' )) {
#'   B[[method]] <- make_hierarchical_mesh_basis(
#'     m,
#'     method = method, alpha = 1.5
#'   )
#' }
#' if (require("ggplot2") && require("patchwork")) {
#'   print(
#'     ggplot(
#'       data =
#'         data.frame(
#'           idx = seq_len(ncol(B$laplace)),
#'           nnz = as.vector(Matrix::colSums(0 != B$laplace)),
#'           laplace = as.vector(Matrix::colSums(B$laplace)),
#'           laplace2 = as.vector(Matrix::colSums(B$laplace2)),
#'           distance = as.vector(Matrix::colSums(B$graphdistance)),
#'           graph = as.vector(Matrix::colSums(B$graphlaplace)),
#'           graph2 = as.vector(Matrix::colSums(B$graphlaplace2))
#'         )
#'     ) +
#'       geom_point(aes(x = idx, y = nnz, color = "nnz")) +
#'       geom_point(aes(x = idx, y = laplace, color = "laplace")) +
#'       geom_point(aes(x = idx, y = laplace2, color = "laplace2")) +
#'       geom_point(aes(x = idx, y = distance, color = "graphdistance")) +
#'       geom_point(aes(x = idx, y = graph, color = "graphlaplace")) +
#'       scale_y_log10() +
#'       scale_x_log10()
#'   )
#'
#'   idx <- seq_len(fm_dof(m))
#'   idx <- seq_len(10)
#'   method <- "laplace"
#'   theta <- qr.solve(B[[method]][, idx, drop = FALSE], m$loc[, 1])
#'   print(
#'     ggplot() +
#'       gg(m,
#'         col = B[[method]][, idx, drop = FALSE] %*% theta - m$loc[, 1],
#'         nx = 60, ny = 60
#'       ) +
#'       geom_fm(data = m, color = ggplot2::alpha("black", 0.1), alpha = 0) +
#'       scale_fill_distiller(palette = "RdBu")
#'   )
#'
#'   ev <- fm_evaluator(m, dims = c(60, 60))
#'   df <- NULL
#'   for (method in c(
#'     "laplace", "laplace2",
#'     "graphdistance",
#'     "graphlaplace", "graphlaplace2"
#'   )) {
#'     fun_cumulative <- 1e-16
#'     for (k in 1:28) {
#'       fun_orig <- as.vector(B[[method]][, k, drop = FALSE])
#'       fun_cumulative <- fun_cumulative + fun_orig
#'       fun_norm <- fun_orig / fun_cumulative
#'       df <- rbind(
#'         df,
#'         data.frame(
#'           x = ev$lattice$loc[, 1],
#'           y = ev$lattice$loc[, 2],
#'           orig = as.vector(fm_evaluate(ev, fun_orig)),
#'           norm = as.vector(fm_evaluate(ev, fun_norm)),
#'           index = k,
#'           fun = paste0(
#'             "fun",
#'             formatC(k, width = 2, format = "d", flag = "0")
#'           ),
#'           method = method
#'         )
#'       )
#'     }
#'   }
#'   df$norm[df$fun == "fun01"] <- NA
#'   df$orig[df$orig == 0] <- NA
#'   df$norm[df$norm == 0] <- NA
#'   print(
#'     ggplot(data = df[df$index <= 4, ]) +
#'       geom_tile(aes(x, y, fill = orig),
#'         na.rm = TRUE
#'       ) +
#'       geom_fm(data = m, color = ggplot2::alpha("black", 0.1), alpha = 0) +
#'       geom_contour(aes(x, y, z = orig),
#'         breaks = seq_len(29 - 15) / (30 - 15), #* 2-0.5,
#'         na.rm = TRUE, col = "black", alpha = 0.5
#'       ) +
#'       geom_contour(aes(x, y, z = norm),
#'         breaks = seq_len(29 - 15) / (30 - 15), #* 2-0.5,
#'         na.rm = TRUE, col = "red", alpha = 0.5
#'       ) +
#'       scale_fill_distiller(palette = "RdBu", limits = c(0, 1)) + #* 2-0.5) +
#'       facet_wrap(vars(fun, method))
#'   )
#' }
#' }
#'
make_hierarchical_mesh_basis <- function(mesh, forward = TRUE, method = NULL,
                                         alpha = 1) {
  if (is.null(alpha)) {
    alpha <- (1 + fmesher::fm_manifold_dim(mesh)) / 2
  }
  stopifnot((alpha >= 1) && (alpha <= 2))
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
      method = method,
      G2 = G2,
      alpha = alpha
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
      method = method,
      G2 = G2,
      alpha = alpha
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

make_basis_fcn <- function(start, inner, D, max_dist, G, Gadj, method = NULL,
                           G2 = NULL, alpha = 1) {
  method <- match.arg(method, c(
    "graphdistance",
    "laplace", "laplace2",
    "graphlaplace", "graphlaplace2"
  ))
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
    "graphdistance" = {
      c(1, (max_dist - D[inner]) / max_dist)
    },
    "laplace" = {
      val <- 0
      if (alpha < 2) {
        b1 <- G[inner, start, drop = FALSE]
        L1 <- G[inner, inner, drop = FALSE]
        val <- val + (2 - alpha) * c(1, -as.vector(Matrix::solve(L1, b1)))
      }
      if (alpha > 1) {
        b2 <- G2[inner, start, drop = FALSE]
        L2 <- G2[inner, inner, drop = FALSE]
        val <- val + (alpha - 1) * c(1, -as.vector(Matrix::solve(L2, b2)))
      }
      val
    },
    "laplace2" = {
      M <- G * (2 - alpha)
      M <- M + G2 * (alpha - 1)
      b <- M[inner, start, drop = FALSE]
      L <- M[inner, inner, drop = FALSE]
      val <- c(1, -as.vector(Matrix::solve(L, b)))
      val
    },
    "graphlaplace" = {
      G1 <- Matrix::Diagonal(nrow(G), Matrix::rowSums(Gadj)) - Gadj
      G2 <- G1 %*% G1
      val <- 0
      if (alpha < 2) {
        b1 <- G1[inner, start, drop = FALSE]
        L1 <- G1[inner, inner, drop = FALSE]
        val <- val + (2 - alpha) * c(1, -as.vector(Matrix::solve(L1, b1)))
      }
      if (alpha > 1) {
        b2 <- G2[inner, start, drop = FALSE]
        L2 <- G2[inner, inner, drop = FALSE]
        val <- val + (alpha - 1) * c(1, -as.vector(Matrix::solve(L2, b2)))
      }
      val
    },
    "graphlaplace2" = {
      G1 <- Matrix::Diagonal(nrow(G), Matrix::rowSums(Gadj)) - Gadj
      G2 <- G1 %*% G1
      M <- (2 - alpha) * G1 + (alpha - 1) * G2
      b <- M[inner, start, drop = FALSE]
      L <- M[inner, inner, drop = FALSE]
      val <- c(1, -as.vector(Matrix::solve(L, b)))
      val
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
