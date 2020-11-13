# Create a stack from a model, data and prediction values
#
# # TODO: make sure data is allowed to be a list with unequal length vectors

make.stack <- function(data,
                       model,
                       y,
                       E = 1,
                       Ntrials = 1,
                       offset = 0,
                       expr = NULL,
                       result = NULL,
                       tag = "BRU.stack",
                       include = NULL,
                       exclude = NULL) {
  included <- parse_inclusion(names(model$effects),
    include = include,
    exclude = exclude
  )
  offsets <- names(model$effects)[vapply(
    model$effects,
    function(x) identical(x$type, "offset"),
    TRUE
  )]
  included <- setdiff(included, offsets)

  # Observations y
  if (length(y) == 1) {
    y <- rep(y, nrow(as.data.frame(data)))
  }

  # Expectation parameter E
  if (length(E) == 1) {
    E <- rep(E, nrow(as.data.frame(data)))
  }

  # Ntrials
  if (length(Ntrials) == 1) {
    Ntrials <- rep(Ntrials, nrow(as.data.frame(data)))
  }


  # Projection matrices (A) and mesh index effects
  A <- amatrix_eval(model$effects[included], data)
  effects <- index_eval(model$effects[included])


  # Taylor approximation
  if (!is.null(expr)) {
    rw <- nlinla.reweight(A, model, data, expr, result)
    A <- rw$A
    taylor.offset <- rw$const
  } else {
    taylor.offset <- 0
  }

  #  # The weirdest workaround ever. Without this, there are convergence problems on ubuntu but not on MacOS ?!?!?!
  #  A <- c(A, list(1))
  #  effects <- c(effects, list(WORKAROUND = runif(dim(A[[1]])[1])))

  # Create and return stack
  stk <- INLA::inla.stack(
    data = list(
      BRU.response = y,
      BRU.E = E,
      BRU.Ntrials = Ntrials,
      BRU.offset = taylor.offset + offset
    ),
    A = A,
    tag = tag,
    effects = effects,
    # Make sure latent components with zero-derivatives aren't removed:
    remove.unused = FALSE
  )
  stk
}
