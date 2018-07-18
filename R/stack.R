# Create a stack from a model, points and prediction alues

make.stack <- function(points,
                       model,
                       y,
                       E = 1,
                       Ntrials = 1,
                       offset = 0,
                       expr = NULL,
                       result = NULL,
                       tag = "BRU.stack"){
  
  
  # Observations y
  if ( length(y) == 1 ) { y = rep(y, nrow(as.data.frame(points))) }
  
  # Expectation parameter E
  if ( length(E) == 1 ) { E = rep(E, nrow(as.data.frame(points))) }
  
  # Ntrials
  if ( length(Ntrials) == 1 ) { Ntrials = rep(Ntrials, nrow(as.data.frame(points))) }
  
  
  # Projection matrices (A) and mesh index effects
  A = lapply(model$effects, amatrix.effect, points)
  effects = lapply(model$effects, index.effect, points)
  
  
  # Taylor approximation
  if ( !is.null(expr) ) {
    rw = nlinla.reweight(A, model, points, expr, result)
    A = rw$A
    taylor.offset = rw$const
  } else {
    taylor.offset = 0
  }

  # The weirdest workaround ever. Without this, there are convergence problems on ubuntu but not on MacOS ?!?!?!
  A = c(A, list(1))
  effects = c(effects, list(WORKAROUND = runif(dim(A[[1]])[1])))

  # Create and return stack
  stk <- INLA::inla.stack(data = list(BRU.response = y,
                                      BRU.E = E,
                                      BRU.Ntrials = Ntrials,
                                      BRU.offset = taylor.offset + offset),
                          A = A,
                          tag = tag,
                          effects = effects)
  stk
}
