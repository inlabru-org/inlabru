library(INLA)
library(inlabru)
library(tibble)

dat1 <- tibble(
  time_to_event = 1:4,
  event_type = rep(c(1, 0), 2),
  cure_covar = matrix(1:8, 4, 2)
)
dat2 <- tibble(
  time_to_event = 1:3 + 10,
  event_type = c(1, 0, 1),
  cure_covar = matrix(1:6, 3, 2) + 10
)

dat1_ <- as_tibble(unclass(with(
  dat1,
  inla.surv(time = time_to_event, event = event_type, cure = cure_covar)
))[-7])
dat2_ <- as_tibble(unclass(with(
  dat2,
  inla.mdata(time_to_event, cure_covar)
)))

stk1 <- inla.stack(dat1_, A = list(diag(4)), effects=list(Intercept=1))
stk2 <- inla.stack(dat2_, A = list(diag(3)), effects=list(Intercept=1))
inla.stack.LHS(stk1)
inla.stack.RHS(stk1)
inla.stack.LHS(stk2)
inla.stack.RHS(stk2)

stk <- inla.stack(stk1, stk2)
inla.stack.LHS(stk)
inla.stack.RHS(stk)


fit <- inla(formula = list(inla.stack.LHS(stk), inla.stack.LHS(stk)) ~ 0 + Intercept,
            family = c("exponentialsurv", "agaussian"),
            data = inla.stack.RHS(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            safe = FALSE)

fit <- inla(formula = list(inla.stack.LHS(stk1), inla.stack.LHS(stk2)) ~ 0 + Intercept,
            family = c("exponentialsurv", "agaussian"),
            data = inla.stack.RHS(stk),
            control.predictor = list(A = inla.stack.A(stk)),
            safe = FALSE)


fitb1 <- bru(
  components = ~ 0 + Intercept(1),
  like(formula = dat1_ ~ .,
       family = "exponentialsurv",
       response_data = list(dat1_ = dat1_)),
  options = list(safe = FALSE)
)
fitb2 <- bru(
  components = ~ 0 + Intercept(1),
  like(formula = dat2_ ~ .,
       family = "agaussian",
       data = dat2),
  options = list(safe = FALSE)
)
fitb <- bru(
  components = ~ 0 + Intercept(1),
  like(formula = as.data.frame(unclass(inla.surv(time = time_to_event, event = event_type, cure = cure_covar))[-7]) ~ .,
       family = "exponentialsurv",
       data = dat1),
  like(formula = inla.mdata(time_to_event, cure_covar) ~ .,
       family = "agaussian",
       data = dat2),
  options = list(safe = FALSE)
)
