library(INLA)
library(inlabru)
library(tibble)
library(ggplot2)

# Linear model, but treated as non-linear, forcing iterations; should be approx
# equivalent to repeated inla.rerun(), which also appears to have the same issue.
# Unclear whether the main cause is premature inner latent field optimisation termination,
# or the hyperparameter optimisation.

set.seed(1234L)
N <- 1000
df <- tibble(
  time = seq_len(N),
  x = exp(time / N * 4) - (time / N)^2 * 50,
  y = 5 + x + rnorm(N, sd = 0.5)
)

ggplot(df) + geom_point(aes(time, y))

bru_options_set(bru_verbose = 3)
# Default is 0.005. Some improvement usually with tolerance = 0.00001
bru_options_set(control.inla = list(tolerance = 0.005))
# Sometimes problem both in parallel and single-threaded, but single usually more stable:
# bru_options_set(num.threads = "1:1")
bru_options_set(inla.mode = "compact")

# The latent field (and sometimes the hyperparameters) "wobble around" between
# iterations, instead of converging to something.
# See "Max deviation from previous" in the verbose output from bru(),
# showing maximum changes of more than 1% of the estimated posterior deviation,
# sometimes as large as 12%.
# The relative change is also shown in the lower left panel of the convergence plot.

fit <- bru(
  components = ~
    Intercept(1) +
    X(time, model = "rw2", constr = TRUE, scale.model = TRUE),
  like(
    formula = y ~ Intercept + X,
    family = "gaussian",
    data = df
  )
)

summary(fit)

bru_convergence_plot(fit)
























set.seed(1234L)
N <- 1000
df <- tibble(
  time = seq_len(N),
  x = exp(time / N * 4) - (time / N)^2 * 50,
  y = 5 + x + rnorm(N, sd = 5)
)

ggplot(df) + geom_point(aes(time, y))

bru_options_set(bru_verbose = 3)
# Default is 0.005. Some improvement usually with tolerance = 0.00001
bru_options_set(control.inla = list(tolerance = 0.005))
# Sometimes problem both in parallel and single-threaded, but single usually more stable:
# bru_options_set(num.threads = "1:1")
bru_options_set(inla.mode = "compact")

# The latent field (and sometimes the hyperparameters) "wobble around" between
# iterations, instead of converging to something.
# See "Max deviation from previous" in the verbose output from bru(),
# showing maximum changes of more than 1% of the estimated posterior deviation,
# sometimes as large as 12%.
# The relative change is also shown in the lower left panel of the convergence plot.

fit2 <- bru(
  components = ~
    Intercept(1) +
    X(time, model = "rw2", constr = TRUE, scale.model = TRUE),
  like(
    formula = y ~ exp(Intercept) + X,
    family = "gaussian",
    data = df
  )
)

summary(fit2)

bru_convergence_plot(fit2)


bru_options_set(control.inla = list(tolerance = 0.0001))

fit3c <- bru(
  components = ~
    Intercept(1) +
    X(time, model = "rw2", constr = TRUE, scale.model = TRUE,
      hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))),
  like(
    formula = y ~ Intercept + exp(X),
    family = "gaussian",
    data = df,
    control.family = list(
      hyper = list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))
    )
  ),
  options = list(bru_max_iter = 40,
                 control.inla = list(optimise.strategy="plain"))
)

summary(fit3)

bru_convergence_plot(fit3b)

