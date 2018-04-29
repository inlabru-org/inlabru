\donttest{
if (require("INLA", quietly = TRUE)) {
  
# The like function's main purpose is to set up models with multiple likelihoods. 
# The following example generates some random covariates which are observed through
# two different random effect models with different likelihoods

# Generate the data 

set.seed(123)

n1 = 200 
n2 = 10

x1 = runif(n1)
x2 = runif(n2)
z2 = runif(n2)

y1 = rnorm(n1, mean = 2 * x1 + 3)
y2 = rpois(n2, lambda = exp(2 * x2 + z2 + 3))

df1 = data.frame(y = y1, x = x1)
df2 = data.frame(y = y2, x = x2, z = z2)

# Single likelihood models and inference using bru are done via

cmp1 = y ~ x + Intercept
fit1 = bru(cmp1, family = "gaussian", data = df1)
summary(fit1)

cmp2 = y ~ x + z + Intercept
fit2 = bru(cmp2, family = "poisson", data = df2)
summary(fit2)

# A joint model has two likelihoods, which are set up using the like function

lik1 = like("gaussian", formula = y ~ x + Intercept, data = df1)
lik2 = like("poisson", formula = y ~ x + z + Intercept, data = df2)

# The union of effects of both models gives the components needed to run bru

jcmp = y ~ x + z + Intercept
jfit = bru(jcmp, lik1, lik2)

# Compare the estimates

p1 = ggplot() + gg(fit1$summary.fixed, bar = TRUE) + ylim(0, 4) + ggtitle("Model 1")
p2 = ggplot() + gg(fit2$summary.fixed, bar = TRUE) + ylim(0, 4) + ggtitle("Model 2")
pj = ggplot() + gg(jfit$summary.fixed, bar = TRUE) + ylim(0, 4) + ggtitle("Joint model")

multiplot(p1, p2, pj)

}
}
