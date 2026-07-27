# Iterative linearised INLA method

## The INLA method for linear predictors

The INLA method is used to compute fast approximative posterior
distribution for Bayesian generalised additive models. The hierarchical
structure of such a model with latent Gaussian components
\boldsymbol{u}, covariance parameters \boldsymbol{\theta}, and measured
response variables \boldsymbol{y}, can be written as \begin{aligned}
\boldsymbol{\theta} &\sim p(\boldsymbol{\theta}) \\
\boldsymbol{u}\|\boldsymbol{\theta} &\sim
\mathcal{N}\\\left(\boldsymbol{\mu}\_u,
\boldsymbol{Q}(\boldsymbol{\theta})^{-1}\right) \\
\boldsymbol{\eta}(\boldsymbol{u}) &= \boldsymbol{A}\boldsymbol{u} \\
\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta} & \sim
p(\boldsymbol{y}\|\boldsymbol{\eta}(\boldsymbol{u}),\boldsymbol{\theta})
\end{aligned} where typically each linear predictor element,
\eta_i(\boldsymbol{u}), is linked to a location parameter of the
distribution for observation y_i, for each i, via a (non-linear) link
function g^{-1}(\cdot). In the R-INLA implementation, the observations
are assumed to be conditionally independent, given \boldsymbol{\eta} and
\boldsymbol{\theta}.

## Approximate INLA for non-linear predictors

The premise for the inlabru method for non-linear predictors is to build
on the existing implementation, and only add a linearisation step. The
properties of the resulting approximation will depend on the nature of
the non-linearity.

Let \widetilde{\boldsymbol{\eta}}(\boldsymbol{u}) be a non-linear
predictor, i.e. a deterministic function of \boldsymbol{u},
\widetilde{\boldsymbol{\eta}} (\boldsymbol{u}) = \textsf{fcn}
(\boldsymbol{u}), and let \overline{\boldsymbol{\eta}}(\boldsymbol{u})
be the 1st order Taylor approximation at \boldsymbol{u}\_0,
\overline{\boldsymbol{\eta}}(\boldsymbol{u}) =
\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}\_0) +
\boldsymbol{B}(\boldsymbol{u} - \boldsymbol{u}\_0) =
\left\[\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}\_0) -
\boldsymbol{B}\boldsymbol{u}\_0\right\] + \boldsymbol{B}\boldsymbol{u} ,
where \boldsymbol{B} is the derivative matrix for the non-linear
predictor, evaluated at \boldsymbol{u}\_0. Hence, we define
\begin{aligned} \boldsymbol{y} \| \boldsymbol{u}, {\boldsymbol{\theta}}
&\overset{d}{=} \boldsymbol{y} \|
\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}), {\boldsymbol{\theta}} \\
&\sim p (\boldsymbol{y} \|
g^{-1}\[\widetilde{\boldsymbol{\eta}}(\boldsymbol{u})\],
{\boldsymbol{\theta}})\\ \end{aligned} The non-linear observation model
p(\boldsymbol{y}\|g^{-1}\[\widetilde{\boldsymbol{\eta}}(\boldsymbol{u})\],\boldsymbol{\theta})
is approximated by replacing the non-linear predictor with its
linearisation, so that the linearised model is defined by

\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) =
p(\boldsymbol{y}\|\overline{\boldsymbol{\eta}}(\boldsymbol{u}),\boldsymbol{\theta})
=
p(\boldsymbol{y}\|g^{-1}\[\overline{\boldsymbol{\eta}}(\boldsymbol{u})\],\boldsymbol{\theta})
\approx
p(\boldsymbol{y}\|g^{-1}\[\widetilde{\boldsymbol{\eta}}(\boldsymbol{u})\],\boldsymbol{\theta})
=
p(\boldsymbol{y}\|\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}),\boldsymbol{\theta})
= \widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) The
non-linear model posterior is factorised as
\widetilde{p}(\boldsymbol{\theta},\boldsymbol{u}\|\boldsymbol{y}) =
\widetilde{p}(\boldsymbol{\theta}\|\boldsymbol{y})\widetilde{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta}),
and the linear model approximation is factorised as
\overline{p}(\boldsymbol{\theta},\boldsymbol{u}\|\boldsymbol{y}) =
\overline{p}(\boldsymbol{\theta}\|\boldsymbol{y})\overline{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta}).

### Fixed point iteration

The remaining step of the approximation is how to choose the
linearisation point \boldsymbol{u}\_\*. For a given linearisation point
\boldsymbol{v}, INLA will compute the posterior mode for
\boldsymbol{\theta}, \widehat{\boldsymbol{\theta}}\_{\boldsymbol{v}} =
\mathop{\mathrm{arg\\max}}\_{\boldsymbol{\theta}}
\overline{p}\_\boldsymbol{v} ( {\boldsymbol{\theta}} \| \boldsymbol{y}
), and the joint conditional posterior mode for \boldsymbol{u},
\widehat{\boldsymbol{u}}\_{\boldsymbol{v}} =
\mathop{\mathrm{arg\\max}}\_{\boldsymbol{u}}
\overline{p}\_\boldsymbol{v} ( \boldsymbol{u} \| \boldsymbol{y},
\widehat{\boldsymbol{\theta}}\_{\boldsymbol{v}} ) .

Define the Bayesian estimation functional[^1]
f(\overline{p}\_{\boldsymbol{v}}) =
(\widehat{\boldsymbol{\theta}}\_{\boldsymbol{v}},\widehat{\boldsymbol{u}}\_{\boldsymbol{v}})
and let f(p)=(\widehat{\boldsymbol{\theta}},\widehat{\boldsymbol{u}})
denote the corresponding posterior modes for the true posterior
distribution, \begin{aligned} \widehat{{\boldsymbol{\theta}}} &=
\mathop{\mathrm{arg\\max}}\_{\boldsymbol{\theta}} p (
{\boldsymbol{\theta}} \| \boldsymbol{y} ), \\ \widehat{\boldsymbol{u}}
&= \mathop{\mathrm{arg\\max}}\_{\boldsymbol{u}} p (\boldsymbol{u} \|
\boldsymbol{y}, \widehat{{\boldsymbol{\theta}}}). \end{aligned}

The fixed point
(\boldsymbol{\theta}\_\*,\boldsymbol{u}\_\*)=f(\overline{p}\_{\boldsymbol{u}\_\*})
should ideally be close to
(\widehat{\boldsymbol{\theta}},\widehat{\boldsymbol{u}}), i.e. close to
the true marginal/conditional posterior mode. We can achieve this for
the conditional latent mode, so that
\boldsymbol{u}\_\*=\mathop{\mathrm{arg\\max}}\_{\boldsymbol{u}} p
(\boldsymbol{u} \| \boldsymbol{y},
\widehat{\boldsymbol{\theta}}\_{\boldsymbol{u}\_\*}).

We therefore seek the latent vector \boldsymbol{u}\_\* that generates
the fixed point of the functional, so that
(\boldsymbol{\theta}\_\*,\boldsymbol{u}\_\*)=f(\overline{p}\_{\boldsymbol{u}\_\*}).

One key to the fixed point iteration is that the observation model is
linked to \boldsymbol{u} only through the non-linear predictor
\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}), since this leads to a
simplified line search method below.

0.  Let \boldsymbol{u}\_0 be an initial linearisation point for the
    latent variables obtained from the initial INLA call. Iterate the
    following steps for k=0,1,2,...
1.  Compute the predictor linearisation at \boldsymbol{u}\_0.
2.  Compute the linearised INLA posterior
    \overline{p}\_{\boldsymbol{u}\_0}(\boldsymbol{\theta}\|\boldsymbol{y}).
3.  Let
    (\boldsymbol{\theta}\_1,\boldsymbol{u}\_1)=(\widehat{\boldsymbol{\theta}}\_{\boldsymbol{u}\_0},\widehat{\boldsymbol{u}}\_{\boldsymbol{u}\_0})=f(\overline{p}\_{\boldsymbol{u}\_0})
    be the initial candidate for new linearisation point.
4.  Let
    \boldsymbol{v}\_\alpha=(1-\alpha)\boldsymbol{u}\_1+\alpha\boldsymbol{u}\_0,
    and find the value \alpha minimises
    \\\widetilde{\eta}(\boldsymbol{v}\_\alpha)-\overline{\eta}(\boldsymbol{u}\_1)\\.
5.  Set the new linearisation point \boldsymbol{u}\_0 equal to
    \boldsymbol{v}\_\alpha and repeat from step 1, unless the iteration
    has converged to a given tolerance.

A potential improvement of step 4 might be to also take into account the
prior distribution for \boldsymbol{u} as a minimisation penalty, to
avoid moving further than would be indicated by a full likelihood
optimisation.

#### Line search

In step 4, we would ideally want \alpha to be  
\mathop{\mathrm{arg\\max}}\_{\alpha} \left\[\ln
p(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta}\_1)\right\]\_{\boldsymbol{u}=\boldsymbol{v}\_\alpha}.
However, since this requires access to the internal likelihood and prior
density evaluation code, we instead use a simpler alternative. We
consider norms of the form
\\\widetilde{\eta}(\boldsymbol{v}\_\alpha)-\overline{\eta}(\boldsymbol{u}\_1)\\
that only depend on the nonlinear and linearised predictor expressions,
and other known quantities, given \boldsymbol{u}\_0, such as the current
INLA estimate of the component wise predictor variances.

Let \sigma_i^2 = \mathrm{Var}\_{\boldsymbol{u}\sim
\overline{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta}\_1)}(\overline{\boldsymbol{\eta}}\_i(\boldsymbol{u}))
be the current estimate of the posterior variance for each predictor
element i. We then define an inner product on the space of predictor
vectors as \langle \boldsymbol{a},\boldsymbol{b} \rangle_V = \sum_i
\frac{a_i b_i}{\sigma_i^2} . The squared norm for the difference between
the predictor vectors
\widetilde{\boldsymbol{\eta}}(\boldsymbol{v}\_\alpha) and
\overline{\boldsymbol{\eta}}(\boldsymbol{u}\_1),with respect to this
inner product, is defined as \\
\widetilde{\boldsymbol{\eta}}(\boldsymbol{v}\_\alpha) -
\overline{\boldsymbol{\eta}}(\boldsymbol{u}\_1)\\^2_V = \sum_i
\frac{\|\widetilde{\boldsymbol{\eta}}\_i(\boldsymbol{v}\_\alpha)-\overline{\boldsymbol{\eta}}\_i(\boldsymbol{u}\_1)\|^2}{\sigma_i^2}
. Using this norm as the target loss function for the line search avoids
many potentially expensive evaluations of the true posterior conditional
log-density. We evaluate
\widetilde{\boldsymbol{\eta}}\_1=\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}\_1)
and make use of the linearised predictor information. Let
\widetilde{\boldsymbol{\eta}}\_\alpha=\widetilde{\boldsymbol{\eta}}(\boldsymbol{v}\_\alpha)
and
\overline{\boldsymbol{\eta}}\_\alpha=\overline{\boldsymbol{\eta}}(\boldsymbol{v}\_\alpha)=(1-\alpha)\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}\_0)+\alpha\overline{\boldsymbol{\eta}}(\boldsymbol{u}\_1).
In other words, \alpha=0 corresponds to the previous linear predictor,
and \alpha=1 is the current estimate from INLA. An exact line search
would minimise
\\\widetilde{\boldsymbol{\eta}}\_\alpha-\overline{\boldsymbol{\eta}}\_1\\.
Instead, we define a quadratic approximation to the non-linear predictor
as a function of \alpha, \breve{\boldsymbol{\eta}}\_\alpha =
\overline{\boldsymbol{\eta}}\_\alpha + \alpha^2
(\widetilde{\boldsymbol{\eta}}\_1 - \overline{\boldsymbol{\eta}}\_1) and
minimise the quartic polynomial in \alpha, \begin{aligned}
\\\breve{\boldsymbol{\eta}}\_\alpha-\overline{\boldsymbol{\eta}}\_1\\^2
&= \\ (\alpha-1)(\overline{\boldsymbol{\eta}}\_1 -
\overline{\boldsymbol{\eta}}\_0) + \alpha^2
(\widetilde{\boldsymbol{\eta}}\_1 - \overline{\boldsymbol{\eta}}\_1)
\\^2 . \end{aligned} If initial expansion and contraction steps are
carried out, leading to an initial guess of \alpha=\gamma^k, where
\gamma\>1 is a scaling factor (see
[`?bru_options`](https://inlabru-org.github.io/inlabru/reference/bru_options.md),
`bru_method$factor`) and k is the (signed) number of expansions and
contractions, the quadratic expression is replaced by \begin{aligned}
\\\breve{\boldsymbol{\eta}}\_\alpha-\overline{\boldsymbol{\eta}}\_1\\^2
&= \\ (\alpha-1)(\overline{\boldsymbol{\eta}}\_1 -
\overline{\boldsymbol{\eta}}\_0) + \frac{\alpha^2}{\gamma^{2k}}
(\widetilde{\boldsymbol{\eta}}\_{\gamma^k} -
\overline{\boldsymbol{\eta}}\_{\gamma^k}) \\^2 , \end{aligned} which is
minimised on the interval \alpha\in\[\gamma^{k-1},\gamma^{k+1}\].

## Posterior non-linearity checks

Whereas the inlabru optimisation method leads to an estimate where \\
\widetilde{\boldsymbol{\eta}} (\boldsymbol{u}\_\*) -
\overline{\boldsymbol{\eta}}(\boldsymbol{u}\_\*)\\=0 for a specific
\boldsymbol{u}\_\*, the overall posterior approximation accuracy depends
on the degree of nonlinearity in the vicinity of \boldsymbol{u}\_\*.
There are two main options for evaluating this nonlinearity, using
sampling from the approximate posterior distribution. The first option
is \begin{aligned} \sum_i \frac{E\_{\boldsymbol{u}\sim
\overline{p}(\boldsymbol{u}\|\boldsymbol{y})}\left\[
\|\overline{\boldsymbol{\eta}}\_i(\boldsymbol{u})-\widetilde{\boldsymbol{\eta}}\_i(\boldsymbol{u})\|^2\right\]}{\mathrm{Var}\_{\boldsymbol{u}\sim
\overline{p}(\boldsymbol{u}\|\boldsymbol{y})}(\overline{\boldsymbol{\eta}}\_i(\boldsymbol{u}))}
, \end{aligned} which is the posterior expectation of the component-wise
variance-normalised squared deviation between the non-linear and
linearised predictor. Note that the normalising variance includes the
variability induced by the posterior uncertainty for
\boldsymbol{\theta}, whereas the \\\cdot\\\_V norm used for the line
search used only the posterior mode. Another option is
E\_{(\boldsymbol{u},\boldsymbol{\theta})\sim
\overline{p}(\boldsymbol{u},\boldsymbol{\theta}\|\boldsymbol{y})}
\left\[\ln \frac{\overline{p}(\boldsymbol{u}
\|\boldsymbol{y},{\boldsymbol{\theta}})}{\widetilde{p}(\boldsymbol{u}\|\boldsymbol{y},{\boldsymbol{\theta}})}\right\]
= E\_{\boldsymbol{\theta}\sim
\overline{p}(\boldsymbol{\theta}\|\boldsymbol{y})} \left\\
E\_{\boldsymbol{u}\sim
\overline{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta})}
\left\[\ln \frac{\overline{p}(\boldsymbol{u}
\|\boldsymbol{y},{\boldsymbol{\theta}})}{\widetilde{p}(\boldsymbol{u}\|\boldsymbol{y},{\boldsymbol{\theta}})}\right\]
\right\\ which is the Kullback–Leibler divergence for the conditional
posterior densities,
\mathsf{KL}\left(\overline{p}\\\middle\\\\\widetilde{p}\right),
integrated over the approximate posterior distribution for
\boldsymbol{\theta}. Implementing this would require access to the
likelihood and prior distribution details. The next section explores
this in more detail.

#### Accuracy

We wish to assess how accurate the approximation is. Thus, we compare
\widetilde{p}(\boldsymbol{u} \| \boldsymbol{y}, \boldsymbol{\theta} )
and \overline{p}(\boldsymbol{u} \|\boldsymbol{y},\boldsymbol{\theta}).
With Bayes’ theorem, \begin{aligned}
p(\boldsymbol{u}\|\boldsymbol{y},{\boldsymbol{\theta}}) &=
\frac{p(\boldsymbol{u},\boldsymbol{y}\|{\boldsymbol{\theta}})}{p(\boldsymbol{y}\|{\boldsymbol{\theta}})}
\\ &= \frac{p(\boldsymbol{y}\|\boldsymbol{u},{\boldsymbol{\theta}})
p(\boldsymbol{u}\|{\boldsymbol{\theta}})}{p(\boldsymbol{y}\|{\boldsymbol{\theta}})},
\end{aligned} where p(\boldsymbol{u}\|\boldsymbol{\theta}) is a Gaussian
density and p(\boldsymbol{y}\|\boldsymbol{\theta}) is a scaling factor
that doesn’t depend on \boldsymbol{u}. We can therefore focus on the
behaviour of \ln p(\boldsymbol{y}\|\boldsymbol{\theta},\boldsymbol{u})
for the exact and linearised observation models.

Recall that the observation likelihood only depends on \boldsymbol{u}
through \boldsymbol{\eta}. Using a Taylor expansion with respect to
\boldsymbol{\eta} and
\boldsymbol{\eta}^\*=\widetilde{\boldsymbol{\eta}}(\boldsymbol{u}\_\*),
\begin{aligned} \ln
p(\boldsymbol{y}\|\boldsymbol{\eta},\boldsymbol{\theta}) &= \ln p
(\boldsymbol{y}\|{\boldsymbol{\theta}},\boldsymbol{\eta}^\*)) \\
&\qquad + \sum_i \left.\frac{\partial}{\partial\eta_i} \ln p
(\boldsymbol{y} \| {\boldsymbol{\theta}}, \boldsymbol{\eta})
\right\|\_{\boldsymbol{\eta}^\*}\cdot (\eta_i - \eta^\*\_i) \\ &\qquad +
\frac{1}{2}\sum\_{i,j}
\left.\frac{\partial^2}{\partial\eta_i\partial\eta_j} \ln p
(\boldsymbol{y} \| {\boldsymbol{\theta}}, \boldsymbol{\eta})
\right\|\_{\boldsymbol{\eta}^\*}\cdot (\eta_i - \eta^\*\_i) (\eta_j -
\eta^\*\_j) + \mathcal{O}(\\\boldsymbol{\eta}-\boldsymbol{\eta}^\*\\^3),
\end{aligned} Similarly, for each component of
\widetilde{\boldsymbol{\eta}}, \begin{aligned}
\widetilde{\eta}\_i(\boldsymbol{u}) &= \eta^\*\_i
+\left\[\left.\nabla\_{u}\widetilde{\eta}\_i(\boldsymbol{u})\right\|\_{\boldsymbol{u}\_\*}\right\]^\top
(\boldsymbol{u} - \boldsymbol{u}\_\*) \\&\quad
+\frac{1}{2}(\boldsymbol{u} -
\boldsymbol{u}\_\*)^\top\left\[\left.\nabla\_{u}\nabla\_{u}^\top\widetilde{\eta}\_i(\boldsymbol{u})\right\|\_{\boldsymbol{u}\_\*}\right\]
(\boldsymbol{u} - \boldsymbol{u}\_\*) +
\mathcal{O}(\\\boldsymbol{u}-\boldsymbol{u}^\*\\^3) \\&= \eta_i^\* +
b_i(\boldsymbol{u}) + h_i(\boldsymbol{u}) +
\mathcal{O}(\\\boldsymbol{u}-\boldsymbol{u}\_\*\\^3) \\&=
\overline{\eta}\_i(\boldsymbol{u}) + h_i(\boldsymbol{u}) +
\mathcal{O}(\\\boldsymbol{u}-\boldsymbol{u}\_\*\\^3) \end{aligned}

where \nabla_u\nabla_u^\top is the Hessian with respect to
\boldsymbol{u}, b_i are linear in \boldsymbol{u}, and h_i are quadratic
in \boldsymbol{u}. Combining the two expansions and taking the
difference between the full and linearised log-likelihoods, we get
\begin{aligned} \ln
\widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) - \ln
\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) &=
\sum_i \left.\frac{\partial}{\partial\eta_i} \ln p (\boldsymbol{y} \|
{\boldsymbol{\theta}}, \boldsymbol{\eta})
\right\|\_{\boldsymbol{\eta}^\*}\cdot h_i(\boldsymbol{u}) +
\mathcal{O}(\\\boldsymbol{u}-\boldsymbol{u}\_\*\\^3) \end{aligned} Note
that the log-likelihood Hessian difference contribution only involves
third order \boldsymbol{u} terms and higher, so the expression above
includes all terms up to second order.

Let g_i^\*=\left.\frac{\partial}{\partial\eta_i} \ln p (\boldsymbol{y}
\| {\boldsymbol{\theta}}, \boldsymbol{\eta})
\right\|\_{\boldsymbol{\eta}^\*} and \boldsymbol{H}^\*\_i =
\left.\nabla\_{u}\nabla\_{u}^\top\widetilde{\eta}\_i(\boldsymbol{u})\right\|\_{\boldsymbol{u}\_\*}
. and form the sum of their products, \boldsymbol{G}=\sum_i
g_i^\*\boldsymbol{H}\_i^\*. Then \begin{aligned} \ln
\widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) - \ln
\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) &=
\frac{1}{2} \sum_i g_i^\* (\boldsymbol{u}-\boldsymbol{u}\_\*)^\top
\boldsymbol{H}\_i^\* (\boldsymbol{u}-\boldsymbol{u}\_\*) +
\mathcal{O}(\\\boldsymbol{u}-\boldsymbol{u}\_\*\\^3) \\&= \frac{1}{2}
(\boldsymbol{u}-\boldsymbol{u}\_\*)^\top \boldsymbol{G}
(\boldsymbol{u}-\boldsymbol{u}\_\*) +
\mathcal{O}(\\\boldsymbol{u}-\boldsymbol{u}\_\*\\^3). \end{aligned} With
\boldsymbol{m}=\mathsf{E}\_\overline{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta})
and
\boldsymbol{Q}^{-1}=\mathsf{Cov}\_\overline{p}(\boldsymbol{u},\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta}),
we obtain \begin{aligned} \mathsf{E}\_{\overline{p}}\left\[
\nabla\_{\boldsymbol{u}} \left\\\ln
\widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) - \ln
\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta})\right\\\right\]
&\approx \boldsymbol{G}(\boldsymbol{m}-\boldsymbol{u}\_\*) , \\
\mathsf{E}\_{\overline{p}}\left\[
\nabla\_{\boldsymbol{u}}\nabla\_{\boldsymbol{u}}^\top \left\\\ln
\widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) - \ln
\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta})\right\\\right\]
&\approx \boldsymbol{G} , \\ \mathsf{E}\_{\overline{p}}\left\[ \ln
\widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta}) - \ln
\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta})\right\]
&\approx \frac{1}{2}
\mathop{\mathrm{tr}}(\boldsymbol{G}\boldsymbol{Q}^{-1}) + \frac{1}{2}
(\boldsymbol{m}-\boldsymbol{u}\_\*)\boldsymbol{G}(\boldsymbol{m}-\boldsymbol{u}\_\*)^\top
. \end{aligned} For each \boldsymbol{\theta} configuration in the INLA
output, we can extract both \boldsymbol{m} and the sparse precision
matrix \boldsymbol{Q} for the Gaussian approximation. The non-sparsity
structure of \boldsymbol{G} is contained in the non-sparsity of
\boldsymbol{Q}, which allows the use of Takahashi recursion
(`inla.qinv(Q)`) to compute the corresponding \boldsymbol{Q}^{-1} values
needed to evaluate the trace
\mathop{\mathrm{tr}}(\boldsymbol{G}\boldsymbol{Q}^{-1}). Thus, to
implement a numerical approximation of this error analysis only needs
special access to the log-likelihood derivatives g_i^\*, as H_i^\* can
in principle be evaluated numerically.

For a given \boldsymbol{\theta}, \begin{aligned}
\mathsf{KL}\left(\overline{p}\\\middle\\\\\widetilde{p}\right) &=
E\_{\overline{p}}\left\[\ln\frac{\overline{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta})}{\widetilde{p}(\boldsymbol{u}\|\boldsymbol{y},\boldsymbol{\theta})}\right\]
\\&= E\_{\overline{p}}\left\[
\ln\frac{\overline{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta})}{\widetilde{p}(\boldsymbol{y}\|\boldsymbol{u},\boldsymbol{\theta})}
\right\] -
\ln\frac{\overline{p}(\boldsymbol{y}\|\boldsymbol{\theta})}{\widetilde{p}(\boldsymbol{y}\|\boldsymbol{\theta})}
. \end{aligned} The first term was approximated above. The second term
can also be approximated using the derived quantities (to be
continued…).

Summary: The form of the observation likelihood discrepancy shows that,
given a linearised posterior
\mathsf{N}(\boldsymbol{m},\boldsymbol{Q}^{-1}), a Gaussian approximation
to the nonlinear model posterior,
\mathsf{N}(\widetilde{\boldsymbol{m}},\widetilde{\boldsymbol{Q}}^{-1}),
can be obtained from
\widetilde{\boldsymbol{Q}}=\boldsymbol{Q}-\boldsymbol{G} and
\widetilde{\boldsymbol{Q}}\widetilde{\boldsymbol{m}}=\boldsymbol{Q}\boldsymbol{m}-\boldsymbol{G}\boldsymbol{u}\_\*.
The K-L divergence becomes \begin{aligned}
\mathsf{KL}\left(\overline{p}\\\middle\\\\\widetilde{p}\right) &\approx
\frac{1}{2} \left\[ \ln\det(\boldsymbol{Q})-
\ln\det(\boldsymbol{Q}-\boldsymbol{G})
-\mathop{\mathrm{tr}}\left(\boldsymbol{G}\boldsymbol{Q}^{-1}\right) +
(\boldsymbol{m}-\boldsymbol{u}\_\*)^\top\boldsymbol{G}(\boldsymbol{Q}-\boldsymbol{G})^{-1}\boldsymbol{G}(\boldsymbol{m}-\boldsymbol{u}\_\*)
\right\] . \end{aligned} When the INLA posterior has mean
\boldsymbol{m}=\boldsymbol{u}\_\*, e.g. for models with additive
Gaussian observation noise, and
\boldsymbol{\theta}=\widehat{\boldsymbol{\theta}}\_{\boldsymbol{u}\_\*},
the last term vanishes.

Note: by implementing the K-L divergence accuracy metric, a by-product
would be improved posterior estimates based on
\widetilde{\boldsymbol{m}} and \widetilde{\boldsymbol{Q}}.

### Well-posedness and initialisation

On a side note, one might be concerned about initialisation at, or
convergence to, a saddle point. Although it is not implemented in
inlabru, we want to talk about the technicality how we define the
initial linearisation point u_0.

Generally speaking, any values of \boldsymbol{u}\_0 work except the case
that the gradient evaluated at \boldsymbol{u}\_0 is \boldsymbol{0}
because the linearisation point will never move away if the prior mean
is also \boldsymbol{0}. In general, this tends to be a saddle point
problem. In some cases the problem can be handled by changing the
predictor parameterisation or just changing the initialisation point
using the `bru_initial` option. However, for true saddle point problems,
it indicates that the predictor parameterisation may lead to a
multimodal posterior distribution or is ill-posed in some other way.
This is a more fundamental problem that cannot be fixed by changing the
initialisation point.

In these examples, where \beta and \boldsymbol{u} are latent Gaussian
components, the predictors 1, 3, and 4 would typically be safe, but
predictor 2 is fundamentally non-identifiable. \begin{aligned}
\boldsymbol{\eta}\_1 &= \boldsymbol{u}, \\ \boldsymbol{\eta}\_2 &= \beta
\boldsymbol{u}, \\ \boldsymbol{\eta}\_3 &= e^\beta \boldsymbol{u}, \\
\boldsymbol{\eta}\_4 &= F\_\beta^{-1} ( \Phi(z\_\beta)) \boldsymbol{u},
\quad z\_{\beta} \sim \mathsf{N}(0,1) . \end{aligned} Note that for
\boldsymbol{\eta}\_3 and \boldsymbol{\eta}\_4, the partial derivatives
with respect to \beta are zero for \boldsymbol{u}=\boldsymbol{0}.
However, the first inlabru iteration will give a non-zero estimate of
\boldsymbol{u}, so that subsequent iteration will involve both \beta and
\boldsymbol{u}.

[^1]: Potential other choices for f(\cdot) include the posterior
    expectation \overline{E}(\boldsymbol{u}\|\boldsymbol{y}) and the
    marginal conditional modes, \left\\\mathop{\mathrm{arg\\max}}\_{u_i}
    \overline{p}\_{\boldsymbol{v}}(u_i\|\boldsymbol{y}),\\i=1,\dots,n\right\\,
    which was used in `inlabru` up to version 2.1.15, which caused
    problems for some nonlinear models. From version 2.2.3, the joint
    conditional mode is used,
    \mathop{\mathrm{arg\\max}}\_{\boldsymbol{u}}
    \overline{p}\_{\boldsymbol{v}}(\boldsymbol{u}\|\boldsymbol{y},\widehat{\boldsymbol{\theta}}\_{\boldsymbol{v}}),
    where
    \widehat{\boldsymbol{\theta}}=\mathop{\mathrm{arg\\max}}\_{\boldsymbol{\theta}}
    \overline{p}\_{\boldsymbol{v}}(\boldsymbol{\theta}\|\boldsymbol{y}).
