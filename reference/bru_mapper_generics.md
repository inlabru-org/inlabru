# Generic methods for bru_mapper objects

A `bru_mapper` sub-class implementation must provide an
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
method. If the model size 'n' and definition values 'values' are stored
in the object itself, default methods
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
and
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
are available. Otherwise the
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
and
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
methods also need to be provided.

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
for constructor methods, and
[bru_get_mapper](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
for hooks to extract mappers from latent model object class objects.

Other mapper methods:
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
