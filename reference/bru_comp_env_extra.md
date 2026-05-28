# Get/set component environment data

Get or set data in the component's `env_extra` and `env` environments.
If `name` is `NULL`, the entire environment is returned.

## Usage

``` r
bru_comp_env_extra(x, name = NULL)

bru_comp_env_extra(x, name) <- value

bru_comp_env(x, name = NULL)

bru_comp_env(x, name) <- value
```

## Arguments

- x:

  A `bru_comp` object

- name:

  The name of the variable to get or set. Names prefixed by `"BRU_"` are
  reserved for internal use and should not be set by users or external
  package authors.

- value:

  The value to set for `name` in `env_extra` (for the setter function).

## Value

For the getter, either the entire environment or the value of `name` in
the environment. For the setters, the modified `bru_comp` object.

## Functions

- `bru_comp_env_extra(x, name) <- value`: Set data in the component's
  `env_extra` environment

- `bru_comp_env()`: Get the component's `env` environment, or an element
  from that environment. Note that in most cases this environment is the
  global R environment.

- `bru_comp_env(x, name) <- value`: Set data in the component's `env`
  environment. Note that in most cases this environment is the global R
  environment, so modifying it is normally not recommended.

## See also

bru_comp

## Examples

``` r
if (bru_safe_inla()) {
  cmp <- bru_comp("x", x)
  as.list(bru_comp_env_extra(cmp))

  cmp <- bru_comp("x", x, model = "fixed")
  as.list(bru_comp_env_extra(cmp))

  bru_comp_env(cmp)
}
#> <environment: 0x624ab5dac8d0>
```
