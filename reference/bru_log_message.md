# Add a log message

Adds a log message.

## Usage

``` r
bru_log_message(
  ...,
  domain = NULL,
  appendLF = TRUE,
  verbosity = 1L,
  allow_verbose = TRUE,
  verbose = NULL,
  verbose_store = NULL,
  x = NULL
)

bru_log_abort(
  msg,
  ...,
  domain = NULL,
  appendLF = TRUE,
  verbosity = 1L,
  allow_verbose = TRUE,
  verbose = FALSE,
  verbose_store = NULL,
  call = rlang::caller_env(),
  .frame = rlang::caller_env()
)

bru_log_warn(
  msg,
  ...,
  domain = NULL,
  appendLF = TRUE,
  verbosity = 1L,
  allow_verbose = TRUE,
  verbose = FALSE,
  verbose_store = NULL,
  call = rlang::caller_env(),
  .frame = rlang::caller_env()
)
```

## Arguments

- ...:

  For `bru_log_message()`, zero or more objects passed on to
  [`base::.makeMessage()`](https://rdrr.io/r/base/message.html). For
  `bru_log_abort()` and `bru_log_warn()`, passed on to
  [`rlang::abort()`](https://rlang.r-lib.org/reference/abort.html) and
  [`rlang::warn()`](https://rlang.r-lib.org/reference/abort.html).

- domain:

  Domain for translations, passed on to
  [`base::.makeMessage()`](https://rdrr.io/r/base/message.html)

- appendLF:

  logical; whether to add a newline to the message. Only used for
  verbose output.

- verbosity:

  numeric value describing the verbosity level of the message

- allow_verbose:

  Whether to allow verbose output. Must be set to FALSE until the
  options object has been initialised.

- verbose:

  logical, numeric, or `NULL`; local override for verbose output. If
  `NULL`, the global option `bru_verbose` or default value is used. If
  `FALSE`, no messages are printed. If `TRUE`, messages with `verbosity`
  \\\le 1\\ are printed. If numeric, messages with `verbosity` \\\le\\
  `verbose` are printed.

- verbose_store:

  Same as `verbose`, but controlling what messages are stored in the
  global log object. Can be controlled via the `bru_verbose_store` with
  [`bru_options_set()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md).

- x:

  A `bru_log` object. If `NULL`, refers to the global `inlabru` log.

- msg:

  character; passed to
  [`base::.makeMessage()`](https://rdrr.io/r/base/message.html)

- call:

  The calling environment.

- .frame:

  The throwing context, for when `.internal` is `TRUE`

## Value

`bru_log_message` returns `invisible(x)`, where `x` is the updated
`bru_log` object, or `NULL`.

## Functions

- `bru_log_abort()`: Store a log message and throw an error.

- `bru_log_warn()`: Store a log message and throw a warning.

## See also

Other inlabru log methods:
[`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md),
[`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md),
[`bru_log_new()`](https://inlabru-org.github.io/inlabru/reference/bru_log_new.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md)

## Examples

``` r
if (interactive()) {
  code_runner <- function() {
    bru_options_set_local(
      # Show messages up to and including level 2 (default 0)
      bru_verbose = 2,
      # Store messages to an including level 3 (default Inf, storing all)
      bru_verbose_store = 3
    )

    bru_log_bookmark("bookmark 1")
    bru_log_message("Test message 1", verbosity = 1)
    bru_log_message("Test message 2", verbosity = 2)
    bru_log_bookmark("bookmark 2")
    bru_log_message("Test message 3", verbosity = 3)
    bru_log_message("Test message 4", verbosity = 4)

    invisible()
  }
  message("Run code")
  code_runner()
  message("Check log from bookmark 1")
  print(bru_log()["bookmark 1"])
  message("Check log from bookmark 2")
  print(bru_log()["bookmark 2"])
}
```
