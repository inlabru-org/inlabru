# Create a `bru_log` object

Create a `bru_log` object, by default empty.

## Usage

``` r
new_bru_log(x = NULL, bookmarks = NULL)
```

## Arguments

- x:

  An optional character vector of log messages, or `data.frame` with
  columns `message`, `timestamp`, and `verbosity`, or a `bru_log`
  object.

- bookmarks:

  An optional `integer` vector of named bookmarks message in `x`.

## Value

A new `bru_log` object.

## See also

Other inlabru log methods:
[`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md),
[`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md),
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md)

## Examples

``` r
x <- new_bru_log()
x <- bru_log_message("Test message", x = x)
print(x)
#> 2026-07-27 22:35:05.074357: Test message
```
