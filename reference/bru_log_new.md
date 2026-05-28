# Create a `bru_log` object

Create a `bru_log` object, by default empty.

## Usage

``` r
bru_log_new(x = NULL, bookmarks = NULL)
```

## Arguments

- x:

  An optional character vector of log messages, or `data.frame` with
  columns `message`, `timestamp`, and `verbosity`.

- bookmarks:

  An optional `integer` vector of named bookmarks message in `x`.

## See also

Other inlabru log methods:
[`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md),
[`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md),
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md)

## Examples

``` r
x <- bru_log_new()
x <- bru_log_message("Test message", x = x)
print(x)
#> 2026-05-26 12:46:05.0878: Test message
```
