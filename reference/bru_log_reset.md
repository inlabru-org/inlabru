# Clear log contents

Clears the log contents up to a given `offset` or `bookmark`. Default:
clear the entire log. When `x` is NULL, the global `inlabru` log is
updated, and `invisible(NULL)` is returned. Otherwise the updated object
is returned (invisibly).

## Usage

``` r
bru_log_reset(x = NULL, bookmark = NULL, offset = NULL)
```

## Arguments

- x:

  A `bru_log` object, or in some cases, and object that can be
  converted/extracted to a `bru_log` object. `NULL` denotes the global
  `inlabru` log object.

- bookmark:

  character; The label for a bookmark with a stored offset.

- offset:

  integer; a position offset in the log, with `0L` pointing at the start
  of the log. If negative, denotes the point `abs(offset)` elements from
  tail of the log. When `bookmark` is non-NULL, the `offset` applies a
  shift (forwards or backwards) to the bookmark list.

## Value

Returns (invisibly) the modified `bru_log` object, or `NULL` (when `x`
is `NULL`)

## See also

Other inlabru log methods:
[`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md),
[`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md),
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md),
[`bru_log_new()`](https://inlabru-org.github.io/inlabru/reference/bru_log_new.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  bru_log_reset()
}
} # }
```
