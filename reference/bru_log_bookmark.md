# Methods for `bru_log` bookmarks

Methods for `bru_log` bookmarks.

## Usage

``` r
bru_log_bookmark(bookmark = "", offset = NULL, x = NULL)

bru_log_bookmarks(x = NULL)
```

## Arguments

- bookmark:

  character; The label for a bookmark with a stored offset.

- offset:

  integer; a position offset in the log, with `0L` pointing at the start
  of the log. If negative, denotes the point `abs(offset)` elements from
  tail of the log. When `bookmark` is non-NULL, the `offset` applies a
  shift (forwards or backwards) to the bookmark list.

- x:

  A `bru_log` object. If `NULL`, the global `inlabru` log is used.

## Value

`bru_log_bookmark()`: Returns the modified `bru_log` object if `x` is
non-NULL.

`bru_log_bookmarks()`: Returns the bookmark vector associated with `x`

## Functions

- `bru_log_bookmark()`: Set a log bookmark. If `offset` is `NULL` (the
  default), the bookmark will point to the current end of the log.

- `bru_log_bookmarks()`: Return a integer vector with named elements
  being bookmarks into the global `inlabru` log with associated log
  position offsets.

## See also

Other inlabru log methods:
[`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md),
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md),
[`new_bru_log()`](https://inlabru-org.github.io/inlabru/reference/new_bru_log.md)
