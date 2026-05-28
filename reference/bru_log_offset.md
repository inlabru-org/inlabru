# Position methods for `bru_log` objects

Position methods for `bru_log` objects.

## Usage

``` r
bru_log_offset(x = NULL, bookmark = NULL, offset = NULL)

bru_log_index(x = NULL, i, verbosity = NULL)
```

## Arguments

- x:

  A `bru_log` object. If `NULL`, the global `inlabru` log is used.

- bookmark:

  character; The label for a bookmark with a stored offset.

- offset:

  integer; a position offset in the log, with `0L` pointing at the start
  of the log. If negative, denotes the point `abs(offset)` elements from
  tail of the log. When `bookmark` is non-NULL, the `offset` applies a
  shift (forwards or backwards) to the bookmark list.

- i:

  indices specifying elements to extract. If `character`, denotes the
  sequence between bookmark `i` and the next bookmark (or the end of the
  log if `i` is the last bookmark)

- verbosity:

  integer value for limiting the highest verbosity level being returned.

## Functions

- `bru_log_offset()`: Utility function for computing log position
  offsets.

- `bru_log_index()`: Utility function for computing index vectors for
  `bru_log` objects.

## See also

Other inlabru log methods:
[`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md),
[`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md),
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md),
[`bru_log_new()`](https://inlabru-org.github.io/inlabru/reference/bru_log_new.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md)
