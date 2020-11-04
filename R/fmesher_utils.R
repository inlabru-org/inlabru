fn_caller_name <- function(which = 0L, override = NULL) {
  if (is.null(override)) {
    which <- -abs(which) - 1L
    if (abs(which) > sys.nframe()) {
      name <- ""
    } else {
      fun <- sys.call(which)
      if (is.null(fun)) {
        name <- ""
      } else {
        name <- as.character(fun)[[1]]
      }
    }
  } else {
    name <- override
  }
  name
}

fn_call_stack <- function(start = 0L, ...) {
  stack <- sys.calls()
  stack <- lapply(
    as.list(stack),
    function(x) as.character(deparse(x))
  )[
    seq_len(max(0, length(stack) - (abs(start) + 1L)))
  ]
  if (length(stack) > 0) {
    msg <-
      paste0(
        "Call stack:\n",
        paste0(
          seq_along(stack),
          ": ",
          lapply(
            stack,
            function(x) {
              paste0(
                vapply(
                  x,
                  function(x) {
                    if (nchar(x) > 80) {
                      paste0(
                        substr(x, 1, 74),
                        " [...]"
                      )
                    } else {
                      x
                    }
                  },
                  ""
                ),
                collapse = paste0("\n   ")
              )
            }
          ),
          collapse = "\n"
        )
      )
  } else {
    msg <- paste0("Call stack: Empty")
  }
  msg
}
