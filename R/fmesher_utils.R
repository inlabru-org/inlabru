# @export
fm_caller_name <- function(which = 0L, override = NULL) {
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

# @export
fm_call_stack <- function(start = 0L, end = 0L, with_numbers = TRUE, ...) {
  stack <- sys.calls()
  stack <- lapply(
    as.list(stack),
    function(x) as.character(deparse(x))
  )[
    start + seq_len(max(0, length(stack) - (abs(end) + 1L) - start))
  ]
  if (length(stack) > 0) {
    msg <-
      paste0(
        if (with_numbers) {
          paste0(seq_along(stack), ": ")
        } else {
          ""
        },
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
        )
      )
  } else {
    msg <- "Empty"
  }
  msg
}


# Inspired by berryFunctions::tryStack
# @export
try_callstack <- function(expr) {
  try_envir <- new.env()
  assign("error_stack", value = NULL, envir = try_envir)
  error_fun <- function(e) {
    # Get whole stack except the handlers
    stack <- fm_call_stack(start = 0, end = 2, with_numbers = FALSE)
    # Remove the try_callstack tryCatch calls part(s),
    # There are 6 of them. First find the try_callstack call (or multiple calls
    # for nested use, which should theoretically (almost) never happen,
    # since the inner call shouldn't fail!)
    self <- which(
      vapply(stack, function(x) grepl("^try_callstack\\(", x), TRUE) |
        vapply(stack, function(x) grepl("^inlabru:::try_callstack\\(", x), TRUE)
    )
    for (idx in rev(self)) {
      stack <- stack[-(idx + seq_len(6))]
      stack[idx] <- "try_callstack(...)"
    }
    stack <- paste0(seq_len(length(stack)), ": ", stack, collapse = "\n")
    assign("error_stack", value = stack, envir = try_envir)
  }
  result <- try(
    withCallingHandlers(
      expr,
      error = error_fun
    ),
    silent = TRUE
  )
  if (inherits(result, "try-error")) {
    result[length(result) + 1] <- paste0(
      try_envir$error_stack,
      collapse = "\n"
    )
  }
  invisible(result)
}
