show_call_stack <- function() {
    stack <- sys.calls()
    stack <- lapply(as.list(stack),
                    function(x) as.character(deparse(x)))[-length(stack)]
    if (length(stack) > 0) {
      message(
        paste0(
          "Call stack:\n",
          paste0(
            seq_along(stack),
            ": ",
            lapply(
              stack, 
              function(x) {
                paste0(
                  vapply(x,
                         function(x) {
                           if (nchar(x) > 80) {
                             paste0(
                               substr(x, 1, 74),
                               " [...]")
                           } else {
                             x
                           }
                         },
                         ""),
                  collapse = paste0("\n   ")
                )
              }
            ),
            collapse = "\n")
        )
      )
    } else {
      message(paste0("Call stack is empty"))
    }
}
