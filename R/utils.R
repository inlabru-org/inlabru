show_call_stack <- function() {
    stack <- sys.calls()
    stack <- lapply(as.list(stack),
                    function(x) as.character(deparse(x)))[-length(stack)]
    if (length(stack) > 0) {
      message(paste0("Call stack:\n",
                     paste0(seq_along(stack),
                            ": ",
                            lapply(stack, 
                                   function(x) substr(x, 1, 70)),
                            collapse = "\n")))
    } else {
      message(paste0("Call stack is empty"))
    }
}

sp_get_crs <- function(...) {
  INLA::inla.sp_get_crs(...)
}

crs_is_null <- function(crs) {
  if (is.null(crs)) {
    TRUE
  } else if (INLA::inla.has_PROJ6()) {
    wkt <- INLA::inla.crs_get_wkt(crs)
    is.null(wkt)
  } else {
    is.na(crs)
  }
}
