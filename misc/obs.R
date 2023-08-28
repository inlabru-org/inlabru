method_classes <- function(f) {
  gsub(
    pattern = paste0("^", f, "\\.([^*]*)\\*?"),
    replacement = "\\1",
    x = format(utils::.S3methods(f))
  )
}


#' Observation model objects
#'
#' Construct an observation model object; new modular alternative to [like()]
#' @param ... Arguments passed on to submethods
#' @export
bru_obs <- function(...) {
  UseMethod("bru_obs")
}

#' @export
#' @rdname bru_obs
bru_obs.character <- function(x, ...) {
  proto_class <- paste0("proto_bru_obs_", x)
  if (proto_class %in% setdiff(method_classes("bru_obs"), "character")) {
    return(do.call(bru_obs, list(structure(list(), class = proto_class), ...)))
  }
  return(structure(list(model = x), class = c("bru_obs_character", "bru_obs")))
}

#' @export
#' @rdname bru_obs
bru_obs.proto_bru_obs_cp <- function(x, E, ...) {
  return(structure(list(model = "possion", E = E), class = c("bru_obs_cp", "bru_obs")))
}
