# @aliases apmaker
# @export
# @param samplers A `[sf]DataFrame` or Spatial[Points/Lines/Polygons]DataFrame`
# object
# TODO 20221109 does domainsssss make more sense?
# @param domain A list/A list of list(s) of named integration definitions, each
# either a vector of factors ,a numeric vector of points given integration
#  weight 1, an `inla.mesh.1d` object, or an `inla.mesh.2d` object. Only those
# domains that are not given in the `samplers` data.frame are used, plus the
# coordinates object, used for the spatial aspect of the `samplers` object.
# @param dnames Names of dimensions
# @param int.args List of arguments passed on to \code{ipoints}
# @return Integration points


apmaker <- function(samplers, domain, dnames,
                    int.args = list(method = "stable", nsub = NULL)) {
  # To allow sf geometry support, should likely change the logic to
  # use the domain specification to determine the type of integration
  # method to call, so that it doesn't need to rely on the domain name.
  # TODO 20221109 For multiple samplers multiple domains, it does have to rely
  # on the domain names. 20221111 They do have to match

  # Mandate both the samplers and domain arguments
  if (is.null(samplers) || is.null(domain)) {
    stop("Samplers or domain argument(s) missing.")
  }

  # Check if both samplers and domain list
  # https://stackoverflow.com/questions/38539654/how-to-know-if-the-data-is-a-list-or-data-frame-in-r
  if (inherits(samplers, "list") && inherits(domain, "list")) {
    is_list <- TRUE
  } else if (inherits(samplers, "list")) {
    multisampler_int <- TRUE
  } else if (inherits(domain, "list")) {
    singlesampler_int <- TRUE
  } else {
    is_list <- FALSE
  }

  # check sf or sp object
  if (is_list) {
    if (sapply(samplers, class) %in% "sf" || sapply(domain, class) %in% "sf") {
      is_sf <- TRUE # fm_as_sfc.inla.mesh
    } else if (sapply(samplers, class) %in% "sp" || sapply(domain, class) %in% "sp") {
        is_sp <- FALSE
      }{
      is_sf <- FALSE
    }
  }

  # How does sf deal with secondary geometry columns?
  # https://r-spatial.github.io/sf/articles/sf6.html
  if (is_sf) {
    cat("The active geometry is", attr(samplers, "sf_column"))
    cat("The active geometry is", attr(domain, "sf_column"))
  }

  # Domain should be more than samplers TODO this is the problem
  if (length(samplers) > length(domain)) {
    stop("There are more samplers items than domain ones.")
  }

  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?
  if (unique(names(samplers)) != unique(names(domain))) {
    warnings("Names of samplers and domain do not match.")
  }

  #####################################

  if ("coordinates" %in% dnames) {
    spatial <- TRUE
  } else {
    spatial <- FALSE
  }

  # Dimensions provided via samplers (except "coordinates")
  samp.dim <- intersect(names(samplers), dnames)

  # Dimensions provided via domain but not via samplers
  nosamp.dim <- setdiff(names(domain), c(samp.dim, "coordinates"))

  # Check if a domain definition is missing
  missing.dims <- setdiff(dnames, c(names(domain), samp.dim))
  if (length(missing.dims > 0)) {
    stop(paste0(
      "Domain definitions missing for dimensions: ",
      paste0(missing.dims, collapse = ", ")
    ))
  }
  extra.dims <- setdiff(names(domain), c(samp.dim, nosamp.dim))
  if (length(missing.dims > 0)) {
    warning(paste0(
      "Unexpected extra domain defintions: ",
      paste0(extra.dims, collapse = ", ")
    ))
  }

  if (spatial) {
    ips <- ipoints(samplers, domain$coordinates,
      group = samp.dim, int.args = int.args
    )
  } else {
    ips <- NULL
  }

  lips <- lapply(nosamp.dim, function(nm) ipoints(NULL, domain[[nm]], name = nm, int.args = int.args))
  ips <- do.call(cprod, c(list(ips), lips))
}
