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
  # on the domain names. or does it contradict?

  # Check if a list of samplers and domains (Double check if the logic is correct)
  list_samplers <- ifelse(length(samplers) > 1, TRUE, FALSE)
  list_domain <- ifelse(length(domain) > 1, TRUE, FALSE)

  # Domain should be more than samplers
  if (length(samplers) > length(domain)) {
    stop("There are more samplers items than domain ones.")
  } else {
    multi_samplers <- TRUE
  }

  # Check if sf object
  if (class(samplers) %in% "sf") {
    is_sf <- TRUE # fm_as_sfc.inla.mesh
  } else {
    is_sf <- FALSE
    }

  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?
  if (unique(names(samplers)) != unique(names(domains))) {
    warnings("Names of samplers and domain do not match.")
  }

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
