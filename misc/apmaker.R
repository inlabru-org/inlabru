#' A S3/S4 generic function to generate stop warning if the object of some class
#'  is not of the specified classes.
#' @param  object An object except NULL
#' @param class A vector of class(es) except NULL
#' @keywords internal
#'
stopifnot_class <- function(object = NULL, class = NULL, ...) {
  # 20221208 Probably better to leave out NULL and missing class in this
  # function, though it deals with it
  class_name <- paste0(sort(class), collapse = ",")
  obj_name <- paste0(as.character(expression(object)))
  UseMethod("stopifnot_class")
}
stopifnot_class.default <- function(object, class) {
  if (!inherits(object, class)) {
    stop(paste0(
      obj_name,
      " of class ",
      class(object),
      "must be of class(es) ",
      class_name
    ))
  }
}
stopifnot_class.list <- function(object, class) {
  if (!all(unlist(lapply(object, function(x) {
    inherits(x, class)
  })))) {
    stop(paste0(
      obj_name,
      " of class ",
      paste0(unlist(lapply(object, class)), collapse = ","),
      "must be of class(es) ",
      class_name
    ))
  }
}

# bru_log_smpdi function ---------------------------------------------------------
# How does sf deal with secondary geometry columns?
# TODO ####
# TODO have to deal with the multidomain samplers
# https://r-spatial.github.io/sf/articles/sf6.html
# TODO #### extra domains we assign the sampler for them
# TODO can use local helper function (with S3 Usemethod() maybe)
#' @param begin the beginning part of the message
bru_log_smpdi <- function(x, start_with = NULL, ...) {
  UseMethod("bru_log_smpdi")
}
bru_log_smpdi.list <- function(x, start_with) {
  for (i in seq_along(x)) {
    bru_log_message(
      paste0(
        start_with,
        x[i],
        " are ",
        attr(x[i], "sf_column"),
        ".\n"
      ),
      verbosity = 2, verbose_store = T
    )
  }
}

#' @aliases apmaker
#' @export
#' @param domain A list/A list of list(s) of named integration definitions, each
#' either character/factor vector, a numeric vector of points given integration
#'  weight 1, an `inla.mesh.1d` object, or an `inla.mesh.2d` object. Only those
#' domains that are not given in the `samplers` data.frame are used, plus the
#' coordinates object, used for the spatial aspect of the `samplers` object.
#' @param samplers A (list of) `[sf]DataFrame` or
#' `Spatial[Points/Lines/Polygons]DataFrame` object(s)
#' @param weight The name of integration weight column in the samplers.
#' TODO 1) how about domain? should not be allowed.
#' 2) allow_names should be bru_weight? 23112022 Only for samplers
#' Default: "weight".
#' It mainly works for line transect at the moment, determining the width of the
#' line to be integrated. See the distance sampling example
#' https://inlabru-org.github.io/inlabru/articles/web/2d_lgcp_distancesampling.html
#' @param int.args List of arguments passed on to \code{ipoints}
#' @return Integration points

# TODO option argument as in bru function with list() bru_int_args
apmaker <- function(domain = NULL, samplers = NULL,
                    weight = "weight",
                    int.args = list(method = "stable", nsub = NULL)) {
  # To allow sf geometry support, should likely change the logic to
  # use the domain specification to determine the type of integration
  # method to call, so that it doesn't need to rely on the domain name.
  # TODO ####
  # TODO handle the weight
  # TODO handle domain not spatial object, can be time points
  # TODO use the unnamed and named columns in list
  # TODO 20221109 For multiple samplers multiple domains, it does have to rely
  # on the domain names. 20221111 They do have to match

  # a template for deprecate
  # 1) dname deprecated
  # 2) samplers should be provided as list
  # 3) ...
  lifecycle::deprecate_soft(
    when = "2.7.1", # TODO next inlabru version
    what = "apmaker(dnames)",
    with = "apmaker(samplers)",
    details = "dnames should be provided in the samplers as the dimension names."
  )

  # Mandate the domain argument to be specified
  if (missing(domain)) {
    stop("Domain argument(s) missing.")
  }
  # Mandate one of domain and samplers to be non-NULL
  if (is.null(domain)) {
    if (is.null(samplers)) {
      stop("Domain and samplers argument(s) NULL.")
    }
    warning("Domain argument(s) NULL. Sampler is used to create the domain")
  }

  # https://stackoverflow.com/questions/38539654/how-to-know-if-the-data-is-a-list-or-data-frame-in-r
  # Turn data frame into a list (standardise the input class)
  if (inherits(domain, "list")) {
    domain_is_list <- TRUE
    # If domain is a list, it must be a named list.
    if (is.null(names(domain))) {
      stop("Domain must be a named list.")
    }
  } else {
    domain_is_list <- FALSE
    domain <- list(domain)
  }

  if (inherits(samplers, "list")) {
    sampler_is_list <- TRUE
  } else {
    sampler_is_list <- FALSE
    samplers <- list(samplers)
  }

  # 20221208 sf sfc Spatial should match inla.mesh, inla.mesh.lattice, raster, SpatRaster
  # character/factor/numeric should match the same classes of domain.
  # numeric is implemented with weight 1 in ipoints. TODO character and factor
  # to be implemented similarly but not changed into numeric.

  # Mandate domain to be the following classes
  # Unused domain is not allowed
  stopifnot_class(
    domain,
    c( # "data.frame", # maybe in the future
      "character",
      "factor",
      "numeric", # here is the split
      "inla.mesh",
      "inla.mesh.1d",
      "inla.mesh.lattice",
      "raster",
      "SpatRaster"
    )
  )
  # Mandate domain to be the following classes
  stopifnot_class(
    samplers,
    c(
      "character", "factor", "numeric", # here is the split
      "data.frame", # 20221212 ipoints allows data.frame samplers
      "sf", "sfc", "Spatial"
    )
  )

  # 20221212 We have to check the matching here, that is
  # 1) certain domain classes  have to match certain samplers classes.
  # TODO 2) in which 1) has to check the domain one by one if it is a list
  if (any(unlist(lapply(domain, function(x) {
    inherits(x, c("character", "factor", "numeric", "inla.mesh.1d"))
  })))) {
    stopifnot_class(samplers, c("character", "factor", "numeric", "data.frame"))
  }
  if (any(unlist(lapply(domain, function(x) {
    inherits(x, cc(
      "inla.mesh",
      "inla.mesh.lattice",
      "raster",
      "SpatRaster"
    ))
  })))) {
    stopifnot_class(samplers, c("sf", "sfc", "Spatial"))
  }

  # Change a mix of sp and sf objects to sf
  sf_samplers <- any(unlist(lapply(samplers, function(x) inherits(x, c("sf", "sfc")))))
  sp_samplers <- any(unlist(lapply(samplers, function(x) inherits(x, "Spatial"))))
  if (sp_samplers && sf_samplers) {
    warning("Both `sf` and `sp` objects in the samplers are detected. Produce `sf` output as desired")
    samplers <- lapply(samplers, sf::st_as_sf)
  }

  # TODO Sort multidomain samplers, single domain samplers,
  # TODO Multidomain samplers happens when a sampler across several domains. How to detect that?
  # TODO we should then do the name check for domain and samplers here
  # TODO remove sampler domains and full domain samplers
  # multidomain sampler the main thing is data dname 23112022 specific column has that meaning
  names_domain <- names(domain)
  names_samplers <- names(samplers)
  names_diff <- setdiff(intersect(names_domain, names_samplers), "weight")
  if(domain_is_list && (length(names_diff)>0)){
    warning("The difference in the names of domain and samplers are ",
            paste0(names_diff,collapse = ","),
            "\n Create samplers for unmathced domains and do nothing on unused samplers.")
  }

  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?
  if (unique(names_samplers) != unique(names_domain)) {
    # TODO have to accommodate weight column in samplers as well, omit this column
    samplers_domain <- intersect(names_samplers, names_domain)
    # TODO for the NULL name case
    bru_log_smpdi(samplers_domain, start_with = "\n The shared samplers and domain: \n")
  }
  # Domain should be more than samplers. However, we can have unused samplers as well.
  else if (length(samplers) > length(domain)) {
    unused_samplers <- setdiff(names_samplers, names_domain)
    bru_log_smpdi(unused_samplers, start_with = "\n The unused samplers: \n")
  } else {
    extra_domain <- setdiff(names_samplers, names_domain)
    bru_log_smpdi(extra_domain, start_with = "\n The extra domain: \n")
  }

  # log the active geometry of the samplers
  if (sf_samplers) {
    bru_log_smpdi(samplers,
      start_with = "The active geometry of the sampler(s)"
    )
  }


  # TODO ####
  # TODO check ibm_values.bru_mapper_factor for character/factor/numeric samplers
  # TODO weight argument to take effect on samplers. The weight should go to
  # the integration part
  ips <- apoints(samplers, domain,
    int.args = int.args, weight_name = weight_name
  )
  # TODO using group_cprod for each domain
  # this is just a draft atm
  lips <- lapply(nosamp.dim, function(nm)
            ipoints(NULL, domain[[nm]], name = nm, int.args = int.args))
  ips <- do.call(group_cprod, c(list(ips), lips, group = samplers_domain))

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
