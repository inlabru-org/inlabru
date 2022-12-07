#' @aliases apmaker
#' @export
#' @param domain A list/A list of list(s) of named integration definitions, each
#' either character/factor vector, a numeric vector of points given integration
#'  weight 1, an `inla.mesh.1d` object, or an `inla.mesh.2d` object. Only those
#' domains that are not given in the `samplers` data.frame are used, plus the
#' coordinates object, used for the spatial aspect of the `samplers` object.
#' @param samplers A (list of) `[sf]DataFrame` or
#' `Spatial[Points/Lines/Polygons]DataFrame` object(s)
#' @param weights The name of integration weights column in the samplers.
#' TODO 1) how about domain? should not be allowed.
#' 2) allow_names should be bru_weight? 23112022 Only for samplers
#' Default: "weights".
#' It mainly works for line transect at the moment, determining the width of the
#' line to be integrated. See the distance sampling example
#' https://inlabru-org.github.io/inlabru/articles/web/2d_lgcp_distancesampling.html
#' @param int.args List of arguments passed on to \code{ipoints}
#' @return Integration points

# TODO option argument as in bru function with list() bru_int_args
apmaker <- function(domain = NULL, samplers = NULL,
                    weights = "weights",
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

  # Mandate both the domain argument not specified or is null
  if (missing(domain) || is.null(domain)) {
    stop("Domain argument(s) missing or NULL.")
  }

  #' A S3/S4 generic function to generate stop warning if the object is not of
  #' the specified classes.
  #' @param  object An object
  #' @param class A vector of class(es)
  #' @keywords internal
  stopifnot_class <- function(object = NULL, class = NULL, ...) {
    # 20221207 Probably better to leave it out not to deal with NULL and missing in this function
    # if (is.null(object) || missing(object)) {
    #   stop(paste0(as.character(expression(object))), " argument missing or NULL.") # TODO It does not work object is NULL
    # }
    # if (is.null(class) || missing(class)) {
    #   stop("Class argument(s) missing or NULL.")
    # }
    class_name <- paste0(sort(class), collapse = ",")
    obj_name <- paste0(as.character(expression(object)), collapse = ",") # TODO checking several objects the same time
    UseMethod("stopifnot_class")
  }
  stopifnot_class.default <- function(object, class) {
    if (!inherits(object, class)) {
      stop(paste0(obj_name), " must be of class(es) ", class_name)
    }
  }
  stopifnot_class.list <- function(object, class) {
    if (!all(unlist(lapply(object, function(x) {
      inherits(x, class)
    })))) {
      stop(paste0(obj_name), " must be of class(es) ", class_name)
    }
  }

  # Mandate domain to be the following classes
  # Unused domain is not allowed
  stopifnot_class(
    domain,
    c(
      # "data.frame", # maybe in the future
      "character",
      "factor",
      "numeric",
      "inla.mesh",
      "inla.mesh.1d",
      "inla.mesh.lattice",
      "raster",
      "spatraster"
    )
  )

  # Mandate domain to be the following classes
  stopifnot_class(samplers,
                  c("numeric", "character", "factor", "sf", "sfc", "Spatial"))

  # https://stackoverflow.com/questions/38539654/how-to-know-if-the-data-is-a-list-or-data-frame-in-r
  # Turn data frame into a list (standardise the input class)
  if (inherits(domain, "list")) {
    is_list_domain <- TRUE
    if (inherits(samplers, "list")) {
      is_list <- TRUE
    }
  } else if (inherits(samplers, "list") && !is_list_domain) {
    single_domain <- TRUE
    domain <- as.list(domain)
  } else {
    is_list <- FALSE
      domain <- as.list(domain)
      samplers <- as.list(samplers)
  }

  # Change a mix of sp and sf objects to sf
  sf_samplers <- lapply(samplers, function(x) inherits(x, c("sf", "sfc")))
  sp_samplers <- lapply(samplers, function(x) inherits(x, "Spatial"))
  if (sp_samplers && sf_samplers) {
    warning("Both sf and sp objects in the samplers are detected. Produce sf output as desired")
    samplers <- lapply(samplers, sf::st_as_sf)
  }

  # TODO Sort multidomain samplers, single domain samplers,
  # TODO Multidomain samplers happens when a sampler across several domains. How to detect that?
  # TODO we should then do the name check for domain and samplers here
  # TODO remove sampler domains and full domain samplers
  # multidomain sampler the main thing is data dname 23112022 specific column has that meaning


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
  # Technically speaking
  # bru_log_smpdi.data.frame <- function(x){
  #   bru_log_message(
  #     paste0(
  #       start_with,
  #       x,
  #       " is ",
  #       attr(x[i], "sf_column"),
  #       ".\n"
  #     ),
  #     verbosity = 2, verbose_store = T
  #   )
  # }

  if (sf_samplers) {
    bru_log_smpdi(samplers,
      start_with = "The active geometry of the sampler"
    )
  }

  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?
  if (unique(names(samplers)) != unique(names(domain))) {
    # TODO have to accommodate weights column in samplers as well, omit this column
    samplers_domain <- intersect(names(samplers), names(domain))
    # TODO for the NULL name case
    bru_log_smpdi(samplers_domain, start_with = "\n The shared samplers and domain: \n")
  }
  # Domain should be more than samplers. However, we can have unused samplers as well.
  else if (length(samplers) > length(domain)) {
    unused_samplers <- setdiff(names(samplers), names(domain))
    bru_log_smpdi(unused_samplers, start_with = "\n The unused samplers: \n")
    bru_log_smpdi(unused_samplers, start_with = "\n The unused samplers: \n")
  } else {
    extra_domain <- setdiff(names(samplers), names(domain))
    bru_log_smpdi(extra_domain, start_with = "\n The extra domain: \n")
  }


  # TODO weights argument to take effect on samplers. The weight should go to
  # the integration part
  ips <- apoints(samplers, domain,
    int.args = int.args, weights_name = weights_name
  )
  # TODO using groupwise cprod for each domain
  # this is just a draft atm
  lips <- lapply(nosamp.dim, function(nm) ipoints(NULL, domain[[nm]], name = nm, int.args = int.args))
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
