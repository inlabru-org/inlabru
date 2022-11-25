# @aliases apmaker
# @export
# @param samplers A (list of) `[sf]DataFrame` or
# Spatial[Points/Lines/Polygons]DataFrame object(s)
# TODO 20221109 does domainsssss make more sense?
# @param domain A list/A list of list(s) of named integration definitions, each
# either a vector of factors ,a numeric vector of points given integration
#  weight 1, an `inla.mesh.1d` object, or an `inla.mesh.2d` object. Only those
# domains that are not given in the `samplers` data.frame are used, plus the
# coordinates object, used for the spatial aspect of the `samplers` object.
# @param weights The name of integration weights column in the samplers.
# TODO 1) how about domain? should not be allowed.  2) allow_names should be bru_weight? 23112022 Only for samplers
# Default: "weights".
# It mainly works for line transect at the moment, determining the width of the
# line to be integrated. See the distance sampling example
# https://inlabru-org.github.io/inlabru/articles/web/2d_lgcp_distancesampling.html
# @param int.args List of arguments passed on to \code{ipoints}
# @return Integration points

# TODO option argument as in bru function with list() bru_int_args
apmaker <- function(domain, samplers,
                    weights = "weights",
                    int.args = list(method = "stable", nsub = NULL)) {
  # To allow sf geometry support, should likely change the logic to
  # use the domain specification to determine the type of integration
  # method to call, so that it doesn't need to rely on the domain name.
  # TODO 20221109 For multiple samplers multiple domains, it does have to rely
  # on the domain names. 20221111 They do have to match

  # Mandate both the domain argument not specified or is null
  if (missing(domain) || is.null(domain)) {
    stop("Domain argument(s) missing or NULL.")
  }

  # Check if both samplers and domain list
  # https://stackoverflow.com/questions/38539654/how-to-know-if-the-data-is-a-list-or-data-frame-in-r
  if (inherits(samplers, "list") && inherits(domain, "list")) {
    is_list <- TRUE
  } else if (inherits(samplers, "list")) {
    singlesampler_int <- TRUE
  } else if (inherits(domain, "list")) {
    multisampler_int <- TRUE
  } else {
    is_list <- FALSE
  }

  # multi domain sampler the main thing is data dname 23112022 specific column has that meaning
  # check sf or sp object, can do a mix of sp and sf objects,
  # TODO if not a list, we should have to check
  # TODO have to convert sp to sf for output if there is a mix of sp and sf
  if (is_list) {
    sf_samplers <- lapply(samplers, function(x) inherits(x, "sf"))
    sf_domain <- lapply(domain, function(x) inherits(x, "sf"))
    sp_samplers <- lapply(samplers, function(x) inherits(x, "Spatial"))
    sp_domain <- lapply(domain, function(x) inherits(x, "Spatial"))
    if (any(sf_domain) && any(sp_samplers)) {
      samplers <- lapply(samplers, sf::st_as_sf)
    }
  }
  # single domain samplers
  if (singlesampler_int) {
    sf_samplers <- lapply(samplers, function(x) inherits(x, "sf"))
    sf_domain <- inherits(domain, "sf")
    sp_samplers <- lapply(samplers, function(x) inherits(x, "Spatial"))
    sp_domain <- inherits(domain, "Spatial")
    if (sf_domain && any(sp_samplers)) {
      samplers <- sf::st_as_sf(samplers)
    }
  }
  # multi domain sampler
  if (multisampler_int) {
    sf_samplers <- inherits(samplers, "sf")
    sf_domain <- lapply(domain, function(x) inherits(x, "sf"))
    sp_samplers <- inherits(samplers, "Spatial")
    sp_domain <- lapply(domain, function(x) inherits(x, "Spatial"))
    if (sf_domain && any(sp_samplers)) {
      samplers <- lapply(samplers, sf::st_as_sf)
    }
  }
  # have to sort this out beforehand then this cannot happen 23112022
  if (!is_list) {
    sf_samplers <- inherits(samplers, "sf")
    sf_domain <- inherits(domain, "sf")
    sp_samplers <- inherits(samplers, "Spatial")
    sp_domain <- inherits(domain, "Spatial")
    if (sf_domain && sp_samplers) {
      samplers <- sf::st_as_sf(samplers)
    }
  }

  # How does sf deal with secondary geometry columns?
  # TODO have to deal with the multidomain samplers
  # https://r-spatial.github.io/sf/articles/sf6.html
  # TODO save to log and verbosity = 2, use seq_along to introduce the indices
  # for each active geometry
  # TODO extra domains we assign the sampler for them
  # TODO can use local helper function (with S3 usemethod() maybe)
  sfsp <- function(x, ...) {
    UseMethod("sfsp")
  }
  sfsp.list <- function(x, begin = NULL){
    bru_log_message(
      paste0(
        begin,
        x[i],
        " are ",
        attr(x[i], "sf_column"),
        ".\n"
      ),
      # TODO to global bru_option_get verbose=logical; if TRUE, print the log message on screen with message(txt). Default: bru_options_get("bru_verbose")
      verbose_store = options$bru_verbose_store,
      verbosity = 2
    )
  }
  sfsp.data.frame


  if (sf_samplers) {
    for (i in seq_along(samplers)) {
      bru_log_message(
        paste0(
          "The active geometry of the sampler",
          samplers[i],
          " is ",
          attr(samplers[i], "sf_column"),
          ".\n"
        ),
        verbose = options$bru_verbose, # TODO to global bru_option_get
        verbose_store = options$bru_verbose_store,
        verbosity = 2
      )
    }
  }

  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?
  if (unique(names(samplers)) != unique(names(domain))) {
    # TODO have to accommodate weights column in samplers as well, omit this column
    samplers_domain <- intersect(names(samplers), names(domain))
    bru_log_message(
      paste0(
        "The shared samplers and domain: \n",
        paste0(samplers_domain, collapse = ", ") # TODO for the NULL case
      ),
      verbose = options$bru_verbose,
      verbose_store = options$bru_verbose_store,
      verbosity = 2
    )

    # Domain should be more than samplers.
    # However, we can have unused samplers as well.
    if (length(samplers) > length(domain)) {
      unused_samplers <- setdiff(names(samplers), names(domain))
      bru_log_message(
        paste0(
          "The unused samplers: \n",
          paste0(unused_samplers, collapse = ", ")
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 2
      )
    } else {
      extra_domain <- setdiff(names(samplers), names(domain))
      bru_log_message(
        paste0(
          "\n  The extra domain: \n",
          paste0(extra_domain, collapse = ", ")
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 2
      )
    }
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
