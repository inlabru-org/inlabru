#' https://github.com/inlabru-org/inlabru/issues/125

# TODO 20220126should not be a S3 but local function and function name should be clearer
# bru_log_list function ---------------------------------------------------------
# How does sf deal with secondary geometry columns?
# TODO ####
# TODO have to deal with the multidomain samplers
# https://r-spatial.github.io/sf/articles/sf6.html
# TODO #### extra domains we assign the sampler for them
# TODO can use local helper function (with S3 Usemethod() maybe)

# @param start_with The message starts with
# @param end_with The message ends with
# @verbosity Default: 2
# @ verbose_store Default: T
bru_log_list <- function(x, attr_names = "sf_column", start_with = NULL, end_with = NULL) {
  for (i in seq_along(x)) {
    bru_log_message(
      paste0(
        start_with,
        i,
        " is/are ",
        attr(x[[i]], attr_names),
        end_with,
        ".\n"
      ),
      verbosity = 2, verbose_store = T
    )
  }
}

# devtools::load_all()
# data(mrsea, package = "inlabru")
# domain = list(coordinates = mrsea$mesh,
#               season = seq_len(4))
# samplers <- mrsea$samplers
# samplers <- list(season = samplers$season, samplers)
# debugonce(inlabru:::apmaker)
# apmaker(samplers=samplers, domain=domain, response="coordinate")

# TODO Extend S3 method with names_list???
# https://stackoverflow.com/questions/18513607/how-to-extend-s3-method-from-another-package-without-loading-the-package
# to deal with samplers names in a dataframe 20230126
# 1) do I need the names of the list or the names within the list(s)?
# 2) Does column names mean sth here?
# 3) What to do with NULL cases, aka unamed list?
#' @name names
#' @export names_list
#' @method extract names within list(s)
#' @title Names within list(s)
#'
names_list <- function(x) {
  lapply(x, function(y) {
    names(y)
  })
}


#' @aliases apmaker
#' @export
#' @param domain A list of named integration definitions, each either
#' character/factor vector, a numeric vector of points given integration
#' weight 1, an `inla.mesh.1d` object, or an `inla.mesh.2d` object. Only those
#' domains that are not given in the `samplers` data.frame are used, plus the
#' coordinates object, used for the spatial aspect of the `samplers` object.
#' TODO * 20221212 For extension, it can be factor or character while keeping
#' the level to map them back to the components
#' @param samplers A (list of unnamed or named element(s) of ) `[sf/sfc]DataFrame` or
#' `Spatial[Points/Lines/Polygons]DataFrame` object(s). Unnamed elements
#' are assumed to be multidomain samperls; named elements are singledomain
#' samplers; domains without corresponding samplers are assumed to be full domain
#' samplers.
#' TODO is response useful here? 20220130
#' @param response Name(s) of the responses. Deprecated: dnames
#' @param weight The name of integration weight column in the samplers. Default: `weight`
#' TODO 1) how about domain? should not be allowed.
#' 2) allow_names should be bru_weight? 23112022 Only for samplers
#' Default: "weight".
#' It mainly works for line transect at the moment, determining the width of the
#' line to be integrated. See the distance sampling example
#' https://inlabru-org.github.io/inlabru/articles/web/2d_lgcp_distancesampling.html
#' @param int.args List of arguments passed on to \code{ipoints}
#' @return Integration points

# TODO option argument as in bru function with list() bru_int_args
apmaker <- function(domain = NULL, samplers = NULL, response = NULL,
                    weight = "weight",
                    int.args = list(method = "stable", nsub = NULL)) {
  # TODO ####
  # TODO To allow sf geometry support, should likely change the logic to
  # use the domain specification to determine the type of integration
  # method to call, so that it doesn't need to rely on the domain name.
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
  # TODO 20221213 leave it for now for backward compatibility
  if (is.null(domain)) {
    if (is.null(samplers)) {
      stop("Domain and samplers argument(s) NULL.")
    }
    warning("Domain argument(s) NULL. Sampler is used to create the domain")
  }

  if (!inherits(domain, "list") ||
    is.null(names(domain))) {
    stop("Domain must be a named list.")
  }

  if (!inherits(samplers, "list")) {
    samplers <- list(samplers)
  }

  #########################################################################################
  # 20221208 sf sfc Spatial should match inla.mesh, inla.mesh.lattice, raster, SpatRaster
  # character/factor/numeric should match the same classes of domain.
  # numeric is implemented with weight 1 in ipoints.
  # TODO character and factor to be implemented similarly but not changed into numeric.

  # Mandate domain to be the following classes
  # Unused domain is not allowed
  # 20221213 not needed, up to another function to deal with it
  # stopifnot_class(
  #   domain,
  #   c( # "data.frame", # maybe in the future
  #     "character",
  #     "factor",
  #     "numeric", # here is the split
  #     "inla.mesh",
  #     "inla.mesh.1d",
  #     "inla.mesh.lattice",
  #     "raster",
  #     "SpatRaster"
  #   )
  # )
  # Mandate domain to be the following classes
  # stopifnot_class(
  #   samplers,
  #   c(
  #     "character", "factor", "numeric", # here is the split
  #     "data.frame", # 20221212 ipoints allows data.frame samplers
  #     "sf", "sfc", "Spatial"
  #   )
  # )
  #####################################################################################
  # 20221213 This function to figure out which samplers the domain integration are
  # class check is not needed here but to be done in S3 method
  # 20221212 We have to check the matching here, that is
  # 1) certain domain classes have to match certain samplers classes.
  # TODO 2) in which 1) has to check the domain one by one if it is a list
  # if (any(unlist(lapply(domain, function(x) {
  #   inherits(x, c("character", "factor", "numeric", "inla.mesh.1d"))
  # })))) {
  #   stopifnot_class(samplers, c("character", "factor", "numeric", "data.frame"))
  # }
  # if (any(unlist(lapply(domain, function(x) {
  #   inherits(x, cc(
  #     "inla.mesh",
  #     "inla.mesh.lattice",
  #     "raster",
  #     "SpatRaster"
  #   ))
  # })))) {
  #   stopifnot_class(samplers, c("sf", "sfc", "Spatial"))
  # }

  # Change a mix of sp and sf objects to sf
  sf_samplers <- unlist(lapply(samplers, function(x) inherits(x, c("sf", "sfc"))))
  sp_samplers <- unlist(lapply(samplers, function(x) inherits(x, "Spatial")))
  if (any(sp_samplers) && any(sf_samplers)) {
    warning("Both `sf` and `sp` objects in the samplers are detected. Produce `sf` output as desired")
    samplers[sp_samplers] <- lapply(samplers[sp_samplers], sf::st_as_sf)
  }

  # TODO 20220126 lapply to extract the names, the current one is not sufficient
  # TODO Sort multidomain samplers, single domain samplers
  # TODO Multidomain samplers happens when a sampler across several domains. How to detect that?
  # TODO we should then do the name check for domain and samplers here
  # TODO remove sampler domains and full domain samplers
  # TODO some thoughts for S3 methods, there should be an extra layer ie function to sort samplers and domain arguments and difine multidomain, singledomain and full domain s3 class
  #######################
  names_domain <- names(domain)
  names_lsamplers <- names(samplers) # so that we know which are named and unnamed
  names_samplers <- unique(unlist(names_list(samplers)))
  names_response <- names(response) # from the formula
  names_reserved <- c(weight) # coordinate and geometry is not required here

  #######################
  # multidomain samplers, ie unnamed element(s) in samplers, for each sampler and then for each domain(lapply)
  # TODO still have to deal with secondary geometry
  if (any("" %in% names_lsamplers)) {
    multidomainsamplers <-
      samplers[unlist(lapply(list(names_lsamplers), function(x) {
        x == ""
      }))]
    lnames <- names_list(multidomainsamplers) # retain the attr(*, "names") while keeping the list structure
    lips_multidomainsamplers <- list()
    for (i in seq_along(multidomainsamplers)) {
      if (!is.null(lnames[[i]])) {
        names_intersect <- setdiff(
          intersect(lnames[[i]], names_domain),
          names_reserved
        )
        # compute ips for each domain with group(block)=names_intersect
        # lapply(list(multidomainsamplers[i], domain, names_intersect), ipoints, )
        lips_multidomainsamplers[i] <- lapply(names_domain, function(nm) ipoints(
          samplers = multidomainsamplers[[i]],
          domain = domain[nm],
          name = nm,
          group = names_intersect, # block=group should be the grouping, say season,
          int.args = int.args
        ))
      } else {
        warning(paste0("The unnamed list", lnames[i], "in the samplers is NULL"))
        lips_multidomainsamplers[i] <- NULL
      }
    }
    # TODO 20230206 cprod step
    lips_multidomainsamplers <- do.call(cprod, lips_multidomainsamplers)
  }

  #######################
  # singledomain samplers, ie named element(s) in samplers
  singledomainsamplers <-
    samplers[unlist(lapply(list(names_lsamplers), function(x) {
      x != ""
    }))]
  if (!is.null(singledomainsamplers)){
    lips_singledomainsamplers <- list()
    for (i in seq_along(singledomainsamplers)) {
      names_intersect <- setdiff(intersect(names(singledomainsamplers),
                                           names_domain), names_reserved)
      lips_singledomainsamplers[i] <- lapply(names_domain, function(nm) ipoints(
        samplers = singledomainsamplers[[i]],
        domain = domain[nm],
        name = nm,
        group = names_intersect, # block=group should be the grouping, say season,
        int.args = int.args
      ))
    }
  }

  # remove sampler domains, then pass them to full domain samplers
  # full domain samplers, i.e. domain with missing samplers
  names_fulldomain <- setdiff(names_domain, c(names_samplers, names_reserved))
  if (!is.null(names_fulldomain)){
    lips_fulldomainsamplers <- list()

  }

  # TODO cprod lips and ips ...
  ips <- do.call(cprod, c(list(ips), lips))


  #############################################################################
  # Warn samplers without domain associated
  # Store the names in domain but not in samplers
  # Warn domains without associated samplers


  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?

  # log the active geometry of the samplers
  # TODO 20220126 when sf_samplers is a vector, if clause does not work
  # 20220130 I think it works now


  # TODO ####
  # TODO check ibm_values.bru_mapper_factor for character/factor/numeric samplers
  # TODO weight argument to take effect on samplers. The weight should go to
  # the integration part

  ##########################################################
}
