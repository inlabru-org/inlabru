#' A S3/S4 generic function to generate stop warning if the object of some class
#'  is not of the specified classes.
#' @param  object An object except NULL
#' @param class A vector of class(es) except NULL
#' @keywords internal
#'
# stopifnot_class <- function(object = NULL, class = NULL, ...) {
#   # 20221208 Probably better to leave out NULL and missing class in this
#   # function, though it deals with it
#   class_name <- paste0(sort(class), collapse = ",")
#   obj_name <- paste0(as.character(expression(object)))
#   UseMethod("stopifnot_class")
# }
# stopifnot_class.default <- function(object, class) {
#   if (!inherits(object, class)) {
#     stop(paste0(
#       obj_name,
#       " of class ",
#       paste0(class(object), collapse = ","),
#       "must be of class(es) ",
#       class_name
#     ))
#   }
# }
# stopifnot_class.list <- function(object, class) {
#   if (!all(unlist(lapply(object, function(x) {
#     inherits(x, class)
#   })))) {
#     stop(paste0(
#       obj_name,
#       " of class ",
#       paste0(unlist(lapply(object, class)), collapse = ","), # TODO more complicated
#       "must be of class(es) ",
#       class_name
#     ))
#   }
# }

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
#' @param samplers A (list of) `[sf/sfc]DataFrame` or
#' `Spatial[Points/Lines/Polygons]DataFrame` object(s)
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
  # To allow sf geometry support, should likely change the logic to
  # use the domain specification to determine the type of integration
  # method to call, so that it doesn't need to rely on the domain name.
  # TODO ####
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
  # For multidomain sampler, the main thing is data aka dname 23112022
  # specific column has that meaning, how to find the response?
  #######################
  # test case
  # a <- data.frame(hi=1:4, hi2=9:10)
  # b <- data.frame(bye=5:8)
  # c <- data.frame(hey=6:9)
  # g <- data.frame(ciao=7:10)
  # e <- data.frame(hej=6:8)
  # domain <- list(a,b,c,e,g) # domain names are unique
  # samp1 <- data.frame(a)
  # samp2 <- data.frame(b,c)
  # samplers <- list(samp1, samp2)
  # samplers <- list(a,list(b,e)) # this is not the case
  # response <- list(g)
  # # names(domain) <- c("a", "b", "c", "b", "a", "a", "e", "c")
  # # names(samplers) <- c("a", "be")
  # # names(response) <- c("g")
  # weight <- "weight"
  #######################

  names_domain <- names(domain)
  names_lsamplers <- names(samplers) # so that we know which are named and unnamed
  names_response <- names(response) # from the formula
  names_reserved <- c(weight) # coordinate and geometry is not required here

  #######################
  # multidomain samplers, ie unnamed element(s) in samplers
  if (any("" %in% names_lsamplers)) {
    multidomainsamplers <-
      samplers[unlist(lapply(list(names_lsamplers), function(x) {
        x == ""
      }))]
    lnames <- names_list(multidomainsamplers)   # retain the attr(*, "names") while keeping the list structure
    if (!is.null(lnames)) {
      lips_multidomainsamplers <- list()
      for (i in seq_along(multidomainsamplers)){
        names_intersect <- setdiff(
          intersect(unlist(lnames[i]), names_domain),
          names_reserved
        )
        # compute ips for each domain with group(block)=names_intersect
        lapply(list(multidomainsamplers[i],domain,names_intersect), ipoints, )
        lips_multidomainsamplers[i] <- ipoints(
          multidomainsamplers,
          domain,
          group = names_intersect, # TODO the block=group isn't sure what to put in, actually what does group mean? Column names of the samplers object (if applicable) for which the integration points are calculated independently and not merged when aggregating to mesh nodes.
          int.args = int.args
        )
      }
    } else {
      warning("The unnamed list in the samplers is NULL")
    }
  }

  #######################
  # singledomain samplers, ie named element(s) in samplers
  singledomain <-
    samplers[unlist(lapply(list(names_lsamplers), function(x) {
      x != ""
    }))]

  #######################
  # remove sampler domains, then pass them to full domain samplers

  #######################
  # full domain samplers, i.e. domain with missing samplers
  names_fulldomain <- setdiff(names_domain, c(names_samplers, names_reserved))
  if (!is.null(names_fulldomain)) {
    lips_fulldomain <- list()
    for (i in seq_along(fulldomain)){
      lips_fulldomain[i] <- ipoints(
        domain = fulldomain[names_intersect[i]],
        group = , # TODO the block=group isn't sure what to put in, actually what does group mean? Column names of the samplers object (if applicable) for which the integration points are calculated independently and not merged when aggregating to mesh nodes.
        int.args = int.args
      )
    }
  }


  #############################################################################
  # Warn samplers without domain associated
  if (!is.null(names_setdiff %in% names_samplers)) {
    bru_log_list(as.list(names_intersect),
      attr_names = "names",
      start_with = "\n The unused samplers: \n"
    )
    warning(paste0("Sampler(s) without associated domain(s): ",
      paste0(names_intersect),
      collapse = ","
    ))
  }
  # Store the names in domain but not in samplers
  # Warn domains without associated samplers
  if (!is.null(names_setdiff %in% names_domain)) {
    extra_domain <- intersect(names_setdiff, names_domains)
    bru_log_list(as.list(extra_domain),
      attr_names = "names",
      start_with = "\n The extra domain: \n"
    )
    warning(paste0("Create sampler(s) for domain(s) without associated sampler(s): ",
      paste0(names_intersect),
      collapse = ","
    ))
  }

  # Check if the names of samplers and domains match. How to establish the link
  # between samplers and domain? If names are not provided, follow the order in
  # list. If names are provided but do not match, what should we do?
  if (unique(names_samplers) != unique(names_domain)) {
    # TODO for the NULL name case
    bru_log_list(samplers_domain, start_with = "\n The shared samplers and domain: \n")
  }

  # log the active geometry of the samplers
  # TODO 20220126 when sf_samplers is a vector, if clause does not work
  # 20220130 I think it works now
  if (any(sf_samplers)) {
    bru_log_list(samplers[sf_samplers],
      start_with = "The active geometry of the sampler(s) "
    )
  }

  # TODO ####
  # TODO check ibm_values.bru_mapper_factor for character/factor/numeric samplers
  # TODO weight argument to take effect on samplers. The weight should go to
  # the integration part
  ips <- apoints(samplers, domain,
    int.args = int.args, weight_name = weight_name
  )
  # TODO using cprod for each domain
  # this is just a draft atm
  lips <- lapply(nosamp.dim, function(nm) {
    ipoints(NULL, domain[[nm]], name = nm, int.args = int.args)
  })
  ips <- do.call(cprod, c(list(ips), lips))
  ##########################################################
}
