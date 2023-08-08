# https://github.com/inlabru-org/inlabru/issues/125

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
# @name names
# @export names_list
# @title Names within list(s)
#
names_list <- function(x) {
  lapply(x, function(y) {
    names(y)
  })
}
