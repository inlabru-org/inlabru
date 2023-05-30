#' @name robins_subset
#' @title robins_subset
#' @docType data
#' @description This is the `robins_subset` dataset, which is a subset of the
#'   full robins data set used to demonstrate a spatially varying trend coefficient
#'   model in Meehan et al. 2019. The dataset includes American Robin counts,
#'   along with time, location, and effort information, from Audubon Christimas Bird
#'   Counts (CBC) conducted in six US states between 1987 and 2016.
#'
#' @usage data(robins_subset)
#'
#' @format The data are a data.frame with variables
#'  \describe{
#'    \item{`circle`:}{
#'    Four-letter code of the CBC circle.
#'    }
#'    \item{`bcr`:}{
#'    Numeric code for the bird conservation region encompassing the
#'    count circle.
#'    }
#'    \item{`state`:}{
#'    US state encompassing the count circle.
#'    }
#'    \item{`year`:}{
#'    calendar year the count was conducted.
#'    }
#'    \item{`std_yr`:}{
#'    transformed year, with 2016 = 0.
#'    }
#'    \item{`count`:}{
#'    number of robins recorded.
#'    }
#'    \item{`log_hrs`:}{
#'    the natural log of party hours.
#'    }
#'    \item{`lon`:}{
#'    longitude of the count circle centroid.
#'    }
#'    \item{`lat`:}{
#'    latitude of the count circle centroid.
#'    }
#'    \item{`obs`:}{
#'    unique record identifier.
#'    }
#'  }
#' @source
#' https://github.com/tmeeha/inlaSVCBC
#'
#'
#' @references
#' Meehan, T.D., Michel, N.L., and Rue, H. 2019. Spatial modeling of Audubon
#' Christmas Bird Counts reveals fine-scale patterns and drivers of relative
#' abundance trends. Ecosphere, 10(4), p.e02707.
#'
#' @examples
#' if (require(ggplot2, quietly = TRUE)) {
#'   data(robins_subset, package = "inlabru") # get the data
#'
#'   # plot the counts for one year of data
#'   ggplot(robins_subset[robins_subset$std_yr == 0, ]) +
#'     geom_point(aes(lon, lat, colour = count + 1)) +
#'     scale_colour_gradient(low = "blue", high = "red", trans = "log")
#' }
NULL
