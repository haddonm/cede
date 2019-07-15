#' @title cede a set of functions to assist working with Fisheries data
#'
#' @description The cede package provides a set of functions to
#'     facilitate the exploration and sketch mapping of fisheries data. Given
#'     a data.frame containing such data histograms of catch, effort, Log(CPUE),
#'     depth, and others by year, or other grouping factor can be plotted for
#'     vsual inspection and summary. Sketch maps of data containinf Lat and Long
#'     data can also be made. Finaly CPUE standardization routines are included
#'     to illustrate alternatives for such analyses. Vignettes (see
#'     browseVignettes("cede") for details) are provided to illustrate each set
#'     of functions. Simulated fisheries data are provided in the package to
#'     simplify the illustrattion of examples.
#'
#' @details This package includes mapping functions, plotting functions, and
#'     functions to facilitate the standardization of CPUE data. It includes
#'     data sets for the coast of Australia and example fisheries data sets.
#' @section Mapping functions:
#' \describe{
#'   \item{addpoints}{literally adds a set of (long,lat) points to a map}
#'   \item{maps}{lists the available mapping functions and their syntax}
#'   \item{plotLand}{add the land mass polygons after fitting data to make
#'       for a tidy map by hiding stray points on land; they are only hidden
#'       you still need to be aware of their presence!}
#'   \item{plotpolys}{plots catch total by a grid over the area plotted}
#'   \item{plotaus}{plots up a schematic map outline of Australia. The only
#'       part visible is defined by leftlong, rightlong, uplat, and downlat}
#' }
#' @section Data Sets:
#' \describe{
#'   \item{aus}{a data frame of 218012 rows x 3 columns depicting
#'       12 different polygons making up a rough map of Australia.
#'       Not exported try head(cede:::aus)}
#'   \item{sps}{a data.frame of 11603 x 10 providing an example dataset. It
#'       contains columns of Year, Month, Vessel (a code), catch_kg, Long, Lat,
#'       Depth, DayNight, Effort, and Zone. While the Lat Long values are
#'       scattered along the west coast of Tasmania this is not data from a
#'       real fishery but is simulated to be like a real fishery. It is for use
#'       when testing and explaining the use of this package. }
#' }
#' @docType package
#' @name cede
NULL

#' @title sps - a 11603 x 10 data.frame for testing cpue and mapping functions
#'
#' @description sps - a 11603 x 10 data.frame for testing CPUE functions
#'    containing simulated trawl shot CPUE data for the years 2003 - 2014,
#'    including details of year, month, day, vessel, catch, longitude,
#'    latitude, depth trawled in meters, a daynight identifier, and effort.
#'    The cpue and log(cpue) can be calculated from this data.
#'
#' @format A data frame with 10603 x 10 variables:
#' \itemize{
#'   \item{Year}{The year in which fishing takes place}
#'   \item{Month}{The month in which fishing took place}
#'   \item{Vessel}{a code uniquely identifying vessels through time}
#'   \item{catch_kg}{literally the catch of the species in the shot in kg}
#'   \item{Long}{the longitude of the start of the trawl shot}
#'   \item{Lat}{the latitude of the start of the trawl shot}
#'   \item{Depth}{the average depth of trawling in meters}
#'   \item{DayNight}{a code denoting the daynight status D = day, N = night,
#'      M = mixed}
#'   \item{Effort}{the hours trawled}
#'   \item{Zone}{The fished area split into three latitudinal zones}
#' }
#' @docType data
#' @name sps
NULL


#' @import graphics
#' @importFrom stats median qnorm qqline qqnorm quantile sd as.formula dnorm lm
#' @importFrom grDevices dev.cur dev.new dev.off png rgb
#' @importFrom utils tail

NULL
