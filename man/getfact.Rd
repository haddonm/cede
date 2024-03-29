% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpue.r
\name{getfact}
\alias{getfact}
\title{getfact extracts parameter estimates for a given factor from a CEout object}
\usage{
getfact(inmat, invar)
}
\arguments{
\item{inmat}{generally this will be a CEout object but it can be
a matrix of coefficients, an lm object, or a gam object}

\item{invar}{the model variable whose parameters are wanted}
}
\value{
a matrix containing the parameters for invar
}
\description{
getfact extracts the parameter estimates for a given factor from a CEout object.
   It does this by searching to rownames of the output parameters of the optimum model.
   It also checks for interaction terms, which for categorical factors is the same as
   determining a trend of the two factors relative to each other. e.g. for Zone:Month
   the outcome is the monthly trend for each zone.
}
