% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpue.r
\name{quants}
\alias{quants}
\title{quants used in apply to count the number > 1 in a vector}
\usage{
quants(invect)
}
\arguments{
\item{invect}{vector of values}
}
\value{
a vector of the c(0.025,0.05,0.5,0.95,0.975) quantiles
}
\description{
quants used in apply to count the number > 1 in a vector
   designed to be used in apply
}
\examples{
x <- rnorm(100,mean=5,sd=1)
quants(x)
y <- matrix(x,nrow=10,ncol=10)
apply(y,2,quants)
}
