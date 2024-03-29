% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ausmaps.R
\name{addpoints}
\alias{addpoints}
\title{addpoints includes a set of points with a set size and symbol}
\usage{
addpoints(
  indat,
  inpch = 20,
  incex = 0.2,
  incol = 2,
  intitle = "",
  jit = F,
  wobble = 1,
  refill = T,
  Long = "Long",
  Lat = "Lat",
  txtout = TRUE,
  namecatch = "catch_kg"
)
}
\arguments{
\item{indat}{a matrix or data.frame containing, at least, columns with the
names 'Long' nad 'Lat'.}

\item{inpch}{default = 20, but can be changed to any valid symbol}

\item{incex}{default 0.2, but can be a constant or variable}

\item{incol}{default 2 (red) but can be any colour or rgb definition}

\item{intitle}{main top center title, the same as in plotaus; this is
mostly useful if the refill variable is set to TRUE}

\item{jit}{default FALSE. If TRUE it will jitter the outcome by the
'wobble' valriable amount. However, for overlapping points it is often
better to use an rgb colour with the last value set to a fraction
e.g. incol=rgb(9)0,0,0,1/5)}

\item{wobble}{the amount to jitter a point id jit = TRUE}

\item{refill}{default TRUE, refill the land after plotting points.}

\item{Long}{used to define the longitude field; defaults to 'Long'}

\item{Lat}{used to define the Latitude field; default to 'Lat'}

\item{txtout}{logical determines whether the number of points added and \
their range will be printed to the screen; defaults to TRUE}

\item{namecatch}{is the name of the catch variable; defaults to 'catch_kg'}
}
\value{
adds the input points to a map; returns nothing.
}
\description{
addpoints given a matrix or data.frame, which must contain
  fields names 'Lat' and 'Long', these value pairs will be added to the plot
  generated by plotaus
}
\examples{
dev.new(height=6.0,width=7.5,noRStudioGD = TRUE)
plotaus()
Long <- c(140,141,142,143)
Lat <- c(-41.1,-41.1,-41.1,-41.1)
indata <- cbind(Long,Lat)
addpoints(indata,incex=1.5)
}
