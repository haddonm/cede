% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ausmaps.R
\name{plotLand}
\alias{plotLand}
\title{plotLand - puts land mass on top of a plot; covers land overlaps}
\usage{
plotLand(incol = "lightgrey", RGB = FALSE, green = NA, blue = NA, alpha = NA)
}
\arguments{
\item{incol}{defaults to "lightgrey" but can be any single number representing
a colour.}

\item{RGB}{should we use an rgb colour; default = FALSE}

\item{green}{the green intensity default = NA}

\item{blue}{the blue intensity}

\item{alpha}{the standard density 1 = full colour 1/10 = one tenth}
}
\value{
Nothing returned but the plot of land is repeated in light blue
}
\description{
plotLand - puts land mass on top of a plot; covers land overlaps
  Used primarily after using plotpolys, which often leaves some oblongs
  sitting over parts of the coast. Using addpoints can also leave points
  on land. However, if you want to illustrate things on land then now you
  can use RGB = TRUE, and if you set alpha to 0.1 this will be transluscent
  andpoints in rivers, etc, will still be visible. If you use this option
  then incol becomes the red intensity and must lie between 0 and 1.
}
\examples{
\dontrun{
  dev.new(height=6.0,width=7.5,noRStudioGD = TRUE)
  plotaus()
  plotLand(incol="red")
  plotLand(incol=1,RGB=TRUE,green=0,blue=0,alpha=(1/10))
}
}
