% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpue.r
\name{turnover}
\alias{turnover}
\title{turnover Estimate turnover of vessels from catch by vessel by year data}
\usage{
turnover(x, minyrs = 1)
}
\arguments{
\item{x}{A matrix of a continuous numeric property by year,
the original usage was to plot catch-by-vessel against year}

\item{minyrs}{limits the analysis to those vessels that remain in the
fishery for at least minyrs years - which would eliminate the occasional
opportunistic fisher who only fishes for one or two years, or whatever
minimum is selected. Vessels with zero catches are not included in case
zeros and NAs are counted as starting and leaving the fishery.}
}
\value{
a matrix of years by Continue, Leave, Start, Total
}
\description{
Estimate turnover of vessels from catch by vessel by year data;
  To specify the minimum number of years that a vessel needs stay in the fishery,
  then give a value to the variable minyrs.
}
\examples{
\dontrun{
library(r4cpue)
data(sps)
cbv <- tapply(sps$catch_kg,list(sps$Vessel,sps$Year),sum,na.rm=TRUE)/1000
dim(cbv)
early <- rowSums(cbv[,1:6],na.rm=TRUE)
late <- rowSums(cbv[,7:14],na.rm=TRUE)
cbv1 <- cbv[order(late,-early),]
plotprep(width=7,height=6)
yearBubble(cbv1,ylabel="Catch by Trawl",vline=2006.5,diam=0.2)
turnover(cbv)
}
}
