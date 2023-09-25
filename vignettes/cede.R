## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	message = FALSE,
	warning = FALSE)
options(knitr.kable.NA = "",
        knitr.table.format = "pandoc")

options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)

library(cede)
suppressPackageStartupMessages(library(knitr))


## ----spscontents, echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data(sps)
head(sps)
properties(sps)

## ----echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# tapply(sps$catch_kg,list(sps$Year,sps$Month),sum,na.rm=TRUE)
round(tapsum(sps,"catch_kg","Year","Month",div=1000),2)

## ----echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tapsum(sps,"catch_kg","Year","Zone",div=1000)

## ----echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
effbyyr <- tapsum(sps,"Effort","Year","Zone",div=1.0)
effbyyr

## ----ploteffort, echo=TRUE, fig.width=7,fig.height=4.5,fig.align="center"---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plotprep(width=7,height=4.5)
ymax <- max(effbyyr,na.rm=TRUE)
label <- colnames(effbyyr)
yrs <- as.numeric(rownames(effbyyr))
par(mfrow=c(1,1),mai=c(0.45,0.45,0.05,0.05)) 
par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7) 
plot(yrs,effbyyr[,label[1]],type="l",lwd=2,col=1,ylim=c(0,ymax),
     ylab="Total Effort (Hours) by Zone per Year",xlab="",
     panel.first=grid())
lines(yrs,effbyyr[,label[2]],lwd=2,col=2)
lines(yrs,effbyyr[,label[3]],lwd=2,col=3)
# the "topright"" can also be "bottomright", "bottom", "bottomleft", "left", 
# # "topleft", "top", "right" and "center". Alternatively one can give an
# x and a y location, see ?legend, which is a function from graphics package
legend("topright",label,col=c(1,2,3),lwd=3,bty="n",cex=1.25)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 # the div argument is redundant because 1000 is the default value, kg -> tonne
tapsum(sps,"catch_kg","Year","DayNight",div=1000)

## ----catchbyvessel, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cbv <- tapsum(sps,"catch_kg","Vessel","Year") # often more vessels than years
total <- rowSums(cbv,na.rm=TRUE)
cbv1 <- cbv[order(total),]   # sort by total catch smallest catches first
round(cbv1,2)


## ----yeasrbubble, fig.width=7,fig.height=5,fig.align="center"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plotprep(width=8,height=6) # not needed in the vignette
to <- turnover(cbv1)
yearBubble(cbv1,ylabel="sqrt(catch-per-vessel)",diam=0.125,txt=c(2,3,4,5),
           hline=TRUE)


## ----echo = TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(to)

## ----includeCPUE, echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sps$CE <- NA     # make explicit space in the data.frame
sps$LnCE <- NA
pick <- which((sps$catch_kg > 0) & (sps$Effort > 0))
sps$CE[pick] <- sps$catch_kg[pick]/sps$Effort[pick]
sps$LnCE[pick] <- log(sps$CE[pick])   # natural log-transformation
# categorize Depth
range(sps$Depth,na.tm=TRUE)  # to aid selection of depth class width
sps$DepCat <- NA
sps$DepCat <- trunc(sps$Depth/25) * 25
table(sps$DepCat)


