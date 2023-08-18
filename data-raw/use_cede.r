
library(cede)
options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)


# Data exploration
data(sps)
properties(sps)
pickcol <- c(3,4,9,11,27,33,35,40,41,47,48,52,53)
sps <- dat[,pickcol]
properties(sps)
dim(sp)
newcol <- c("Lat","Long","grid","Csq","FishID","Fishery","Method",
            "Year","Month","Effort","Effunit","catch",
            "whwt")
colnames(sps) <- newcol
properties(sps)
# No variables have NA values, teh catch_kg and effort variabales are all > 0
pickcol <- c(1,2,9,10,12,13)
nvar <- length(pickcol)
for (i in 1:nvar) {
sps[,pickcol[i]] <- as.numeric(sps[,pickcol[i]])
}
properties(sps)
yearNA(sps)
pick <- which((sps$Fishery == "Ocean Trap and Line"))
length(pick)
sps1 <- droplevels(sps[pick,])
pick <- which((sps$Year == "2008/09   "))
length(pick)
sps1 <- droplevels(sps1[-pick,])


properties(sps1)

yearNA(sps1)



tapsum(sps1,"catch","Year","Month")
tapsum(sps1,"catch","Year","Fishery")
table(sps1$Year,sps1$Effunit)

sps1$CE <- NA 
properties(sps1)
pickC <- which(sps1$catch > 0)
sps1$CE[pickC] <- sps1$catch[pickC]
sps1$LnCE <- NA
sps1$LnCE[pickC] <- log(sps1$CE[pickC])

pickMeth <- which(sps1$Effunit %in% c("Number of hooks","Number of lures"))
length(pickMeth)
sps1 <- droplevels(sps1[pickMeth,])



#  Hook Method ---------
pickM <- which(sps1$Effunit == "Number of hooks")
sps2 <- droplevels(sps1[pickM,])
range(sps2$LnCE,na.rm=TRUE)
plotprep(width=9,height=7)
outH <- histyear(sps2,Lbound=-1.25,Rbound = 8.75,inc=0.250,
                 pickvar="LnCE",varlabel="Log(CE) HookNUM",
                 plots=c(3,3),vline=4)

plotprep(width=6,height=8)
leftlong <- 149;  rightlong <- 155 
uplat <- -27;  downlat <- -39 
plotaus(leftlong,rightlong,uplat,downlat) 
addpoints(sps2,incex=0.75)
addpoints(sps3,incol="blue",incex=0.75)


#  Lure Method ------------
pickM <- which(sps1$Effunit == "Number of lures")
sps3 <- droplevels(sps1[pickM,])
range(sps3$LnCE,na.rm=TRUE)
plotprep(width=9,height=7)
outH <- histyear(sps3,Lbound=-1.25,Rbound = 8.75,inc=0.250,
                 pickvar="LnCE",varlabel="Log(CE) LureNUM",
                 plots=c(3,3),vline=4)

tmp <- table(sps1$FishID,sps1$Effunit)
tmp[order(tmp[,1]),]


cbv <- tapsum(sps1,"catch","FishID","Year")
total <- rowSums(cbv,na.rm=TRUE)
sorttot <- sort(total)
ans <- cumsum(sorttot)
plot(ans)


pick <- which(total > 2)
cbv1 <- cbv[pick,]
cbv2 <- cbv1[order(total[pick]),]
dim(cbv2)
head(cbv2)
colnames(cbv2) <- 2009:2017
plotprep(width=8,height=5)
yearBubble(cbv2,ylabel="sqrt(catch-per-vessel)",diam=0.25,
           txt=c(8, 16,18,20),
           hline=TRUE)
to <- turnover(cbv2)


# Add CPUE
sps$CE <- NA
sps$LnCE <- NA
pick <- which((sps$catch_kg > 0) & (sps$Effort > 0))
sps$CE[pick] <- sps$catch_kg[pick]/sps$Effort[pick]
sps$LnCE[pick] <- log(sps$CE)
# categorize Depth
sps$DepCat <- NA
sps$DepCat <- trunc(sps$Depth/25) * 25
table(sps$DepCat)

range(sps$LnCE,na.rm=TRUE)

plotprep(width=8,height=5)
histyear(sps,Lbound=-1.75,Rbound=8.5,inc=0.25,pickvar="LnCE",years="Year",
         varlabel="log(CPUE",plots=c(3,4))

range(sps$Depth)
plotprep(width=8,height=5)
outH <- histyear(sps,Lbound=0,Rbound=375,inc=12.5,pickvar="Depth",years="Year",
                 varlabel="Depth (m)",plots=c(3,4),vline=120)

range(sps$Effort)
plotprep(width=8,height=5)
outH <- histyear(sps,Lbound=0,Rbound=10,inc=0.25,pickvar="Effort",years="Year",
                 varlabel="Effort (Hrs)",plots=c(3,4),vline=NA)

plotprep(width=8,height=5)
plot(sps$Effort,sps$catch_kg,type="p",pch=16,col=rgb(1,0,0,1/5),ylim=c(0,1500))
abline(h=0.0,col="grey")



plotprep(width=7,height=5.5)
leftlong <- 143;  rightlong <- 150
uplat <- -40;  downlat <- -44.6
plotaus(leftlong,rightlong,uplat,downlat,gridon=1.0)
addpoints(sps)
plotLand(incol="blue")

plotprep(width=7,height=5.5)
leftlong <- 143;  rightlong <- 150
uplat <- -40;  downlat <- -44.6
plotaus(leftlong,rightlong,uplat,downlat,gridon=1.0)
plotpolys(sps,leftlong,rightlong,uplat,downlat,gridon=0.1,leg="left",
          mincount=3)
plotLand(incol=rgb(1,0.5,0.5,1))


# Standardization

labelM <- c("Year","grid","FishID","Month","Effunit")
sps4 <- makecategorical(labelM,sps1)
mod <- makeonemodel(labelM)
out <- dosingle(mod,sps4)
str(out,max.level=1)
summary(out$optModel)
anova(out$optModel)
cbind(out$Results,out$StErr)
yrs <- facttonum(out$years); yrs

geom <-geomean(sps1$CE)
catch <- tapsum(sps,"catch_kg","Year")

plotprep(width=7,height=4.5)
plotstand(out,bars=TRUE,geo=geom)


invect <- sps1$CE


# Conduct a GAM
library(nlme)
library(mgcv)
library(gamm4)

modelGam <- gam(LnCE ~ s(Long,Lat) + Year + grid + FishID + Month + Effunit, data = sps1)
anova(modelGam)
modelGam <- gam(LnCE ~ s(Long,Lat) + Year + FishID + Month + Effunit, data = sps1)
anova(modelGam)



answer <- getfact(modelGam,"Year")
opti <- answer[,"Scaled"]
answer

anova(modelGam)
summary(modelGam)


plotprep(width=4.5,height=7)
plot(modelGam,ylim=c(-44.5,-40),xlim=c(143.5,146.5))
title(ylab=list("Latitude", cex=1.0, font=7),
      xlab=list("Longitude", cex=1.0, font=7))
plotLand("pink")

plotprep(width=7,height=4.5)
plotstand(out,bars=TRUE)
lines(facttonum(out$years),opti,col=4,lwd=2)
legend("bottomleft",c("GLM","GAM"),col=c(1,4),lwd=3,bty="n",cex=1.2)


plotprep(width=7,height=6)
plotstand(out,bars=TRUE,catch=catch)

columns <- c("adjR2","incR2","RSS","MSS","Npar","nobs","AIC")
nummod <- length(labelM)
results <- as.data.frame(matrix(0,nrow=nummod,ncol=length(columns),
                         dimnames=list(labelM,columns)))
for (i in 1:nummod) {
   mod <- makeonemodel(labelM[1:i])
   out <- dosingle(mod,sps1)
   outsum <- summary(out$optModel)
   aov <- anova(out$optModel)
   RSS <- tail(aov$"Sum Sq",1)
   df <- aov$Df
   nobs <- sum(df) + 1
   numfact <- length(df) - 1
   npars <- sum(df[1:numfact]) + 1
   AIC <- nobs * log(RSS/nobs) + (2 * npars)
   results[i,] <- c(outsum$adj.r.squared,NA,RSS,sum(aov$"Sum Sq")-RSS,npars,nobs,AIC)
}
results[2:nummod,"incR2"] <- results[2:nummod,"adjR2"]-results[1:(nummod-1),"adjR2"]
round(results,4)















