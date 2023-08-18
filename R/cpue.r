
#' @title addnorm - adds a normal distribution to a histogram of a data set.
#'
#' @description  addnorm - adds a normal distribution to a histogram of a data
#'    set. This is generally to be used to illustrate whether log-transformation
#'    normalizes a set of catch or cpue data.
#' @param inhist - is the output from a call to 'hist' (see examples)
#' @param xdata -  is the data that is being plotted in the histogram.
#' @param inc - defaults to a value of 0.01; is the fine grain increment used to
#'    define the normal curve. The histogram will be coarse grained relative to
#'    this.
#' @return a list with a vector of 'x' values and a vector of 'y' values (to be
#'    used to plot the fitted normal probability density function), and a vector
#'    used two called 'stats' containing the mean and sandard deviation of the
#'    input data
#' @export addnorm
#' @examples
#' x <- rnorm(1000,mean=5,sd=1)
#' dev.new(height=6,width=4,noRStudioGD = TRUE)
#' par(mfrow= c(1,1),mai=c(0.5,0.5,0.3,0.05))
#' par(cex=0.85, mgp=c(1.5,0.35,0), font.axis=7)
#' outH <- hist(x,breaks=25,col=3,main="")
#' nline <- addnorm(outH,x)
#' lines(nline$x,nline$y,lwd=3,col=2)
#' print(nline$stats)
addnorm <- function(inhist,xdata,inc=0.01) {
   lower <- inhist$breaks[1]
   upper <- tail(inhist$breaks,1)
   cw <- inhist$breaks[2]-inhist$breaks[1]
   x <- seq(lower,upper, inc) #+ (cw/2)
   avCE <- mean(xdata,na.rm=TRUE)
   sdCE <- sd(xdata,na.rm=TRUE)
   N <- length(xdata)
   ans <- list(x=x,y=(N*cw)*dnorm(x,avCE,sdCE),stats=c(avCE,sdCE,N))
   return(ans)
} # end of addnorm

#' @title countgtzero used in apply to count the number of zeros in a vector
#'
#' @description countgtzero used in apply to count the number of zeros in a vector
#' @param invect vector of values
#' @return A single value of zero or the number of zeros
#' @export
#' @examples
#' x <- rnorm(30,mean=1,sd=1)
#' countgtzero(x)
#' cat(30 - countgtzero(x),"out of ",30, " are less than 0  \n")
#' y <- matrix(x,nrow=6,ncol=5)
#' apply(y,1,countgtzero)
countgtzero <- function(invect) {
   pick <- which(invect > 0)
   return(length(pick))
}


#' @title dosingle conducts a standardization of indat using the inmodel
#'
#' @description dosingle conducts a standardization of indat using the
#'    inmodel.
#'
#' @param inmodel the formula used in the analysis; usually from the
#'    function makeonemodel.
#' @param indat the data.frame containing the data to be analysed.
#'
#' @return a list with a similar structure to the out object, so not
#'    a CEout class member but can be used with plotstand
#' @export  dosingle
dosingle <- function(inmodel,indat) {  # inmodel=mod; indat=sps1
   ans <- lm(inmodel,data=indat)
   bits <- unlist(strsplit(as.character(inmodel)," "))
   modcoef <- summary(ans)$coefficients
   years <- getfact(modcoef,bits[3])
   yrs <- sort(unique(indat[,bits[3]]))

   geo <- paste0(bits[2]," ~ ",bits[3])
   ans2 <- lm(as.formula(geo),indat)
   modcoefG <- summary(ans2)$coefficients
   yearsG <- getfact(modcoefG,bits[3])
   Results <- cbind("Year"=yearsG[,"Scaled"],"optimum"=years[,"Scaled"])
   rownames(Results) <- yrs
   StErr <- Results
   StErr[,1] <- yearsG[,"SE"]
   StErr[,2] <- years[,"SE"]
   optimum <- 2
   result <- list(Results=Results,StErr=StErr, Optimum=optimum,
                  modelcoef=modcoef,optModel=ans,modelG=ans2,years=yrs)
   return(result)
}  # end of dosingle

#' @title facttonum converts a vector of numeric factors into numbers
#'
#' @description facttonum converts a vector of numeric factors into numbers.
#'     If the factors are not numeric then the outcome will be a series of NA.
#'     It is up to you to apply this function only to numeric factors. A warning
#'     will be thrown if the resulting output vector contains NAs
#'
#' @param invect the vector of numeric factors to be converted back to numbers
#'
#' @return an output vector of numbers instead of the input factors
#' @export
#'
#' @examples
#' \dontrun{
#'  DepCat <- as.factor(rep(seq(100,600,100),2)); DepCat
#'  5 * DepCat[3]
#'  as.numeric(levels(DepCat))  # #only converts the levels not the replicates
#'  DepCat <- facttonum(DepCat)
#'  5 * DepCat[3]
#'  x <- factor(letters)
#'  facttonum(x)
#' }
facttonum <- function(invect){
   if (class(invect) == "factor") {
      outvect <- suppressWarnings(as.numeric(levels(invect))[invect])
   }
   if (class(invect) == "numeric") outvect <- invect
   if (any(is.na(outvect)))
      warning("NAs produced, your input vector may have non-numbers present \n")
   return(outvect)
} # end of facttonum

#' @title geomean: the geometric mean of a vector corrected for log-normal bias
#'
#' @description Calculate the geometric mean of a vector corrected for log-normal
#'   bias. NAs and zeros are removed from consideration.
#' @param invect is a vector of numbers in linear space.
#' @return The bias-corrected geometric mean of the vector
#' @export geomean
#' @examples
#' x <- c(1,2,3,4,5,6,7,8,9)
#' geomean(x)
## Calculate the geometric mean of a vector corected for log-normal bias
geomean <- function(invect) {
   pick <- which((invect <= 0.0))
   if (length(pick) == 0) {
      avCE <- mean(log(invect),na.rm=TRUE)
      stdev <- sd(log(invect),na.rm=TRUE)
   } else {
      avCE <- mean(log(invect[-pick]),na.rm=TRUE)
      stdev <- sd(log(invect[-pick]),na.rm=TRUE)
   }
   gmean <- exp(avCE + (stdev^2)/2)
   return(gmean)
}  # end of geomean

#' @title getfact extracts parameter estimates for a given factor from a CEout object
#'
#' @description getfact extracts the parameter estimates for a given factor from a CEout object.
#'    It does this by searching to rownames of the output parameters of the optimum model.
#'    It also checks for interaction terms, which for categorical factors is the same as
#'    determining a trend of the two factors relative to each other. e.g. for Zone:Month
#'    the outcome is the monthly trend for each zone.
#'
#' @param inmat generally this will be a CEout object but it can be
#'    a matrix of coefficients, an lm object, or a gam object
#' @param invar the model variable whose parameters are wanted
#'
#' @return a matrix containing the parameters for invar
#' @export getfact
getfact <- function(inmat,invar) {  # inmat=modcoef; invar=bits[3]
   allowable <- c("matrix","CEout","lm","gam")
   if (inherits(inmat,allowable)) {
   whatclass <- class(inmat)[1]
   if ("gam" %in% whatclass) {
         whatclass <- "gam"
    } else {
     if ("lm" %in% whatclass) {
        whatclass <- "lm"
     }
    }
   }
   if (whatclass %in% allowable) {
      if (whatclass == "matrix") pardat <- inmat
      if (whatclass == "CEout") pardat <- inmat$Parameters$coefficients
      if (whatclass == "lm") pardat <- summary(inmat)$coefficients
      if (whatclass == "gam") pardat <- summary(inmat)$p.table
   } else {
      stop("Input matrix is not, in fact, a matrix of coefficients")
   }
   pick <- grep(invar,rownames(pardat))
   if (length(pick) == 0) stop("The selected factor is not in the selected model!  \n")
   startmat <- pardat[pick,]               # isolate rows containing variable in question
   pickI <- grep(":",rownames(startmat))  # check for interaction terms and split off
   if (length(pickI) > 0) {               # if present
      intermat <- startmat[pickI,]
      startmat <- startmat[-pickI,]
   }
   lnce <- startmat[,"Estimate"]         # first do the non-interaction terms
   se <- startmat[,"Std. Error"]
   tval <- startmat[,"t value"]
   Prob <- startmat[,"Pr(>|t|)"]
   backtran <- exp(lnce + (se * se)/2)
   ans <- cbind(c(1.0,backtran),c(0,se),c(0,lnce),
                scaleCE(c(1.0,backtran),avCE=1.0),
                c(NA,tval),c(NA,Prob)
   )
   colnames(ans) <- c("Coeff","SE","LogCE","Scaled","t value","Prob")
   rownames(ans)[1] <- invar
   norigvar <- dim(ans)[1]
   if (length(pickI) > 0) {              # do interaction terms if they exist
      terms <- unlist(strsplit(rownames(intermat),":"))
      nrow <- dim(intermat)[1]
      firstvar <- grep(invar,terms)
      nfirst <- length(unique(terms[firstvar]))
      nsecond <- length(unique(terms[-firstvar]))
      if ((nfirst * nsecond) != nrow)
         stop(paste0("something wrong with the interactions terms for ",invar," \n",
                     "the first and second variables are not balanced"))
      start <- 1
      for (i in 1:(nrow/nfirst)) {  # i <- 1
         finish <- (start+nfirst-1)
         tmp <- intermat[start:finish,]
         lnce <- tmp[,"Estimate"]         # first do the non-interaction terms
         se <- tmp[,"Std. Error"]
         backtran <- exp(lnce + (se * se)/2)
         tmpans <- cbind(c(1.0,backtran),c(0,se),c(0,lnce),scaleCE(c(1.0,backtran),avCE=1.0))
         ans <- rbind(ans,tmpans)
         start <- finish + 1
      }
   }
   return(ans)
} # end of getfact


#' @title histyear plots a histogram of a given variable for each year available
#'
#' @description histyear plots a histogram of a given variable for each year
#'     available
#'
#' @param x the data.frame of data with at least a 'Year' and pickvar present
#' @param Lbound leftbound on all histograms, defaults to -3.5
#' @param Rbound right bound on all histograms, defaults to 12.25
#' @param inc  the class width of the histogram, defaults to 0.25
#' @param pickvar which variable to plot each year default = 'LnCE'
#' @param years which variable name identifies the yaer column, default='Year'
#' @param varlabel what label to use on x-axis, default = 'log(CPUE)'
#' @param vline an optional vertical line to aid interpretation. If it is
#'     numeric it will be added to each plot
#' @param plots how many plots to generate, default = c(3,3)
#'
#' @return a matrix of the year, mean value, stdev, and N number of
#'     observations. It also plots a histogram for each year and fits a
#'     normal distribution to each one.
#' @export
#'
#' @examples
#' \dontrun{
#' print("still to be developed")
#' }
histyear <- function(x,Lbound=-3.5,Rbound=12.25,inc=0.25,
                     pickvar="LnCE",years="Year",varlabel="log(CPUE)",
                     vline=NA,plots=c(3,3)) {
   yrs <- sort(unique(x[,years]))
   nyr <- length(yrs)
   columns <- c("Year","maxcount","Mean","StDev","N","Min","Max")
   results <- matrix(0,nrow=nyr,ncol=length(columns),dimnames=list(yrs,columns))
   par(mfcol=plots,mai=c(0.25,0.25,0.05,0.05),oma=c(1.2,1.0,0.0,0.0))
   par(cex=0.75, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
   for (yr in 1:nyr) {
      pick <- which(x[,years] == yrs[yr])
      outh <- hist(x[pick,pickvar],breaks=seq(Lbound,Rbound,inc),col=2,main="",xlab="",ylab="")
      mtext(paste0("  ",yrs[yr]),side=3,outer=F,line=-2,font=7,cex=0.9,adj=0)
      mtext(paste0("  ",length(pick)),side=3,outer=F,line=-3,font=7,cex=0.9,adj=0)
      if (is.numeric(vline)) abline(v=vline,col=4,lwd=2)
      if (pickvar != "catch_kg") {
         pickmax <- which.max(outh$counts)
         ans <- addnorm(outh,x[pick,pickvar])
         lines(ans$x,ans$y,col=3,lwd=2)
         results[yr,] <- c(yrs[yr],outh$mids[pickmax],ans$stats,
                           range(x[pick,pickvar],na.rm=TRUE))
      }
   }
   mtext("Frequency",side=2,outer=T,line=0.0,font=7,cex=1.0)
   mtext(varlabel,side=1,outer=T,line=0.0,font=7,cex=1.0)
   return(results)
} # end of histyear

#' @title lefthist draws a histogram up the y-axis
#'
#' @description lefthist translates a histogram from along the x-axis to
#'     flow along the y-axis - it transposes a histogram.
#'
#' @param x a vector of the data to be plotted
#' @param bins the breaks from the histogram, can be a single number of a
#'     sequence of values; defaults to 25
#' @param mult the multiplier for the maximum count in the histogram. Becomes
#'     the upper limit of teh x-axis.
#' @param col the colour for the histogram polygons; default = 2
#' @param lwd the line width for each polygon; default = 1
#' @param width the width of each bar in teh histogram; default = 0.9
#' @param border the colour for the border line; default = 1 = black
#' @param xinc the step size for the x-axis (counts) labels; default= NA,
#'     which means the increment will equal the bin width.
#' @param yinc the step size for the y-axis (breaks) labels; default= 1.
#' @param title the title for the left-histogram; defaults to ""
#' @param xlabel the xlab; defaults to ""
#' @param ylabel the ylab; defaults to "Frequency"
#' @param cex the size of text in teh plot. defaults = 1.0
#' @param textout prints input data range to console; default = FALSE
#' @param hline if this has a value a horizontal line will be plotted;
#'     default = NA
#'
#' @return the output from hist but done so invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- rnorm(1000,mean=5,sd=1)
#' dev.new(width=6,height=4,noRStudioGD = TRUE)
#' par(mai=c(0.45,0.45,0.05,0.05))
#' lefthist(dat)
#' lefthist(dat,textout=TRUE,width=0.8,border=3)
#' }
lefthist <- function(x,bins=25,mult=1.025,col=2,lwd=1,width=0.9,border=1,
                     xinc=1,yinc=NA,title="",xlabel="Frequency",ylabel="",
                     cex=1.0,textout=FALSE,hline=NA) {
   outh <- hist(x,breaks=bins,plot=FALSE)
   cw <- outh$breaks[2]-outh$breaks[1]
   newcount <- c(outh$counts,0)
   ymax <- max(newcount,na.rm=TRUE) * mult
   nvalues <- length(newcount)
   values <- outh$breaks
   if (is.na(yinc)) yinc <- values[2] - values[1]
   xlabs <- seq(0,(ymax+(2 * xinc)),xinc)
   xlabs <- xlabs[xlabs < ymax]
   plot(seq(0,ymax,length=nvalues),values,type="n",xlim=c(0,ymax),
        xaxs="i",ylim=c(range(values)),yaxs="r",xlab="",ylab="",xaxt="n",yaxt="n")
   grid()
   axis(side=1,at=xlabs,labels=xlabs)
   title(ylab=list(ylabel,cex=cex),xlab=list(xlabel,cex=cex))
   values1 <- seq(values[1],values[nvalues],yinc)
   axis(side=2,at=values1,labels=values1,cex.lab=cex)
   for (i in 1:nvalues) {  # i <- 1
      y1 <- values[i]
      y2 <- values[i] + (cw * width)
      yplot <- c(y1,y1,y2,y2,y1)
      xplot <- c(0,newcount[i],newcount[i],0,0)
      if (is.null(border)) border <- col
      polygon(xplot,yplot,col=col,border=border,lwd=lwd)
   }
   if (!is.na(hline)) abline(h=hline,col=(col+2))

   if (textout) cat("  input data range: ",range(x,na.rm=TRUE),"\n\n")
   return(invisible(outh))
}  # end of lefthist

#' @title makecategorical - converts given variables into categorical factors
#'
#' @description Given a list of variables, as character strings, this function
#'   converts each variable into a factor after checking no mistake has been
#'   made with the spelling of the variable name. A copy of the data is made
#'   so that the non factored variables can be retained.
#' @param labelModel is the set of variables from the data.frame that are to
#'   be converted into categorical factors.
#' @param indat is the data.frame that is to be analysed in the standardization
#' @return a copy of the database with the selected variables converted into
#'     factors
#' @export makecategorical
#' @examples
#' \dontrun{
#' data(sps)
#' labelM <- c("Year","Vessel","Month")
#' sps1 <- makecategorical(labelM,sps)
#' }
makecategorical <- function(labelModel,indat) {
   Interact <- grep(":",labelModel)
   nInteract <- length(Interact)
   numvars <- length(labelModel) - nInteract
   for (fac in 1:numvars) {
      if (length(indat[,labelModel[fac]]) > 0) {
         indat[,labelModel[fac]] <- factor(indat[,labelModel[fac]])
      } else { warning(paste0("Factor name ",labelModel[fac],
                              "does not appear in data.frame"))
      }
   }
   return(indat)
} # end of makecategorical

#' @title makeonemodel generates a single model for use in dosingle
#'
#' @description makeonemodel puts together the formula needed to
#'     run a statistical model, but, different to makemodels, it
#'     only generates a single model.
#'
#' @param labelModel a vector of labels for each factor to be included
#'    in the analysis
#' @param dependent the name of the dependent variable; defaults to
#'    LnCE
#'
#' @return a formula of 'dependent ~ labelModel components
#' @export makeonemodel
#'
#' @examples
#' labelM <- c("Year","Vessel","DepCat","Zone:Month")
#' makeonemodel(labelM)
#' makeonemodel(labelM,dependent="LnCE")
makeonemodel <- function(labelModel,dependent = "LnCE") {
   numvars <- length(labelModel)
   interterms <- grep(":", labelModel)
   ninter <- length(interterms)
   form <- paste0(dependent, " ~ ", labelModel[1])
   if (numvars > 1) {
      for (i in 2:(numvars - ninter))
         form <- paste0(form, " + ", labelModel[i])
      if (ninter > 0) for (i in interterms) {
         form <- paste0(form, " + ", labelModel[i])
      }
   }
   form <- as.formula(form)
   return(form)
} # end of makeonemodel


#' @title plotprep: sets up a window and the par values for a single plot
#'
#' @description plotprep: sets up a window and the par values for a single plot.
#'   it checks to see if a graphics device is open and opens a new one if not.
#'   This is simply a utility function to save typing the standard syntax.
#'   Some of the defaults can be changed. Typing the name without () will
#'   provide a template for modification. If 'windows' is called repeatedly this
#'   will generate a new active graphics device each time leaving the older ones
#'   inactive but present. For quick exploratory plots this behaviour is not
#'   wanted, hence the check if an active device exists already or not.
#' @param width defaults to 6 inches = 15.24cm - width of plot
#' @param height defaults to 3 inches = 7.62cm - height of plot
#' @param plots defaults to c(1,1), but arranges multiple plots. If used it may
#'    be necessary to print out this code and adjust the mai and oma variables
#' @param usefont default is 7 (bold Times); 1 = sans serif, 2 = sans serif bold
#' @param cex default is 0.85, the size of font used for text within the plots
#' @param xmtext default is TRUE; if plots is not c(1,1) this alters the mai and
#'    oma variables for the x-axis to allow for mtexting and avoid the x title
#' @param ymtext default is TRUE; if plots is not c(1,1) this alters the mai and
#'    oma variables for the y-axis to allow for mtexting and avoid the y title
#' @param newdev reuse a previously defined graphics device or make a new one;
#'    defaults to TRUE
#' @param rows defaults to TRUE, determines whether to use mfrow or mfcol
#' @param filename defaults to "" = do not save to a filename. If a filename is
#' @return Checks for and sets up a graphics device and sets the default plotting
#'   par values. This changes the current plotting options!
#' @export plotprep
#' @examples
#' x <- rnorm(1000,mean=0,sd=1.0)
#' plotprep()
#' hist(x,breaks=30,main="",col=2)
plotprep <- function(width=6,height=3.6,plots=c(1,1),usefont=7,cex=0.85,
                     xmtext=TRUE,ymtext=TRUE,
                     newdev=TRUE,rows=TRUE,filename="") {
   if  ((names(dev.cur()) != "null device") & (newdev)) suppressWarnings(dev.off())
   lenfile <- nchar(filename)
   if (lenfile > 3) {
      end <- substr(filename,(lenfile-3),lenfile)
      if (end != ".png") filename <- paste0(filename,".png")
      png(filename=filename,width=width,height=height,units="in",res=300)
   } else {
      if (names(dev.cur()) %in% c("null device","RStudioGD"))
         dev.new(width=width,height=height,noRStudioGD = TRUE)
   }
   firstmai <- 0.45; secondmai <- 0.45; thirdmai <- 0.1
   firstoma <- 0.0; secondoma <- 0.0; thirdoma <- 0.0
   if (sum(plots) != 2) {
      if (xmtext) {
         firstmai <- 0.25
         thirdmai <- 0.05
         firstoma <- 1.0
         thirdoma <- 0.1
      }
      if (ymtext) {
         secondmai <- 0.25
         secondoma <- 1.0
      }
   }
   maival <- c(firstmai,secondmai,thirdmai,0.05)
   omaval <- c(firstoma,secondoma,thirdoma,0.0)
   if (rows) {
      par(mfrow = plots,mai=maival,oma=omaval)
   } else {
      par(mfcol = plots,mai=maival,oma=omaval)
   }
   par(cex=cex, mgp=c(1.35,0.35,0), font.axis=usefont,font=usefont,font.lab=usefont)
   if (lenfile > 0) cat("\n Remember to place 'graphics.off()' after the plot \n")
} # end of plotprep


#' @title plotstand - plot optimum model from standLM  vs Geometric mean
#'
#' @description plot optimum model from standLM  vs Geometric mean.
#'   Has options that allow for log-normls P% intervals around each time
#'   period's parameter estimate. Also can rescale the graph to have an average
#'   the same as the geometric mean average of the original time series of data.
#' @param stnd is the list output from standLM
#' @param bars is a logical T or F determining whether to put confidence bounds
#'   around each estimate; defaults to FALSE
#' @param geo is an estimate of the original geometric mean catch rate across
#'   all years. If this is > 0.0 it is used to rescale the graph to the
#'   nominal scale, otherwise the mean of each time-series will be 1.0, which
#'   simplifies visual comparisons. geo defaults to 0.0.
#' @param P is the percentile used for the log-normal confidence bounds, if
#'   they are plotted; defaults to 95.
#' @param catch if it is desired to plot the catch as well as the CPUE
#'   then a vector of catches needs to be input here
#' @param usefont enables the font used in the plot to be modified. Most
#'   publications appear to prefer usefont=1; defaults to 7 - Times bold
#' @return a plot of the model with the smallest AIC (solid line) and the
#'   geometric mean (model 1, always = LnCE ~ Year, the dashed line). 'Year'
#'   could be some other time step.
#' @export plotstand
#' @examples
#' \dontrun{
#' data(sps)
#' splabel = "SpeciesName"
#' labelM <- c("Year","Vessel","Month")
#' sps1 <- makecategorical(labelM[1:3],sps)
#' mods <- makemodels(labelM)
#' out <- standLM(mods,sps1,splabel)
#' plotprep()
#' plotstand(out, bars=TRUE, P=90,geo=100.0,usefont=1)
#' plotstand(out)
#' }
plotstand <- function(stnd,bars=FALSE,geo=0.0,P=95,catch=NA,usefont=7) {
   result <- stnd$Results
   if (geo > 0.0) result <- stnd$Results*geo
   sterr <- stnd$StErr
   whichM <- stnd$WhichM
   optimum <- stnd$Optimum
   years <- rownames(result)
   fishyr <- FALSE
   if (nchar(years[1]) > 4) {
      yrs <- 1:length(years)
      fishyr <- TRUE
   } else {
      yrs <- as.numeric(rownames(result))
   }
   laby <- paste(stnd$Label," CPUE",sep="")
   if (bars) {
      Zmult <- -qnorm((1-(P/100))/2.0)
      lower <- result[,optimum] * exp(-Zmult*sterr[,optimum])
      upper <- result[,optimum] * exp(Zmult*sterr[,optimum])
      ymax <- max(result[,1],result[,optimum],upper,na.rm=TRUE)*1.025
   } else {
      ymax <- max(result[,1],result[,optimum],na.rm=TRUE)*1.025
   }
   if (length(catch) > 1) par(mfrow= c(2,1)) else par(mfrow= c(1,1))
   par(mai=c(0.4,0.5,0,0), oma=c(0,0,0.25,0.25))
   par(cex=0.85, mgp=c(1.5,0.3,0), font.axis=usefont)
   plot(yrs,result[,1],type="l",lty=2,lwd=2,ylim=c(0,ymax),yaxs="i",ylab="",
        xlab="",xaxs="r")
   lines(yrs,result[,optimum],lwd=3)
   if (bars) {
      arrows(x0=yrs[-1],y0=lower[-1],x1=yrs[-1],y1=upper[-1],
             length=0.035,angle=90,col=2,lwd=2,code=3)
   }
   title(ylab=list(laby, cex=1.0, col=1, font=usefont))
   if (geo > 0.0) {
      abline(h=geo,col="grey")
   } else {
      abline(h=1.0,col="grey")
   }
   if (length(catch) > 1) {
      if (length(catch) != length(yrs)) {
         warning("input catch data has incorrect number of years")
      }
      ymax <- max(catch,na.rm=TRUE) * 1.05
      plot(yrs,catch,type="b",pch=16,cex=0.8,lwd=2,ylim=c(0,ymax),yaxs="i",ylab="",
           xlab="",xaxs="r")
      grid()
      title(ylab=list("Catch", cex=1.0, col=1, font=usefont))
   }
} # End of plotstand


#' @title properties - used to check a data.frame before standardization
#'
#' @description properties - used to check a data.frame before
#'     standardization
#' @param indat the data.frame containing the data fields to be used
#'     in the subsequent standardization. It tabulates the number of
#'     NAs and the number of unique values for each variable and finds
#'     the minimum and maximum of the numeric variables
#' @param dimout determines whether or noth the dimensions of the data.frame
#'     are printed to the screen or not; defaults to FALSE
#' @return a data.frame with the rows being each variable from the input
#'     input data.frame and the columns being the number of NAs, the
#'     number of unique values, and minimum and maximum (where possible).
#' @export properties
#' @examples
#' \dontrun{
#' data(ab)
#' properties(ab)
#' }
properties <- function(indat,dimout=FALSE) {
   if(dimout) print(dim(indat))
   isna <- sapply(indat,function(x) sum(is.na(x)))
   uniques <- sapply(indat,function(x) length(unique(x)))
   clas <- sapply(indat,class)
   numbers <- c("integer","numeric")
   pick <- which(clas %in% numbers)
   minimum <- numeric(length(uniques))
   maximum <- minimum
   minimum[pick] <- sapply(indat[,pick],min,na.rm=TRUE)
   maximum[pick] <- sapply(indat[,pick],max,na.rm=TRUE)
   index <- 1:length(isna)
   props <- as.data.frame(cbind(index,isna,uniques,clas,minimum,
                                maximum,t(indat[1,])))
   colnames(props) <- c("Index","isNA","Unique","Class","Min",
                        "Max","Example")
   for (i in 5:6) props[,i] <- as.numeric(props[,i])
   return(props)
} # end of properties


#' @title qqplotout plots up a single qqplot for a lm model
#'
#' @description qqplotout generates a single qqplot in isolation from the
#'     plot of a model's diagnostics. It is used with lefthist to
#'     illustrate how well a model matches a normal distribution
#'
#' @param inmodel the optimum model from standLM or dosingle
#' @param title a title for the plot, defaults to 'Normal Q-Q Plot'
#' @param cex the size of the font used, defaults to 0.9
#' @param ylow the lower limit of the residuals
#' @param yhigh he upper limit of the residuals
#' @param plotrug a logical value determinning whether a rug is included
#'
#' @return currently nothing, but it does generate a qqplot to the current
#'     device
#' @export
#'
#' @examples
#' y <- rep(1:100,2)
#' x <- rnorm(200,mean=10,sd=1)
#' model <- lm(y ~ x)
#' dev.new(width=6,height=3.5,noRStudioGD = TRUE)
#' par(mai=c(0.45,0.45,0.15,0.05),font.axis=7)
#' qqplotout(model,ylow=-50,yhigh=50)
qqplotout <- function(inmodel, title="Normal Q-Q Plot", cex=0.9,
                      ylow=-5,yhigh=5,plotrug=FALSE)  {
   resids <- inmodel$residuals
   labs <- cex
   qqnorm(resids, ylab=list("Standardized Residuals", cex=labs, font=7),
          xlab=list("Theoretical Quantiles", cex=labs, font=7),
          main=list(title,cex=labs,font=7),ylim=c(ylow,yhigh))
   qqline(resids, col=2,lwd=2)
   grid()
   if (plotrug) rug(resids)
   abline(v=c(-2.0,2.0),col="grey")
}  # end of qqplotout

#' @title qqdiag generates a qqplot with a histogram of residuals
#'
#' @description qqdiag generates a qqplot with a complementary histogram of
#'     the residuals to illustrate the proportion of all residuals along the
#'     qqline. If the qqline deviates from the expected straigt line, which
#'     is red i colour to make for simpler comparisons, then the histogram
#'     enables one to estiamte what proportion of records deviate from
#'     normality. The zero point is identified with a line, as are the
#'     approximate 5% and 95% percentiles. In both cases > 5% is above or
#'     below the blue lines, with < 90% in between depending on the
#'     proportions in each class. To get a more precise estimate use the
#'     invisibly returned histogram values.
#'
#' @param inmodel the optimum model being considered
#' @param plotrug a logical term determining whether a rug is plotted on the
#'     qqplot.
#' @param bins defaults to NA, but can be set to a given series
#' @param hline Include some horizontal lines on the histogram. defaults to 0.
#' @param xinc the increment for tick marks on the xaxis of the histogram
#' @param yinc the increment for tick marks on the y-axis of the histogram
#' @param ylab the y-axis label for the histogram, defaults to 'residuals'
#'
#' @return plots a graph and invisibly returns the output from the histogram
#' @export
#'
#' @examples
#'
#' y <- rep(1:100,2)
#' x <- rnorm(200,mean=10,sd=1)
#' model <- lm(y ~ x)
#' dev.new(width=6,height=3.5,noRStudioGD = TRUE)
#' par(mai=c(0.45,0.45,0.15,0.05),font.axis=7)
#' qqdiag(model,xinc=1,yinc=10,bins=seq(-55,50,2.5))
qqdiag <- function(inmodel,plotrug=FALSE,bins=NA,hline=0.0,
                   xinc=100,yinc=1,ylab="residuals") {
   layout(matrix(c(1,2),ncol=2),widths=c(5,2.5))
   par(mai=c(0.45,0.45,0.15,0.05),oma=c(0.0,0,0.0,0.0))
   par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
   resids <- inmodel$residuals
   qs <- quantile(resids,probs=c(0.025,0.05,0.95,0.975))
   if (!is.numeric(bins)) {
      loy <- min(resids); hiy <- max(resids)
      scale <- trunc(100*(hiy - loy)/35) / 100
      loy <- round(loy - (scale/2),2); hiy <- round(hiy + scale,2)
      bins <- seq(loy,hiy,scale)
   } else {
      loy <- min(bins); hiy <- max(bins)
   }
   qqplotout(inmodel,plotrug=plotrug,ylow=loy,yhigh=hiy)
   abline(h=qs,lwd=c(1,2,2,1),col=4)
   outL <- lefthist(resids,bins=bins,hline=0.0,yinc=yinc,xinc=xinc,
                    ylabel=ylab,width=0.9,border=1)
   abline(h=qs,lwd=c(1,2,2,1),col=4)
   ans <- addnorm(outL,resids)
   lines(ans$y,ans$x,lwd=2,col=3)
   return(invisible(outL))
}  # end of qqdiag

#' @title quants used in apply to count the number > 1 in a vector
#'
#' @description quants used in apply to count the number > 1 in a vector
#'    designed to be used in apply
#' @param invect vector of values
#' @return a vector of the c(0.025,0.05,0.5,0.95,0.975) quantiles
#' @export
#' @examples
#' x <- rnorm(100,mean=5,sd=1)
#' quants(x)
#' y <- matrix(x,nrow=10,ncol=10)
#' apply(y,2,quants)
quants <- function(invect) {
   ans <- quantile(invect,probs = c(0.025,0.05,0.5,0.95,0.975),na.rm=T)
   return(ans)
}



#' @title scaleCE - Function to scale a vector of CPUE to a mean of one of avCE
#'
#' @description Function to scale a vector of CPUE to a mean of one or avCE.
#'   The use of a mean of one means that visual comparisons between different
#'   time-series becomes visually simplified. The avCE option could be used
#'   to scale the CPUE to the average geometric mean - so as to put it on
#'   the nominal scale
#' @param invect a vector of linear scale CPUE
#' @param avCE defaults to one but can be set to a particular value
#' @return the time-series of CPUE re-scaled to a mean of one or avCE
#' @export scaleCE
#' @examples
#' ce <- c(0.4667187,1.2628564,0.8442146,0.9813531, 0.5554076,0.7426321)
#' scaleCE(ce)
#' scaleCE(ce,100.0)
scaleCE <- function(invect,avCE=1.0) {
   average <- mean(invect,na.rm=TRUE)
   ans <- avCE * (invect/average)
   return(ans)
} # end of scaleCE

#' @title tapsum simplifies the use of tapply for summarizing variables
#'
#' @description data exploration commonly uses the tapply function and tapsum
#'     simplifies its use when obtaining the sum of any variable relative to
#'     other variables. For example it is common to want the total catch by
#'     year and, for example, Month, DepCat, Zone, etc.
#'
#' @param indat the data.frame containing the raw fishery data
#' @param first the variable name (in quotes) being summed
#' @param second the first grouping variable
#' @param third the second grouping variable, defaults to NA
#' @param div defaults to 1000 to change Kg to tonnes. set t0 1.0 or NA to
#'     avoid its influence
#'
#' @return a vector or matrix of sums of the pickvar by the frist and optionally
#'      the second grouping variable
#' @export
#'
#' @examples
#' \dontrun{
#'   data(sps)
#'   tapsum(sps,"catch_kg","Year","Month")
#' }
tapsum <- function(indat,first,second,third=NA,div=1000) {
   if (is.na(third)) {
      result <- tapply(indat[,first],indat[,second],sum,na.rm=TRUE)
   } else {
      result <- tapply(indat[,first],list(indat[,second],indat[,third]),
                       sum,na.rm=TRUE)
   }
   if (is.numeric(div)) result <- result/div
   return(result)
} # end of tapsum

#' @title turnover Estimate turnover of vessels from catch by vessel by year data
#'
#' @description Estimate turnover of vessels from catch by vessel by year data;
#'   To specify the minimum number of years that a vessel needs stay in the fishery,
#'   then give a value to the variable minyrs.
#' @param x A matrix of a continuous numeric property by year,
#'    the original usage was to plot catch-by-vessel against year
#' @param minyrs limits the analysis to those vessels that remain in the
#'    fishery for at least minyrs years - which would eliminate the occasional
#'    opportunistic fisher who only fishes for one or two years, or whatever
#'    minimum is selected. Vessels with zero catches are not included in case
#'    zeros and NAs are counted as starting and leaving the fishery.
#' @return a matrix of years by Continue, Leave, Start, Total
#' @export turnover
#' @examples
#' \dontrun{
#' library(r4cpue)
#' data(sps)
#' cbv <- tapply(sps$catch_kg,list(sps$Vessel,sps$Year),sum,na.rm=TRUE)/1000
#' dim(cbv)
#' early <- rowSums(cbv[,1:6],na.rm=TRUE)
#' late <- rowSums(cbv[,7:14],na.rm=TRUE)
#' cbv1 <- cbv[order(late,-early),]
#' plotprep(width=7,height=6)
#' yearBubble(cbv1,ylabel="Catch by Trawl",vline=2006.5,diam=0.2)
#' turnover(cbv)
#' }
turnover <- function(x,minyrs=1) {
   years <- as.numeric(colnames(x))
   ny <- length(years)
   count <- apply(x,1,countgtzero)
   pick <- which(count == 0)
   if (length(pick) > 0) {
      warning("Some Rows only have zeros or NAs")
      x <- x[-pick,]
   }
   count <- apply(x,1,countgtzero)
   pick <- which(count > (minyrs - 1))
   if (length(pick) > 0) x <- x[pick,]
   columns <- c("Continue","Leave","Start","Total")
   turnover <- matrix(0,nrow=ny,ncol=length(columns),dimnames=list(years,
                                                                   columns))
   turnover[1,1] <- length(which(x[,1] > 0))
   for (yr in 2:ny) {
      pair <- x[,(yr-1):yr]
      pickC <- which((pair[,1] > 0) & (pair[,2] > 0))
      pickL <- which((pair[,1] > 0) & ((is.na(pair[,2])) |
                                          (pair[,2] < 0.001)))
      pickS <- which(((is.na(pair[,1])) | (pair[,2] < 0.001)) &
                        (pair[,2] > 0))
      turnover[yr,1:3] <- c(length(pickC),length(pickL),length(pickS))
   }
   turnover[,4] <- turnover[,1] + turnover[,3]
   return(turnover)
} # end of turnover

#' @title yearBubble: Generates a bubbleplot of x against Year.
#'
#' @description yearBubble: Generates a bubbleplot of x against Year.
#' @param x - a matrix of variable * Year; although it needn't be year
#' @param xlabel - defaults to nothing but allows a custom x-axis label
#' @param ylabel - defaults to nothing but allows a custom y-axis label
#' @param diam - defaults to 0.1, is a scaling factor to adjust bubble size
#' @param vline - defaults to NA but allows vertical ablines to higlight regions
#' @param txt - defaults are lines to vessel numbers, catches, catches, maximumY
#' @param Fyear - defaults to FALSE, if TRUE generates a fishing year x-axis
#' @param xaxis - defaults to TRUE, allows for a custom x-axis if desired by
#'    using something like axis(1,at=years,labels=years).
#' @param yaxis - defaults to TRUE, allows for a custom y-axis if desired by
#'    using something like axis(side=2,at=years,labels=years).
#' @param hline - defaults to FALSE
#' @param nozero - defaults to FALSE, if TRUE replaces all zeros with NA so they
#'    do not appear in the plot
#' @return - invisible, vectors of catch and vessels by year, and radii matrix
#' @export yearBubble
#' @examples
#' \dontrun{
#' data(sps)
#' cbv <- tapply(sps$catch_kg,list(sps$Vessel,sps$Year),sum,na.rm=TRUE)/1000
#' dim(cbv)
#' early <- rowSums(cbv[,1:6],na.rm=TRUE)
#' late <- rowSums(cbv[,7:14],na.rm=TRUE)
#' cbv1 <- cbv[order(late,-early),]
#' plotprep(width=7,height=6)
#' yearBubble(cbv1,ylabel="Catch by Trawl",vline=2006.5,diam=0.2)
#' }
yearBubble <- function(x,xlabel="",ylabel="",diam=0.1,vline=NA,txt=c(4,6,9,11),
                       Fyear=FALSE,xaxis=TRUE,yaxis=TRUE,hline=FALSE,nozero=FALSE) {
   nyrs <- dim(x)[2]
   if (Fyear) {
      tyrs <- colnames(x)  # assumes a yyyy/yyyy format
      if (nchar(tyrs[1]) != 9) warning("Wrong fishing year format for yearBubble \n")
      years <- as.numeric(substr(tyrs,1,4))
   } else { years <- as.numeric(colnames(x)) # assumes columns are years
   }
   nves <- length(rownames(x))
   yvar <- seq(1,nves,1)
   if (nozero) {
      pick <- which(x == 0)
      x[pick] <- NA
   }
   radii <- sqrt(x)
   biggest <- max(radii,na.rm=TRUE)
   catch <- colSums(x,na.rm=TRUE)   # total annual catches
   numves <- apply(x,2,function(x1) length(which(x1 > 0))) # num vess x year
   answer <- list(catch,numves,radii) # generate output
   names(answer) <- c("Catch","Vessels","Radii")
   xspace <- 0.3
   if (nchar(xlabel) > 0) xspace <- 0.45
   par(mfrow= c(1,1))
   par(mai=c(xspace,0.45,0.1,0.1), oma=c(0.0,0.0,0.0,0.0))
   par(cex=0.85, mgp=c(1.5,0.3,0), font.axis=7,font=7)
   xt <- "s"
   yt <- "s"
   if (!xaxis) xt <- "n"
   if (!yaxis) yt <- "n"
   plot(years,years,type="n",xlab="",ylab="",ylim=c(0,(nves+txt[4])),yaxs="r",
        yaxt=yt,xaxt=xt,xaxs="r")
   if (hline) abline(h=yvar,col="grey")
   for (x in 1:nyrs) {
      yr <- years[x]
      odd.even<-x%%2
      if (odd.even == 0) text(yr,nves+txt[3],round(catch[x],0),cex=0.65,font=7)
      else text(yr,nves+txt[2],round(catch[x],0),cex=0.65,font=7)
      text(yr,nves+txt[1],numves[x],cex=0.8,font=7)
      mult <- max(radii[,x],na.rm=TRUE)/biggest
      symbols(rep(yr,nves),yvar,circles=radii[,x],inches=diam*mult,
              bg=rgb(1, 0, 0, 0.5), fg = "black",xlab="",ylab="",add=TRUE)
   }

   if (length(vline) > 0) abline(v=c(vline),col="grey")
   title(ylab=list(ylabel, cex=1.0, col=1, font=7))
   return(invisible(answer))
} # end of YearBubble



#' @title yearNA - counts NAs per year in each numeric field in a data.frame
#'
#' @description yearNA - counts the number of NAs in each year of each numeric
#'    field in a data.frame and outputs the results as a matrix
#' @param indat the data.frame whose numeric fields are to be considered
#' @param years identifies the name of the "Year" field in the data.frame
#' @param empty logical default=FALSE, determines whether columns which have
#'    no NA present are printed or not
#' @return a matrix of years x numeric fields with the number of records per
#'    field per year as reference
#' @export yearNA
#' @examples
#' year <- sort(rep(1990:1994,5))
#' columns <- c("year","Var1","Var2","Var3")
#' dat <- matrix(runif(100),nrow=25,ncol=4,dimnames=list(year,columns))
#' dat[trunc(100*runif(20))] <- NA
#' dat[,1] <- year
#' print(dat)
#' yearNA(as.data.frame(dat),years="year")
yearNA <- function(indat,years="Year",empty=FALSE) {
   records <- table(indat[,years])
   nna <- function(x) sum(is.na(x))
   columns <- colnames(indat)
   clas <- sapply(indat, class)
   numbers <- c("integer", "numeric")
   pick <- which(clas %in% numbers)
   num <- length(pick)
   ans <- NULL
   for (i in 1:num) {
      ans <- cbind(ans,tapply(indat[,columns[pick[i]]],indat[,years],nna))
   }
   ans <- cbind(ans,records)
   colnames(ans) <- c(columns[pick],"Records")
   if (!empty) {
      pick <- which(colSums(ans,na.rm=TRUE) > 0)
      if (length(pick) > 0)  ans <- ans[,pick]
   }
   return(ans)
}  #end of yearNA
