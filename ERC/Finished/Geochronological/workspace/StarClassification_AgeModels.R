# R code to implement a start classification system for classical (non-Bayesian) age-depth modelling.
# Maarten Blaauw <maarten.blaauw@qub.ac.uk>
# see for more info on clam: accompanying manual.html and Blaauw 2010 (Quaternary Geochronology 5: 512-518)
# This script was originally developed for the European Pollen Database in Giesecke et al 2012 and adapted for the Latin American Pollen Database in Flantua, Blaauw, Hooghiemstra, 2016.

# Note: the default values can be changed permanently within this file
# -----------------------------------------------------------------------------------

#source("clam.R")

correct.files <- function()
  {
    for(i in list.files("Cores/"))
      {
	i <<- i
	base.name <- paste("Cores/", i, "/", i, sep="")

	if(file.exists(paste(base.name, "_ages.csv", sep="")))
	  file.rename(paste(base.name, "_ages.csv", sep=""), paste(base.name, ".csv", sep=""))

	tmp <- read.csv(paste(base.name, ".csv", sep=""))
 	if(suppressWarnings(min(tmp[,4])) < 1) # no errors of 0 yr!!!
 	    stop(i, " has <1 yr errors")
	write.csv(tmp, paste(base.name, ".csv", sep=""), quote=F, row.names=F)

      	tmp <- read.csv(paste(base.name, "_depth.csv", sep=""))
	tmp <- ordered(unique(tmp[,2]))
	write.table(tmp, paste(base.name, "_depths.txt", sep=""), row.names=F, col.names=F, quote=F)
    }
  }

runname <- c("_interpolated", "_polyn_regr", "_cubic_spline", "_smooth_spline", "_loess")

EPD.stars <- function(core, type, straight=4, straight.var=.2, bad.extrapol=2000, good.extrapol=1000, best.extrapol=500, plot.points=TRUE, youngest=c(), cc=1, hiatus=c(), outliers=c(), slump=c(), postbomb=1, ...)
  {
    dfile <- paste("Cores/", core, "/", core, "_depths.txt", sep="")
    if(file.exists(dfile))
      depths.file <- TRUE else depths.file <- FALSE

    clam(core, type, storedat=TRUE, cc=cc, hiatus=hiatus, outliers=outliers, slump=slump, postbomb=postbomb, depths.file=depths.file, depthseq=c())

    if(!file.exists(dfile))
      write.table(min(dat$depth) : max(dat$depth), dfile, row.names=FALSE, col.names=FALSE)

    dd <- read.table(dfile)[,1]
    dd <<- dd
    whichd <- which(calrange[,1] %in% dd) # produce graphs using many depths, but later analyse only the requested depths
    if(min(dat$depth) < min(dd)) # if there are dated levels above the depths with proxy information...
      whichd <- c(1, whichd) # then add the first dated depth to the required depths
    d <- calrange[whichd,1]; age.min <- min(calrange[,2], dat$model); age.max <- max(calrange[,3], dat$model); age <- calrange[whichd,4]
    section <- c(); constant <- c(); good <- c(); far <- c(); best <- c()

    d <<- d
    if(length(dat$model) > straight) # enough dates to look for straight sections?
      for(i in (straight+1):length(dat$model)) # find approximately straight sections
	{
	  acc <- min(which(d >= dat$depth[i-straight])) : max(which(d <= dat$depth[i]))
	  acc <- diff(age[acc]) / diff(d[acc])
	  acc <- acc/mean(acc) # normalise
	  if(length(acc)> 1 && (var(acc) < straight.var)) #  is acc.var within acceptable range?
	    section <- c(section, i)
        }

    if(length(section) > 0)
      for(i in section) # which dates are part of approximately straight sections?
	constant <- sort(unique(c(constant,
	  min(which(d >= dat$depth[i-straight]))[1] : max(which(d <= dat$depth[i]))[1])))

    for(i in dat$model) # which sections extrapolate within "bad.extrapol"?
      {
	temp <- which(age <= (i+bad.extrapol))
	far <- sort(unique(c(far, which(age[temp] >= (i-bad.extrapol)))))
      }

    for(i in dat$model) # which sections extrapolate within "good.extrapol"?
      {
	temp <- which(age <= i+good.extrapol)
	good <- sort(unique(c(good, which(age[temp] >= (i-good.extrapol)))))
      }

    for(i in dat$model) # which sections extrapolate within "best.extrapol"?
      {
	temp <- which(age <= i+best.extrapol)
	best <- sort(unique(c(best, which(age[temp] >= (i-good.extrapol)))))
      }

#-------------------
    star <- array(0, dim=c(length(d), 6))
    star[,1] <- d
    star[constant,2] <- 1
    star[far,3] <- 1
    star[good,4] <- 1
    star[best,5] <- 1
# Commented out the condition for reversals
#    if(min(diff(age)) < 0)
#      star[,c(2,3,4,5)] <- 0 # no stars at all if core has age reversals
    star[,6] <- rowSums(star[,2:5])

    if(plot.points)
    # if(min(diff(age)) < 0)
    # points(max(age.max), min(d), pch=4, cex=3, col=2, font=2) else
	{
	  points(rep(age.max-0.00*(age.max-age.min), length(constant)), d[constant], pch=10, col=2, cex=.2)
	  points(rep(age.max-0.01*(age.max-age.min), length(far)), d[far], pch=10, col=3, cex=.2)
	  points(rep(age.max-0.02*(age.max-age.min), length(good)), d[good], pch=10, col=4, cex=.2)
	  points(rep(age.max-0.03*(age.max-age.min), length(best)), d[best], pch=10, col=5, cex=.2)
	}
    if(type==1)
      dev.copy2pdf(file=paste("Cores/", core, "/", core, "_interpolated.pdf", sep=""))
    if(type==4)
      dev.copy2pdf(file=paste("Cores/", core, "/", core, "_smoothspline.pdf", sep=""))

    colnames(star) <- c("depth", "constant", "bad.extra", "good.extra", "best.extra", "stars")
    write.table(star, paste(dat$coredir, core, runname[type], "_stars.txt", sep=""), row.names=FALSE, sep="\t")
  }

Agemodel.stars <- function(core, cc=c(), postbomb=1, outliers=c(), slump=c(), hiatus=c(), youngest=c())
  {
    tmp <- read.csv(paste("Cores/", core, "/", core, ".csv", sep=""))
    n.dates <- length(unique(tmp[,6]))

    youngest <- c()
    if(tmp[1,8] == "TOP_A" || tmp[1,8] == "TOP_B" || tmp[1,8] == "TOP_C")
    if(!is.na(tmp[1,3]))
       youngest <- tmp[1,3]
    if(length(outliers) == 0) Outliers <- 0 else Outliers <- length(outliers)
    if((n.dates - Outliers) > 1) # linear interpolation if >1 date
      EPD.stars(core, 1, plotpdf=FALSE, youngest=youngest, cc=cc, postbomb=postbomb, outliers=outliers, slump=slump, hiatus=hiatus)
    if((n.dates - Outliers) > 3) # smooth spline if more than 4 dates
      EPD.stars(core, 4, plotpdf=FALSE, youngest=youngest, cc=cc, postbomb=postbomb, outliers=outliers, slump=slump, hiatus=hiatus)
  }

#source("script_per_agemodel.R")
