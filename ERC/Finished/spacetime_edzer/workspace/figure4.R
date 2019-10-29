library(spacetime)
if (!require("plm")) install.packages("plm", repos = "http://cloud.r-project.org")
if (!require("diveMove")) install.packages("diveMove", repos = "http://cloud.r-project.org")
if (!require("trip")) install.packages("trip", repos = "http://cloud.r-project.org")
if (!require("adehabitatLT")) install.packages("adehabitatLT", repos = "http://cloud.r-project.org")
if (!require("cshapes")) install.packages("cshapes", repos = "http://cloud.r-project.org")

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
set.seed(1331)

###################################################
sp = cbind(x = c(0,0,1), y = c(0,1,1))
row.names(sp) = paste("point", 1:nrow(sp), sep="")
library("sp")
sp = SpatialPoints(sp)

###################################################
library("gstat")
data("wind")
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84"

###################################################
wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
wind$jday = as.numeric(format(wind$time, '%j'))

###################################################
stations = 4:15
windsqrt = sqrt(0.5148 * as.matrix(wind[stations])) # knots -> m/s
Jday = 1:366
windsqrt = windsqrt - mean(windsqrt)
daymeans = sapply(split(windsqrt, wind$jday), mean)
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocities = apply(windsqrt, 2, function(x) { x - meanwind })

###################################################
wind.loc = wind.loc[match(names(wind[4:15]), wind.loc$Code),]
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
rownames(pts) = wind.loc$Station
pts = SpatialPoints(pts, CRS("+proj=longlat +datum=WGS84"))

###################################################
library("rgdal")
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84")
pts = spTransform(pts, utm29)


wind.data = stConstruct(velocities, space = list(values = 1:ncol(velocities)), 
                        time = wind$time, SpatialObj = pts, interval = TRUE)

library("lattice")
library(xts)

yrs = 1970:1986
time = as.POSIXct(paste(yrs, "-01-01", sep=""), tz = "GMT")

data("Produc")

library("RColorBrewer")
library("maptools")
m = map2SpatialLines(
  map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84"
m = spTransform(m, utm29)

grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),
                    proj4string = proj4string(m))

wind.data = wind.data[, "1961-04"]

n = 10
tgrd = xts(1:n, seq(min(index(wind.data)), max(index(wind.data)), length=n))
pred.grd = STF(grd, tgrd)

model="Gau"
v = vgmST("separable", space = vgm(0.6, method, 750000), time = vgm(1, method, 1.0 * 3600 * 24), sill=0.6)
wind.ST = krigeST(values ~ 1, wind.data, pred.grd, v)
colnames(wind.ST@data) <- "sqrt_speed"

layout = list(list("sp.lines", m, col='grey'),
              list("sp.points", pts, first=F, cex=.5))
stplot(wind.ST, col.regions=brewer.pal(11, "RdBu")[-c(10,11)],
       at=seq(-1.375,1,by=.25),
       par.strip.text = list(cex=.7), sp.layout = layout)
