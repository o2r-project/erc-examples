#' @get /figure1
#' @png (width = 700) 
function(newValue0) { 
startAnalysis <- Sys.time() 
newValue0=as.numeric(newValue0) 
library(ggplot2)
library("dplyr")
library("latticeExtra")
load("mjoBestUse.RData")
Intensity = newValue0
Tracks.df = filter(Tracks.df, Int >= Intensity) 
fullData.df = load("best.use.2014.Rdata")
begin = 1975; end = 2014
library("rgdal")
ll = "+proj=longlat +ellps=WGS84"

Tracks.sdf = Tracks.df
coordinates(Tracks.sdf) = c("lon", "lat")
proj4string(Tracks.sdf) = CRS(ll)
library("raster")
r = raster(ncol = 10, nrow = 5, 
           xmn = -100, xmx = -20, 
           ymn = 10, ymx = 50)

Tracks.grid = rasterize(Tracks.sdf, r,
                        field = 'DIntDt',
                        fun = mean)
sid1=1675
Example.sdf = subset(Tracks.sdf, 
                     Sid == sid1 | 
                     Sid == 1677 |
                     Sid == 1683)
Example.grid = rasterize(Example.sdf, r, 
                         field = 'DIntDt',
                         fun = mean)
library("RColorBrewer")
library("rasterVis")
range(values(Example.grid), na.rm = TRUE)
rng = seq(0, 1.5, 1.5)
cr = brewer.pal(9, "Purples")
cr = cr[-(1:3)]
vals = levelplot(Example.grid, margin = FALSE, 
          xlab = NULL, ylab = NULL, 
          col.regions = cr, at = rng, 
          colorkey = NULL,
          border = "white", border.lwd = 2, pretty=TRUE,
          par.settings = list(fontsize = list(text = 15)))
library(mapproj)
library(maptools)

outlines = as.data.frame(map("world", xlim = c(-100, -20), 
                             ylim = c(10, 50), 
                             plot = FALSE)[c("x", "y")],
                             color = "gray")
map = geom_path(aes(x, y), inherit.aes = FALSE, data = outlines, 
                alpha = .8, show_guide = FALSE, color = "blue")
ext = as.vector(extent(r))
boundaries = map("world", fill = TRUE, xlim = ext[1:2], 
                 ylim = ext[3:4], plot = FALSE)
IDs = sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <<- map2SpatialPolygons(boundaries, IDs = IDs,
                              proj4string = CRS(projection(r)))
require(gridExtra)
require(sp)

Int2a.df = filter(Tracks.df, Sid == sid1)
Int2b.df = filter(Tracks.df, Sid == 1677)
Int2c.df = filter(Tracks.df, Sid == 1683)
holdera = data.frame(Int2a.df$lon, Int2a.df$lat)
holderb = data.frame(Int2b.df$lon, Int2b.df$lat)
holderc = data.frame(Int2c.df$lon, Int2c.df$lat)
linesa <<- SpatialLines(list(Lines(list(Line(holdera)),
                                 ID = 'Int2a.df$Sid')))
linesb <<- SpatialLines(list(Lines(list(Line(holderb)),
                                 ID = 'Int2b.df$Sid')))
linesc <<- SpatialLines(list(Lines(list(Line(holderc)),
                                 ID = 'Int2c.df$Sid')))


tracking = vals + 
  latticeExtra::layer(sp.polygons(bPols, col = gray(.8))) +
  latticeExtra::layer(sp.lines(linesa, lwd = 2, col = gray(.4))) + 
  latticeExtra::layer(sp.lines(linesb, lwd = 2, col = gray(.4))) + 
  latticeExtra::layer(sp.lines(linesc, lwd = 2, col = gray(.4))) +
  latticeExtra::layer(panel.text(-63.5, 27.5, 'Gustav', col = gray(.4))) +
  latticeExtra::layer(panel.text(-51, 16.5, 'Ike', col = gray(.4))) + 
  latticeExtra::layer(panel.text(-81.5, 16.25, 'Omar', col = gray(.4)))
plot(tracking)

endAnalysis <- Sys.time() 
totaltime <- difftime(endAnalysis, startAnalysis, units="secs") 
message(paste("Total time is: ", totaltime)) 
}