#` # Reclassification and rasterisation(discretisation) of HaeanCover data
#` [doi:10.1594/PANGAEA.823677](doi:10.1594/PANGAEA.823677)
#` Corresponding to Bumsuk Seo ([Bumsuk.Seo@uni-bayreuth.de]Bumsuk.Seo@uni-bayreuth.de)

path.Pangaea <- "." # Set your Pangaea data directory here 
setwd(path.Pangaea)

lccs.filename <- "HaeanCover_Legend.xls" 
igbp.filename <- "HaeanCover_IGBP_lookup.xls"

library(rgdal)
library(xlsx)
library(RColorBrewer)
library(rgeos)
library(parallel)
library(raster)
library(lattice)

n.threads <- 1

plotting <- TRUE # To plot and save the figures. Plotting could take a few minutes.  
saving <- TRUE   # To save the reclassified data into ESRI Shape file format.


reclassLULC <- function(LULC, rcl, columns) {
    
    tar <- org <- LULC@data
    
    # tar[] <- NA
    
    rcl <- rcl[!is.na(rcl[,1]),]
    if (ncol(rcl) > 2) { 
        
        additionalInfo <- matrix(NA, nrow=nrow(tar), ncol=12)
        
        
        colnames(additionalInfo) <- c("Mgmt2009", "Mgmt2010", "Mgmt2011", "DblCrp2009",  "DblCrp2010",  "DblCrp2011", "Mixed2009", "Mixed2010","Mixed2011", "Note2009", "Note2010", "Note2011")
        
        columns.qa <- c("qa2009", "qa2010", "qa2011")
        
        #     str(tar[,columns.qa])
    }
    
    for (i in 1:length(columns)) { 
        
        for (j in 1:nrow(rcl)) {
            idx1 <- org[,columns[i]] == as.character(rcl$original.types[j])
            idx1[is.na(idx1)] <- FALSE
            
            tar[idx1, columns[i]] <- as.character(rcl$replacement.types[j])
            
            # additionalInfo[idx1, (i-1) + c(1,4,7)] <- as.character(rcl[j, 3:5]) 
            # Not working in this way 
            
            if (ncol(rcl) > 2) { 
                
                additionalInfo[idx1, (i-1) + 1] <- as.character(rcl[j, 3]) 
                additionalInfo[idx1, (i-1) + 4] <- as.character(rcl[j, 4]) 
                additionalInfo[idx1, (i-1) + 7] <- as.character(rcl[j, 5]) 
                additionalInfo[idx1, (i-1) + 10] <- as.character(rcl[j, 6]) 
                
                # Update QA info
                
                if (!is.na(rcl[j, 7])) { 
                    tar[idx1,columns.qa[i]] <- as.character(rcl[j, 7])
                    # print(tar[idx1,columns.qa[i]])
                    # stop()
                }
            }
        }
        
        # No data & NA
        idx.na <- is.na(org[,columns[i]])
        tar[idx.na, columns[i]] <- "NA" # Character 'NA'
        
        if (ncol(rcl) > 2) { 
            
            tar[idx.na, columns.qa[i]] <- as.character(rcl[rcl[,1]=="no data", 7]) # QA 
            
            
            additionalInfo[idx.na, (i-1) + 1] <- as.character(rcl[rcl[,1]=="no data", 3]) 
            additionalInfo[idx.na, (i-1) + 4] <- as.character(rcl[rcl[,1]=="no data", 4]) 
            additionalInfo[idx.na, (i-1) + 7] <- as.character(rcl[rcl[,1]=="no data", 5]) 
            additionalInfo[idx.na, (i-1) + 10] <- as.character(rcl[rcl[,1]=="no data", 6]) 
        }
        
    }
    if (ncol(rcl) > 2) { 
        
        LULC@data <- cbind(tar, additionalInfo, stringsAsFactors=F)
        
    } else {
        
        LULC@data <- cbind(tar, stringsAsFactors=F)
        
    }
    return(LULC)
}



# Calculate intersected areas between polygons and a raster grid

polygonArea <- function(sp1, sp2, ras, cols.df, n.threads=1) { 
    
    
    cols.occur <- names(table(sp1@data))
    # cols.occur <- as.character(sort(as.numeric(cols.occur)))
    
    polygonAreaInner <- function(i) { 
        
        frac <- rep(NA, length(cols.occur))
        pix <- sp2[i,]
        
        L.pix <- gIntersection(sp1, pix, byid=T, id= sp1[[1]])
        
        a <- gArea(L.pix, byid=T)
        names(a) <- names(L.pix) # repeated types
        prop <- by(a, FUN=sum, INDICES=names(L.pix)) / sum(a)
        frac[ match( unlist(names(prop)), cols.occur)] <- as.numeric(prop)        
        
        return(frac)
    }
    
    intersected <- gIntersects(sp1, sp2, byid = T, prepared = T)
    
    idx <- which(rowSums(intersected) != 0)
    
    result <- vector("list", length = length(idx))
    
    if (n.threads ==1) { 
        
        cat("pixel id ")
        for (i in 1:length(idx)) {
            cat(idx[i], ", ")
            result[[i]] <- polygonAreaInner(idx[i])
        }
        
    } else {
        result <- mclapply(idx, mc.cores=n.threads, mc.preschedule=F, FUN= polygonAreaInner)
    }
    
    stopifnot(length(result)==length(idx)) # mc.preschedule should be set to FALSE! Otherwise a problem. 
    
    result.m <- matrix(data = NA, nrow = length(sp2), ncol = length(cols.occur))
    result.m[idx,] <- do.call(rbind, result)
    
    frac.rs <- vector("list", length(cols.occur))
    
    for (i in 1:length(cols.occur)) { 
#         print(i)
        frac.rs[[i]] <- setValues(ras, result.m[,i])
    }
    
    frac.rs <- stack(frac.rs)
    names(frac.rs) <- cols.df[,2][match(cols.occur, cols.df[,1])]
    
    return(frac.rs)
    
}


# Decide a raster value based on fractions
maxType <- function(x.vec, NA.code, place=1, code=NULL, tol=1E-1) {
    
    if (diff(range(x.vec, na.rm=T)) < tol) { 
        print("No dominant type")
        res <- NA.code
    } else {
        
        idx <- order(x.vec, decreasing = T)[place]
        
        if(x.vec[idx] > tol) {
            
            if (is.null(code)) {
                res <- idx
            } else {
                res <- code[idx]
            }
            
        } else {
            print("No dominant type")
            res <- NA.code
        }
        
    }
    
    return(res)
}

 
#` proj4strings
proj4.LL <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" # WGS84 EPSG:4326 
proj4.UTM52N <- "+proj=utm +zone=52 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" # UTM52N (WGS84) EPSG:32652
proj4.MODIS <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" 


#` Read the original data
HaeanCover <- readOGR(dsn=".", layer="HaeanCover", verbose=T, dropNULLGeometries=T, stringsAsFactors=FALSE)

#` Boundaries 
boundary.UTM <- SpatialPolygonsDataFrame(gUnaryUnion(HaeanCover), data=data.frame(LULC="boundary"), match.ID=F)
boundary.SINUS <- spTransform(boundary.UTM, CRS(proj4.MODIS))


#` ## Haean Cover (S1)
#` Convert mixed LULC types to singular LULC types
preconversion.table <- read.xlsx(file=lccs.filename, sheetName = "correction")
colnames(preconversion.table) <- c("original.types", "replacement.types")
HaeanCover.pre <- reclassLULC(HaeanCover, rcl=preconversion.table, columns= c("LULC2009", "LULC2010", "LULC2011"))

haean.pre.levels <-sort( unique(c(HaeanCover.pre$LULC2009, HaeanCover.pre$LULC2010, HaeanCover.pre$LULC2010)))
haean.pre.levels <- c(haean.pre.levels[haean.pre.levels!="NA"], "NA")



#` Reclassification 
conversion.table.all <- read.xlsx(file=lccs.filename, sheetIndex="types")
conversion.table <- conversion.table.all[,c("LULC.class", "Presence.of.Vegetation", "Edaphic.Condition", "Artificiality.of.Cover", "Life.form", "Crop.Type", "Life.cycle", "Group")]


#` ## S2 Groups
conversion.table.S2 <- subset(conversion.table, select=c("LULC.class", "Group"))
colnames(conversion.table.S2) <- c("original.types", "replacement.types")

HaeanCover.S2 <- reclassLULC(HaeanCover.pre, rcl=conversion.table.S2, columns= c("LULC2009", "LULC2010", "LULC2011"))

HaeanCover.S2@data$LULC2009 <- factor(HaeanCover.S2@data$LULC2009)
HaeanCover.S2@data$LULC2010 <- factor(HaeanCover.S2@data$LULC2010)
HaeanCover.S2@data$LULC2011 <- factor(HaeanCover.S2@data$LULC2011)

#` ## FAO-LCCS 
conversion.table.LCCS <- subset(conversion.table, select=c("LULC.class", "Artificiality.of.Cover"))
colnames(conversion.table.LCCS) <- c("original.types", "replacement.types")

HaeanCover.LCCS <- reclassLULC(HaeanCover.pre, rcl=conversion.table.LCCS, columns= c("LULC2009", "LULC2010", "LULC2011"))

HaeanCover.LCCS@data$LULC2009 <- factor(HaeanCover.LCCS@data$LULC2009)
HaeanCover.LCCS@data$LULC2010 <- factor(HaeanCover.LCCS@data$LULC2010)
HaeanCover.LCCS@data$LULC2011 <- factor(HaeanCover.LCCS@data$LULC2011)

# [1] "Artificial Surfaces and Associated Area"
# [2] "Artificial Waterbodies, Snow and Ice"
# [3] "Bare Area"
# [4] "Cultivated and Managed Terrestrial Area"
# [5] "Cultivated Aquatic or Regularly Flooded Areas "
# [6] "Natural and Semi-Natural Aquatic or Regularly Flooded Vegetation"
# [7] "Natural and Semi-Natural Terrestrial Vegetation"
# [8] "Natural Waterbodies, Snow and Ice"

### Step 4: IGBP 17 classes
# IGBP

### 17 classes
conversion.tableIGBP <- read.xlsx(file=igbp.filename, sheetName="IGBPlookup")
colnames(conversion.tableIGBP) <- c("original.types", "replacement.types")

HaeanCover.IGBP <- reclassLULC(HaeanCover.pre, rcl=conversion.tableIGBP, columns= c("LULC2009", "LULC2010", "LULC2011"))

# unique(c(unlist(HaeanCover.IGBP[[1]]), HaeanCover.IGBP[[2]], HaeanCover.IGBP[[3]]))

igbp.desc <- read.xlsx(file= igbp.filename, sheetName="IGBP_desc", stringAsFactors=F)
igbp.colors <- rgb(igbp.desc[,-c(1,2)]/255) # color map http://duckwater.bu.edu/lc/igbp.txt
igbp.colors[16] <- "dimgrey" # Missing data
igbp.colors[19] <- "white" # Missing data

igbp.levels <- as.character(igbp.desc$IGBP.class)
igbp.labels <- igbp.desc$MODIS.12Q1.Class

HaeanCover.IGBP@data$LULC2009 <-  factor(HaeanCover.IGBP@data$LULC2009, levels=igbp.levels, labels= igbp.labels)
HaeanCover.IGBP@data$LULC2010 <-  factor(HaeanCover.IGBP@data$LULC2010, levels=igbp.levels, labels= igbp.labels)
HaeanCover.IGBP@data$LULC2011 <-  factor(HaeanCover.IGBP@data$LULC2011, levels=igbp.levels, labels= igbp.labels)

# MODIS12Q1Class    Code    Name
# 0     WAT    Water Bodies
# 1     ENF	Evergreen Needleleaf Forests 
# 2	    EBF	Evergreen Broadleaf Forests 
# 3	    DNF	Deciduous Needleleaf Forests
# 4	    DBF	Deciduous Broadleaf Forests
# 5	    MF	Mixed Forests
# 6	    CSH	Closed Shrublands
# 7	    OSH	Open Shrublands
# 8	    WSA	Woody Savannas
# 9	    SAV	Savannas
# 10	GRA	Grasslands
# 11	WET	Permanent Wetlands
# 12	CRO	Croplands
# 13	URB	Urban and Built-Up Lands
# 14	CVM	Cropland/Natural Vegetation Mosaics
# 15    SNO	Snow and Ice
# 16	BRN	Barren or Sparsely Vegetated
# 99	INT	Interrupted Areas
# 100	MD Missing Data

if (saving)  {
    writeOGR(HaeanCover.IGBP, dsn=".", "HaeanCover_IGBP", driver="ESRI Shapefile", overwrite_layer = T, verbose = T)                     
}



### MODIS comparison
 
Haean_MODISgrid <- raster(nrows=28, ncols=38, xmn= 11174639, xmx=11192245, ymn= 4249968, ymx=4262940, crs=proj4.MODIS, resolution = c(463.3127,463.3127), vals=100)

#' Mask pixels outside of the catchment
MODIS.r.500 <- mask(Haean_MODISgrid, boundary.SINUS)

MODIS.p.500 <- rasterToPolygons(MODIS.r.500, na.rm = F)
MODIS.p.500.UTM <- spTransform(MODIS.p.500, CRSobj = CRS(proj4.UTM52N))
 
 

HaeanCover.IGBP.SINUS <-spTransform(HaeanCover.IGBP, CRS(proj4.MODIS))

HaeanCover.IGBP.SINUS$LULC2009 <- as.character(HaeanCover.IGBP.SINUS$LULC2009)
HaeanCover.IGBP.SINUS$LULC2010 <- as.character(HaeanCover.IGBP.SINUS$LULC2010)
HaeanCover.IGBP.SINUS$LULC2011 <- as.character(HaeanCover.IGBP.SINUS$LULC2011)


igbp.levels <- as.character(igbp.desc$IGBP.class)
igbp.labels <- igbp.desc$MODIS.12Q1.Class
igbp.df <- cbind(igbp.labels, labs=igbp.labels)
 

# Calculate fractional cover
system.time(HaeanCover.IGBP.fraction.2009 <- polygonArea(HaeanCover.IGBP.SINUS[,"LULC2009"], MODIS.p.500, MODIS.r.500, cols.df =igbp.df, n.threads=n.threads))  # 500m/17 types took 177 secs with 7 cores 
  
system.time(HaeanCover.IGBP.fraction.2010 <- polygonArea(HaeanCover.IGBP.SINUS[,"LULC2010"], MODIS.p.500, MODIS.r.500, cols.df = igbp.df, n.threads=n.threads)) 

system.time(HaeanCover.IGBP.fraction.2011 <- polygonArea(HaeanCover.IGBP.SINUS[,"LULC2011"], MODIS.p.500, MODIS.r.500, cols.df = igbp.df, n.threads=n.threads))   
 
res.temp <- lapply(list(HaeanCover.IGBP.fraction.2009, HaeanCover.IGBP.fraction.2010, HaeanCover.IGBP.fraction.2011), FUN=function (x) { 
    xv <- getValues(x) 
    xv[is.na(x)@data@values] <- 0
    x <- setValues(x, xv)
    # rowSums(getValues(x))
    # x <- mask(x, boundary.SINUS)
})

LULC.IGBP.2009.fraction.m <- res.temp[[1]]@data@values
LULC.IGBP.2010.fraction.m <- res.temp[[2]]@data@values
LULC.IGBP.2011.fraction.m <- res.temp[[3]]@data@values

cols.occur <- lapply(list(LULC.IGBP.2009.fraction.m, LULC.IGBP.2010.fraction.m, LULC.IGBP.2011.fraction.m), FUN = function(x) substr(colnames(x), start = 2, stop = length(colnames(x))))

na.code <- 100

LULC.IGBP.2009.m <- as.numeric(apply(LULC.IGBP.2009.fraction.m, MARGIN=1, FUN=maxType, tol=0.1, NA.code = na.code, code=cols.occur[[1]]) )
LULC.IGBP.2010.m <- as.numeric(apply(LULC.IGBP.2010.fraction.m, MARGIN=1, FUN=maxType, tol=0.1, NA.code = na.code, code=cols.occur[[2]]) )
LULC.IGBP.2011.m <- as.numeric(apply(LULC.IGBP.2011.fraction.m, MARGIN=1, FUN=maxType,tol=0.1, NA.code = na.code, code= cols.occur[[3]]) )

table(LULC.IGBP.2009.m)
table(LULC.IGBP.2010.m)
table(LULC.IGBP.2011.m)
  
 
LULC.IGBP.poly <- MODIS.p.500.UTM

LULC.IGBP.poly@data <- data.frame(LULC2009=factor(LULC.IGBP.2009.m, levels=igbp.labels, labels= igbp.levels), LULC2010 = factor(LULC.IGBP.2010.m, levels=igbp.labels, labels= igbp.levels), LULC2011 = factor(LULC.IGBP.2011.m, levels=igbp.labels, labels= igbp.levels)) 

 
LULC.IGBP.poly.masked <- crop(LULC.IGBP.poly, boundary.UTM, byid=T)

 

plot.haean.igbp.vector.2009 <- spplot(HaeanCover.IGBP, "LULC2009", col.regions=igbp.colors, lty=0, edge.col = "transparent", colorkey=F)
plot.haean.igbp.vector.2010 <- spplot(HaeanCover.IGBP, "LULC2010", col.regions=igbp.colors, lty=0, edge.col = "transparent", colorkey=F)
plot.haean.igbp.vector.2011 <- spplot(HaeanCover.IGBP, "LULC2011", col.regions=igbp.colors, lty=0, edge.col = "transparent", colorkey=F)

plot.haean.igbp.sinus.2009 <- spplot(LULC.IGBP.poly.masked, "LULC2009", col.regions= igbp.colors, lty=0, edge.col = "transparent", colorkey=F) #, sub="(a) MODIS Land cover (2010)")
plot.haean.igbp.sinus.2010 <- spplot(LULC.IGBP.poly.masked, "LULC2010", col.regions= igbp.colors, lty=0, edge.col = "transparent", colorkey=F) #, sub="(a) MODIS Land cover (2010)")
plot.haean.igbp.sinus.2011 <- spplot(LULC.IGBP.poly.masked, "LULC2011", col.regions= igbp.colors, lty=0, edge.col = "transparent", colorkey=F) #, sub="(a) MODIS Land cover (2010)")

 
if (plotting) { 
    cairo_pdf(filename="Haean_IGBP17_HaeanCover.pdf", width=12, height=8, antialias="none")
    
    plot.new()
    trellis.par.set(axis.line=list(col=NA), )
    
    igbp.levels[15]<- "Cropland/Natural Vegetation\nMosaics"
    legend(x=0.8, y= 1, legend=igbp.levels, fill=igbp.colors, cex=1, box.lty=0, xpd=NA) # , text.width=0.2)
    
    text("(a)", x= -0.1, y=1.1, xpd=NA)
    text("(b)", x= 0.2, y=1.1,  xpd=NA)
    text("(c)", x= 0.5, y=1.1,  xpd=NA)
    text("(d)", x= -0.1, y=0.4, xpd=NA)
    text("(e)", x= 0.2, y=0.4,  xpd=NA)
    text("(f)", x= 0.5, y=0.4,  xpd=NA)
    
    print(plot.haean.igbp.vector.2009, split=c(1, 1, 4, 2), newpage=FALSE)
    print(plot.haean.igbp.vector.2010, split=c(2, 1, 4, 2), newpage=FALSE)
    print(plot.haean.igbp.vector.2011, split=c(3, 1, 4, 2), newpage=FALSE)
    
    print(plot.haean.igbp.sinus.2009, split=c(1, 2, 4, 2), newpage=FALSE)
    print(plot.haean.igbp.sinus.2010, split=c(2, 2, 4, 2), newpage=FALSE)
    print(plot.haean.igbp.sinus.2011, split=c(3, 2, 4, 2), newpage=FALSE)
    
    dev.off()
    
}


if (saving)  {
    writeOGR(HaeanCover.LCCS, dsn=".", "HaeanCover_LCCS", driver="ESRI Shapefile", overwrite_layer = T, verbose = T)              
    writeOGR(HaeanCover.S2, dsn=".", "HaeanCover_S2", driver="ESRI Shapefile", overwrite_layer = T, verbose = T)
    writeOGR(HaeanCover.IGBP, dsn=".", "HaeanCover_IGBP", driver="ESRI Shapefile", overwrite_layer = T, verbose = T)         
    writeOGR(LULC.IGBP.poly.masked, dsn=".", "HaeanCover_IGBP_discretized", driver="ESRI Shapefile", overwrite_layer = T, verbose = T)                
}