###################################################################
###################################################################
#-----Landslides length and width estimation-v-01------------------
###################################################################
###################################################################
#
#
#    The self-explanatory R script performs the analysis of 
#    landslide flow direction using 2 steps:
# 1. Generate the minimum-area oriented bounding box for every 
#    landslide polygon, and the corner points and midpoints for 
#    the bbox sides; add altitudes for the midpoints and compute
#    the square of the difference raised to the power of two.
# 2. Crop the DEM and compute flow length; extract the maximum
#    slope length. 
#    The landslide width and length are then estimated based on
#    the oriented along the flow direction bounding box.
#    Please cite the article bellow if you use this script.
#    Niculita Mihai, 2016, Automatic landslide length and width 
#    estimation based on the geometric processing of the bounding
#    box and the geomorphometric analysis of DEMs, Natural Hazards 
#    and Earth System Science, 16, 1-10, doi:10.5194/nhess-16-1-2016
#     !!! ATENTION !!! slope length is a computational demanding
#    operation so if you can please insert the slope length as a 
#    column in the landslide inventory shapefile (slMFD field).
#    *please also add an id field.
#    *see the attached shapefiles for a reference on the script
#    variables and use them to format your data.
#    *modify accordingly for any change in the input data.
#
###################################################################
#------Set working directory---------------------------------------
#------Please edit the code for you working directory--------------
###################################################################
 

###################################################################
#------Install necessary packages----------------------------------
###################################################################
pkg=c("sp","rgdal","maptools","rgeos","raster","shotGroups","RSAGA",
      "reshape","data.table","stringr","pROC", "vcd", "caret")
install.packages(pkg, dependencies=T)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
library(shotGroups)
library(reshape)
library(data.table)
library(RSAGA)
library(stringr)
library(pROC)
library(vcd)
library(caret)
getOption("pROCProgress")
options("pROCProgress" = list(name = "text", width = NA, style = 3, char = "Step"))
getOption("pROCProgress")
###################################################################
#------see and set set the memory limit----------------------------
###################################################################
memory.limit()
#memory.limit(size=24493)
###################################################################
#------set SAGA environment----------------------------------------
###################################################################
#--ussually is best to move a SAGA binary package to RSAGA folder from the library
#myenv=rsaga.env(workspace="e:/Dropbox/Dropbox (Licente_Disertatii)/2016_NHESS",
#                path="C:/Users/MIHAI/Documents/R/win-library/3.2/RSAGA/SAGA-GIS",
#                modules="C:/Users/MIHAI/Documents/R/win-library/3.2/RSAGA/SAGA-GIS/modules",
#                version="2.2.3", cores=8)
###################################################################
#------read the shapefile------------------------------------------
###################################################################
shape <- readShapeSpatial("landslide_inventory", 
                          proj4string = CRS("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 +y_0=500000 +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,-0.053,0.844,-0.84 +units=m +no_defs"), 
                          repair=TRUE, force_ring=T,
                          verbose=TRUE)
###################################################################
#------create the list storing the minimum bounding boxes
###################################################################
coord.list=list()
mbbox.list=list()
bbox.list=list()
polygon.bbox.list=list()
polygons.bbox.list=list()
bbox.sorted.list=list()
ids=1:4
for (i in 1:length(shape))
{
  coord.list[[i]]=shape@polygons[[i]]@Polygons[[1]]@coords
  mbbox.list[[i]]=getMinBBox(coord.list[[i]])
  bbox.list[[i]]=data.frame(mbbox.list[[i]]$pts)
}
###################################################################
#------populating the polygon shapefile with the details of
#------the bounding boxes corners and midpoints
###################################################################
ps <- lapply(bbox.list, Polygon)
for (i in 1:length(shape))
{
  bbox.list[[i]]$ids=ids
  bbox.list[[i]]$poly_id=shape@data$id[[i]]
  shape@data$width[[i]]=mbbox.list[[i]]$width
  shape@data$height[[i]]=mbbox.list[[i]]$height
  shape@data$angle[[i]]=mbbox.list[[i]]$angle
  shape@data$p1x1y1[[i]]=bbox.list[[i]][1,1]
  shape@data$p1x1y2[[i]]=bbox.list[[i]][1,2]
  shape@data$p2x1y1[[i]]=bbox.list[[i]][2,1]
  shape@data$p2x1y2[[i]]=bbox.list[[i]][2,2]
  shape@data$p3x1y1[[i]]=bbox.list[[i]][3,1]
  shape@data$p3x1y2[[i]]=bbox.list[[i]][3,2]
  shape@data$p4x1y1[[i]]=bbox.list[[i]][4,1]
  shape@data$p4x1y2[[i]]=bbox.list[[i]][4,2]
  #
  shape@data$mdp1x1y1[[i]]=(bbox.list[[i]][1,1]+bbox.list[[i]][2,1])/2
  shape@data$mdp1x1y2[[i]]=(bbox.list[[i]][1,2]+bbox.list[[i]][2,2])/2
  shape@data$mdp2x1y1[[i]]=(bbox.list[[i]][2,1]+bbox.list[[i]][3,1])/2
  shape@data$mdp2x1y2[[i]]=(bbox.list[[i]][2,2]+bbox.list[[i]][3,2])/2
  shape@data$mdp3x1y1[[i]]=(bbox.list[[i]][3,1]+bbox.list[[i]][4,1])/2
  shape@data$mdp3x1y2[[i]]=(bbox.list[[i]][3,2]+bbox.list[[i]][4,2])/2
  shape@data$mdp4x1y1[[i]]=(bbox.list[[i]][4,1]+bbox.list[[i]][1,1])/2
  shape@data$mdp4x1y2[[i]]=(bbox.list[[i]][4,2]+bbox.list[[i]][1,2])/2
}
p1 <- lapply(seq_along(ps), function(i) Polygons(list(ps[[i]]), ID = shape@polygons[[i]]@ID[[1]]))
sp1=SpatialPolygons(p1, proj4string = CRS("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 
                                          +y_0=500000 +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,
                                          -0.053,0.844,-0.84 +units=m +no_defs"))
df=shape@data
SPDF <- SpatialPolygonsDataFrame(sp1, df) 
writeOGR(SPDF, ".", "landslide_inventory_mbbox", driver="ESRI Shapefile")
###################################################################
#------create points from bounding box corners---------------------
###################################################################
points=do.call("rbind", lapply(bbox.list, as.data.frame))
xy <- points[,c(1,2)]
spdf <- SpatialPointsDataFrame(coords = xy, data = points,
            proj4string = CRS("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 +y_0=500000 
                              +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,-0.053,0.844,-0.84 +units=m +no_defs"))
writeOGR(spdf, ".", "landslide_inventory_mbbox_points", driver="ESRI Shapefile")
###################################################################
#------create midpoints from bounding box sides--------------------
###################################################################
midpoints.list=list()
for (i in 1:length(shape))
{
  midpoints.list[[i]]=data.frame(((bbox.list[[i]][1,1]+bbox.list[[i]][2,1])/2),
                                 ((bbox.list[[i]][1,2]+bbox.list[[i]][2,2])/2))
  colnames(midpoints.list[[i]])=c("x","y")
  row2=c(((bbox.list[[i]][2,1]+bbox.list[[i]][3,1])/2),(bbox.list[[i]][2,2]+bbox.list[[i]][3,2])/2)
  row3=c(((bbox.list[[i]][3,1]+bbox.list[[i]][4,1])/2),(bbox.list[[i]][3,2]+bbox.list[[i]][4,2])/2)
  row4=c(((bbox.list[[i]][4,1]+bbox.list[[i]][1,1])/2),(bbox.list[[i]][4,2]+bbox.list[[i]][1,2])/2)
  midpoints.list[[i]]=rbind(midpoints.list[[i]], row2)
  midpoints.list[[i]]=rbind(midpoints.list[[i]], row3)
  midpoints.list[[i]]=rbind(midpoints.list[[i]], row4)
  midpoints.list[[i]]$ids=ids
  midpoints.list[[i]]$poly_id=shape@data$id[[i]]
}
#
midpoints=do.call("rbind", lapply(midpoints.list, as.data.frame))
xxyy <- midpoints[,c(1,2)]
midspdf <- SpatialPointsDataFrame(coords = xxyy, data = midpoints,
                               proj4string = CRS("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 +y_0=500000 
                              +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,-0.053,0.844,-0.84 +units=m +no_defs"))
writeOGR(midspdf, ".", "landslide_inventory_midpoints", driver="ESRI Shapefile")
###################################################################
#------add DEM values to corner points and midpoints--------------
###################################################################
#---try with SRTM
#dem=raster("dem_miletin_srtm.asc", crs=("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 
#                                        +y_0=500000 +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,
#                                   -0.053,0.844,-0.84 +units=m +no_defs"))
#---or with high resolution DEM
dem=raster("dem_miletin_srtm.asc", crs=("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 
                                        +y_0=500000 +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,
                                        -0.053,0.844,-0.84 +units=m +no_defs"))
zpoints=extract(dem, spdf, df=T, sp=T)
writeOGR(zpoints, ".", "landslide_inventory_mbbox_points_dem", driver="ESRI Shapefile")
midzpoints=extract(dem, midspdf, df=T, sp=T)
writeOGR(midzpoints, ".", "landslide_inventory_midpoints_dem", driver="ESRI Shapefile")
#
z=zpoints@data
shapef=cast(z, poly_id~ids)
zz=data.frame(shapef)
shape1=merge(shape, zz, by.x="id", by.y="poly_id", df=T, sp=T)
colnames(shape1@data)[which(names(shape1@data) == "X1")] <- "Zcnr1"
colnames(shape1@data)[which(names(shape1@data) == "X2")] <- "Zcnr2"
colnames(shape1@data)[which(names(shape1@data) == "X3")] <- "Zcnr3"
colnames(shape1@data)[which(names(shape1@data) == "X4")] <- "Zcnr4"
#
mz=midzpoints@data
mshapef=cast(mz, poly_id~ids)
zzz=data.frame(mshapef)
shape1=merge(shape1, zzz, by.x="id", by.y="poly_id", df=T, sp=T)
colnames(shape1@data)[which(names(shape1@data) == "X1")] <- "Zmidp1"
colnames(shape1@data)[which(names(shape1@data) == "X2")] <- "Zmidp2"
colnames(shape1@data)[which(names(shape1@data) == "X3")] <- "Zmidp3"
colnames(shape1@data)[which(names(shape1@data) == "X4")] <- "Zmidp4"
###################################################################
#------Compute the needed values-----------------------------------
###################################################################
shape1@data$dz12=sqrt((shape1@data$Zcnr1-shape1@data$Zcnr2)^2)
shape1@data$dz23=sqrt((shape1@data$Zcnr2-shape1@data$Zcnr3)^2)
shape1@data$dz34=sqrt((shape1@data$Zcnr3-shape1@data$Zcnr4)^2)
shape1@data$dz41=sqrt((shape1@data$Zcnr4-shape1@data$Zcnr1)^2)
shape1@data$L=sqrt((shape1@data$p1x1y1-shape1@data$p2x1y1)^2+(shape1@data$p1x1y2-shape1@data$p2x1y2)^2)
shape1@data$W=sqrt((shape1@data$p2x1y1-shape1@data$p3x1y1)^2+(shape1@data$p2x1y2-shape1@data$p3x1y2)^2)
shape1@data$elongated=ifelse(shape1@data$L>=shape1@data$W,"Yes","No")
shape1@data$p1x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$p1x1y1,shape1@data$p2x1y1)
shape1@data$p1x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$p1x1y2,shape1@data$p2x1y2)
shape1@data$p2x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$p2x1y1,shape1@data$p3x1y1)
shape1@data$p2x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$p2x1y2,shape1@data$p3x1y2)
shape1@data$p3x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$p3x1y1,shape1@data$p4x1y1)
shape1@data$p3x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$p3x1y2,shape1@data$p4x1y2)
shape1@data$p4x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$p4x1y1,shape1@data$p1x1y1)
shape1@data$p4x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$p4x1y2,shape1@data$p1x1y2)
shape1@data$oL=sqrt((shape1@data$p1x1y1-shape1@data$p2x1y1)^2+(shape1@data$p1x1y2-shape1@data$p2x1y2)^2)
shape1@data$oW=sqrt((shape1@data$p2x1y1-shape1@data$p3x1y1)^2+(shape1@data$p2x1y2-shape1@data$p3x1y2)^2)
shape1@data$mdp1x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp1x1y1,shape1@data$mdp2x1y1)
shape1@data$mdp1x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp1x1y2,shape1@data$mdp2x1y2)
shape1@data$mdp2x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp2x1y1,shape1@data$mdp3x1y1)
shape1@data$mdp2x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp2x1y2,shape1@data$mdp3x1y2)
shape1@data$mdp3x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp3x1y1,shape1@data$mdp4x1y1)
shape1@data$mdp3x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp3x1y2,shape1@data$mdp4x1y2)
shape1@data$mdp4x1y1=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp4x1y1,shape1@data$mdp1x1y1)
shape1@data$mdp4x1y2=ifelse(shape1@data$elongated=="Yes",shape1@data$mdp4x1y2,shape1@data$mdp1x1y2)
shape1@data$mL=sqrt((shape1@data$mdp2x1y1-shape1@data$mdp4x1y1)^2+(shape1@data$mdp2x1y2-shape1@data$mdp4x1y2)^2)
shape1@data$mW=sqrt((shape1@data$mdp1x1y1-shape1@data$mdp3x1y1)^2+(shape1@data$mdp1x1y2-shape1@data$mdp3x1y2)^2)
shape1@data$Zcnr1=ifelse(shape1@data$elongated=="Yes",shape1@data$Zcnr1,shape1@data$Zcnr2)
shape1@data$Zcnr2=ifelse(shape1@data$elongated=="Yes",shape1@data$Zcnr2,shape1@data$Zcnr3)
shape1@data$Zcnr3=ifelse(shape1@data$elongated=="Yes",shape1@data$Zcnr3,shape1@data$Zcnr4)
shape1@data$Zcnr4=ifelse(shape1@data$elongated=="Yes",shape1@data$Zcnr4,shape1@data$Zcnr1)
shape1@data$Zmidp1=ifelse(shape1@data$elongated=="Yes",shape1@data$Zmidp1,shape1@data$Zmidp2)
shape1@data$Zmidp2=ifelse(shape1@data$elongated=="Yes",shape1@data$Zmidp2,shape1@data$Zmidp3)
shape1@data$Zmidp3=ifelse(shape1@data$elongated=="Yes",shape1@data$Zmidp3,shape1@data$Zmidp4)
shape1@data$Zmidp4=ifelse(shape1@data$elongated=="Yes",shape1@data$Zmidp4,shape1@data$Zmidp1)
shape1@data$dZcnr12=sqrt((shape1@data$Zcnr1-shape1@data$Zcnr2)^2)
shape1@data$dZcnr43=sqrt((shape1@data$Zcnr4-shape1@data$Zcnr3)^2)
shape1@data$dZcnr14=sqrt((shape1@data$Zcnr1-shape1@data$Zcnr4)^2)
shape1@data$dZcnr23=sqrt((shape1@data$Zcnr2-shape1@data$Zcnr3)^2)
shape1@data$dZmdip13=sqrt((shape1@data$Zmidp1-shape1@data$Zmidp3)^2)
shape1@data$dZmdip24=sqrt((shape1@data$Zmidp2-shape1@data$Zmidp4)^2)
###################################################################
#--------------Step 1---Use midpoint altitude differences----------
###################################################################
#step1=1 long landslides; step1=0 wide landslides
shape1@data$step1=ifelse(shape1@data$dZmdip24>shape1@data$dZmdip13, "Long", "Wide")#if YES then LONG, if NOT then WIDE
shape1@data$step1=as.factor(shape1@data$step1)
summary(shape1@data$step1)
###################################################################
#--------------Step 2---Use slope length---------------------------
###################################################################
#--If the computer does not have enough resources or if the script is taking too long
#--there is the option to use the slope length which can be computed and attached as 
#--a value in SAGA/GRASS or other software; in our case the slope length is given in 
#--the base inventory as column slD8 and slMFD
#------Compute Slope length----------------------------------------
#--a vector having ids and slope length for any landslide polygon
#---extract the landslide polygons as separate shapes
#rsaga.geoprocessor("shapes_tools",9,env=myenv,list(					
#  SHAPES="landslides_miletin_s70.shp",
#  FIELD="id",
#  LIST="landslide",
#  NAMING=1))
#---cut the filled dem for every landslide polygon
#for (i in 1:length(shape))
#{
#  rsaga.geoprocessor("shapes_grid",7,env=myenv,list(					
#    OUTPUT=paste("dem_ldl_", sprintf("%04.0f", i), ".sgrd", sep=""),
#    INPUT="miletin_lidar5m_filled.sgrd",
#    POLYGONS=paste("landslide_", sprintf("%04.0f", i), ".shp", sep=""),
#    NODATA=0))
#}
#---compute slope length for every dem
#for (i in 1:length(shape))
#{
#  rsaga.geoprocessor("ta_hydrology",6,env=myenv,list(					
#    ELEVATION=paste("dem_ldl_", sprintf("%04.0f", i), ".sgrd", sep=""),
#    LENGTH=paste("dem_ldl_slD8_", sprintf("%04.0f", i), ".sgrd", sep=""),
#    METHOD=0))
#}
#---merge polygons
#---!!!!!! ATENTION !!!!!!
#---this does not work in RSAGA so it does need to be manually done
#---!!!!!! ATENTION !!!!!!
#sd="e:/Dropbox/Dropbox (Licente_Disertatii)/2016_NHESS/shp"
#shp.files=list.files(sd, ".shp", full.names=F)
#shp.files=paste(shp.files, sep=",", collapse=",")
#shp.files=str_replace_all(shp.files, ".shp,", ",\")
#shp.files=paste("/", shp.files, sep="")
#
#rsaga.geoprocessor("shapes_tools", 2,env=myenv,list(					
#    INPUT=shp.files,
#    MERGED="landslide_database_merged.shp",
#    SRCINFO=0))
#--read the shape with shape lengths
#---!!!!!! ATENTION !!!!!!
#--this is slower in RSAGA so it can be does separatelly in SAGA
#---!!!!!! ATENTION !!!!!!
#slength=raster("slope_length_ldl.asc", crs=("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 
#                                           +x_0=500000 +y_0=500000 +ellps=krass +towgs84=33.4,
#                                          -146.6,-76.3,-0.359,-0.053,0.844,-0.84 +units=m +no_defs"))
#shape2=extract(slength, shape1, fun=max, df=T, sp=T)
#shape1@data$sl=NULL
#shape_length <- readShapeSpatial("miletin_landslides_s70_slope_length.shp", 
#                          proj4string = CRS("+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 +y_0=500000 +ellps=krass +towgs84=33.4,-146.6,-76.3,-0.359,-0.053,0.844,-0.84 +units=m +no_defs"), 
#                          repair=TRUE, force_ring=T,
#                          verbose=TRUE)
#shape1=merge(shape1, shape_length, by.x="id",df=T, sp=T)
#shape1@data[is.na(shape1@data)]=0
#---D8
#shape1@data$step2=ifelse(shape1@data$step1=="Long", ifelse(((shape1@data$mL)*1)>=shape1@data$slD8, "Wide", "Long"), "Wide")
#shape1@data$step2=as.factor(shape1@data$step2)
#summary(shape1@data$step2)
#writeOGR(shape1, ".", "miletin_ldl_database_D8", driver="ESRI Shapefile")
#---MFD
shape1@data$step2=ifelse(shape1@data$step1=="Long", ifelse(((shape1@data$mL)*1)>=shape1@data$slMFD, "Wide", "Long"), "Wide")
shape1@data$step2=as.factor(shape1@data$step2)
summary(shape1@data$step2)
writeOGR(shape1, ".", "landslide_inventory_long_wide_result.shp", driver="ESRI Shapefile")
#----------------------------SCRIPT END----------------------------
###################################################################
#--------------Validation------------------------------------------
###################################################################
validation1=data.frame(shape1@data$step1,shape1@data$real_type)
validation2=data.frame(shape1@data$step2,shape1@data$real_type)
ff1=ftable(validation1)
ff2=ftable(validation2)
pdf("fourfoldplot_step1.pdf")
fourfoldplot(ff1,col=c("blue","red"),conf.level = 0.95)
dev.off()
pdf("fourfoldplot_step2.pdf")
fourfoldplot(ff2,col=c("blue","red"),conf.level = 0.95)
dev.off()
#
shape1@data$step1n=ifelse(shape1@data$step1=="Long",1,2)
shape1@data$step2n=ifelse(shape1@data$step2=="Long",1,2)
shape1@data$real_typen=ifelse(shape1@data$real_type=="Long",1,2)
shape1@data$check1=ifelse(shape1@data$step1==shape1@data$real_type, "Right", "Wrong")
shape1@data$check1_type=ifelse(shape1@data$check1=="Wrong",shape1@data$step1,"OK")
shape1@data$check2=ifelse(shape1@data$step2==shape1@data$real_type, "Right", "Wrong")
shape1@data$check2_type=ifelse(shape1@data$check2=="Wrong",shape1@data$step2,"OK")
shape1@data$real_typef=shape1@data$real_type
#
pdf("ROC_step1.pdf")
rocobj1 <- plot.roc(shape1@data$step1n, shape1@data$real_typen, percent=TRUE,  ci=TRUE, of="se", 
                   specificities=seq(0, 100, 5), ci.type="shape", ci.col='#1c61b6AA')
dev.off()
pdf("ROC_step2.pdf")
rocobj2 <- plot.roc(shape1@data$step2n, shape1@data$real_typen, percent=TRUE,  ci=TRUE, of="se", 
                   specificities=seq(0, 100, 5), ci.type="shape", ci.col='#1c61b6AA')
dev.off()
#
tb1=table(validation1)
tb2=table(validation2)
cm1=confusionMatrix(tb1)
cm2=confusionMatrix(tb2)
cm1
cm2
