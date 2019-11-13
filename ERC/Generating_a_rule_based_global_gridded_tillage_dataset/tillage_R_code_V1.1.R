############################################################################################################
#(C) 2018, Potsdam Institute for Climate Impact Research (PIK)
#Contact: vera.porwollik@pik-potsdam.de
#Authors, and contributors <see METADATA>
#Licensed under MIT <see README.txt file>
#Citation: Porwollik, Vera; Rolinski, Susanne; Müller, Christoph (2019): A global gridded data set on tillage - R-code (V. 1.1). GFZ Data Services. http://doi.org/10.5880/PIK.2019.010
#Accompanying tillage dataset: Porwollik, Vera; Rolinski, Susanne; Müller, Christoph (2019): A global gridded data set on tillage. V. 1.1. GFZ Data Services. http://doi.org/10.5880/PIK.2019.009
#Supplement to article: Porwollik, V., Rolinski, S., Heinke, J., and Müller, C.: Generating a rule-based global gridded tillage dataset, Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2018-152, in review, 2018
############################################################################################################
#Input data sets used:
#a) physical cropland (ha) from SPAM2005; (GeoTIFF - file format)
#International Food Policy Research Institute (IFPRI); International Institute For Applied Systems Analysis (IIASA). (2017). <spam2005V3r1_global_phys_area.geotiff.zip> [Data set]. Harvard Dataverse. https://doi.org/10.7910/dvn/dhxbjx/k5hvukand 
#b) grid cell allocation key to countries from SPAM2005; (.txt- file format)
#International Food Policy Research Institute (IFPRI); International Institute For Applied Systems Analysis (IIASA). (2017). cell5m_allockey_xy.dbf.zip [Data set]. Harvard Dataverse. https://doi.org/10.7910/dvn/dhxbjx/lvrjlf
#c) Conservation Agriculture area for all countries, around year 2005; (.csv - file format)
#FAO (2016). Conservation Agriculture. AQUASTAT. Main Database Food and Agriculture Organization of the United Nations (FAO), URL: http://www.fao.org/nr/water/aquastat/data/query/index.html
#d) Aridity index (1961-1990); (.adf - file format)
#FAO: FAO GEONETWORK (2016) Global map of aridity - 10 arc minutes (GeoLayer). FAO, Rome, Italy. URL: http://www.fao.org/geonetwork/srv/en/main.home?uuid=221072ae-2090-48a1-be6f-5a88f061431a
#e) Field sizes; (.img - file format)
#Fritz, S., See, L., McCallum, I., You, L., Bun, A., Moltchanova, E., . Obersteiner, M. (2015). Mapping global cropland and field size. Global Change Biology, 21(5), 1980-1992. doi:10.1111/gcb.12838
#f) Soilsgrid - Absolute depth to bedrock (in cm); (GeoTIFF - file format)
#Hengl, T., de Jesus, J. M., MacMillan, R. A., Batjes, N. H., Heuvelink, G. B. M., Ribeiro, ... Gonzalez, M. R. (2014). SoilGrids1km - Global soil information based on automated mapping , College of Global Change and Earth System Science, Beijing Normal University/ISRIC - World Soil Information, PLoS ONE, 9, e105992, https://doi.org/10.1371/journal.pone.0105992.
#g) GLADIS - Water erosion; (.adf - file format)
#Nachtergaele, F. O., Petri, M., Biancalani, R., van Lynden, G., and van Velthuizen, H. (2011). Global Land Degradation Information System (GLADIS). An information database for land degradation assessment at global level. Technical report of the LADA FAO/UNEP Project. http://www.fao.org/fileadmin/templates/solaw/files/thematic_reports/SOLAW_thematic_report_3_land_degradation.pdf
#h) World Bank - Level of income as GNI per country, year 2005; (.xls - file format)
#World Bank (2017). World Development indicators- historical classification by income URL: https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups
#for sample input calculation use sections 1 to 10 for comparison to sample results per tillage system area (sample_calc)
##################################################################################################

####1.Pre-calculation steps####
#create directory and path to folder for input data
#path_input <- ".../input_data/"                                 
#download input data sets stated above and adjust data formats (e.g. convert csv or xls to txt file format)

#create directory and path to ouput data
#path_calc <- ".../sample_output/"

#load packages
require(raster)
require(fields)
require(ncdf4)

####2.Harmonization of input datasets####
sample_calc <- T                     #set to TRUE for sample input and output calculation else set to FALSE

if (sample_calc){
  final.ext <- extent(-180,180,-45,90)
  r <- raster(res=45,ext=final.ext)
  origin(r) <- 0
} else {
  final.ext <- extent(-180,180,-56,84) #geopraphic coverage of dataset
  r <- raster(res=1/12,ext=final.ext)  #resolution to 0.083333...°
  r[] <- NA
}

###SPAM2005 physical cropland in (ha)###
if (sample_calc) {
  brick_crop <- brick(nc=8,nr=3,xmn=-180,xmx=180,ymn=-45,ymx=90,nl=42,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  #sample SPAM_TA
  numb <- rep(1:24,times=42)
  SPAM_TA <- brick_crop
  SPAM_TA[] <- numb
  #sample SPAM_TR
  numb_tr <- numb-1/3
  SPAM_TR <- brick_crop
  SPAM_TR[] <- numb_tr
} else {
  spam_prod_level <- c("TA","TI","TR","TH","TL","TS")
  file_fin <- c("A","I","R","H","L","S")
  crop_vec <- c("WHEA","RICE","MAIZ","BARL","REST","OOIL","TOBA","TEAS","COCO","RCOF","ACOF","OFIB","COTT","SUGB","SUGC","OILP","VEGE","TEMF","TROF","PLNT","BANA","CNUT","GROU","OTRS","CASS","YAMS","SWPO","POTA","SESA","RAPE","SUNF","SOYB","OPUL","LENT","PIGE","COWP","CHIC","BEAN","OCER","SORG","SMIL","PMIL")

  for(ll in 1:length(spam_prod_level)){
    for(cc in 1:length(crop_vec)){
      cropx <- raster(paste0(path_input,"SPAM2005V3r1_global_A_",spam_prod_level[ll],"_",crop_vec[cc],"_",file_fin[ll],".tif"))
      if(cc==1){
      SPAM_X <- cropx
      } else {
      SPAM_X <- stack(SPAM_X,cropx)
      }
    }
  origin(SPAM_X) <- 0
  SPAM_X <- crop(SPAM_X,final.ext)
  SPAM_X <- assign(paste0("SPAM_",spam_prod_level[ll]),crop(SPAM_X,final.ext))
  readAll(SPAM_X)                                               
  save(list=paste0("SPAM_",spam_prod_level[ll]),file=paste0(path_calc,"SPAM_",c_spam_prod_level[ll],".Rdata"))
  rm(list=paste0("SPAM_",spam_prod_level[ll])) 
  }
}
#load(paste0(path_calc,"SPAM_",spam_prod_level[ll]),".Rdata")

#for calculation with other physical cropland (landuse and spatial crop area distribution), adjust here the dimensions of the resulting raster stack

###SoilGrids Absolute depth to bedrock (in cm)###
if (sample_calc){
  #sample depth to bedrock
  num <- c(1:24)
  soilgrid_res <- r
  soilgrid_res[] <- num
} else {
  soilgrid_rast <- raster(paste0(path_input,"r_BDTICM_M_10km_ll.tif")) #here resolution and extend need change
  soilgrid_res <- resample(soilgrid_rast,r,method="bilinear")
}
origin(soilgrid_res) <- 0
save(soilgrid_res,file=paste0(path_calc,"soilgrid_res.Rdata"))
#load(paste0(path_calc,"soilgrid_res.Rdata"))

flat_area <- soilgrid_res
flat <- flat_area<15                                             #every grid cell with depth less than 15 cm as too flat for deeper tillage
flat_area[!flat] <- NA
origin(flat_area) <- 0
save(flat_area,file=paste0(path_calc,"flat_area.Rdata")) 

if (sample_calc){
  #sample fields
  b <- c(10:33)
  fields_interpol <- r
  fields_interpol[] <- b
} else {
  fiel_rast <- raster(paste0(path_input,"field_size_10_40_cropland.img"))
  fiel_rast_aggr <- aggregate(fiel_rast,fact=10,fun=modal,na.rm=T) #decrease resolution and find most frequent value of ordinal variable
  fields <- extend(sizes.in_1,final.ext)
  origin(fields) <- 0
  #save(fields,file=paste0(path_calc,"Fritz_fields.Rdata"))

  #interpolate missing field size for grid cells where SPAM2005 reports cropland
  spam_ta_no_field_info <- mask(SPAM_TA,fields,inverse=T)         #42 stack without field size
  spam_ta_no_field_info_sum <- sum(spam_ta_no_field_info,na.rm=T) #0, 8586.4  (min, max)
  origin(spam_ta_no_field_info_sum) <- 0

  spam_ta_no_field_info_sum[spam_ta_no_field_info_sum==0] <- NA  
  fields[is.na(fields)] <- 0
  fie <- mask(fields,spam_ta_no_field_info_sum,inverse=T) 
  #length(r[is.na(fie)])                                           #43.711 of spam_sum cells are having no field size info 

  fill.na <- function(x,i=18241) {                                #for 191x191 moving window
    if(is.na(x)[i]) {             
      return(round(mean(x[x!=0],na.rm=TRUE),0))
    } else return(x[i])
  }  
  #for the following step: "Pass the fill.na function to raster::focal and check results. The pad argument creates virtual rows  /columns of NA values to keep the vector length constant along the edges of the raster. This is why we can always expect, e.g.   the fifth value of the vector to be the focal value in a 3x3 window thus, the index i=5 in the fill.na function."
  fie_1 <- focal(fie,w=matrix(1,191,191),fun=fill.na,pad=TRUE,na.rm=FALSE)
  fie_2 <- focal(fie_1,w=matrix(1,191,191),fun=fill.na,pad=TRUE,na.rm=FALSE)
  fields_interpol <- focal(fie_2,w=matrix(1,191,191),fun=fill.na,pad=TRUE,na.rm=FALSE)
  fields_interpol[is.na(fields_interpol)] <- 20  #set logit-neutral value to grid cell of Hawaian Islands, because of too large distance to other reporting cells
  fields_interpol[fields_interpol==0] <- NA
}
origin(fields_interpol) <- 0
save(fields_interpol,file=paste0(path_calc,"fields_interpol.Rdata"))

###GLADIS: Water erosion###
if (sample_calc) {
  #sample gladis
  nu <- c(1:24)
  gladis_hdr <- r
  gladis_hdr[] <- nu
} else {
  gladis_adf <- raster(paste0(path_input,"hdr.adf")) #0, 12110.03 (min,max)
  crs(gladis_adf) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  origin(gladis_adf) <- 0
  gladis_adf_ext <- extend(gladis_adf,final.ext)   
  gladis_crop <- crop(gladis_adf_ext,final.ext)
  gladis_ext <- setExtent(gladis_crop,final.ext,keepres=TRUE) 
  gladis_hdr <- setExtent(gladis_ext,final.ext)
}
origin(gladis_hdr) <- 0
save(gladis_hdr,file=paste0(path_calc,"gladis_hdr.Rdata"))

###Aridity index### 
#classification original dataset: Hyperarid AI < 0.05, Arid 0.05 < AI < 0.20, Semi-arid 0.20 < AI < 0.50, Sub-humid 0.50 < AI < 0.65
if (sample_calc){
  #sample aridity
  n <- c(1:23,NA)
  aridity_res <- r
  aridity_res[] <- n
} else {
  ori_aridity <- raster(paste0(path_input,"w001001.adf"))               #is 10 arc min
  aridity_crop <- crop(ori_aridity,final.ext)
  disag_ori_aridity <- disaggregate(aridity_crop,fact= 0.166667/(1/12)) #here resolution results in 0.0833335 only
  aridity_res <- resample(disag_ori_aridity,r,method="bilinear")        #0, 10.478  (min, max)
}
origin(aridity_res) <- 0
save(aridity_res,file=paste0(path_calc,"aridity_res.Rdata"))       #has few zeros, NAs for ocean

###allocation of grid cells to countries### 
if (sample_calc) {
  #sample alloc_rast - in case of sample calculation, take entire cropland as one country
  alloc_rast <- r
  alloc_rast[] <- 240
} else {  
  allokey <- read.delim(file=paste0(path_input,"lut_cell5m_iso3_allockey.txt"))
  allo_iso <- allokey[,c(3,4,2)]           #lon, lat, ISO3 letter code for countries
  allo_ord <- allo_iso[order(allo_iso$ISO3),]
  levels(allo_ord$ISO3) <- as.numeric(c(1:239)) #here levels are stored as factors but need to be numeric for rasterize, levels are automatically sorted alphabetically 
  d <- as.numeric(new$ISO3)
  new <- cbind(new[,1:2],d)
  alloc_extent <- extent(-179.9583,179.9583,-89.95833,83.62500)
  rr <- raster(res=1/12,ext=alloc_extent)
  alloc_raster <- rasterize(new[,1:2],rr,field=new[,3]) 
  alloc_crop <- crop(alloc_raster,alloc_extent)
  alloc_rast <- extend(alloc_crop,final.ext)
}
origin(alloc_rast) <- 0
save(alloc_rast,file=paste0(path_calc,"spam_alloc_rast.Rdata"))  #1:239 ISO3 codes from grid cell allocation file

###From allocation raster to country income level mask###
if (sample_calc){
  #sample income
  high_raster <- alloc_rast  #all grid cells get high income level
} else {
  alloc_levels <- read.delim(file=paste0(path_input,"spam_alloc_country_code_raster.txt"),header=T) #file generated based on World Bank income levels for 2005, 239 observations of three variables: 1.) ISO3 code, 2.) value in alloc_rast c(1:239), and 3.) sorted acc. to income levels 1-4 as low, lower-middle, upper-middle, and high income 
  high <- alloc_levels[c(113:205),2]
  high_raster <- alloc_rast
  high_raster[high_raster%in%high] <- 400
  high_raster[!high_raster%in%c(400)] <- NA 
}
save(high_raster,file=paste0(path_calc,"high_raster.Rdata"))

####3.Logit_model for calculating CA likelihood with four spatial predictor variables####
###3a) Crop_mix, as ratio of 22 CA suitable crop area over total sum of cropland per grid cell###
#Get area suitable for downscaling per country per grid cell: for combi with logit model result
SPAM_TR_flat_out <- mask(SPAM_TR,flat_area,inverse=T)
fields_interpol[fields_interpol>=20] <- NA        #set NA where fields are larger than 2 ha
small <- mask(SPAM_TR_flat_out,fields_interpol)   #42 crops,  extract crop specific ha on small fields
high_income_small_area <- mask(small,high_raster)
origin(high_income_small_area) <- 0
save(high_income_small_area,file=paste0(path_calc,"high_income_small_area.Rdata")) #42 stack of small field sizes cropland area in high and upper-middle income countries
#gc()
annuals <- subset(SPAM_TR_flat_out,c(1,2,3,4,5,7,13,14,17,c(24:42),23)) #29 stack
load(paste0(path_calc,"fields_interpol.Rdata")) 
large <- fields_interpol>=20                                            #here all cells with fields larger 2 ha 
fields_interpol[!large] <- NA 
annuals_large <- mask(annuals,fields_interpol)
annuals_small_high_income <- subset(high_income_small_area,c(1,2,3,4,5,7,13,14,17,c(24:42),23)) #29 stack
for (i in 1:29) {
  print(paste("building annuals_all_stack:",i))
  annuals_stack <- stack(annuals_large[[i]],annuals_small_high_income[[i]])
  annuals_sum <- sum(annuals_stack,na.rm=T)
  if (i > 1) annuals_all_stack <- stack(annuals_all_stack,annuals_sum) else annuals_all_stack <- annuals_sum
}
names(annuals_all_stack) <- c("whea","rice","maiz","barl","rest","toba","cott","sugb","vege","orts","cass","yams","swpo","pota","sesa","rape","sunf","soyb","opul","lent","pige","cowp","chic","bean","ocer","sorg","smil","pmil","grou")
grains_stack <- subset(annuals_all_stack,c(1,3,4,5,6,7,9,15:29)) #22 (excluding 7 crop types:rice,sugb,cass,yams,swpo,pota,orts)
grains_sum <- sum(grains_stack,na.rm=T)
grains_sum[grains_sum==0] <- NA
sum_spam_ta <- sum(SPAM_TA,na.rm=T)
sum_spam_ta[sum_spam_ta==0] <- NA
crop_mix <- grains_sum/sum_spam_ta
crop_mix[!is.finite(crop_mix)] <- NA  #check values for NANs
origin(crop_mix) <- 0
save(crop_mix,file=paste0(path_calc,"crop_mix.Rdata")) 
#gc()

###3b) interpolate erosion for the correct amount of values###
gladis <- mask(gladis_hdr,crop_mix)  #set NA to erosion where no cropland is reported
b <- mask(crop_mix,gladis,inverse=T) 
#length(b[is.na(b)])                  #7.252.586
gladis[is.na(gladis)] <- 999
erosion <- mask(gladis,b,inverse=T)  #set NA to cells, where crop_mix has values but gladis not
#length(erosion[is.na(erosion)])      #5014 to fill with a value
erosion[is.na(erosion)] <- 12        #setting missing erosion values to 12, where cropland is reported but NA in original gladis
erosion[erosion==999] <- NA          #length 0s: 4320, length is.na=6.815.668
save(erosion,file=paste0(path_calc,"erosion.Rdata")) 

###3c) interpolate aridity###
if (sample_calc){
  aridity_res[is.na(aridity_res)] <- 0.65
  aridity <- aridity_res
} else {
  m_aridity <- mask(aridity_res,crop_mix)
  b <- mask(crop_mix,m_aridity,inverse=T) 
  m_aridity[is.na(m_aridity)] <- 999
  ari <- mask(m_aridity,b,inverse=T) 
  #length(ari[is.na(ari)])               #1230 of crop_mix cells are having no aridity info
  fill.na <- function(x,i=1301) {      
    if(is.na(x)[i]) {           
      return(round(mean(x[x!=999],na.rm=TRUE),10))
    } else return(x[i])
  }  
  aridity <- focal(ari,w=matrix(1,51,51),fun=fill.na,pad=TRUE,na.rm=FALSE)
  aridity[is.na(aridity)] <- 0.65   #as a neutral value for islands too far away for interpolation
  aridity[aridity==999] <- NA
  #length(aridity[is.na(aridity)])   #6.815.668 same as target amount of NAs from crop_mix
}
save(aridity,file=paste0(path_calc,"aridity.Rdata"))

###3d) field_size###
field_size <- mask(fields_interpol,crop_mix) 
save(field_size,file=paste0(path_calc,"field_size.Rdata"))

#here calculate Spearmann's rank correlation coefficient (r) between input variables
eros_vec <- as.vector(getValues(erosion))      
fields_vec <- as.vector(getValues(field_size))
aridity_vec <- as.vector(getValues(aridity))
crop_mix_vec <- as.vector(getValues(crop_mix))

field_eros <- cor(fields_vec,eros_vec,method='spearman',use='complete.obs') #complete.obs: discard the entire row if an NA occurs
field_arid <- cor(fields_vec,aridity_vec,method='spearman',use='complete.obs')
field_crop_mix <- cor(fields_vec,crop_mix_vec,method='spearman',use='complete.obs')
eros_arid <- cor(eros_vec,aridity_vec,method='spearman',use='complete.obs')
eros_crop_mix <- cor(eros_vec,crop_mix_vec,method='spearman',use='complete.obs')
arid_crop_mix <- cor(aridity_vec,crop_mix_vec,method='spearman',use='complete.obs')

cor_vec <- as.matrix(c(field_eros,field_arid,field_crop_mix,eros_arid,eros_crop_mix,arid_crop_mix))
rownames(cor_vec) <- c('field_eros','field_arid','field_crop_mix','eros_arid','eros_crop_mix','arid_crop_mix')
colnames(cor_vec) <- 'correlation coefficient'

#fn <- paste0(path_calc,"tab_cor_input_vars.txt")
#write.table(x=cor_vec,file=fn,sep="\t")  #generates Table 3

###3e) Logit model and sensitivity test###
#CA area values of 54 countries around year 2005
if (sample_calc){
  #sample CA area
  ca_FAO <- data.frame("sample_calc",NA,10000,240) 
} else {
  ca_FAO <- read.delim(file=paste0(path_input,"max_CA_countries_around_2005_based_on_all_countries.txt"),header=T) 
}

ca_countries <- ca_FAO[,4]      #numeric ISO3 country codes in raster
ca_country_names <- ca_FAO[,1]  #names of countries

xmid <- c(20,12,0.65,0.5)       #thresholds based on literature findings used as mid-points of functions
kvalue <- c(1/4,1/60,-5,10)     #reference slopes 
sens_range <- c(1,99,999,2,0.5) #modifying factors to slopes of reference functions and variable combinations
sens_names <- c("ref","drop","only","slope_100_more","slope_50_less") #varying slope and variable combinations
k_names <- c("field_k","erosion_k","aridity_k","crop_mix_k")

mat_cor <- matrix(NA,nrow=length(kvalue),ncol=length(sens_range)-1)
rownames(mat_cor) <- k_names
colnames(mat_cor) <- sens_names[-1]

arr_country_cor <- array(data=NA,dim=c(length(kvalue),length(sens_range)-1,length(ca_countries)),dimnames=list(k_names,sens_names[-1],ca_country_names))

for(se in 1:length(sens_range)){
  print(paste("sensitivity test:",se))
  for(kk in 1:length(kvalue)){
    act_k <- kvalue
    act_xmid <- xmid
    if(se==1 & kk!=1){
      next                                      #skip to next kk
    } else {
      if(se==2){
        act_k[kk] <- 0                          #each k_run gets a zero for dropping one variable
        act_xmid[kk] <- 0  
      } else if(se==3) {                        #each run sets three vars to zero and keeps only one input layer
        act_k[act_k!=act_k[kk]] <- 0
        act_xmid[act_xmid!=act_xmid[kk]] <- 0
      } else {
        act_k[kk] <- kvalue[kk]*sens_range[se]  #here loop over all sensitivity slope variations
      }
    }
    fac_1 <- act_k[1]*field_size
    fac_2 <- act_k[2]*erosion 
    fac_3 <- act_k[3]*aridity
    fac_4 <- act_k[4]*crop_mix
    b_stack <- stack(fac_1,fac_2,fac_3,fac_4)
    b <- sum(b_stack,na.rm=F)        
    b[is.infinite(b)] <- NA                     
    f <- -act_xmid[1]*act_k[1]-act_xmid[2]*act_k[2]-act_xmid[3]*act_k[3]-act_xmid[4]*act_k[4]   
    d <- b+f  
    e <- d*(-1)    
    logit <- 1/(1+exp(e))  
    logit[!is.finite(logit)] <- NA    
    if(se==1 & kk==1){
      assign(paste0("logit_",sens_names[se]),logit)
      save(list=paste0("logit_",sens_names[se]),file=paste0(path_calc,"logit_",sens_names[se],".Rdata"))
      logit_old_c <- as.vector(getValues(logit))
    } else {
      assign(paste0("logit_",sens_names[se],"_",k_names[kk]),logit)
      save(list=paste0("logit_",sens_names[se],"_",k_names[kk]),file=paste0(path_calc,"logit_",sens_names[se],"_",k_names[kk],".Rdata"))
      rm(list=paste0("logit_",sens_names[se],"_",k_names[kk]))
      logit_new_c <- as.vector(getValues(logit))
      mat_cor[kk,se-1] <- round(cor(logit_new_c,logit_old_c,method='spearman',use='complete.obs'),digits=7)
      for (ii in 1:length(ca_countries)) {    
        print(paste("country:",ii))
        alloc <- alloc_rast==ca_countries[ii]   #e.g. Argentinia cells length=40468
        alloc[alloc==0] <- NA                   #set to NA for following mask command 
        ref_cut <- mask(get(paste0("logit_",sens_names[1])),alloc)   #ref_logit...
        sens_cut <- mask(logit,alloc) 
        ref_vec <- as.vector(getValues(ref_cut))
        sens_vec <- as.vector(getValues(sens_cut))
        arr_country_cor[kk,se-1,ii] <- round(cor(ref_vec,sens_vec,method='spearman',use='complete.obs'),digits=7)
        mat_1 <- arr_country_cor[,,ii]
        if (ii>1) mat_all_countries <- as.matrix(rbind(mat_all_countries,mat_1)) else mat_all_countries <- mat_1
      }
    }
    gc()
  }  
}
rownames(mat_all_countries) <- rep(ca_country_names,each=4)

#fn <- paste0(path_calc,"tab_cor_logit_sensitivities.txt")  #Table 4 in article
#write.table(x=mat_cor,file=fn,sep="\t")

#fn <- paste0(path_calc,"tab_cor_logit_sens_countries.txt") #Table S6 in SI of article
#write.table(x=mat_all_countries,file=fn,sep="\t") 

####4.Conservation Agriculture####
#downscale algorithm of reported national CA area values closest to year 2005 for 54 countries
#load(paste0(path_calc,logit_ref,".Rdata"))
ca_countries <- ca_FAO[,4]                        #ISO3 country code
area_FAO <- ca_FAO[,3]                            #reported CA area per country

pot_ca_area <- mask(grains_stack,logit_ref)       #stack of 22, same subset as for crop_mix above under 3.
pot_ca_area[pot_ca_area==0] <- NA
origin(pot_ca_area) <- 0
names(pot_ca_area) <- c("whea","maiz","barl","rest","toba","cott","vege","sesa","rape","sunf","soyb","opul","lent","pige","cowp","chic","bean","ocer","sorg","smil","pmil","grou")
#readAll(pot_ca_area)
#save(pot_ca_area,file=paste0(path_calc,"pot_ca_area.Rdata")) 

new_area <- sum(pot_ca_area,na.rm=T)    
new_b <- stack(logit_ref,new_area)      
df_all_countries <- NULL

for (i in 1:length(ca_countries)) {               #loop over countries
  gc()
  print(paste("country:",i))
  alloc <- alloc_rast==ca_countries[i]            
  alloc[alloc==0] <- NA                           #set rest of world to NA for following mask command 
  d <- mask(new_b,alloc)                          #creates 2 layer stack
  e <- rasterToPoints(d)                          #creates list object
  n <- as.data.frame(e)                           #df of 4 variables per grid cell: longitude, latitude, logit_probability, CA area
  colnames(n)[c(3,4)] <- c("prob","area")
  both_sets_order <- n[order(-n[,3]),]            #order probability index decreasing 
  thresholdValue <- area_FAO[i]                   #set target downscale ha value for each of the counties from ca_FAO
  print(paste("tresholdvalue:",thresholdValue))
  if(sum(both_sets_order[,4])>thresholdValue){
    io <- length(which(cumsum(both_sets_order[,4])<=thresholdValue)) #length of cell index until target value is reached in area_column
    jo <- sum(both_sets_order[,4][1:io],na.rm=T)  #indicates cum_sum of area_col until treshold is reached
    print(paste("sum jo:",jo))
    ao <- length(1:(io+1))                        #evaluating potential lower deviation from original treshold with one further cell with cropland fraction
    joa <- sum(both_sets_order[,4][1:ao],na.rm=T) 
    a <- abs(jo-thresholdValue)                   #result of difference can be '=' or '>' 0 
    v <- abs(joa-thresholdValue) 
    both_sets_order$CA <- -999
    if (a<v){
      both_sets_order$CA[1:io] <- 999
      print(paste("CA area downscaled is slightly smaller than tresholdvalue"))
    } else { 
      both_sets_order$CA[1:ao] <- 999
      print(paste("CA area downscaled is slightly larger than tresholdvalue"))
    }
    ca_crops <- both_sets_order[which(both_sets_order[,5]==999),c(1,2,4)] #subset df with 3 columns as x,y,area
  } else {
    sum_remaining_area <- sum(both_sets_order[,4],na.rm=T) 
    print(paste("area remaining:",sum_remaining_area))                    #states which amount is in fact available
    both_sets_order$CA  <- -999
    b <- both_sets_order$area>0
    both_sets_order$CA[which(b)] <- 999
    ca_crops <- both_sets_order[which(both_sets_order[,5]==999),c(1,2,4)] 
  }
  df_all_countries <- rbind(df_all_countries,ca_crops) 
}

ca_mask <- rasterize(df_all_countries[,1:2],r,field=df_all_countries[,3]) #one layer
origin(ca_mask) <- 0
Conservation_Agriculture <- mask(pot_ca_area,ca_mask) #puts NA to Conservation_Agriculture where ca_mask has NAs, stack of 22 crops
names(Conservation_Agriculture) <- c("whea","maiz","barl","rest","toba","cott","vege","sesa","rape","sunf","soyb","opul", "lent","pige","cowp","chic","bean","ocer","sorg","smil","pmil","grou")
origin(Conservation_Agriculture) <- 0
save(Conservation_Agriculture,file=paste0(path_calc,"Conservation_Agriculture.Rdata"))    #22 crop-specific area stack
sum_Conservation_Agriculture <- sum(Conservation_Agriculture,na.rm=T)                     #sum over 22 crops to one layer
sum_Conservation_Agriculture_sum <- cellStats(sum_Conservation_Agriculture,'sum',na.rm=T) #110.189.921 ha for ~2005, 54 countries

tab_ca_crops <- cbind(as.data.frame(cellStats(pot_ca_area,'sum',na.rm=T)),as.data.frame(cellStats(Conservation_Agriculture,'sum',na.rm=T)))
tab_ca_crops$share_down_on_pot <- round(tab_ca_crops[,2]/tab_ca_crops[,1],digits=2)*100

#fn <- paste0(path_calc,"tab_ca_crops.txt")     
#write.table(x=tab_ca_crops,file=fn,sep="\t")  #Table 5 in article

#finding out how much area is in fact downscaled per country
#total of 110.449.988 ha sum reported, e.g. downscaled area for Korea only: 2.390,1 instead of 23.000 and for New Zealand  and 78.517,8 instead of 162.000 ha available
ca_countries <- ca_FAO[,4] #ISO code of ountries
d_df <- NULL

for (i in 1:length(ca_countries)) {
  print(paste("country:",i))
  alloc <- alloc_rast==ca_countries[i]            #e.g. Argentinia raster length=40468
  alloc[alloc==0] <- NA                           #set to NA for following mask command 
  d <- mask(sum_Conservation_Agriculture,alloc)   #creates stack with 2 layer 
  d_sum <- cellStats(d,'sum',na.rm=T)
  d_df <- rbind(d_df,d_sum)
}
rownames(d_df) <- ca_FAO[,1]
#fn <- paste0(path_calc,"tab_ca_downscaled.txt")
#write.table(x=d_df,file=fn,sep="\t")             #Table S7 in SI of article

####5.Traditional annual tillage####
#rainfed and irrigated annual crops on small fields and of low or lower-middle income countries
SPAM_TA_flat_out <- mask(SPAM_TA,flat_area,inverse=T)
load(paste0(path_calc,"fields_interpol.Rdata")) 
fields_interpol[fields_interpol>=20] <- NA          #mask out areas where fields are equal or larger than 2 ha
small <- mask(SPAM_TA_flat_out,fields_interpol) #42 crops,extract crop-specific area on small fields
small_low_income <- mask(small,high_raster,inverse=T)                     #mask out areas with small fields in high income countries which presumably also are mechanized or some other sort of higher input production
traditional_annual_tillage <- subset(small_low_income,c(1,2,3,4,5,7,13,14,17,c(24:42),23)) #29 stack subsetting for annuals 
origin(traditional_annual_tillage) <- 0
readAll(traditional_annual_tillage) 
save(traditional_annual_tillage,file=paste0(path_calc,"traditional_annual_tillage.Rdata")) 
sum_traditional_annual_tillage <- sum(traditional_annual_tillage,na.rm=T)                       #one layer, sum over 29 annual crops
sum_traditional_annual_tillage_sum <- cellStats(sum_traditional_annual_tillage,'sum',na.rm=T)   #401.527.885 ha, for Table 6 in article

####6.Traditional rotational tillage####
#rainfed and irrigated perennial crops on small fields and in low or lower-middle income countries
traditional_rotational_tillage <- subset(small_low_income,c(6,8,9,10,11,12,15,16,18,19,20,21,22)) #13 perm crops rainfed and irrigated
origin(traditional_rotational_tillage) <- 0
traditional_rotational_tillage[traditional_rotational_tillage==0] <- NA
names(traditional_rotational_tillage) <- c("ooil","teas","coco","rcof","acof","ofib","sugc","oilp","temf","trof","plnt","bana","cnut")
save(traditional_rotational_tillage,file=paste0(path_calc,"traditional_rotational_tillage.Rdata")) 
sum_traditional_rotational_tillage <- sum(traditional_rotational_tillage,na.rm=T) 
sum_traditional_rotational_tillage_sum <- cellStats(sum_traditional_rotational_tillage,'sum',na.rm=T)   #65.050.925 ha, for table 6 in paper

####7.Reduced tillage####
#on cropland, where soils generally are shallower than 15 cm depth to bedrock
pre_reduced_till <- mask(SPAM_TA,flat_area)                                   #set NA to cropland where flat_area has NA
names(pre_reduced_till) <- names(SPAM_TA)
origin(pre_reduced_till) <- 0
save(pre_reduced_till,file=paste0(path_calc,"pre_reduced_till.Rdata"))
#sum_redu_till <- cellStats(pre_reduced_till,'sum',na.rm=T)
#sum_redu_till_sum <- sum(sum_redu_till,na.rm=T)                           

####8.Rotational tillage#### 
#for rainfed and irrigated perennials on large fields (mechanized)
#gc()
SPAM_TA_flat_out <- mask(SPAM_TA,flat_area,inverse=T)
perms <- subset(SPAM_TA_flat_out,c(6,8,9,10,11,12,15,16,18,19,20,21,22))  #13 stack
load(paste0(path_calc,"fields_interpol.Rdata"))
large <- fields_interpol>=20 
fields_interpol[!large] <- NA                                             #cells, where fields larger than 2 ha
perm_large <- mask(perms,fields_interpol)

#irrigated and rainfed perennials on small fields but in high income countries
#gc()
load(paste0(path_calc,"fields_interpol.Rdata"))      #reload 
small <- fields_interpol<20 
fields_interpol[!small] <- NA                           #cells, where fields smaller than 2 ha
perm_small <- mask(perms,fields_interpol)
perm_small_high_income <- mask(perm_small,high_raster)  #subset small field area in high income countries

for (i in 1:13) {                                       #adding rainfed and irrigated cropland in cells with large fields and the small field fraction of high and upper-middle income countries
  print(paste("From 13 adding permanent crop number:",i))
  crop <- stack(perm_large[[i]],perm_small_high_income[[i]])
  b <- sum(crop,na.rm=T)
  if (i > 1) pre_rot_till <- stack(pre_rot_till,b) else pre_rot_till <- b
}

load(paste0(path_calc,"soilgrid_res.Rdata"))
soilgrid_res[!soilgrid_res<20] <- NA
flat_red_rot <- soilgrid_res
origin(flat_red_rot) <- 0    #15, 19.99583  (min, max)

rotational_tillage <- mask(pre_rot_till,flat_red_rot,inverse=T)
rotational_tillage[rotational_tillage==0] <- NA
names(rotational_tillage) <- c("ooil","teas","coco","rcof","acof","ofib","sugc","oilp","temf","trof","plnt","bana","cnut")
origin(rotational_tillage) <- 0
save(rotational_tillage,file=paste0(path_calc,"rotational_tillage.Rdata")) 
sum_rotational_tillage <- sum(rotational_tillage,na.rm=TRUE)  
sum_rotational_tillage_sum <- cellStats(sum_rotational_tillage,'sum',na.rm=T) #74.179.844 ha

####9.Conventional annual tillage####
#rm(list=ls(all=T));gc()
#filling up stack and sort order of crop-specific layers
sum_spam_ta <- sum(SPAM_TA,na.rm=T)
sum_spam_ta_sum <- cellStats(sum_spam_ta,'sum',na.rm=T) #1.131.438.612 ha total cropland or 12.600 for sample dataset
gc()

r[] <- NA
#sum perennials
for (i in 1:13) {
  print(paste("From 13 total crops adding perennial crop:",i))
  perm_pre_stack <- stack(rotational_tillage[[i]],traditional_rotational_tillage[[i]])  
  perm_pre_stack_sum <- sum(perm_pre_stack,na.rm=T)
  if (i>1) perms_stack <- stack(perms_stack,perm_pre_stack_sum) else perms_stack <- perm_pre_stack_sum
}
names(perms_stack) <- names(rotational_tillage)

perms_sum_ext <- stack(r,r,r,r,r,perms_stack[[1]],r,perms_stack[[2]],perms_stack[[3]],perms_stack[[4]],perms_stack[[5]],perms_stack[[6]],r,r,perms_stack[[7]],perms_stack[[8]],r,perms_stack[[9]],perms_stack[[10]],perms_stack[[11]],perms_stack[[12]],perms_stack[[13]],r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r) #now 42

#gc()
traditional_annual_tillage_ext <- stack(traditional_annual_tillage[[1]],traditional_annual_tillage[[2]],traditional_annual_tillage[[3]],traditional_annual_tillage[[4]],traditional_annual_tillage[[5]],r,traditional_annual_tillage[[6]],r,r,r,r,r,traditional_annual_tillage[[7]],traditional_annual_tillage[[8]],r,r,traditional_annual_tillage[[9]],r,r,r,r,r,traditional_annual_tillage[[29]],traditional_annual_tillage[[10]],traditional_annual_tillage[[11]],traditional_annual_tillage[[12]],traditional_annual_tillage[[13]],traditional_annual_tillage[[14]],traditional_annual_tillage[[15]],traditional_annual_tillage[[16]],traditional_annual_tillage[[17]],traditional_annual_tillage[[18]],traditional_annual_tillage[[19]],traditional_annual_tillage[[20]],traditional_annual_tillage[[21]],traditional_annual_tillage[[22]],traditional_annual_tillage[[23]],traditional_annual_tillage[[24]],traditional_annual_tillage[[25]],traditional_annual_tillage[[26]],traditional_annual_tillage[[27]],traditional_annual_tillage[[28]]) #42

ca_down <- stack(Conservation_Agriculture[[1]],r,Conservation_Agriculture[[2]],Conservation_Agriculture[[3]],Conservation_Agriculture[[4]],r,Conservation_Agriculture[[5]],r,r,r,r,r,Conservation_Agriculture[[6]],r,r,r,Conservation_Agriculture[[7]],r,r,r,r,r,Conservation_Agriculture[[22]],r,r,r,r,r,Conservation_Agriculture[[8]],Conservation_Agriculture[[9]],Conservation_Agriculture[[10]],Conservation_Agriculture[[11]],Conservation_Agriculture[[12]],Conservation_Agriculture[[13]],Conservation_Agriculture[[14]],Conservation_Agriculture[[15]],Conservation_Agriculture[[16]],Conservation_Agriculture[[17]],Conservation_Agriculture[[18]],Conservation_Agriculture[[19]],Conservation_Agriculture[[20]],Conservation_Agriculture[[21]]) #now 42

#gc()
perms_sum_1 <- -1*perms_sum_ext
redu_till_1 <- -1*pre_reduced_till
traditional_annual_tillage_1 <- -1*traditional_annual_tillage_ext
ca_down_1 <- -1*ca_down

origin(SPAM_TA) <- 0
origin(perms_sum_1) <- 0
origin(redu_till_1) <- 0
origin(traditional_annual_tillage_1) <- 0
origin(ca_down_1) <- 0      

gc()
for (i in 1:42) {
  gc()
  print(paste("From 42 total crops subtracting for conventional tillage crop:",i))
  pre_stack <- stack(SPAM_TA[[i]],perms_sum_1[[i]],redu_till_1[[i]],traditional_annual_tillage_1[[i]],ca_down_1[[i]])
  sum_pre_stack <- sum(pre_stack,na.rm=T)
  sum_pre_stack[sum_pre_stack==0] <- NA
  if (i>1) pre_conv_till <- stack(pre_conv_till,sum_pre_stack) else pre_conv_till <- sum_pre_stack
}
#names(pre_conv_till) <- names(SPAM_TA)
#sum_pre_conv_till <- sum(pre_conv_till,na.rm=T)
#sum_pre_conv_till_sum <- cellStats(sum_pre_conv_till,'sum',na.rm=T) #465.630.785 ha

#Separating area of conventional tillage, which has depth to bedrock as deep as 15 but shallower than 20 cm
load(paste0(path_calc,"soilgrid_res.Rdata"))
still_flat <- soilgrid_res
flat <- still_flat >=15 & still_flat<20 
still_flat[!flat] <- NA
save(still_flat,file=paste0(path_calc,"still_flat.Rdata"))

mech_shallow <- mask(pre_conv_till,still_flat) #cropland is too flat for 20 cm tillage so is passed to reduced tillage
save(mech_shallow,file=paste0(path_calc,"mech_shallow.Rdata"))
#sum_mech_shallow <- sum(mech_shallow,na.rm=T)
#length(sum_mech_shallow[sum_mech_shallow>0])                      #365 cells have values larger than 0
#sum_mech_shallow_sum <- cellStats(sum_mech_shallow,'sum',na.rm=T) #592.923,5 ha  

conventional_annual_tillage <- mask(pre_conv_till,still_flat,inverse=T)
#names(conventional_annual_tillage) <- names(SPAM_TA)
#conventional_annual_tillage[conventional_annual_tillage==0] <- NA
#readAll(conventional_annual_tillage)
save(conventional_annual_tillage,file=paste0(path_calc,"conventional_annual_tillage.Rdata"))                  
sum_conventional_annual_tillage <- sum(conventional_annual_tillage,na.rm=T)     
sum_conventional_annual_tillage_sum <- cellStats(sum_conventional_annual_tillage,'sum',na.rm=T) #465.037.862 ha

####Adding last portions to reduced tillage####
redu_rot <- mask(pre_rot_till,flat_red_rot)                #portion of 15-20 cm depth from perennials to reduced tillage
#save(redu_rot,file=paste0(path_calc,"redu_rot.Rdata"))
sum_redu_rot <- sum(redu_rot,na.rm=T)
sum_redu_rot_sum <- cellStats(sum_redu_rot,'sum',na.rm=T)
sum_redu_rot_sum   #32.418,5

redu_rot_ext <- stack(r,r,r,r,r,redu_rot[[1]],r,redu_rot[[2]],redu_rot[[3]],redu_rot[[4]],redu_rot[[5]],redu_rot[[6]],r,r,redu_rot[[7]],redu_rot[[8]],r,redu_rot[[9]],redu_rot[[10]],redu_rot[[11]],redu_rot[[12]],redu_rot[[13]],r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r)

for (i in 1:42){
  print(paste("stacking and summing crop:",i))
  d <- stack(pre_reduced_till[[i]],mech_shallow[[i]],redu_rot_ext[[i]])
  d_sum <- sum(d,na.rm=T)
  d_sum[d_sum==0] <- NA
  if (i>1) reduced_tillage <- stack(reduced_tillage,d_sum) else reduced_tillage <- d_sum
}
#names(reduced_tillage) <- names(SPAM_TA)
origin(reduced_tillage) <- 0
save(reduced_tillage,file=paste0(path_calc,"reduced_tillage.Rdata"))
sum_reduced_tillage <- sum(reduced_tillage,na.rm=T)
sum_reduced_tillage_sum <- cellStats(sum_reduced_tillage,'sum',na.rm=T) #15.440.283 ha

####9.Scenario Conservation Agriculture area####
#get large fields size and CA suitable cropland under reduced tillage
pre_reduced_till[pre_reduced_till==0] <- NA
#getting large fields in low income
redu_ca_suit <- subset(pre_reduced_till,c(1,3,4,5,7,13,17,29,30,31,32,33,34,35,36,37,38,39,40,41,42,23)) #subset of 22 CA annual crops
low_redu_ca_suit <- mask(redu_ca_suit,high_raster,inverse=T)      #stack of 22; get area of low or lower-middle income countries with large field sizes
load(paste0(path_calc,"fields_interpol.Rdata"))
large <- fields_interpol >= 20 
fields_interpol[!large] <- NA 
low_income_large_redu <- mask(low_redu_ca_suit,fields_interpol)   #sum: 636.797,1 ha
high_income_redu <- mask(redu_ca_suit,high_raster) #for all cropland area with all field sizes in high income countries; 482.652,3 ha

for (i in 1:22){
  print(paste("stacking crop:",i))
  redu_stack <- stack(high_income_redu[[i]],low_income_large_redu[[i]])
  redu_stack_sum <- sum(redu_stack,na.rm=T)
  if (i>1) pot_ca_redu <- stack(pot_ca_redu,redu_stack_sum) else pot_ca_redu <- redu_stack_sum
}

for (i in 1:22) {             #adding the CA suitable area per crop from reduced tillage to the stack of potential CA
  print(paste("building stack:",i))
  stack_redu_pot <- stack(pot_ca_redu[[i]],pot_ca_area[[i]])
  stack_redu_pot_sum <- sum(stack_redu_pot,na.rm=T)                 
  if (i>1) scenario_ca_area <- stack(scenario_ca_area,stack_redu_pot_sum) else scenario_ca_area <- stack_redu_pot_sum
}
#names(scenario_ca_area) <- names(Conservation_Agriculture)
save(scenario_ca_area,file=paste0(path_calc,"scenario_ca_area.Rdata")) 
sum_scenario_ca_area <- cellStats(scenario_ca_area,'sum',na.rm=T)
sum_scenario_ca_area_sum <- sum(sum_scenario_ca_area,na.rm=T)          #466.996.460 ha

####10.Total sums of tillage system areas####
mat_till_types_area <- as.matrix(c(sum_conventional_annual_tillage_sum,sum_rotational_tillage_sum,sum_reduced_tillage_sum,sum_traditional_rotational_tillage_sum,sum_traditional_annual_tillage_sum,sum_Conservation_Agriculture_sum,sum_scenario_ca_area_sum))
rownames(mat_till_types_area) <- c('conventional_annual_tillage','rotational_tillage','reduced_tillage','traditional_rotational_tillage','traditional_annual_tillage','CA_downscaled','Scenario_CA_area')
mat_till_types_area

#fn <- paste0(path_calc,"tab_tillage_types_area.txt")
#write.table(x=mat_till_types_area,file=fn,sep="\t")  #Table 6 in article

#for comparison of sample data total output results per tillage area in hectares (ha):
#conventional_annual_tillage:     806.6667
#rotational_tillage:             1430.0000
#reduced_tillage:                7251.6667
#traditional_rotational_tillage:    0.0000
#traditional_annual_tillage:        0.0000
#CA_downscaled:                  4216.6667
#Scenario_CA_area:                    6526.6667

#for tillage system areas per country 
allokey <- read.delim(file=paste0(path_input,"lut_cell5m_iso3_allockey.txt"))
allo_iso <- allokey[,c(2,5)]
b <- unique(allo_iso)
coun_names <- b

df_all_countries <- NULL

for (i in 1:239){
  print(paste("doing country:",i))
  gc()
  coun_alloc <- alloc_rast==i #depicts area of country 
  coun_alloc[coun_alloc==0] <- NA 
  cro_mask <- mask(sum_spam_ta,coun_alloc)
  cropland_sum <- cellStats(cro_mask,'sum',na.rm=T)  
  co_till <- mask(sum_conventional_annual_tillage,coun_alloc)
  co_till_sum <- cellStats(co_till,'sum',na.rm=T)
  re_till <- mask(sum_reduced_tillage,coun_alloc)
  re_till_sum <- cellStats(re_till,'sum',na.rm=T)
  n_till <- mask(sum_Conservation_Agriculture,coun_alloc)
  n_till_sum <- cellStats(n_till,'sum',na.rm=T)
  ro_till <- mask(sum_rotational_tillage,coun_alloc)
  ro_till_sum <- cellStats(ro_till,'sum',na.rm=T)
  tra_ro_till <- mask(sum_traditional_rotational_tillage,coun_alloc)
  tra_ro_till_sum <- cellStats(tra_ro_till,'sum',na.rm=T)
  trad_ann_till <- mask(sum_traditional_annual_tillage,coun_alloc)
  trad_ann_till_sum <- cellStats(trad_ann_till,'sum',na.rm=T)
  scen_no_till <- mask(sum_scenario_ca_area,coun_alloc)
  scen_no_till_sum <- cellStats(scen_no_till,'sum',na.rm=T)
  till_area <- cbind(cropland_sum,co_till_sum,n_till_sum,re_till_sum,ro_till_sum,tra_ro_till_sum,trad_ann_till_sum,scen_no_till_sum)
  df_all_countries <- rbind(df_all_countries,till_area)
}
tillage_countries <- as.data.frame(df_all_countries)
rownames(tillage_countries) <- coun_names
#fn <- paste0(path_calc,"tillage_per_country.txt")
#write.table(x=tillage_countries,file=fn,sep="\t")  #Table S9 in SI of article

####11.Area weighted means of input variables used in logit for each of the tillage types####
#gc()
sum_pot_ca_area <- sum(pot_ca_area,na.rm=T) #reduce to one layer only

vari_stack <- stack(sum_spam_ta,sum_pot_ca_area,sum_Conservation_Agriculture,sum_traditional_annual_tillage,sum_traditional_rotational_tillage,sum_rotational_tillage,sum_reduced_tillage,sum_conventional_annual_tillage,sum_scenario_ca_area)

load(paste0(path_calc,"aridity_res.Rdata"))
load(paste0(path_calc,"fields_interpol.Rdata"))
load(paste0(path_calc,"crop_mix.Rdata"))
load(paste0(path_calc,"gladis_hdr.Rdata")) 

#getting crop_mix of 22 annuals on large fields in low income or all field size in high income of reduced tillage system area
pot_ca_redu_sum <- sum(pot_ca_redu,na.rm=T) #get from reduced tillage above
pot_ca_redu_sum[pot_ca_redu_sum==0] <- NA
spam_ta_sum <- sum(SPAM_TA,na.rm=T)
spam_ta_sum[spam_ta_sum==0] <- NA 
crop_mix_redu <- round(pot_ca_redu_sum/spam_ta_sum,digits=2)  #0.02, 1  (min, max)
crop_mix_redu[!is.finite(crop_mix_redu)] <- NA 

#getting crop_mix for scenario layer combining crop_mix (cells > 0: 441.574) and crop_mix_redu 
crop_mix_scen <- merge(crop_mix,crop_mix_redu) 

#loop over tillage area subsets
for (i in 1:9){
  e <- cellStats(((mask(aridity,vari_stack[i]))*vari_stack[i]),'sum',na.rm=T)/cellStats(vari_stack[i],'sum',na.rm=T)
  f <- cellStats(((mask(fields_interpol,vari_stack[i]))*vari_stack[i]),'sum',na.rm=T)/cellStats(vari_stack[i],'sum',na.rm=T)
  if (i=9) g <- cellStats(((mask(crop_mix_scen,vari_stack[i]))*vari_stack[i]),'sum',na.rm=T)/cellStats(vari_stack[i],'sum',na.rm=T)
  else g <- cellStats(((mask(crop_mix,vari_stack[i]))*vari_stack[i]),'sum',na.rm=T)/cellStats(vari_stack[i],'sum',na.rm=T)
  h <- cellStats(((mask(gladis_hdr,vari_stack[i]))*vari_stack[i]),'sum',na.rm=T)/cellStats(vari_stack[i],'sum',na.rm=T)
  mat <- cbind(e,f,g,h)
  if (i>1) fin_mat <- rbind(fin_mat,mat) else fin_mat <- mat
}

#fn <- paste0(path_calc,"area_weighted_means_input_vars_tillage_system_areas.txt")
#write.table(x=fin_mat,file=fn,sep="\t")  #Table S8 of article

####12. Comparison of derived tillage systems to area of two tillage intensities by Erb et al, 2016####
#unit: hectares
low_int <- (sum_reduced_tillage_sum+sum_Conservation_Agriculture_sum+sum_traditional_rotational_tillage_sum+sum_rotational_tillage_sum)
high_int <- (sum_traditional_annual_tillage_sum+sum_conventional_annual_tillage_sum)

low_int_share <- (low_int/sum_spam_ta_sum_km2)*100      #shares in percent
high_int_share <- (high_int/sum_spam_ta_sum_km2)*100

siebert_fallow <- 440000000  #ha

low_int_fal <- low_int+siebert_fallow
glob_cropland <- sum_spam_ta_sum_km2+siebert_fallow

low_int_fal_share <- (low_int_fal/glob_cropland)*100
high_int_fal_share <- (high_int/glob_cropland)*100

erb_low_int <- 473000000     #values from Erb et al, 2016 their Table S2
erb_high_int <- 743000000
erb_cropland <- 1216000000

erb_low_int_share <- (erb_low_int/erb_cropland)*100
erb_high_int_share <- (erb_high_int/erb_cropland)*100

comp_low_int <- cbind(low_int,low_int_share,low_int_fal,low_int_fal_share,erb_low_int,erb_low_int_share)
comp_high_int <- cbind(high_int,low_int_share,high_int,high_int_fal_share,erb_high_int,erb_high_int_share)

comp_high_low <- rbind(comp_low_int,comp_high_int)
colnames(comp_high_low) <- c("Tillage area this study (ha)","Tillage area this study (%)","Tillage area this study + fallow (ha)","Tillage area this study + fallow (%)","Tillage area (ha) (Erb et al., 2016)","Tillage area (%) (Erb et al., 2016)")

#fn <- paste0(path_calc,"tillage_intensities_comparison.txt")
#write.table(x=comp_high_low,file=fn,sep="\t")  #Table 7 in article

####raster to netcdf####
if (sample_calc){
  final.ext <- extent(-180,180,-45,90)
  r <- raster(res=45,ext=final.ext)
  origin(r) <- 0
} else {
  final.ext <- extent(-180,180,-56,84) #geopraphic coverage of dataset
  r <- raster(res=1/12,ext=final.ext)  #resolution to 0.083333...°
  r[] <- NA
}

rotational_tillage_ext <- stack(r,r,r,r,r,rotational_tillage[[1]],r,rotational_tillage[[2]],rotational_tillage[[3]],rotational_tillage[[4]],rotational_tillage[[5]],rotational_tillage[[6]],r,r,rotational_tillage[[7]],rotational_tillage[[8]],r,rotational_tillage[[9]],rotational_tillage[[10]],rotational_tillage[[11]],rotational_tillage[[12]],rotational_tillage[[13]],r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r)  #42

traditional_rotational_tillage_ext <- stack(r,r,r,r,r,traditional_rotational_tillage[[1]],r,traditional_rotational_tillage[[2]],traditional_rotational_tillage[[3]],traditional_rotational_tillage[[4]],traditional_rotational_tillage[[5]],traditional_rotational_tillage[[6]],r,r,traditional_rotational_tillage[[7]],traditional_rotational_tillage[[8]],r,traditional_rotational_tillage[[9]],traditional_rotational_tillage[[10]],traditional_rotational_tillage[[11]],traditional_rotational_tillage[[12]],traditional_rotational_tillage[[13]],r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r)  #42

ca_down <- stack(Conservation_Agriculture[[1]],r,Conservation_Agriculture[[2]],Conservation_Agriculture[[3]],Conservation_Agriculture[[4]],r,Conservation_Agriculture[[5]],r,r,r,r,r,Conservation_Agriculture[[6]],r,r,r,Conservation_Agriculture[[7]],r,r,r,r,r,Conservation_Agriculture[[22]],r,r,r,r,r,Conservation_Agriculture[[8]],Conservation_Agriculture[[9]],Conservation_Agriculture[[10]],Conservation_Agriculture[[11]],Conservation_Agriculture[[12]],Conservation_Agriculture[[13]],Conservation_Agriculture[[14]],Conservation_Agriculture[[15]],Conservation_Agriculture[[16]],Conservation_Agriculture[[17]],Conservation_Agriculture[[18]],Conservation_Agriculture[[19]],Conservation_Agriculture[[20]],Conservation_Agriculture[[21]]) #now 42

scenario_down <- stack(scenario_ca_area[[1]],r,scenario_ca_area[[2]],scenario_ca_area[[3]],scenario_ca_area[[4]],r,scenario_ca_area[[5]],r,r,r,r,r,scenario_ca_area[[6]],r,r,r,scenario_ca_area[[7]],r,r,r,r,r,scenario_ca_area[[22]],r,r,r,r,r,scenario_ca_area[[8]],scenario_ca_area[[9]],scenario_ca_area[[10]],scenario_ca_area[[11]],scenario_ca_area[[12]],scenario_ca_area[[13]],scenario_ca_area[[14]],scenario_ca_area[[15]],scenario_ca_area[[16]],scenario_ca_area[[17]],scenario_ca_area[[18]],scenario_ca_area[[19]],scenario_ca_area[[20]],scenario_ca_area[[21]]) #now 42

traditional_annual_tillage_ext <- stack(traditional_annual_tillage[[1]],traditional_annual_tillage[[2]],traditional_annual_tillage[[3]],traditional_annual_tillage[[4]],traditional_annual_tillage[[5]],r,traditional_annual_tillage[[6]],r,r,r,r,r,traditional_annual_tillage[[7]],traditional_annual_tillage[[8]],r,r,traditional_annual_tillage[[9]],r,r,r,r,r,traditional_annual_tillage[[29]],traditional_annual_tillage[[10]],traditional_annual_tillage[[11]],traditional_annual_tillage[[12]],traditional_annual_tillage[[13]],traditional_annual_tillage[[14]],traditional_annual_tillage[[15]],traditional_annual_tillage[[16]],traditional_annual_tillage[[17]],traditional_annual_tillage[[18]],traditional_annual_tillage[[19]],traditional_annual_tillage[[20]],traditional_annual_tillage[[21]],traditional_annual_tillage[[22]],traditional_annual_tillage[[23]],traditional_annual_tillage[[24]],traditional_annual_tillage[[25]],traditional_annual_tillage[[26]],traditional_annual_tillage[[27]],traditional_annual_tillage[[28]]) #stack, 42

ca_cro_nam_vec <- c("whea","maiz","barl","rest","toba","cott","vege","sesa","rape","sunf","soyb","opul","lent","pige","cowp" ,"chic" ,"bean","ocer","sorg","smil","pmil","grou") #22
ca_cro_nr_vec <- c(1,3,4,5,7,13,17,23,29,30,31,32,33,34,35,36,37,38,39,40,41,42)
tub_ric_nam_vec <- c("rice","sugb","orts","cass","yams","swpo","pota")
tub_ric_nr_vec <- c(2,14,24,25,26,27,28)
pere_nam_vec <- c("tea","cocoa","robusta_coffee","arabica_coffee","other_fibre_crops","sugarcane","oilpalm","temperate_fruit","tropical_fruit","plantain","banana","coconut","other_oil_crops")
pere_nr_vec <- c(6,8,9,10,11,12,15,16,18,19,20,21,22) #13

for (i in ca_cro_nr_vec) {
  print(paste("From 22 stacking ca_crop:",i)) 
  b <- raster(res=1/12,ext=final.ext)
  b[] <- 0 #produces a layer with same extent and resolution as SPAM_stuff for proper calc
  conventional_annual_tillage_fac <- conventional_annual_tillage[[i]]
  b[conventional_annual_tillage_fac>0] <- 1
  traditional_annual_tillage_ext_fac <- traditional_annual_tillage_ext[[i]]
  b[traditional_annual_tillage_ext_fac>0] <- 2
  reduced_tillage_fac <- reduced_tillage[[i]]
  b[reduced_tillage_fac>0] <- 3
  ca_down_fac <- ca_down[[i]]
  b[ca_down_fac>0] <- 4
  b[b==0] <-NA
  if (i>1) b_stack <- stack(b_stack,b) else b_stack <- b  
}

for (j in tub_ric_nr_vec){
  print(paste("From 7 stacking tuber or rice crop:",j))
  d <- raster(res=1/12,ext=final.ext)
  d[] <- 0 #produces a layer with same extent and resolution as SPAM_stuff for proper calc
  conventional_annual_tillage_fac <- conventional_annual_tillage[[j]]
  d[conventional_annual_tillage_fac>0] <- 1
  traditional_annual_tillage_ext_fac <- traditional_annual_tillage_ext[[j]]
  d[traditional_annual_tillage_ext_fac>0] <- 2
  reduced_tillage_fac <- reduced_tillage[[j]]
  d[reduced_tillage_fac>0] <- 3
  d[d==0] <- NA
  if (j>2) d_stack <- stack(d_stack,d) else d_stack <- d 
}

for (k in pere_nr_vec){
  print(paste("From 13 stacking perennial crop:",k))  
  g <- raster(res=1/12,ext=final.ext)
  g[] <- 0 
  rotational_tillage_fac<- rotational_tillage_ext[[k]]
  g[rotational_tillage_fac>0] <- 5
  traditional_rotational_tillage_fac <- traditional_rotational_tillage_ext[[k]]
  g[traditional_rotational_tillage_fac>0] <- 6
  g[g==0] <- NA
  if (k>6) g_stack <- stack(g_stack,g) else g_stack <- g
}

scenario_down_sum <- sum(scenario_ca_area,na.rm=T)
scenario_down_sum[scenario_down_sum>0] <- 7
scenario_down_sum[scenario_down_sum==0] <- NA

tillage_43 <- stack(b_stack[[1]],d_stack[[1]],b_stack[[2]],b_stack[[3]],b_stack[[4]],g_stack[[1]],b_stack[[5]],g_stack[[2]],g_stack[[3]],g_stack[[4]],g_stack[[5]],g_stack[[6]],b_stack[[6]],d_stack[[2]],g_stack[[7]],g_stack[[8]],b_stack[[7]],g_stack[[9]],g_stack[[10]],g_stack[[11]],g_stack[[12]],g_stack[[13]],b_stack[[8]],d_stack[[3]],d_stack[[4]],d_stack[[5]],d_stack[[6]],d_stack[[7]],b_stack[[9]],b_stack[[10]],b_stack[[11]],b_stack[[12]],b_stack[[13]],b_stack[[14]],b_stack[[15]],b_stack[[16]],b_stack[[17]],b_stack[[18]],b_stack[[19]],b_stack[[20]],b_stack[[21]],b_stack[[22]],scenario_down_sum)

till_brick <- brick(tillage_43,values=TRUE)

data_array <- as.array(till_brick)       
array_perm <- aperm(data_array,c(2,1,3))    #exchange dims of lat and lon 

filename <- "tillage.nc4"
path_fn <- paste0(path_calc,filename)
missval <- NaN

londim <- ncdim_def("lon","degrees_east",seq(-180+((1/12)*1/2),180-((1/12)*1/2),1/12),longname="Longitude") 
latdim <- ncdim_def("lat","degrees_north",seq(84-((1/12)*1/2),-56+((1/12)*1/2),-1/12),longname="Latitude")

var1d <- ncvar_def("whea_till","category 1-4",list(londim,latdim),missval,longname="Wheat tillage",compression=9)
var2d <- ncvar_def("rice_till","category 1-3",list(londim,latdim),missval,longname="Rice tillage",compression=9)
var3d <- ncvar_def("maiz_till","category 1-4",list(londim,latdim),missval,longname="Maize tillage",compression=9)
var4d <- ncvar_def("barl_till","category 1-4",list(londim,latdim),missval,longname="Barley tillage",compression=9)
var5d <- ncvar_def("rest_till","category 1-4",list(londim,latdim),missval,longname="Rest tillage",compression=9)
var6d <- ncvar_def("ooil_till","category 5-6",list(londim,latdim),missval,longname="Other oilcrops tillage",compression=9)
var7d <- ncvar_def("toba_till","category 1-4",list(londim,latdim),missval,longname="Tobacco tillage",compression=9 )
var8d <- ncvar_def("teas_till","category 5-6",list(londim,latdim),missval,longname="Teas tillage",compression=9)
var9d <- ncvar_def("coco_till","category 5-6",list(londim,latdim),missval,longname="Cocoa tillage",compression=9)
var10d <- ncvar_def("rcof_till","category 5-6",list(londim,latdim),missval,longname="Robusta coffee tillage",compression=9)
var11d <- ncvar_def("acof_till","category 5-6",list(londim,latdim),missval,longname="Arabica coffee tillage",compression=9)
var12d <- ncvar_def("ofib_till","category 5-6",list(londim,latdim),missval,longname="Other fibre crops tillage",compression=9)
var13d <- ncvar_def("cott_till","category 1-4",list(londim,latdim),missval,longname="Cotton tillage",compression=9)
var14d <- ncvar_def("sugb_till","category 1-3",list(londim,latdim),missval,longname="Sugarbeet tillage",compression=9)
var15d <- ncvar_def("sugc_till","category 5-6",list(londim,latdim),missval,longname="Sugarcane tillage",compression=9)
var16d <- ncvar_def("oilp_till","category 5-6",list(londim,latdim),missval,longname="Oilpalm tillage",compression=9)
var17d <- ncvar_def("vege_till","category 1-4",list(londim,latdim),missval,longname="Vegetables tillage",compression=9)
var18d <- ncvar_def("temf_till","category 5-6",list(londim,latdim),missval,longname="Temperate fruit tillage",compression=9)
var19d <- ncvar_def("trof_till","category 5-6",list(londim,latdim),missval,longname="Tropical fruit tillage",compression=9)
var20d <- ncvar_def("plnt_till","category 5-6",list(londim,latdim),missval,longname="Plantain tillage",compression=9)
var21d <- ncvar_def("bana_till","category 5-6",list(londim,latdim),missval,longname="Banana tillage",compression=9)
var22d <- ncvar_def("cnut_till","category 5-6",list(londim,latdim),missval,longname="Coconut tillage",compression=9)
var23d <- ncvar_def("grou_till","category 1-4",list(londim,latdim),missval,longname="Groundnut tillage",compression=9)
var24d <- ncvar_def("orts_till","category 1-3",list(londim,latdim),missval,longname="Other roots tillage",compression=9)
var25d <- ncvar_def("cass_till","category 1-3",list(londim,latdim),missval,longname="Cassava tillage",compression=9)
var26d <- ncvar_def("yams_till","category 1-3",list(londim,latdim),missval,longname="Yams tillage",compression=9)
var27d <- ncvar_def("swpo_till","category 1-3",list(londim,latdim),missval,longname="Sweet potato tillage",compression=9)
var28d <- ncvar_def("pota_till","category 1-3",list(londim,latdim),missval,longname="Potato tillage",compression=9)
var29d <- ncvar_def("sesa_till","category 1-4",list(londim,latdim),missval,longname="Sesameseed tillage",compression=9)
var30d <- ncvar_def("rape_till","category 1-4",list(londim,latdim),missval,longname="Rapeseed tillage",compression=9)
var31d <- ncvar_def("sunf_till","category 1-4",list(londim,latdim),missval,longname="Sunflower tillage",compression=9)
var32d <- ncvar_def("soyb_till","category 1-4",list(londim,latdim),missval,longname="Soybean tillage",compression=9)
var33d <- ncvar_def("opul_till","category 1-4",list(londim,latdim),missval,longname="Other pulses tillage",compression=9)
var34d <- ncvar_def("lent_till","category 1-4",list(londim,latdim),missval,longname="Lentils tillage",compression=9)
var35d <- ncvar_def("pige_till","category 1-4",list(londim,latdim),missval,longname="Pigeopea tillage",compression=9)
var36d <- ncvar_def("cowp_till","category 1-4",list(londim,latdim),missval,longname="Cowpea tillage",compression=9)
var37d <- ncvar_def("chic_till","category 1-4",list(londim,latdim),missval,longname="Chicpea tillage",compression=9)
var38d <- ncvar_def("bean_till","category 1-4",list(londim,latdim),missval,longname="Bean tillage",compression=9)
var39d <- ncvar_def("ocer_till","category 1-4",list(londim,latdim),missval,longname="Other cereals tillage",compression=9)
var40d <- ncvar_def("sorg_till","category 1-4",list(londim,latdim),missval,longname="Sorghum tillage",compression=9)
var41d <- ncvar_def("smil_till","category 1-4",list(londim,latdim),missval,longname="Small millet tillage",compression=9)
var42d <- ncvar_def("pmil_till","category 1-4",list(londim,latdim),missval,longname="Pearl millet tillage",compression=9)
var43d <- ncvar_def("scenario_ca_area","category 7",list(londim,latdim),missval,longname="Scenario Conservation Agriculture area",compression=9)

nc <- nc_create(path_fn,list(var1d,var2d,var3d,var4d,var5d,var6d,var7d,var8d,var9d,var10d,var11d,var12d,var13d,var14d,var15d,var16d,var17d,var18d,var19d,var20d,var21d,var22d,var23d,var24d,var25d,var26d,var27d,var28d,var29d,var30d,var31d,var32d,var33d,var34d,var35d,var36d,var37d,var38d,var39d,var40d,var41d,var42d,var43d),force_v4=T)

ncvar_put(nc,var1d,array_perm[,,1])
ncvar_put(nc,var2d,array_perm[,,2])
ncvar_put(nc,var3d,array_perm[,,3])
ncvar_put(nc,var4d,array_perm[,,4])
ncvar_put(nc,var5d,array_perm[,,5])
ncvar_put(nc,var6d,array_perm[,,6])
ncvar_put(nc,var7d,array_perm[,,7])
ncvar_put(nc,var8d,array_perm[,,8])
ncvar_put(nc,var9d,array_perm[,,9])
ncvar_put(nc,var10d,array_perm[,,10])
ncvar_put(nc,var11d,array_perm[,,11])
ncvar_put(nc,var12d,array_perm[,,12])
ncvar_put(nc,var13d,array_perm[,,13])
ncvar_put(nc,var14d,array_perm[,,14])
ncvar_put(nc,var15d,array_perm[,,15])
ncvar_put(nc,var16d,array_perm[,,16])
ncvar_put(nc,var17d,array_perm[,,17])
ncvar_put(nc,var18d,array_perm[,,18])
ncvar_put(nc,var19d,array_perm[,,19])
ncvar_put(nc,var20d,array_perm[,,20])
ncvar_put(nc,var21d,array_perm[,,21])
ncvar_put(nc,var22d,array_perm[,,22])
ncvar_put(nc,var23d,array_perm[,,23])
ncvar_put(nc,var24d,array_perm[,,24])
ncvar_put(nc,var25d,array_perm[,,25])
ncvar_put(nc,var26d,array_perm[,,26])
ncvar_put(nc,var27d,array_perm[,,27])
ncvar_put(nc,var28d,array_perm[,,28])
ncvar_put(nc,var29d,array_perm[,,29])
ncvar_put(nc,var30d,array_perm[,,30])
ncvar_put(nc,var31d,array_perm[,,31])
ncvar_put(nc,var32d,array_perm[,,32])
ncvar_put(nc,var33d,array_perm[,,33])
ncvar_put(nc,var34d,array_perm[,,34])
ncvar_put(nc,var35d,array_perm[,,35])
ncvar_put(nc,var36d,array_perm[,,36])
ncvar_put(nc,var37d,array_perm[,,37])
ncvar_put(nc,var38d,array_perm[,,38])
ncvar_put(nc,var39d,array_perm[,,39])
ncvar_put(nc,var40d,array_perm[,,40])
ncvar_put(nc,var41d,array_perm[,,41])
ncvar_put(nc,var42d,array_perm[,,42])
ncvar_put(nc,var43d,array_perm[,,43])

#varid=0 for global attributes of output dataset
ncatt_put(nc,varid=0,"Licensed under Open Data Commons Open Database License (ODbL):","<see README.pdf>")
ncatt_put(nc,varid=0,"Accompanying R-script:","Porwollik, Vera; Rolinski, Susanne; Müller, Christoph (2019): A global gridded data set on tillage - R-code. V. 1.1. GFZ Data Services. http://doi.org/10.5880/PIK.2019.010")
ncatt_put(nc,varid=0,"Supplement to article:","Generating a global gridded tillage dataset, Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2018-152, in review, 2018")
ncatt_put(nc,varid=0,"Citation:","Porwollik, Vera; Rolinski, Susanne; Müller, Christoph (2019): A global gridded data set on tillage V. 1.1. GFZ Data Services. http://doi.org/10.5880/PIK.2019.009")
ncatt_put(nc,varid=0,"Contact:","vera.porwollik@pik-potsdam.de")
ncatt_put(nc,varid=0,"Title:","A global gridded tillage area dataset")     
ncatt_put(nc,varid=0,"Institution:","(C) 2018 Potsdam Institute for Climate Impact Research (PIK)")
nc_close(nc)
