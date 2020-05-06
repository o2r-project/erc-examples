#!/usr/bin/env Rscript
#Use central value of various statistics, like CentralParameters.R
#except for N2O: Fig SI.1. 
#calculated N2O damage ratios for Fig SI.3
#Also calculate optimal timescale for various HFCs and N2O (Table 2)

#Code expects files for 4 RCP scenarios, and a GDP projection file. 

library(ggplot2)  #for graphing
library(stats)    #for uniroot 
library(reshape2) #for melt command, which is used for graphing data format
library(caTools)  #for trapz, as it is closer to integrate than a straight sum.
library(gridExtra)
library(grid)
library(extrafont)
library(Cairo)

#import functions
source("src/metricsfunctions.R")

#PARAMETER LISTS
scenariolist <- c("RCP6")
CSlist <- c(3.92) #climate sensitivity
damagelist <- c(2) #damage exponent
forcimblist <- c(0.84) #Medhaug et al. forcing imbalance
tempofflist <- c(0.6) #0.6 is GISStemp2011 relative to 1950-1980. 0.8 is NCA 2014 estimate from 1880-2012.
GDPlist <- c("GDPref")
#discountlist <- c(1.03) #only 3% discount rate
discountlist <- c(seq(1.005,1.15, by = 0.001)) 
#hfclifetimelist <- c(seq(5,100, by = 5))
#common fluorinated GHGs & lifetimes:
#HFC-134a: 13.4 years
#PFC-14: 50,000
#HFC-23: 222
#HFC-125: 28.2
#NF3: 500
#HCFC-22: 11.9 (not interesting as it is so similar to HFC-134a)
#HCFC-122: 1 (not common, but convenient)
hfclifetimelist <- c(1, 13.4, 28.2, 222, 50000)


#NOTES ON PARAMETERS
#Higher emissions have offsetting effects: higher future forcing makes CO2 more valuable
#higher future CO2 concentration makes CH4 more valuable
#higher climate sensitivity makes CO2 more valuable
#Larger exponent makes CO2 more valuable
#A larger forcing imbalance likely makes CO2 relatively more valuable
#A larger baseline temperature makes CH4 relatively more valuable
#larger GDP growth makes CO2 more valuable (assuming no Ramsey approach)

#OTHER CONSTANTS
damage1 <- 0.00236 #from DICE2016R Excel spreadsheet. 
#damage1 won't matter until there's a second damage coefficient

#equivalentfactor <- 24.8 #24.8 is ratio of CO2exp to CH4exp with discount=3
#equivalentfactor <- 312 #for N2O, discount 3
equivalentfactor <- 312 #for general damage ratio calculations
#calculating pulse masses
ch4emission <- 28.3 #28.3 Mt should be 10 ppb
co2emission <- ch4emission*12/44/1000 * equivalentfactor
#1000 for Mt -> Gt
#12/44 conversion to get from GtCO2 -> GtC
#needed because we use 2.12 factor in CO2 function which is based on GtC
hfcemission <- ch4emission
n2oemission <- ch4emission
#hfclifetime <- 13.4

#READFILES

#RCP concentration & forcing data from http://www.pik-potsdam.de/~mmalte/rcps/
#For 2011, CH4 = 1803, CO2 = 391, N2O = 324, according to IPCC Table 8.2
#However, the RCP concentrations don't exactly match. Therefore, a calculated
#GWP using the RCPs would be slightly different than using the IPCC concentrations.

RCP3PD <- read.csv("data/rcp3PDgasconcs.csv", header=TRUE, stringsAsFactors = FALSE)
RCP45 <- read.csv("data/rcp45gasconcs.csv", header=TRUE, stringsAsFactors = FALSE)
RCP6 <- read.csv("data/rcp6gasconcs.csv", header=TRUE, stringsAsFactors = FALSE)
RCP85 <- read.csv("data/rcp85gasconcs.csv", header=TRUE, stringsAsFactors = FALSE) 

dyears = length(RCP3PD$co2conc)-1  #calculate length of dataset

#GDP from Nordhaus DICE 2016 (http://www.econ.yale.edu/~nordhaus/homepage/DICEmodels09302016.htm)
gdp <- read.csv("data/gdpovertime.csv", header=TRUE, stringsAsFactors = FALSE)

#slowcons = half the growth rate
#fastcons = one and a half times the growth rate
#need to interpolate between 5 year intervals from Nordhaus:
#GDP in trillion dollars
GDPref <- approx(x = gdp$year, y = gdp$consumption, xout = seq(2011, dyears+2011), method = "linear")
GDPlow <- approx(x = gdp$year, y = gdp$slowcons, xout = seq(2011, dyears+2011), method = "linear")
GDPhigh <- approx(x = gdp$year, y = gdp$fastcons, xout = seq(2011, dyears+2011), method = "linear")

#create dataframe for final data
discountnames <- paste('discount',discountlist, sep="")
mycolnames <- c('scenario','gdp','CS','damage','forcimb','tempoff',discountnames)
totallength <- length(scenariolist)*length(GDPlist)*length(CSlist)*length(damagelist)*
  length(forcimblist)*length(tempofflist)
sensitivitymatrix <- matrix(ncol=length(mycolnames), nrow=totallength)
sensitivityframe <- data.frame(sensitivitymatrix)
sensitivityframe <- setNames(sensitivityframe, mycolnames)
iterationnum <- 1

discountmatrix <- matrix(ncol = length(hfclifetimelist)+2, nrow=length(discountlist))
discountframe <- data.frame(discountmatrix)
discountframenames <- c(hfclifetimelist, "n2o", "discount")
discountframe <- setNames(discountframe, discountframenames)

hfctimescalevector <- vector(length = length(hfclifetimelist))
hfclifetimeiter <- 1

for(hfclifetime in hfclifetimelist)
{
#start looping over all variables
#which, in this case, is only 1 variable each
for(scenarioiter in scenariolist)
{
  concentrations <- get(scenarioiter)
  
  for (gdpiter in GDPlist)  
  {
    gdpseq <- get(gdpiter)
    
    #Adding GDP data to concentration data. 
    fulldata <- cbind(concentrations, gdpseq$y)
    names(fulldata)[6] <- "gdp"
    
    for (CSiter in CSlist)
      for (damageiter in damagelist)
        for (forcimbiter in forcimblist)
          for (tempoffiter in tempofflist)
          {
            #calculate radiative forcing increment from CH4 & CO2 emissions:
            fulldata$hfcrad <- hfcforc(hfcemission, 
                                       fulldata$year-fulldata$year[1], hfclifetime) -
              hfcforc(0,  fulldata$year-fulldata$year[1], hfclifetime)
            fulldata$co2rad <- co2forc(co2emission, fulldata$co2conc, 
                                       fulldata$year-fulldata$year[1]) -
              co2forc(0, fulldata$co2conc, 
                      fulldata$year-fulldata$year[1])
            
            #calculate how much larger totforc is in year 1 than the 
            #forcing imbalance estimated from Medhaug et al.
            forcimbalance <- forcimbiter - fulldata$totforc[1]
            
            #calculate temperature in BAU and in two pulse scenarios for each year      
            fulldata$backtemp <- tempcalc(fulldata$year-fulldata$year[1], fulldata$totforc, forcimbalance, CSiter)
            forcvector <- fulldata$totforc+fulldata$hfcrad
            fulldata$hfctemp <- tempcalc(fulldata$year-fulldata$year[1], forcvector, forcimbalance, CSiter)
            forcvector <- fulldata$totforc+fulldata$co2rad
            fulldata$co2temp <- tempcalc(fulldata$year-fulldata$year[1], forcvector, forcimbalance, CSiter)
            
            #set negative net temperature to zero so it doesn't cause errors with exponents
            fulldata$hfctemp <- ifelse(fulldata$hfctemp < -tempoffiter,-tempoffiter,fulldata$hfctemp)
            fulldata$co2temp <- ifelse(fulldata$co2temp < -tempoffiter, -tempoffiter,fulldata$co2temp)
            fulldata$backtemp <- ifelse(fulldata$backtemp < -tempoffiter, -tempoffiter,fulldata$backtemp)
            
            #calculate the "damage" for each year for background & two pulses
            fulldata$hfcexp <- damage1 * (fulldata$hfctemp + tempoffiter)^damageiter * gdpseq$y
            fulldata$co2exp <- damage1 * (fulldata$co2temp + tempoffiter)^damageiter *gdpseq$y
            fulldata$backexp <- damage1 * (fulldata$backtemp + tempoffiter)^damageiter * gdpseq$y
            
 
            fulldata$discounthfc <- (fulldata$hfcexp-fulldata$backexp)/(1.03^(fulldata$year-fulldata$year[1]))
            fulldata$discountco2 <- (fulldata$co2exp-fulldata$backexp)/(1.03^(fulldata$year-fulldata$year[1]))
            
            dataforgraph <- with(fulldata,as.data.frame(cbind(year,hfcrad,co2rad,hfctemp-backtemp,co2temp-backtemp,hfcexp-backexp,co2exp-backexp,discounthfc,discountco2)))
            datanames <- c("year","hfcrad","co2rad","hfctemp","co2temp","hfcexp","co2exp","discounthfc","discountco2")
            colnames(dataforgraph) <- datanames
            melteddata <- melt(dataforgraph, id = "year", measured=c("hfcrad","co2rad","hfctemp","co2temp","hfcexp","co2exp","discounthfc","discountco2"))

            #discount difference between pulse scenario and BAU scenario 
            #back to present using several discount rates
            #calculate the ratio of the discounted net present value of the two pulses
            #and save in discountvector. 
            discountvector <- vector() #maybe initialize with length of discountlist?
            discountnum <- 1
            for (discountiter in discountlist)
            {
              netexphfc <- trapsum2(fulldata$hfcexp, discountiter)
              netexpco2 <- trapsum2(fulldata$co2exp, discountiter)
              netexpback <- trapsum2(fulldata$backexp, discountiter)
              discountvector[discountnum] <- (netexphfc-netexpback)/(netexpco2-netexpback)
              discountnum <- discountnum + 1
            }

            discountframe[,toString(hfclifetime)] <- discountvector

          }
  }
}

constantco2 <- rep(concentrations$co2conc[1], times = length(concentrations$co2conc))
constantch4 <- rep(concentrations$ch4conc[1], times = length(concentrations$ch4conc))
constantn2o <- rep(concentrations$n2oconc[1], times = length(concentrations$n2oconc))
netexphfc <- trapsum2(fulldata$hfcexp, 1.03)
netexpco2 <- trapsum2(fulldata$co2exp, 1.03)
netexpback <- trapsum2(fulldata$backexp, 1.03)
ratioofimpact <- (netexphfc-netexpback)/(netexpco2-netexpback)
hfctimescalevector[hfclifetimeiter] <- ratioofimpact

hfclifetimeiter <- hfclifetimeiter + 1
}


########################
##AND AGAIN FOR N2O ###
########################


for(scenarioiter in scenariolist)
{
  concentrations <- get(scenarioiter)
  
  for (gdpiter in GDPlist)  
  {
    gdpseq <- get(gdpiter)
    
    #Adding GDP data to concentration data. 
    fulldata <- cbind(concentrations, gdpseq$y)
    names(fulldata)[6] <- "gdp"
    
    for (CSiter in CSlist)
      for (damageiter in damagelist)
        for (forcimbiter in forcimblist)
          for (tempoffiter in tempofflist)
          {
            #calculate radiative forcing increment from CH4 & CO2 emissions:
            fulldata$n2orad <- n2oforc(n2oemission, fulldata$ch4conc, fulldata$n2oconc,
                                       fulldata$year-fulldata$year[1]) -
              n2oforc(0,  fulldata$ch4conc, fulldata$n2oconc, 
                                      fulldata$year-fulldata$year[1])
            fulldata$co2rad <- co2forc(co2emission, fulldata$co2conc, 
                                       fulldata$year-fulldata$year[1]) -
              co2forc(0, fulldata$co2conc, 
                      fulldata$year-fulldata$year[1])
            
            #calculate how much larger totforc is in year 1 than the 
            #forcing imbalance estimated from Medhaug et al.
            forcimbalance <- forcimbiter - fulldata$totforc[1]
            
            #calculate temperature in BAU and in two pulse scenarios for each year      
            fulldata$backtemp <- tempcalc(fulldata$year-fulldata$year[1], fulldata$totforc, forcimbalance, CSiter)
            forcvector <- fulldata$totforc+fulldata$n2orad
            fulldata$n2otemp <- tempcalc(fulldata$year-fulldata$year[1], forcvector, forcimbalance, CSiter)
            forcvector <- fulldata$totforc+fulldata$co2rad
            fulldata$co2temp <- tempcalc(fulldata$year-fulldata$year[1], forcvector, forcimbalance, CSiter)
            
            #set negative net temperature to zero so it doesn't cause errors with exponents
            fulldata$n2otemp <- ifelse(fulldata$n2otemp < -tempoffiter,-tempoffiter,fulldata$n2otemp)
            fulldata$co2temp <- ifelse(fulldata$co2temp < -tempoffiter, -tempoffiter,fulldata$co2temp)
            fulldata$backtemp <- ifelse(fulldata$backtemp < -tempoffiter, -tempoffiter,fulldata$backtemp)
            
            #calculate the "damage" for each year for background & two pulses
            fulldata$n2oexp <- damage1 * (fulldata$n2otemp + tempoffiter)^damageiter * gdpseq$y
            fulldata$co2exp <- damage1 * (fulldata$co2temp + tempoffiter)^damageiter *gdpseq$y
            fulldata$backexp <- damage1 * (fulldata$backtemp + tempoffiter)^damageiter * gdpseq$y
            
  
            fulldata$discountn2o <- (fulldata$n2oexp-fulldata$backexp)/(1.03^(fulldata$year-fulldata$year[1]))
            fulldata$discountco2 <- (fulldata$co2exp-fulldata$backexp)/(1.03^(fulldata$year-fulldata$year[1]))
            
            dataforgraph <- with(fulldata,as.data.frame(cbind(year,n2orad,co2rad,n2otemp-backtemp,co2temp-backtemp,n2oexp-backexp,co2exp-backexp,discountn2o,discountco2)))
            datanames <- c("year","n2orad","co2rad","n2otemp","co2temp","n2oexp","co2exp","discountn2o","discountco2")
            colnames(dataforgraph) <- datanames
            melteddata <- melt(dataforgraph, id = "year", measured=c("n2orad","co2rad","n2otemp","co2temp","n2oexp","co2exp","discountn2o","discountco2"))
            
            #Plot Figure SI.1 from melteddata
            ytitle <- expression("mW m"^-2)
            n2otitle <- expression(paste("N"[2], "O", sep=""))
            co2title <- expression(CO[2])
            labeltitle <- c(n2otitle, co2title)
            radplot <- ggplot(subset(melteddata, variable %in% c("n2orad","co2rad")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point() + 
              labs(title="(a) Radiative forcing", x=element_blank(), y=ytitle)+
              theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(1,1),
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+scale_y_continuous(labels=function(x)x*1000)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            tempplot <- ggplot(subset(melteddata, variable %in% c("n2otemp","co2temp")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point()+ 
              labs(title="(b) Temperature", x=element_blank(), y="mK")+
              theme(legend.title=element_blank(), legend.justification=c(1,1), legend.position=c(1,1), 
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+scale_y_continuous(labels=function(x)x*1000)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            damageplot <- ggplot(subset(melteddata, variable %in% c("n2oexp","co2exp")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point()+ 
              labs(title="(c) Damages",x="Year", y="Billion $")+
              theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(1,1), 
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            discountplot <- ggplot(subset(melteddata, variable %in% c("discountn2o","discountco2")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point()+ 
              labs(title="(d) Discounted",x="Year", y="Billion $")+
              theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(1,1), 
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            
            gridplot <- grid.arrange(radplot, tempplot, damageplot, discountplot,ncol = 2)
            ggsave("results/n2ogridplot.pdf", plot = gridplot, device = cairo_pdf, width = 6, height = 6, units = "in", dpi = 600)
            
            #discount difference between pulse scenario and BAU scenario 
            #back to present using several discount rates
            #calculate the ratio of the discounted net present value of the two pulses
            #and save in discountvector. 
            discountvector <- vector() 
            discountnum <- 1
            for (discountiter in discountlist)
            {
              netexpn2o <- trapsum2(fulldata$n2oexp, discountiter)
              netexpco2 <- trapsum2(fulldata$co2exp, discountiter)
              netexpback <- trapsum2(fulldata$backexp, discountiter)
              discountvector[discountnum] <- (netexpn2o-netexpback)/(netexpco2-netexpback)
              discountnum <- discountnum + 1
            }
            
            discountframe[,"n2o"] <- discountvector
            discountframe[,"discount"] <- 100*(discountlist-1)
            
          }
  }
}

netexpn2o <- trapsum2(fulldata$n2oexp, 1.03)
netexpco2 <- trapsum2(fulldata$co2exp, 1.03)
netexpback <- trapsum2(fulldata$backexp, 1.03)
ratioofimpact <- (netexpn2o-netexpback)/(netexpco2-netexpback)

#hfctimescalevector for hfc ratio of impacts
#ratioofimpact is n2o ratio of impact

#Figure SI.3, N2O. 
n2otitle <- expression(paste("N"[2], "O damage ratio", sep=""))
n2oaxis <- "Damage ratio" 
n2oplot <- ggplot(discountframe, aes(x=discount, y=equivalentfactor*n2o)) + geom_point() + 
  labs(title=n2otitle, x="Discount rate", y=n2oaxis)+
  theme(text=element_text(family="Arial", size = 16), plot.title=element_text(size=16))

hfctitle <- "HFC-134a damage ratio"
hfcaxis <- "Damage ratio" 
hfcplot <- ggplot(discountframe, aes(x=discount, y=equivalentfactor*discountframe[2])) + geom_point() + 
  labs(title=hfctitle, x="Discount rate", y=hfcaxis)+
  theme(text=element_text(family="Arial", size = 16), plot.title=element_text(size=16))

ggsave("results/N2Odamageratio.pdf", plot = n2oplot, device = cairo_pdf, width = 3, height = 3, units = "in", dpi = 600)

#GWP100/damageratio for Table 2:
HFC134aGWP100ratio <- calcHFCGWP(100, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[2])/(hfctimescalevector[2]*equivalentfactor)
HFC134aGWP20ratio <- calcHFCGWP(20, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[2])/(hfctimescalevector[2]*equivalentfactor)
HFC23GWP100ratio <- calcHFCGWP(100, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[4])/(hfctimescalevector[4]*equivalentfactor)
HFC23GWP20ratio <- calcHFCGWP(20, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[4])/(hfctimescalevector[4]*equivalentfactor)
PFC14GWP100ratio <- calcHFCGWP(100, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[5])/(hfctimescalevector[5]*equivalentfactor)
PFC14GWP20ratio <- calcHFCGWP(20, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[5])/(hfctimescalevector[5]*equivalentfactor)

N2OGWP100ratio <- calcN2OGWP(100, co2emission/equivalentfactor, n2oemission, constantco2, constantch4, constantn2o)/(ratioofimpact*equivalentfactor)
N2OGWP20ratio <- calcN2OGWP(20, co2emission/equivalentfactor, n2oemission, constantco2, constantch4, constantn2o)/(ratioofimpact*equivalentfactor)

hfc134timescale <- calcHFCGWProot(hfctimescalevector[2]*equivalentfactor, co2emission/equivalentfactor, hfcemission, constantco2, hfclifetimelist[2])
#Optimal timescales for N2O, HFC-23, and PFC-14 were calculated manually
#for N2O and HFC-23 by taking the timescale where the damage ratio peaks
#for PFC-14 recognizing that it may exist, but goes beyond the timescale of this calculation.

print("reached end of code")
