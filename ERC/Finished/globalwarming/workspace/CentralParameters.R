#!/usr/bin/env Rscript
#Use central value of all parameters to show relative CH4 to CO2 impacts over time
#Graph 1 (set of 4): 
#show forcing, temperature, damages over time, discounted damages

#produces Figure 1, methane values for Table 2, and methane plot for Figure SI.3

#Code expects files for 4 RCP scenarios, and a GDP projection file. 

library(ggplot2)  #for graphing
library(stats)    #for uniroot 
library(reshape2) #for melt command, which is used for graphing data format
library(caTools)  #for trapz, as it is closer to integrate than a straight sum.
library(gridExtra) #for plots
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

#equivalentfactor <- 24.95 #24.95 is ratio of CO2exp to CH4exp with discount=3
equivalentfactor <- 24.95

#calculating pulse masses
ch4emission <- 28.3 #28.3 Mt should be 10 ppb
co2emission <- ch4emission*12/44/1000 * equivalentfactor
#1000 for Mt -> Gt
#12/44 conversion to get from GtCO2 -> GtC
#needed because we use 2.12 factor in CO2 function which is based on GtC

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
#need to interpolate between 5 year intervals from Nordhaus
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

discountmatrix <- matrix(ncol = 2, nrow=length(discountlist))
discountframe <- data.frame(discountmatrix)
discountframenames <- c("ch4", "discount")
discountframe <- setNames(discountframe, discountframenames)

iterationnum <- 1

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
            fulldata$ch4rad <- ch4forc(ch4emission, fulldata$ch4conc, fulldata$n2oconc, 
                                       fulldata$year-fulldata$year[1]) -
              ch4forc(0, fulldata$ch4conc, fulldata$n2oconc, 
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
            forcvector <- fulldata$totforc+fulldata$ch4rad
            fulldata$ch4temp <- tempcalc(fulldata$year-fulldata$year[1], forcvector, forcimbalance, CSiter)
            forcvector <- fulldata$totforc+fulldata$co2rad
            fulldata$co2temp <- tempcalc(fulldata$year-fulldata$year[1], forcvector, forcimbalance, CSiter)
            
            #set negative net temperature to zero so it doesn't cause errors with exponents
            fulldata$ch4temp <- ifelse(fulldata$ch4temp < -tempoffiter,-tempoffiter,fulldata$ch4temp)
            fulldata$co2temp <- ifelse(fulldata$co2temp < -tempoffiter, -tempoffiter,fulldata$co2temp)
            fulldata$backtemp <- ifelse(fulldata$backtemp < -tempoffiter, -tempoffiter,fulldata$backtemp)
            
            #calculate the "damage" for each year for background & two pulses
            fulldata$ch4exp <- 1000*damage1 * (fulldata$ch4temp + tempoffiter)^damageiter * gdpseq$y
            fulldata$co2exp <- 1000*damage1 * (fulldata$co2temp + tempoffiter)^damageiter *gdpseq$y
            fulldata$backexp <- 1000*damage1 * (fulldata$backtemp + tempoffiter)^damageiter * gdpseq$y
            
            fulldata$discountch4 <- (fulldata$ch4exp-fulldata$backexp)/(1.03^(fulldata$year-fulldata$year[1]))
            fulldata$discountco2 <- (fulldata$co2exp-fulldata$backexp)/(1.03^(fulldata$year-fulldata$year[1]))
            
            #can calculate what percent occurs in first N years:
            #trapsumshort(fulldata$ch4exp-fulldata$backexp,1.03,87)/trapsumshort(fulldata$ch4exp-fulldata$backexp,1.03,480)

            dataforgraph <- with(fulldata,as.data.frame(cbind(year,ch4rad,co2rad,ch4temp-backtemp,co2temp-backtemp,ch4exp-backexp,co2exp-backexp,discountch4,discountco2)))
            datanames <- c("year","ch4rad","co2rad","ch4temp","co2temp","ch4exp","co2exp","discountch4","discountco2")
            colnames(dataforgraph) <- datanames
            melteddata <- melt(dataforgraph, id = "year", measured=c("ch4rad","co2rad","ch4temp","co2temp","ch4exp","co2exp","discountch4","discountco2"))
            
            #Figure 1 from Paper
            #Plot Graph 1 set from melteddata
            #ytitle <- expression(mW/m^2)
            ytitle <- expression("mW m"^-2)
            ch4title <- expression(CH[4])
            co2title <- expression(CO[2])
            labeltitle <- c(ch4title, co2title)
            radplot <- ggplot(subset(melteddata, variable %in% c("ch4rad","co2rad")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point() + 
              labs(title="(a) Radiative forcing", x=element_blank(), y=ytitle)+
              theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(1,1),
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+scale_y_continuous(labels=function(x)x*1000)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            tempplot <- ggplot(subset(melteddata, variable %in% c("ch4temp","co2temp")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point()+ 
              labs(title="(b) Temperature", x=element_blank(), y="mK")+
              theme(legend.title=element_blank(), legend.justification=c(1,1), legend.position=c(1,1), 
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+scale_y_continuous(labels=function(x)x*1000)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            damageplot <- ggplot(subset(melteddata, variable %in% c("ch4exp","co2exp")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point()+ 
              labs(title="(c) Damages",x="Year", y="Billion $")+
              theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(1,1), 
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(ylim=c(0,30),xlim=c(2005,2295))
            discountplot <- ggplot(subset(melteddata, variable %in% c("discountch4","discountco2")), aes(x=year, y=value, color=variable, shape=variable)) + geom_point()+ 
              labs(title="(d) Discounted",x="Year", y="Billion $")+
              theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(1,1), 
                    text=element_text(family="Arial", size = 15), plot.title=element_text(size=15),
                    axis.text.x = element_text(hjust=0.7))+
              scale_colour_discrete(labels=labeltitle)+
              scale_shape_discrete(labels=labeltitle)+coord_cartesian(xlim=c(2005,2295))
            
            gridplot <- grid.arrange(radplot, tempplot, damageplot, discountplot,ncol = 2)
            ggsave("results/Figure1.gridplot.pdf", plot = gridplot, device = cairo_pdf, width = 12.7, height = 12.7, units = "cm", dpi = 600)
            
            #discount difference between pulse scenario and BAU scenario 
            #back to present using several discount rates
            #calculate the ratio of the discounted net present value of the two pulses
            #and save in discountvector. 
            discountvector <- vector() 
            discountnum <- 1
            for (discountiter in discountlist)
            {
              netexpch4 <- trapsum2(fulldata$ch4exp, discountiter)
              netexpco2 <- trapsum2(fulldata$co2exp, discountiter)
              netexpback <- trapsum2(fulldata$backexp, discountiter)
              discountvector[discountnum] <- (netexpch4-netexpback)/(netexpco2-netexpback)
              discountnum <- discountnum + 1
            }
            
            discountframe[,"ch4"] <- discountvector
            discountframe[,"discount"] <- 100*(discountlist-1)

            #Figure SI.3, CH4 panel. 
            ch4title <- expression(paste("CH"[4], " damage ratio", sep=""))
            ch4axis <- "Damage ratio" 
            ch4plot <- ggplot(discountframe, aes(x=discount, y=equivalentfactor*ch4)) + geom_point() + 
              labs(title=ch4title, x="Discount rate", y=ch4axis)+
              theme(text=element_text(family="Arial", size = 16), plot.title=element_text(size=16))
            
            ggsave("results/CH4damageratio.pdf", plot = ch4plot, device = cairo_pdf, width = 3, height = 3, units = "in", dpi = 600)
            
            
            #Can also calculate what % of discounted damage comes in 
            #first N years.  
            #for discount rate 1.03, 90% of CO2 comes in 157 years, 95% in 189 years
            #for CH4: 87 years! 
            #Might be interesting to compare to Hector or MAGICC
          }
  }
}

constantco2 <- rep(concentrations$co2conc[1], times = length(concentrations$co2conc))
constantch4 <- rep(concentrations$ch4conc[1], times = length(concentrations$ch4conc))
constantn2o <- rep(concentrations$n2oconc[1], times = length(concentrations$n2oconc))
netexpch4 <- trapsum2(fulldata$ch4exp, 1.03)
netexpco2 <- trapsum2(fulldata$co2exp, 1.03)
netexpback <- trapsum2(fulldata$backexp, 1.03)
ratioofimpact <- (netexpch4-netexpback)/(netexpco2-netexpback)
#For Table2:
CH4GWP100ratio <- calcGWP(100,co2emission/equivalentfactor, ch4emission, constantco2, constantch4, constantn2o)/(ratioofimpact*equivalentfactor)
CH4GWP20ratio <- calcGWP(20,co2emission/equivalentfactor, ch4emission, constantco2, constantch4, constantn2o)/(ratioofimpact*equivalentfactor)
CH4optimaltimescale <- calcGWProot(ratioofimpact*equivalentfactor, co2emission/equivalentfactor, ch4emission, constantco2, constantch4, constantn2o)



#For examining other time lengths, use trapsumshort: e.g.
#(trapsumshort(fulldata$ch4exp,1.015,490)-trapsumshort(fulldata$backexp,1.015,490))/(trapsumshort(fulldata$co2exp,1.015,490)-trapsumshort(fulldata$backexp,1.015,490))
#(trapsumshort(fulldata$ch4exp,1.015,290)-trapsumshort(fulldata$backexp,1.015,290))/(trapsumshort(fulldata$co2exp,1.015,290)-trapsumshort(fulldata$backexp,1.015,290))
# which shows that at a 1.5% discount rate, 490 years yields a damage ratio of 11.02,
#compared to 12.0 for 290 years.
#Or, for percent of damage within X years:
#(trapsumshort(fulldata$co2exp, 1.02, 287)-trapsumshort(fulldata$backexp, 1.02, 287))/(trapsum(fulldata$co2exp, 1.02)-trapsum(fulldata$backexp, 1.02))
#showing that 95% of CO2 damage from an emissions pulse occurs in the first 287 years at a 2% discount rate. 

print("reached end of code")
