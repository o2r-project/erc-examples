#!/usr/bin/env Rscript
#code supporting Equivalent GWP timescale paper for submission to ESD
#A sensitivity analysis calculating the relative damage of a pulse of CO2
#to a pulse of CH4, calculated over uncertain parameters and discount rates.
#GWP equivalent timescales are also calculated. 

#Produces Figure 2 and data for Table 1. 

#Code expects files for 4 RCP scenarios, and a GDP projection file. 

library(ggplot2)  #for graphing
library(stats)    #for uniroot, in metricsfunctions
library(caTools)  #for trapz, in metricsfunctions
library(extrafont) #for plots
library(Cairo)     #for plots

#import functions
source("src/metricsfunctions.R")

#PARAMETER LISTS
scenariolist <- c("RCP3PD", "RCP45", "RCP6", "RCP85")
CSlist <- c(1.5, 3.92, 4.5) #climate sensitivity. 3.92 is chosen to match IPCC default parameters.
damagelist <- c(1.5, 2, 3) #damage exponent
forcimblist <- c(0.75, 0.84, 0.93) #Medhaug et al. forcing imbalance range
tempofflist <- c(0, 0.6, 0.8) #0.6 is GISStemp2011 relative to 1950-1980. 0.8 is NCA 2014 estimate from 1880-2012.
GDPlist <- c("GDPlow","GDPref","GDPhigh") #GDPlow is 0.5 times the GDPref growth rate
#and GDPhigh is 1.5 times. GDPref from DICE 2016. 
discountlist <- c(seq(1.005,1.15, by = 0.001)) #discount rates from .5% to 15%, by 0.1%.

#NOTES ON PARAMETERS
#Higher emissions have offsetting effects: 
#higher total future forcing (and therefore temperature) makes CO2 more valuable
#because of non-linear damage function
#but higher future CO2 concentration makes CH4 more valuable because CO2 forcing
#is logarithmic (similarly, higher CH4 concentrations make CO2 more valuable)
#higher climate sensitivity makes CO2 more valuable
#Larger exponent makes CO2 more valuable
#A larger forcing imbalance likely makes CO2 relatively more valuable
#A larger baseline temperature change makes CH4 relatively more valuable
#larger GDP growth makes CO2 more valuable (assuming no Ramsey approach)

#OTHER CONSTANTS
damage1 <- 0.00236 #from DICE2016R Excel spreadsheet. 
#this constant has no effect on CH4/CO2 ratio in current formulation
#if the damage function was to be of the form ax^2 + bx^6 then it would matter.

#calculating pulse masses
ch4emission <- 28.3 #28.3 Mt should be 10 ppb
co2emission <- ch4emission*12/44/1000 
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
#need to interpolate between 5 year intervals from Nordhaus:
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

#start looping over all variables
for(scenarioiter in scenariolist) #scenarioiter <- scenariolist[[1]]
{
  concentrations <- get(scenarioiter)
  fulldata <- concentrations
  fulldata$gdp <- NA_real_
  
  #calculate radiative forcing increment from CH4 & CO2 emissions:
  fulldata$ch4rad <- ch4forc(ch4emission, fulldata$ch4conc, fulldata$n2oconc, 
                             fulldata$year-fulldata$year[1]) -
    ch4forc(0, fulldata$ch4conc, fulldata$n2oconc, 
            fulldata$year-fulldata$year[1])
  fulldata$co2rad <- co2forc(co2emission, fulldata$co2conc, 
                             fulldata$year-fulldata$year[1]) -
    co2forc(0, fulldata$co2conc, 
            fulldata$year-fulldata$year[1])
  
  #loop over CS, tempoffset, and forcing imbalance:
  for (CSiter in CSlist) # CSiter <- CSlist[[1]]
    for (tempoffiter in tempofflist) # tempoffiter <- tempofflist[[1]]
      for (forcimbiter in forcimblist) # forcimbiter <- forcimblist[[1]]
      { 
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
        
        #loop over GDPs and damage exponents
        for (gdpiter in GDPlist) #gdpiter <- GDPlist[[1]]  
          for (damageiter in damagelist) # damageiter <- damagelist[[1]]
          {
            gdpseq <- get(gdpiter)
            
            #Adding GDP data to concentration data. 
            fulldata$gdp <- gdpseq$y
            
            #calculate the "damage" for each year for background & two pulses
            fulldata$ch4exp <- damage1 * (fulldata$ch4temp + tempoffiter)^damageiter * gdpseq$y
            fulldata$co2exp <- damage1 * (fulldata$co2temp + tempoffiter)^damageiter *gdpseq$y
            fulldata$backexp <- damage1 * (fulldata$backtemp + tempoffiter)^damageiter * gdpseq$y
            
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
            
            #add non-numerical parameters to data frame
            sensitivityframe[iterationnum, 1:2] <- c(scenarioiter, gdpiter)
            #put all numerical parameters in a row, plus calculated discount values
            newrow <- c(CSiter, damageiter, forcimbiter, tempoffiter, discountvector)
            sensitivityframe[iterationnum,3:ncol(sensitivityframe)] <- newrow
            #to see how far along the calculation is.
            if (iterationnum %% 10 == 1) print(c("iteration: ",iterationnum))
            iterationnum <- iterationnum + 1
          }
      }
}

write.csv(sensitivityframe, "results/finalframeESD.csv")

discountquant <- matrix(ncol=7, nrow=length(discountnames))
colnames(discountquant) <- c("GWPmin", "GWP10", "GWP25", "GWPmed", "GWP75", "GWP90", "GWPmax")
rownames(discountquant) <- discountnames
for(iter in 1:length(discountnames))
{
  discountholder <- sensitivityframe[discountnames[iter]]
  discountquant[iter,] <- quantile(discountholder[,], probs = c(0,0.1,0.25,0.5,0.75,0.9,1.0))
}
write.csv(discountquant, "results/discountquantESD.csv")

constantco2 <- rep(concentrations$co2conc[1], times = length(concentrations$co2conc))
constantch4 <- rep(concentrations$ch4conc[1], times = length(concentrations$ch4conc))
constantn2o <- rep(concentrations$n2oconc[1], times = length(concentrations$n2oconc))
GWPquant <- matrix(ncol=7, nrow=length(discountnames))
colnames(GWPquant) <- c("GWPmin", "GWP10", "GWP25", "GWPmed", "GWP75", "GWP90", "GWPmax")
rownames(GWPquant) <- discountnames
for (discountrow in 1:length(discountnames))
{
  for (discountcol in 1:7)
  {
  GWPquant[discountrow,discountcol] <- calcGWProot(discountquant[discountrow, discountcol], co2emission, ch4emission, constantco2, constantch4, constantn2o)
  }
}
write.csv(GWPquant, "results/GWPquantESD.csv")

#graphs median, 90/10, 25/75 lines: could be more aesthetic. Could add max/min.
#consider updating "discount" to be "(discount-1)*100". 
GWPplusone <- as.data.frame(cbind(GWPquant, discountlist))
colnames(GWPplusone) <- c("GWPmin", "GWP10", "GWP25", "GWPmed", "GWP75", "GWP90", "GWPmax","discount")
#Figure 2 in paper (except redrawn using other software, and including max & min)
fig2 <- ggplot(GWPplusone)+geom_line(aes(x=discount,y=GWP25))+geom_line(aes(x=discount,y=GWPmed))+
  geom_line(aes(x=discount,y=GWP75))+geom_line(aes(x=discount,y=GWP10))+geom_line(aes(x=discount,y=GWP90))+
  coord_cartesian(ylim=c(0,200))+labs(title="equivalent GWP by discount rate",x="Discount Rate", y="equivalent GWP year length")
ggsave("results/fig2dummy.pdf", plot = fig2, device = cairo_pdf, width = 6.4, height = 6.4, units = "cm", dpi = 600)
#flipping axes, for GWP as the :
ggplot(GWPplusone)+geom_line(aes(x=discount,y=GWP25))+geom_line(aes(x=discount,y=GWPmed))+
  geom_line(aes(x=discount,y=GWP75))+coord_flip(ylim=c(0,200))+
  labs(title="discount rate by GWP",y="equivalent GWP year length",x="Discount Rate")

#to look at the max/min/median with one given parameter option, use:
#quantile(sensitivityframe[which(sensitivityframe$CS == 4.5),"discount1.03"], probs = c(0,0.25,0.5,0.75,1.0))
#max(sensitivityframe[which(finalframe$CS == 4.5),"discount1.03"])

#Table 1 calculations below: the inverse is then taken of the offset and imbalance ratios
GDPratio <- median(sensitivityframe[which(sensitivityframe$gdp == "GDPlow"),"discount1.03"])/median(sensitivityframe[which(sensitivityframe$gdp == "GDPhigh"),"discount1.03"])
Damageratio <- median(sensitivityframe[which(sensitivityframe$damage == min(damagelist)),"discount1.03"])/median(sensitivityframe[which(sensitivityframe$damage == max(damagelist)),"discount1.03"])
Scenarioratio <- median(sensitivityframe[which(sensitivityframe$scenario == "RCP3PD"),"discount1.03"])/median(sensitivityframe[which(sensitivityframe$scenario == "RCP85"),"discount1.03"])
Offsetratio <- median(sensitivityframe[which(sensitivityframe$tempoff == min(tempofflist)),"discount1.03"])/median(sensitivityframe[which(sensitivityframe$tempoff == max(tempofflist)),"discount1.03"])
CSratio <- median(sensitivityframe[which(sensitivityframe$CS == min(CSlist)),"discount1.03"])/median(sensitivityframe[which(sensitivityframe$CS == max(CSlist)),"discount1.03"])
imbalanceratio <- median(sensitivityframe[which(sensitivityframe$forcimb == min(forcimblist)),"discount1.03"])/median(sensitivityframe[which(sensitivityframe$forcimb == max(forcimblist)),"discount1.03"])


#GRAPH DEVELOPMENT
#sec.axis works by transforming the first axis, so need to design a function that takes the CH4/CO2 ratio
#as an input, and outputs the number of years. Uniroot would be ideal, but doesn't work,
#so instead fit a polynomial to the breakpoints and create a function using the
#polynomial coefficients. See below. 

GWProotvector <- vector()
ratiovector <- vector()
tempvar <- 1
#solve for GWP year lengths equal to the breakpoints for the secondary axis:
breaklist <- c(200,100,60,40,30,20)
for (iter in breaklist)
{
 ratiovector[tempvar] <- calcGWP(iter, co2emission, ch4emission, constantco2, constantch4, constantn2o)
 GWProotvector[tempvar] <- iter
 tempvar <- tempvar + 1
}
#fifth order polynomial fit yields exact values for 6 breakpoints.
fit4 <- lm(GWProotvector~poly(ratiovector,5,raw=TRUE))

 GWPfit <- function(ratioinput)
{
coefs <- coef(fit4)
result <- coefs[1]+(coefs[2]*ratioinput)+
  (coefs[3]*ratioinput^2)+(coefs[4]*ratioinput^3)+
 (coefs[5]*ratioinput^4)+(coefs[6]*ratioinput^5)
 return(result)
}

ratioplusone <- as.data.frame(cbind(discountquant, discountlist))
colnames(ratioplusone) <- c("GWPmin", "GWP10", "GWP25", "GWPmed", "GWP75", "GWP90", "GWPmax","discount")

ytitle <- expression(CH[4]~to~CO[2]~ratio)
myplot <- ggplot(ratioplusone)+geom_line(aes(x=discount, y = GWP25))+geom_line(aes(x=discount,y=GWPmed))+
  geom_line(aes(x=discount,y=GWP75))+geom_line(aes(x=discount,y=GWP10))+geom_line(aes(x=discount,y=GWP90))+
 labs(title="GWP and damage ratio",x="Discount rate", y=ytitle)+
theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(0,1))+
scale_colour_discrete(labels=c("integrated radiative forcing","integrated temperature","tempsquare"))+
scale_shape_discrete(labels=c("integrated radiative forcing","integrated temperature","tempsquare"))+
coord_cartesian(ylim=c(5,95))

#need to use a GWPfit because the program can't handle using calcGWProot directly.
secaxisplot <- myplot +  scale_y_continuous(breaks=c(10,30,50,70,90), sec.axis =
                     sec_axis(~GWPfit(.),name="equivalent GWP year length",breaks=breaklist))

ggsave("results/doubleaxis.pdf", plot = secaxisplot, device = cairo_pdf, width = 6.4, height = 6.4, units = "cm", dpi = 600)


print("reached end of code")

