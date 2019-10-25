#FUNCTION DEFINITIONS

# libraries used in this file:
# library(stats)    #for uniroot 
# library(caTools)  #for trapz, as it is closer to integrate than a straight sum.

#CO2 functions

#CO2 mass from a pulse emission based on Table 8.SM.10 from IPCC AR5
co2mass <- function(yearsincepulse)
{
  return(0.2173+0.2240*exp(-yearsincepulse/394.4)+0.2824*exp(-yearsincepulse/36.54)+0.2763*exp(-yearsincepulse/4.304))
}

#CO2 forcing based on Table 8.SM.1, IPCC AR5
#Mass to ppm conversion based on CO2-C: 12/29*5.13 = 2.12 Gt/ppm
co2forc <- function(pulse, co2conc, yearsincepulse)
{
  return(5.35*log((co2conc[yearsincepulse+1]+pulse*co2mass(yearsincepulse)/2.12)/co2conc[yearsincepulse+1])) 
}

#CH4 functions

#forcing functions based on Table 8.SM.1
ftotcalc <- function(ch4new, ch4old, n2o)
{
  return(0.036*(sqrt(ch4new)-sqrt(ch4old))-(fcalc(ch4new,n2o)-fcalc(ch4old,n2o)))
}

fcalc <- function(ch4, n2o)
{
  return(0.47*log(1+0.0000201*(ch4*n2o)^0.75+0.00000000000000531*ch4*(ch4*n2o)^1.52))
}

#12.4 year lifetime based on Table 8.A.1
ch4mass <- function(yearsincepulse)
{
  return(exp(-yearsincepulse/12.4)) 
}

#molecular weight of air: 29 g/mole
#mass of atmosphere: 5.13 x 10^18 kg
#Mass to ppb conversion: 16/29*5.13 = 2.83 Mt/ppb
#include factor of 1.65 for strat water vapor & tropospheric ozone
ch4forc <- function(pulse, ch4conc, n2oconc, yearsincepulse)
{
  return(ftotcalc(ch4conc[yearsincepulse+1]+pulse*ch4mass(yearsincepulse)/2.83, ch4conc[1], n2oconc[yearsincepulse+1])*1.65)
}

#HFC functions

#hfc134a lifetime from Table 8.A.1 is 13.4 years
#but function is designed for arbitrary lifetime to emulate any hfc/pfc/sf6
hfcmass <- function(yearsincepulse, lifetime)
{
  return(exp(-yearsincepulse/lifetime))
}

#hfc134a conversion: 102/29*5.13 = 18.04. So mass/18.04 = ppb
#0.16 is radiative efficiency per ppb from Table 8.A.1
#however, mass & radiative efficiency don't matter for GWP timescale purposes
#even though they matter for the actual damage ratio. 
hfcforc <- function(pulse, yearsincepulse, lifetime)
{
  return((pulse*hfcmass(yearsincepulse, lifetime)/18.04)*0.16)
}


#N2O functions
#forcing functions based on Table 8.SM.1
n2oftotcalc <- function(ch4new, n2onew, n2oold)
{
  return(0.12*(sqrt(n2onew)-sqrt(n2oold))-(fcalc(ch4new,n2onew)-fcalc(ch4new,n2oold)))
}

#121 year lifetime based on Table 8.A.1
n2omass <- function(yearsincepulse)
{
  return(exp(-yearsincepulse/121)) 
}

#molecular weight of air: 29 g/mole
#mass of atmosphere: 5.13 x 10^18 kg
#Mass to ppb conversion: 44/29*5.13 = 7.78 Mt/ppb
#account for CH4 indirect effect: 1-0.36*(1.65)*RE(CH4)/RE(N2O)
#assume RE(CH4) and RE(N2O) are constant and from Appendix 8.A:  =0.363/3
#so factor is approx 0.928
n2oforc <- function(pulse, ch4conc, n2oconc, yearsincepulse)
{
  return(n2oftotcalc(ch4conc[yearsincepulse+1], n2oconc[yearsincepulse+1]+pulse*n2omass(yearsincepulse)/7.78, n2oconc[1])*0.928)
}



#Different options to calculate net present value of a vector over time given a discount rate
#Integrate itself only worked in limited cases.
#trapsum more closely emulated integrate than discountsum
#for constant concentration scenarios, particularly for temperature integration. 
#In the future, possibly look at simpson's rule as an alternative

discountsum <- function(disvector, discountrate)
{
  dissum <- 0
  for(i in 0:(length(disvector)-1))
  {
    dissum <- dissum + disvector[i+1]/(discountrate^i)
  }
  return(dissum)
}

trapsum <- function(disvector, discountrate)
{ 
  dissum <- vector()
  for(i in 1:(length(disvector)))
  {
    dissum[i] <- disvector[i]/(discountrate^(i-1))
  } 
  return(trapz(1:length(dissum),dissum[1:length(dissum)]))
  #alternate exact integral:
  #return(sintegral(1:length(dissum), dissum[1:length(dissum)]))
} 


#trapsum2 is faster than trapsum, yields the same results
trapsum2 <- function(disvector, discountrate) {
  last_i <- length(disvector)
  i <- 1:last_i
  discount_factor <- 1 / (discountrate^(i-1))
  dissum <- disvector * discount_factor
  
  # trapezoid rule with special case of every year
  sum(dissum) - (dissum[1] + dissum[last_i])/2 
}



#temperature calculations: see 8.SM.11.2 from AR5
#forcing imbalance should be the difference between the existing imbalance in 2011,
#and the RCP forcing in 2011 (before any perturbation). Literature estimates:
#0.62 W/m2 (http://onlinelibrary.wiley.com/doi/10.1002/2014GL060962/full)
#0.58 W/m2 (https://www.giss.nasa.gov/research/briefs/hansen_16/)
#http://onlinelibrary.wiley.com/doi/10.1002/2015JD023264/full:
#0.5 to 1W/m2 based on Hansen et al., 2011; Loeb et al., 2012; Trenberth et al., 2014; Wild et al., 2015
#http://www.nature.com/nclimate/journal/v6/n2/full/nclimate2876.html:
#Also, 0.5 to 1 w/m2 based on Hansen, Loeb, Trenberth, Allan
#Medhaug et al. https://www.nature.com/nature/journal/v545/n7652/full/nature22315.html
#0.75 to 0.93 W/m2 
#climate sensitivity = (sum of coefficients)*(W/m2 for doubling CO2 = 3.7)
#IPCC default is therefore (0.631 + 0.429)*3.7 = 3.92
#Scaling to 1.5 would be: 0.241 + 0.164 
#Scaling to 4.5 woudl be: 0.728 + .492


tempcalc <- function(yearsincepulse, forc, forcingimbalance, CSiter) #update to include heat uptake?
{
  CSfactor <- CSiter/3.92
  
  totalTemp <- numeric(length(yearsincepulse)) #set totalTemp to be a vector of zeros
  for (iter in 0:length(yearsincepulse))
    totalTemp <- totalTemp +  ifelse(iter > yearsincepulse, 0, (forc[iter+1] + forcingimbalance)*
                                       (0.631*CSfactor/8.4*exp((iter-yearsincepulse)/8.4)+0.429*CSfactor/409.5*exp((iter-yearsincepulse)/409.5)))
  
  return(totalTemp)
}

#calculating the CH4 GWP in one call
#assumes constant concentrations because of the integrate function
calcGWP <- function(gwpyears, co2emi, ch4emi, co2concin, ch4concin, n2oconcin)
{
  co2agwp <- integrate(co2forc, lower=0, upper=gwpyears, pulse = co2emi, co2conc = co2concin)
  ch4agwp <- integrate(ch4forc, lower = 0, upper=gwpyears, pulse = ch4emi, ch4conc = ch4concin, n2oconc = n2oconcin)
  return(ch4agwp$val/co2agwp$val)
}

#calculating the HFC GWP in one call
#assumes constant concentrations because of the integrate function
calcHFCGWP <- function(gwpyears, co2emi, hfcemi, co2concin, lifetime)
{
  co2agwp <- integrate(co2forc, lower=0, upper=gwpyears, pulse = co2emi, co2conc = co2concin)
  hfcagwp <- integrate(hfcforc, lower = 0, upper=gwpyears, pulse = hfcemi, lifetime = lifetime)
  return(hfcagwp$val/co2agwp$val)
}

#calculating the N2O GWP in one call
#assumes constant concentrations because of the integrate function
calcN2OGWP <- function(gwpyears, co2emi, n2oemi, co2concin, ch4concin, n2oconcin)
{
  co2agwp <- integrate(co2forc, lower=0, upper=gwpyears, pulse = co2emi, co2conc = co2concin)
  n2oagwp <- integrate(n2oforc, lower = 0, upper=gwpyears, pulse = n2oemi, ch4conc = ch4concin, n2oconc = n2oconcin)
  return(n2oagwp$val/co2agwp$val)
}

#calculating the number of years one would need for a given CH4/CO2 ratio
#NOTE: DOES NOT SOLVE FOR ch4co2ratio < 13 (GWP years of 274), or ch4co2ratio>95 (GWP years of 14.5)
#Second Note: if the gas lifetime is different, the ratio limits will be different!
#the latter limit makes sense: as the ch4co2 ratio approaches the relative radiative forcing,
#solutions will be hard. 
#while this is defined by limits, function doesn't seem to work much beyond these points...
#may only work for constant concentrations
calcGWProot <- function(ch4co2ratio, co2emi, ch4emi, co2concin, ch4concin, n2oconcin)
{
  rootvalue <- ifelse(ch4co2ratio<13,275, ifelse(ch4co2ratio>95, 14, 
                                       uniroot(function(y) calcGWP(y, co2emi, ch4emi, co2concin, ch4concin, n2oconcin)-ch4co2ratio, lower=14, upper=280)))
  return(unlist(rootvalue))
}

#SOLVE LIMITS:
#works down to a small timescale of 1 (or smaller) and up to 490
#challenge: if lifetime increases, function can become non-monotonic!
calcHFCGWProot <- function(hfcco2ratio, co2emi, hfcemi, co2concin, lifetime)
{
  smalltimescale <- 20
  bigtimescale <- 250
  lowerlimit <- calcHFCGWP(smalltimescale, co2emi, hfcemi, co2concin, lifetime)
  upperlimit <- calcHFCGWP(bigtimescale, co2emi, hfcemi, co2concin, lifetime)
  rootvalue <- ifelse(hfcco2ratio<upperlimit,bigtimescale, ifelse(hfcco2ratio>lowerlimit, smalltimescale, 
                                                 uniroot(function(y) calcHFCGWP(y, co2emi, hfcemi, co2concin, lifetime)-hfcco2ratio, lower=smalltimescale, upper=bigtimescale)))
  return(unlist(rootvalue))
}
