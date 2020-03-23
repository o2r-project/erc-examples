#!/usr/bin/env Rscript
#plot GWP value vs. timescale for CH4 and for HFCs of different lifetimes
#produces Figure SI.2 for the Sarofim & Giordano ESD paper

#Code expects files for RCP85 scenario

library(ggplot2)  #for graphing
library(stats)    #for uniroot 
library(reshape2) #for melt command, which is used for graphing data format
library(caTools)  #for trapz, as it is closer to integrate than a straight sum.
library(gridExtra)
library(grid)

#import functions
source("src/metricsfunctions.R")

concentrations <- read.csv("data/rcp85gasconcs.csv", header=TRUE, stringsAsFactors = FALSE) 
constantco2 <- rep(concentrations$co2conc[1], times = length(concentrations$co2conc))
constantch4 <- rep(concentrations$ch4conc[1], times = length(concentrations$ch4conc))
constantn2o <- rep(concentrations$n2oconc[1], times = length(concentrations$n2oconc))

ch4emission <- 28.3 #28.3 Mt should be 10 ppb
co2emission <- ch4emission*12/44/1000 
#1000 for Mt -> Gt
#12/44 conversion to get from GtCO2 -> GtC
#needed because we use 2.12 factor in CO2 function which is based on GtC
hfcemission <- ch4emission
n2oemission <- ch4emission

hfclifetime <- 13.4
timescale <- c(60,120,180)

calcGWP(100, co2emission, ch4emission, constantco2, constantch4, constantn2o)
calcHFCGWP(100, co2emission, hfcemission, constantco2, hfclifetime)

timescalelist <- c(seq(1,300, by = 1))
CH4resultmatrix <- matrix(ncol=2, nrow = length(timescalelist))
CH4resultframe <- data.frame(CH4resultmatrix)
CH4resultframe[,1] <- timescalelist
for (timescale in timescalelist)
{
  CH4resultframe[timescale,2] <- calcGWP(CH4resultframe[timescale,1], co2emission, ch4emission, constantco2, constantch4, constantn2o)
}

N2Oresultmatrix <- matrix(ncol=2, nrow = length(timescalelist))
N2Oresultframe <- data.frame(N2Oresultmatrix)
N2Oresultframe[,1] <- timescalelist
for (timescale in timescalelist)
{
  N2Oresultframe[timescale,2] <- calcN2OGWP(N2Oresultframe[timescale,1], co2emission, n2oemission, constantco2, constantch4, constantn2o)
}

#Figure SI.2
n2oaxis <- expression(paste("N"[2], "O GWP", sep=""))
n2otitle <- expression(paste("N"[2], "O GWP", sep=""))
n2oplot <- ggplot(N2Oresultframe, aes(x=X1, y=X2)) + geom_point() + 
  labs(title=n2otitle, x="Timescale", y=n2oaxis)+
  theme(text=element_text(family="Arial", size = 16), plot.title=element_text(size=16))

ch4axis <- expression(paste("CH"[4], " GWP", sep=""))
ch4title <- expression(paste("CH"[4], " GWP", sep=""))
ch4plot <- ggplot(CH4resultframe, aes(x=X1, y=X2)) + geom_point() + 
  labs(title=ch4title, x="Timescale", y=ch4axis)+
  theme(text=element_text(family="Arial", size = 16), plot.title=element_text(size=16))

ggsave("results/N2OGWP.pdf", plot = n2oplot, device = cairo_pdf, width = 6.4, height = 6.4, units = "cm", dpi = 600)
ggsave("results/CH4GWP.pdf", plot = ch4plot, device = cairo_pdf, width = 6.4, height = 6.4, units = "cm", dpi = 600)


