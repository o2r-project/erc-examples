#Launching Reservoir Operations model

#Specify the working directory
setwd("P:/Research/Modeling/SunshineCity")

#Install shiny package if not previously installed
install.packages("shiny")

#load package
library(shiny)

#launch model
runApp("ReservoirOperations")
