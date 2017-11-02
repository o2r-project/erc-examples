# server.R

library(ggplot2)
library(plyr)
library(reshape2)
Streamflow <- read.csv("SyntheticStreamflow.csv")
#Q <- Streamflow[,"Q1"]

shinyServer(
  function(input, output, clientData, session) {
    source("shFun.R", local = TRUE)
    source("multiplot.R", local = TRUE)
    
    #Set simulation time period
    SU <- 5 
    TP <- 50 + SU             #Time period
    
    #Demand Parameters
    beta <- 1/1000           #Background efficiency increase [.]
    alpha <- 0.15            #Fractional adoption rate [.]
    Dmin <- 200/1000^3       #Minimum possible demand  [km^3/yr]
    
    #Shortage Memory Parameters
    muS <- 0.05              #Forgetting rate [.]
    Tp <- 0.4
    
    #Population Parameters
    deltaB <- 0.04           #Birth rate [.]    
    deltaD <- 0.03           #Death rate [.]
    deltaI <- 0.05           #Immigration rate [.]
    deltaE <- 0.03           #Emmigration rate [.]
    
    #Hydrologic Parameters
    muH <- 2.0               #Mean annual streamflow [L^3]
    sigmaH <- 0.50           #Standard deviation of streamfow [L^3]
    rhoH <- 0.60             #Lag 1 autocorrelation [.]
    etaH <- 0.001            #Evaporation rate [L]
    sigmaT <- 0.1            #Average slope of reservoir []
    Vmax <- 2.0              #Reservoir capacity [L^3]
    Kp <- 1                  #Hedging slope
    
    #Intial Conditions
    V1 <- 0.5                #Initial storage
    P1 <- 1000000            #Initial population
    M1 <- 0                  #Initial shortage memory
    D1 <- 400/1000^3         #Initial per cap. demand
    A1 <- 10                 #Reservoir base area [L^2]
    
    output$Plot <- renderPlot({
      Q <- Streamflow[,input$StreamSeq]
      
      ResSOP <- shFun(Q, TP, SU, V1, P1, D1, M1, A1, 
                      sigmaT, 1, etaH, deltaB, deltaD, deltaI, 
                      deltaE, muS, alpha, beta, Vmax, Dmin, Tp)
      ResSOP <- subset(ResSOP, Year > 0)
      
      ResHP <- shFun(Q, TP, SU, V1, P1, D1, M1, A1, 
                     sigmaT, input$Kp, etaH, deltaB, deltaD, deltaI, 
                     deltaE, muS, alpha, beta, Vmax, Dmin, Tp)
      ResHP <- subset(ResHP, Year > 0)
      
      colnames(ResHP) <- c("Year","StreamHP","VolHP","MemHP","DemHP","PopHP","TotHP")
      Res <- merge(ResSOP,ResHP,by="Year")
      
      Qdata <- Res[,c("Year","Stream")]
      p1 <- ggplot(data=Qdata, aes(x=Year, y=Stream)) +
        geom_line(size=1) +
        ggtitle(paste("Annual Streamflow")) +
        ylab("Streamflow (km3/yr)") +
        scale_y_continuous(limits=c(0, 3.5)) + 
        theme(legend.position="bottom")
      
      Voldata <- Res[,c("Year","Vol","VolHP")]
      Voldata <- melt(Voldata, id=c("Year")) 
      p2 <- ggplot(data=Voldata, aes(x=Year, y=value, group = variable)) +
        geom_line(aes(linetype=variable), size=1) +
        ggtitle(paste("Storage Volume")) +
        ylab("Volume (km3)") +
        scale_y_continuous(limits=c(0, 2)) + 
        scale_linetype_manual(values=c("solid", "dotted"), name="Operating\nPolicy", labels=c("SOP", "HP"))+
        theme(legend.position="bottom")
      
      Memdata <- Res[,c("Year","Mem","MemHP")]
      Memdata <- melt(Memdata, id=c("Year")) 
      p3 <- ggplot(data=Memdata, aes(x=Year, y=value, group = variable)) +
        geom_line(aes(linetype=variable), size=1) +
        ggtitle(paste("Shortage Awareness")) +
        ylab("Shortage Awareness") +
        scale_y_continuous(limits=c(0, 1)) + 
        scale_linetype_manual(values=c("solid", "dotted"), name="Operating\nPolicy", labels=c("SOP", "HP"))+
        theme(legend.position="bottom")
      
      Demdata <- Res[,c("Year","Dem","DemHP")]
      Demdata <- mutate(Demdata, Dem = Dem*1000^3, DemHP = DemHP*1000^3)
      Demdata <- melt(Demdata, id=c("Year")) 
      p4 <- ggplot(data=Demdata, aes(x=Year, y=value, group = variable)) +
        geom_line(aes(linetype=variable), size=1) +
        ggtitle(paste("Per Capita Demand")) +
        ylab("Per Cap. Demand (m3/yr)") +
        scale_y_continuous(limits=c(200, 400)) + 
        scale_linetype_manual(values=c("solid", "dotted"), name="Operating\nPolicy", labels=c("SOP", "HP"))+
        theme(legend.position="bottom")
      
      Popdata <- Res[,c("Year","Pop","PopHP")]
      Popdata <- mutate(Popdata, Pop = Pop/1000000, PopHP = PopHP/1000000)
      Popdata <- melt(Popdata, id=c("Year")) 
      p5 <- ggplot(data=Popdata, aes(x=Year, y=value, group = variable)) +
        geom_line(aes(linetype=variable), size=1) +
        ggtitle(paste("Population")) +
        ylab("Population (Millions)") +
        scale_linetype_manual(values=c("solid", "dotted"), name="Operating\nPolicy", labels=c("SOP", "HP"))+
        theme(legend.position="bottom")
      
      Totdata <- Res[,c("Year","Tot","TotHP")]
      Totdata <- melt(Totdata, id=c("Year")) 
      p6 <- ggplot(data=Totdata, aes(x=Year, y=value, group = variable)) +
        geom_line(aes(linetype=variable), size=1) +
        ggtitle(paste("Total Demand")) +
        ylab("Total Demand (km3/yr)") +
        scale_linetype_manual(values=c("solid", "dotted"), name="Operating\nPolicy", labels=c("SOP", "HP"))+
        theme(legend.position="bottom")
      
      Plot <- multiplot(p1, p2, p3, p4, p5, p6, cols=3)
      
    })
    
  })