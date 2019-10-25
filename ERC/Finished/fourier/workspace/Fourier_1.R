##########################################################################
#
#  Fast Fourier Transform of hourly time-series temperature data to 
#    1) detect key spectra which largely compose the variability
#    2) decompose it to mean, seasonal and diel cycles, and noise
#
# Input:   an hourly recorded temperature data (text file)
#           1st column for "Time" and 2nd column for "Temperature"
#           "example_data1.txt" is analyzed.
# parameter: threshold value for noise removal (threshold, ?C)
# Output: figure for summary and decomposed temperature as csv
#
#                    2016.07    Masahiro Ryo
#                         masahiroryo@gmail.com
#
#
##########################################################################

###############################################################
###     PREPARATION PART
###############################################################

#setwd("C:/")
data<-read.table("example_data1.txt",sep="\t", header=TRUE)
Time <-as.POSIXct(strptime(data$Time,format="%Y.%m.%d %H:%M"))
  # default format for time: YYYY.MM.DD HH:MM
  # e.g. 2013.05.17 23:00
N <- nrow(data)

##########################################################################
# FFT (Fast Fourier Transform) to capture seasonal and dial periodicity
###########################################################################
### Fourier Transform
FFT_Temperature <- fft(data$Temperature)/N  
FFT_max <- max(Mod(FFT_Temperature[2:(N/2)]))


### Extracting only key spectrum for modeling
 ################################################################
 # Tips: what should be interpreted?                            #
 #   What are the key seasonalities in data?                    #
 #    -> check the length of key periodic cycles (specs_day)    #
 #   How strong are they?                                       #
 #    -> check the absolute value of corresponding coefficient  #
 #       (abs(FFT_Temperature))                             #
 #   How much do they compose data variability?                 #
 #    -> check the correlation and visualized plot              #
 #       after modeling used only the key seasonalities         #
 ################################################################
### threshold selection #######################
threshold <- 0.2 # this can be measurement-specific ####
# Note: threshold can be arbitrary, try different values and visually evaluate from figures
# a recommendation is to set it firstly to the value of minimum unit of measurement  

### selecting spectra whose intensity is higher than the threshold
### assuming they dominantly constitute seasonal and diel cycles
FFT_Tmain <- which(Mod(FFT_Temperature) > threshold) 
specs <- 1/(FFT_Tmain/N)   # period lengths of spectra selected by the threshold (T=1/f) [hour]
nspecs <- length(specs)

FFT_seasonal <- FFT_Tmain[which(specs>30*24 & specs<367*24)] # extracting seasonal component
FFT_seasonal <- append(FFT_seasonal, FFT_Tmain[nspecs-which(is.element(FFT_Tmain,FFT_seasonal))+2])  # also corresponding ones
FFT_diel <- FFT_Tmain[which(specs>11 & specs<25)]              #  diel component
FFT_diel <- append(FFT_diel, FFT_Tmain[nspecs-which(is.element(FFT_Tmain,FFT_diel))+2])

### removing noises whose spectral intensity is lower than the threshold
FFT_Temperature_selected <- replace(FFT_Temperature,c(-FFT_Tmain),0)    
FFT_Temperature_seasonal <- replace(FFT_Temperature,c(-FFT_seasonal),0)
FFT_Temperature_diel     <- replace(FFT_Temperature,c(-FFT_diel),0)

### transforming from frequency field to time field
Temperature_selected  <- Re(fft(FFT_Temperature_selected, inverse=TRUE))  
Temperature_seasonal <- Re(fft(FFT_Temperature_seasonal, inverse=TRUE))
Temperature_diel     <- Re(fft(FFT_Temperature_diel,     inverse=TRUE))
ext.factor <- data$Temperature - abs(Temperature_selected) # noise (external factor)
# Note: Only Real part of the inversed FFT is needed in time series


### Plotting and comparing the FFT estimation and the observed data
par(mfrow=c(3,2))
par(mar=c(2,5,1,1))
plot(Mod(FFT_Temperature[1:(N/2)]), ylab="spectral intensity")
plot(Time, data$Temperature, type="l", col="black", ylab="Observed")
plot(Time, Temperature_selected, type="l", col="red", ylab="Mean + Seasonal + Diel" )
plot(Time, Temperature_seasonal, type="l", col="red", ylab="Seasonal")
plot(Time, Temperature_diel, type="l", col="red", ylab="Diel")
plot(Time, ext.factor, type="l",col="gray", ylab="external factor (noise)")

### exporting results as csv:
result <- data.frame(cbind(as.character(data$Time), data$Temperature, Temperature_selected,
                           Temperature_seasonal,Temperature_diel, ext.factor))
colnames(result) <- c("Time", "Observed", "Modeled", "Seasonal", "Diel", "Noise")
write.csv(result, "FFT_estimation_result.csv")
