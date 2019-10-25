#' @get /figure2b
#' @png (width = 700) 
function(newValue0) { 
startAnalysis <- Sys.time() 
newValue0=as.numeric(newValue0) 
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
threshold = newValue0 # this can be measurement-specific ####
FFT_Tmain <- which(Mod(FFT_Temperature) > threshold) 
specs <- 1/(FFT_Tmain/N)   # period lengths of spectra selected by the threshold (T=1/f) [hour]
nspecs <- length(specs)
FFT_diel <- FFT_Tmain[which(specs>11 & specs<25)]              #  diel component
FFT_diel <- append(FFT_diel, FFT_Tmain[nspecs-which(is.element(FFT_Tmain,FFT_diel))+2])
FFT_Temperature_diel     <- replace(FFT_Temperature,c(-FFT_diel),0)
Temperature_diel     <- Re(fft(FFT_Temperature_diel,     inverse=TRUE))
plot(Time, Temperature_diel, type="l", col="red", ylab="Diel")

endAnalysis <- Sys.time() 
totaltime <- difftime(endAnalysis, startAnalysis, units="secs") 
message(paste("Total time is: ", totaltime)) 
}