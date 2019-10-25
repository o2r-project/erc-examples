#' @get /figure2a
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
# Note: threshold can be arbitrary, try different values and visually evaluate from figures
# a recommendation is to set it firstly to the value of minimum unit of measurement  

### selecting spectra whose intensity is higher than the threshold
### assuming they dominantly constitute seasonal and diel cycles
FFT_Tmain <- which(Mod(FFT_Temperature) > threshold) 
specs <- 1/(FFT_Tmain/N)   # period lengths of spectra selected by the threshold (T=1/f) [hour]
nspecs <- length(specs)

FFT_seasonal <- FFT_Tmain[which(specs>30*24 & specs<367*24)] # extracting seasonal component
FFT_seasonal <- append(FFT_seasonal, FFT_Tmain[nspecs-which(is.element(FFT_Tmain,FFT_seasonal))+2])  # also corresponding ones
FFT_Temperature_seasonal <- replace(FFT_Temperature,c(-FFT_seasonal),0)
Temperature_seasonal <- Re(fft(FFT_Temperature_seasonal, inverse=TRUE))
plot(Time, Temperature_seasonal, type="l", col="red", ylab="Seasonal")

endAnalysis <- Sys.time() 
totaltime <- difftime(endAnalysis, startAnalysis, units="secs") 
message(paste("Total time is: ", totaltime)) 
}