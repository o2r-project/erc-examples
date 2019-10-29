#' @get /figure2
#' @png 
function(newValue){ 
  startAnalysis <- Sys.time() 
  library("forecast")
  library("expsmooth") # required for the data
  approx = as.logical(newValue)
  aafit1 <- auto.arima(bonds,max.P=0,max.Q=0,D=0,approximation=approx)
  arimafit <- auto.arima(usnetelec)
  aafit3 <- auto.arima(ukcars,approximation=approx)
  aafit4 <- auto.arima(log(visitors),approximation=approx)
  par(mfrow=c(2,2))
  plot(forecast(aafit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")
  plot(forecast(arimafit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")
  plot(forecast(aafit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")
  plot(forecast(aafit4),lambda=0,xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")
  
  endAnalysis <- Sys.time() 
  totaltime <- difftime(endAnalysis, startAnalysis, units="secs") 
  message(paste("Total time is: ", totaltime)) 
}