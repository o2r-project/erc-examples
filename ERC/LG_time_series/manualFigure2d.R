library("forecast")
library("expsmooth") # required for the data
aafit4 <- auto.arima(log(visitors),approximation=FALSE)
plot(forecast(aafit4),lambda=0,xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")