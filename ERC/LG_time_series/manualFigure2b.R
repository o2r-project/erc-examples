library("forecast")
library("expsmooth") # required for the data
arimafit <- auto.arima(usnetelec)
plot(forecast(arimafit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")