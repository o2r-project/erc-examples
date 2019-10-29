library("forecast")
library("expsmooth") # required for the data
aafit1 <- auto.arima(bonds,max.P=10,max.Q=0,D=0,approximation=FALSE)
plot(forecast(aafit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")