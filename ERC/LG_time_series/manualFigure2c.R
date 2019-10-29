library("forecast")
library("expsmooth") # required for the data
aafit3 <- auto.arima(ukcars,approximation=FALSE)
plot(forecast(aafit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")