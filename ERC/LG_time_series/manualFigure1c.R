library("forecast")
library("expsmooth") # required for the data
model = "ZZZ"
fit3 <- ets(ukcars, model)
plot(forecast(fit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")