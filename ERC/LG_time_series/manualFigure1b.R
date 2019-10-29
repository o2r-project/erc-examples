library("forecast")
library("expsmooth") # required for the data
model = "ZZZ"
etsfit <- ets(usnetelec, model)
plot(forecast(etsfit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")