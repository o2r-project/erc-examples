library("forecast")
library("expsmooth") # required for the data
model = "ZZZ"
fit4 <- ets(visitors, model)
plot(forecast(fit4),xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")