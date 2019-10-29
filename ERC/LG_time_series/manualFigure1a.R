library("forecast")
library("expsmooth") # required for the data
model = "ZZZ"
fit1 <- ets(bonds, model)
plot(forecast(fit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")