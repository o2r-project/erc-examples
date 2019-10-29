## CODE FOR PAPER: Hyndman and Khandakar (JSS, 2008)
#install.packages(c("forecast", "expsmooth"))
library("forecast")
library("expsmooth") # required for the data

## EXPONENTIAL SMOOTHING

fit1 <- ets(bonds)
etsfit <- ets(usnetelec)
fit3 <- ets(ukcars)
fit4 <- ets(visitors)

# Figure 1
par(mfrow=c(2,2))
plot(forecast(fit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")
plot(forecast(etsfit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")
plot(forecast(fit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")
plot(forecast(fit4),xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")

# Print the model for US net electricity
etsfit
# In-sample accuracy of model
accuracy(etsfit)
# Compute forecasts from the model
fcast <- forecast(etsfit)
# Plot the forecasts with prediction intervals
plot(fcast)
# Print the forecasts
fcast
# Fit an ETS model to the first 45 months
fit <- ets(usnetelec[1:45])
# Apply the same model to the last 10 months without re-estimating the parameters
test <- ets(usnetelec[46:55],model=fit)
# Look at the in-sample accuracy
accuracy(test)
# Look at the out-of-sample accuracy
accuracy(forecast(fit,10), usnetelec[46:55])


## ARIMA MODELLING

aafit1 <- auto.arima(bonds,max.P=0,max.Q=0,D=0,approximation=FALSE)
arimafit <- auto.arima(usnetelec)
aafit3 <- auto.arima(ukcars,approximation=FALSE)
aafit4 <- auto.arima(log(visitors),approximation=FALSE)

# Forecast from the ARIMA model
fcast <- forecast(arimafit)
# Plot the forecasts
plot(fcast)
# Summary of the forecasts
summary(fcast)

# Figure 2
par(mfrow=c(2,2))
plot(forecast(aafit1),xlab="Year",ylab="Percentage per annum",main="(a) US 10-year bonds yield")
plot(forecast(arimafit),xlab="Year",ylab="Billion kwh",main="(b) US net electricity generation")
plot(forecast(aafit3),xlab="Year",ylab="Thousands of cars",main="(c) UK passenger motor vehicle production")
plot(forecast(aafit4),lambda=0,xlab="Year",ylab="Thousands of people",main="(d) Overseas visitors to Australia")

# COMPARISONS
summary(forecast(etsfit))
summary(forecast(arimafit))
