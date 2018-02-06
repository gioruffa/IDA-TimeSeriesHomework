## IDA HW.3
library(forecast)
library(xlsx)
library(ggplot)
library(ggplot2)
library(astsa)
##Load the data from the file
NPO_Time_Series=read.xlsx("~/UniData/IDA/TimeSeries/data_g11.xlsx", sheetIndex = 1)
summary(NPO_Time_Series$Tasa)
class(NPO_Time_Series)
##Change the data to time series and frequency= 12 for month
##Start and end the begining and ending of the time data set
NPO_Time_Series = ts(data = NPO_Time_Series$Tasa, frequency = 12, start=c(2000,1) , end=c(2004,12), names=c("total credit") )

autoplot(NPO_Time_Series) # the first thing that is naively observed is that rather than one single
#trend we can see two trends, a descenting one until the beginning of 2003 and a rising one after

##DA STAMPARE
tsdisplay(NPO_Time_Series, plot.type="partial") #the ACF goes to zero very slowly, we need to normalize ths series
#the series is clearly non stationary but the signature may loook like an AR series
#it seems very difficult to dect seasonality or ciclicity
tsdisplay(NPO_Time_Series, plot.type="scatter")
#the variance does not seem to increase, so the log transformation shouldn't be very helpful
tsdisplay(NPO_Time_Series, plot.type="spectrum")

BoxCox.lambda(NPO_Time_Series,lower=0, upper=2) # value very close to zero, lets go for the log
tsdisplay(log(NPO_Time_Series), plot.type="scatter") #in fact they look pretty much the same

ggmonthplot(NPO_Time_Series) #huge variations but the mean is consistent 
ggseasonplot(NPO_Time_Series) #does not show any year season
#but we should get rid of the trend

#we will stay with additive decomposition since the log as no use
NPO.decomposed = decompose(NPO_Time_Series, type= "additive")
plot(NPO.decomposed) 

stl1.NPO <- stl(NPO_Time_Series, t.window=5,s.window="periodic", robust=TRUE)
plot(stl1.NPO)
stl2.NPO <- stl(NPO_Time_Series, t.window=15,s.window="periodic", robust=TRUE)
plot(stl2.NPO) #I like this one more ->  

stl.NPO <- stl(NPO_Time_Series, s.window="periodic", robust=TRUE)
stl.NPO #seasonal, trend, remainder
plot(stl.NPO)##DA STAMPARE

NPO.notrend = NPO_Time_Series - stl2.NPO$time.series[,"trend"]
par(mfrow=c(2,1))
plot(NPO_Time_Series)
plot(NPO.notrend)
plot(NPO.notrend - stl2.NPO$time.series[,"seasonal"])

par(mfrow=c(3,1))
plot(NPO_Time_Series)
plot(NPO_Time_Series - stl2.NPO$time.series[,"seasonal"])
#questo e' interessante-> se rimouviamo il trend la serie comincia a sembrare stazionaria
tsdisplay(NPO.notrend) 

stl2.NPO

spec.pgram(NPO_Time_Series) 

library("tseries")

adf.test(NPO.notrend)
##
##  stl forecasting...
##
#plotting the original series and the seasonal adjusted one, to see the trend
par(mfrow=c(2,1))
plot(NPO_Time_Series)
plot(seasadj(stl2.NPO))

fit1 <- ses(NPO_Time_Series, alpha=0.2, initial="simple", h=3)
fit2 <- ses(NPO_Time_Series, alpha=0.6, initial="simple", h=3)
fit3 <- ses(NPO_Time_Series, h=3)

fit1 <- holt(NPO_Time_Series, alpha=0.8, beta=0.2, initial="simple", h=5) 
fit2 <- holt(NPO_Time_Series, alpha=0.8, beta=0.2, initial="simple", exponential=TRUE, h=5) 
# Results for first model:
fit1$model$state
fitted(fit1)
fit1$mean

fit3 <- holt(NPO_Time_Series, alpha=0.8, beta=0.2, damped=TRUE, initial="simple", h=5) 
par(mfrow=c(1,1))
plot(fit2, type="o", ylab="Air passengers in Australia (millions)", xlab="Year", 
     fcol="white", plot.conf=FALSE)
lines(fitted(fit1), col="blue") 
lines(fitted(fit2), col="red")
lines(fitted(fit3), col="green")
lines(fit1$mean, col="blue", type="o") 
lines(fit2$mean, col="red", type="o")
lines(fit3$mean, col="green", type="o")
legend("top", lty=1, col=c("black","blue","red","green"), 
       c("Data","Holt's linear trend","Exponential trend","Additive damped trend"))

##
## Holt Winters -> good for seasonality so it may not work here
##
NPO.HW<- HoltWinters(NPO_Time_Series)
NPO.HW
names(NPO.HW)
NPO.HW$coefficients
head(NPO.HW$fitted)
plot(NPO.HW) #or
plot(NPO.HW, ylab="Credit", type="l")
lines(NPO.HW$fitted[,1], col="red")

pred.NPO.HW <- predict(NPO.HW, n.ahead=24, prediction.interval=TRUE)
plot(NPO.HW, pred.NPO.HW, ylab="Credit", xlab="", main="")
labs <- c("Observed values", "Fitted values", "Predicted intervals")
legend("bottomleft", lty=rep(1,3), col=c("black", "red", "blue"), legend=labs)

wplot(pred.NPO.HW)

##
##  box jenkins
##

#lets make the series stationary.
acf2(NPO_Time_Series)
tsdisplay(diff(NPO_Time_Series))
tsdisplay(diff(diff(NPO_Time_Series)))
sd(diff(NPO_Time_Series))
sd(diff(diff(NPO_Time_Series)))

adf.test(diff(NPO_Time_Series)) #ok
pp.test(NPO_Time_Series) #ok
pp.test(diff(NPO_Time_Series)) #ok
ndiffs(NPO_Time_Series)
#Ok we are sure that differentiating to 1 is ok to obtain a stationary serie

#happy with first order, still we do not have a seasonal component
NPO.firstdiff = diff(NPO_Time_Series)
tsdisplay(NPO.firstdiff) #as expected differenciation gets rid of the trend, so the degree of differencing should be 1
tsdisplay(NPO.firstdiff, plot.type = "spectrum")

acf2(NPO_Time_Series) #The ACF looks like an AR with dumped sine, the PACF cutoff seems to be 1 and it is positive
#So we can go for AR(1) or AR(2)
#I honestly do not see anything in the PACF plot that suggest an MA component
NPO.ar1=Arima(NPO_Time_Series, order=c(1,1,0))
NPO.ar2=Arima(NPO_Time_Series, order=c(2,1,0))

summary(NPO.ar1)
summary(NPO.ar2) #slightly better MAE and RMSE

#we must check the residuals to be "white noise", lets seee
par(mfrow=c(2,1))
plot(NPO.ar1$residuals)
plot(NPO.ar2$residuals)
#they basically look the same -> it is clear that is a low (-3) value as an outlier

acf2(NPO.ar1$residuals)
acf2(NPO.ar2$residuals)
#they are both below zero but in ar2 they seem more "random". The Box-lyung test is done for that
Box.test(NPO.ar1$residuals, lag=12, fitdf = 1)
Box.test(NPO.ar2$residuals, lag=12, fitdf = 2) #degree of fredom suggested depending on p+q, see help
jarque.bera.test(NPO.ar1$residuals) #normality failed
jarque.bera.test(NPO.ar2$residuals) #normality failed

#it may be for some outliers
which.max(abs(NPO.ar1$residuals))
which.max(abs(NPO.ar1$residuals[-6]))

hist(NPO.ar1$residuals,breaks = 20)
hist(NPO.ar1$residuals[-6],breaks = 20)

boxplot(NPO.ar1$residuals)
library(car)
qqPlot(NPO.ar1$residuals)
library(outliers)
#get only outliers

NPO_Time_Series_clean_ar1 = NPO_Time_Series[!outlier(NPO.ar1$residuals, logical=TRUE)]
NPO.clean.ar1 = Arima(NPO_Time_Series_clean_ar1, order=c(1,1,0))
jarque.bera.test(NPO.clean.ar1$residuals) #not fantastic but ok

NPO_Time_Series_clean_ar2 = NPO_Time_Series[!outlier(NPO.ar2$residuals, logical=TRUE)]
NPO.clean.ar2 = Arima(NPO.clean.ar2, order=c(2,1,0))
jarque.bera.test(NPO.clean.ar2$residuals) #not fantastic but ok

#
# Let's try with an auto model
#

best.NPO=auto.arima(NPO_Time_Series,d=1,max.p=10,max.q=10)
jarque.bera.test(best.NPO$residuals)

which.max(abs(best.NPO$residuals))
best.clean.NPO=auto.arima(NPO_Time_Series[-6],d=1,max.p=10,max.q=10)
jarque.bera.test(best.clean.NPO$residuals)
qqPlot(best.clean.NPO$residuals)

summary(best.clean.NPO)
summary(NPO.clean.ar1)
summary(NPO.clean.ar2)

acf2(best.clean.NPO$residuals)
acf2(NPO.clean.ar1$residuals)
acf2(NPO.clean.ar2$residuals) #has best autocorrelation and best jarque brera non automatic
qqPlot(NPO.clean.ar2$residuals)


#
# ARIMA (AR in this case, forecasting)
#

plot(forecast(NPO.clean.ar2, h=20))

getrmse <- function(x,h,...)
{
  train.end <- time(x)[length(x)-h]   #train data end
  test.start <- time(x)[length(x)-h+1]  #test data start
  train <- window(x,end=train.end) #extract train data
  test <- window(x,start=test.start)  #extract test data
  fit <- Arima(train,...) # fit model with train data
  fc <- forecast(fit,h=h) # forecast with model
  return(accuracy(fc,test)[2,"RMSE"]) #compare forecast with test data, extract the rmse
}

getrmse(NPO_Time_Series[-6],h=20,order=c(0,1,0))
getrmse(NPO_Time_Series[-6],h=20,order=c(1,1,0))
getrmse(NPO_Time_Series[-6],h=20,order=c(2,1,0)) #still the best
getrmse(NPO_Time_Series[-6],h=20,order=c(0,1,1))
getrmse(NPO_Time_Series[-6],h=20,order=c(2,1,1)) #the actual best one


#and the winner is -> ar(2)
#justify that we pick it even if its more complicated

NPO_Time_Series
