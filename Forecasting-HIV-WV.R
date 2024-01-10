#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                   R SCRIPT                                  # 
#       Forecasting HIV Cases in Western Visayas Using Google Trends Data     #
#        By: Daniel David M. Pamplona, UP Visayas, dmpamplona@up.edu.ph       #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#### PART 1. EXPLORATORY DATA ####

# Load Data
library(lubridate)
setwd("")  # Specify file location here
masterdata <- read.csv(".\\HIV-GTD-data.csv")
masterdata$Time <- ym(masterdata$Time)

# Plot of HIV Data and Google Trends Data
ggplot(data = masterdata, aes(x = Time)) +
  geom_point(aes(y = hiv, colour = "HIV"), shape=17, alpha= 0.5) + 
  geom_smooth(aes(y = hiv, colour = "HIV"), span = 0.2) +
  geom_point(aes(y = gtd, colour = "GTD"), shape=16, alpha= 0.5) + 
  geom_smooth(aes(y = gtd, colour = "GTD"),linetype="dashed", span = 0.2) +
  scale_colour_manual("", breaks = c("HIV", "GTD"), values = c("#d32c3d","#316690")) +
  scale_x_date(date_minor_breaks = "1 month")  +
  xlab(NULL) +
  ylab("Monthly HIV Cases and Google Trends Data") +
  theme_grey()

# Correlation between HIV and GTD
gtd <- ts(masterdata$gtd, start=c(2012, 1), end=c(2022, 12), frequency=12)
hiv <- ts( masterdata$hiv , start=c(2012, 1), end=c(2022, 12), frequency=12)
period <- ts(cbind(hiv,gtd))
gtd.l <- cbind(
  AdLag0 = gtd,
  AdLag1 = stats::lag(gtd,-1), 
  AdLag2 = stats::lag(gtd,-2),
  AdLag3 = stats::lag(gtd,-3)) %>% 
  head(NROW(period))

cor.test(hiv, gtd.l[,"AdLag0"], method="kendall")
cor.test(hiv, gtd.l[,"AdLag1"], method="kendall")
cor.test(hiv, gtd.l[,"AdLag2"], method="kendall")
cor.test(hiv, gtd.l[,"AdLag3"], method="kendall")

# Correlation between HIV and GTD Before Pandemic
gtd <- ts(masterdata[1:97,]$gtd, start=c(2012, 1), end=c(2020,1), frequency=12)
hiv <- ts( masterdata[1:97,]$hiv , start=c(2012, 1), end=c(2020, 1), frequency=12)
cor.test(hiv, gtd, method="kendall")


#### PART 2: DATA CORRECTION ####

# Method 1: Outlier Replacement via Linear Interpolation
library(forecast)
lambda <- BoxCox.lambda(masterdata$hiv + 0.5)
outliers <- tsoutliers(masterdata$hiv,lambda=lambda)
masterdata2 <- masterdata
masterdata2$hiv[outliers$index] <- outliers$replacements

# Method 2: Outlier Replacement via MA Imputation
library(imputeTS)
masterdata3 <- masterdata             
masterdata3$hiv[outliers$index] <- NA
masterdata3$hiv <- na_ma(masterdata3$hiv, weighting = "simple")

# Plot the Results of Outlier Replacement
library(ggplot2)
out <- masterdata$Time[outliers$index]
# Original Data
p1 <- ggplot(data=masterdata, aes(y=hiv, x=Time)) +
  geom_line(linewidth = 0.6, color="#8E1E29") +
  scale_x_date(date_minor_breaks = "1 month") +
  geom_vline(xintercept=out, linetype="dashed", color="#316690",
             alpha = 0.8, linewidth=0.3) + ylim(-10,235) +
  xlab("") + ylab("HIV (Uncorrected)") +
  theme_grey()
# Linear Interpolation
p2 <- ggplot(data=masterdata2, aes(y=hiv, x=Time)) +
  geom_line(linewidth = 0.6, color="#8E1E29") +
  scale_x_date(date_minor_breaks = "1 month") +
  geom_vline(xintercept=out, linetype="dashed", color="#316690",
             alpha = 0.8, linewidth=0.3) + ylim(-10,235) +
  xlab("") + ylab("HIV (Linear Interpolation)") +
  theme_grey()
# MA Imputation
p3 <- ggplot(data=masterdata3, aes(y=hiv, x=Time)) +
  geom_line(linewidth = 0.6, color="#8E1E29") +
  scale_x_date(date_minor_breaks = "1 month") +
  geom_vline(xintercept=out, linetype="dashed", color="#316690",
             alpha = 0.8, linewidth=0.3) + ylim(-10,235) +
  xlab("") + ylab("HIV (MA Imputation)") +
  theme_grey()

# Plot Data Distributions
p4 <- ggplot(data=masterdata, aes(x=hiv)) +
  geom_histogram(aes(y=..density..), colour=1, fill="white") +
  geom_density(fill="#316690", color="#316690", alpha=0.4) +
  ylab("Density") + xlab("") + xlim(-50, 275) +
  ylim(0,0.015) + theme_grey()
p5 <- ggplot(data=masterdata2, aes(x=hiv)) +
  geom_histogram(aes(y=..density..), colour=1, fill="white") +
  geom_density(fill="#316690", color="#316690", alpha=0.4) +
  ylab("Density") + xlab("") + xlim(-50, 275) +
  ylim(0,0.015) + theme_grey()
p6 <- ggplot(data=masterdata3, aes(x=hiv)) +
  geom_histogram(aes(y=..density..), colour=1, fill="white") +
  geom_density(fill="#316690", color="#316690", alpha=0.4) +
  ylab("Density") + xlab("") + xlim(-50, 275) +
  ylim(0,0.015) + theme_grey()

#### PART 3. BUILDING FORECAST MODELS ####
# Use masterdata2 to call linearly interpolated data
# Use masterdata3 to call MA-imputated data
gtd <- ts(masterdata3$gtd, start=c(2012, 1), end=c(2022, 12), frequency=12)
hiv <- ts( masterdata3$hiv , start=c(2012, 1), end=c(2022, 12), frequency=12)
period <- ts(cbind(hiv,gtd))
gtd.l <- cbind(
  AdLag0 = gtd,
  AdLag1 = stats::lag(gtd,-1), 
  AdLag2 = stats::lag(gtd,-2),
  AdLag3 = stats::lag(gtd,-3)) %>% 
  head(NROW(period))

# Model 1 (SARIMA Without GTD)
fit.m1 <- auto.arima(hiv, trace = TRUE, lambda=lambda, approximation = FALSE )
fit.m1
checkresiduals(fit.m1)
fcast <- function(y, h){forecast(fit.m1, h=h)}
e <- tsCV(hiv, fcast, h=2, window=73)

RMSE <- sqrt(mean(e^2, na.rm=TRUE)); RMSE
MAE <- mean(abs(e), na.rm=TRUE); MAE
MAPE <- mean(abs(e/hiv), na.rm=TRUE)*100; MAPE

# Model 2 (SARIMA With GTD)
fit.m2 <- auto.arima(hiv, trace = TRUE, xreg=gtd, lambda=lambda, approximation = FALSE )
fit.m2; checkresiduals(fit.m2)
fcast <- function(y,h,xreg,newxreg) {
  forecast(fit.m2, xreg=newxreg, h=h)
}
e <- tsCV(hiv, fcast, h=2, window=73, xreg=gtd)  # Repeat RMSE, MAE, and MAPE codes

# Model 3 (SARIMA With GTD Lag 0,1,2,3)
fit.m3 <- auto.arima(hiv, trace = TRUE, xreg=gtd.l, lambda=lambda, approximation = FALSE )
fit.m3; checkresiduals(fit.m3)
fcast <- function(y,h,xreg,newxreg) {
  forecast(fit.m3, xreg=newxreg, h=h)
}
e <- tsCV(hiv, fcast, h=2, window=73, xreg=gtd.l) # Repeat RMSE, MAE, and MAPE codes

# Model 4 (NNAR Without GTD)
fit.m4 <- nnetar(hiv, lambda=lambda)
fit.m4;
fcast <- function(y, h){forecast(fit.m4, h=h)}
e <- tsCV(hiv, fcast, h=2, window=73)
resid <- fit.m4$residuals
Box.test(resid, type="Ljung-Box")

# Model 5 (NNAR With GTD)
fit.m5 <- nnetar(hiv, xreg=gtd, lambda=lambda)
fit.m5
fcast <- function(y,h,xreg,newxreg) {
  forecast(fit.m5, xreg=newxreg, h=h)
}
e <- tsCV(hiv, fcast, h=2, window=73, xreg=gtd)
resid <- fit.m5$residuals

# Model 6 (NNAR With GTD Lag 0,1,2,3)
fit.m6 <- nnetar(hiv, xreg=gtd.l, lambda=lambda)
fit.m6
fcast <- function(y,h,xreg,newxreg) {
  forecast(fit.m6, xreg=newxreg, h=h)
}
e <- tsCV(hiv, fcast, h=2, window=73, xreg=gtd.l)
resid <- fit.m6$residuals

#### PART 4. FORECASTS for 2023 Using GTD Data ####

# Generate Forecasts Using the GTD 2023 Data on HIV/AIDS
new.gtd <- ts(c(gtd,43,42,46,42,45,29,31,32,49,52), start=c(2012,1), end=c(2023,10), frequency=12)
new.gtd.l <- cbind(
  AdLag0 = new.gtd,
  AdLag1 = stats::lag(new.gtd,-1), 
  AdLag2 = stats::lag(new.gtd,-2),
  AdLag3 = stats::lag(new.gtd,-3))
new.gtd.l <- ts(new.gtd.l, start=c(2012,1), end=c(2023,10), frequency=12)
new.gtd.l <- tail(new.gtd.l, 10)
new.gtd <- tail(new.gtd, 10)

fcast.2023.m1 <- forecast(fit.m1, 10) 
p1 <- autoplot(fcast.2023.m1, fcol="#8E1E29") +  xlab("") + 
  ylab("") + ylim(0,250)

fcast.2023.m2 <- forecast(fit.m2, xreg=new.gtd) 
p2 <- autoplot(fcast.2023.m2, fcol="#8E1E29") +  xlab("") + 
  ylab("") + ylim(0,250)

fcast.2023.m3 <- forecast(fit.m3, xreg=new.gtd.l) 
p3 <- autoplot(fcast.2023.m3, fcol="#8E1E29") +  xlab("") + 
  ylab("") + ylim(0,250)

fcast.2023.m4 <- forecast(fit.m4, PI=TRUE,  10) 
p4 <- autoplot(fcast.2023.m4, fcol="#8E1E29") +  xlab("") + 
  ylab("") + ylim(0,250)

fcast.2023.m5 <- forecast(fit.m5, PI=TRUE,  xreg=new.gtd) 
p5 <- autoplot(fcast.2023.m5, fcol="#8E1E29") +  xlab("") + 
  ylab("") + ylim(0,250)

fcast.2023.m6 <- forecast(fit.m6, PI=TRUE, xreg=new.gtd.l)
p6 <- autoplot(fcast.2023.m6, fcol="#8E1E29") + xlab("") +
  ylab("") + ylim(0,250)


#### NOTHING FOLLOWS
