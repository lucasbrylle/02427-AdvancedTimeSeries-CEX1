# Following this tutorial: https://rpubs.com/JSHAH/481706
#install.packages("tseries")
#install.packages("forecast")

library(tseries)
library(forecast)


df = read.csv("~/Documents/Uni/Advanced TSA/comp_ex_1_scrips_2018/DataPart5.csv")
df
df$x
 

# Get overview over data
summary(df)
start(df$x)
end(df$x)
time(df$x)
plot(df$x[0:50], type = "l", main = "Time Series Data", xlab = "Time", ylab = "Value")


# Perform the ADF test
adf.test(df$x)

dev.off() 
# Adjusting margins and layout to reduce white space and keep title visible
par(mfrow = c(2, 1),    # 2 rows, 1 column
    oma = c(2, 2, 2, 1), # Outer margins: bottom, left, top, right
    mar = c(4, 3, 1, 1)) # Set inner margins: bottom, left, top, right

# Plot ACF
acf(df$x, lag.max = 25, xaxt = "n", main = "", xlab="", ylab="")  # No individual title
title("ACF & PACF", outer = TRUE, line = 0.5)   # Add a common title for both plots
axis(1, at = seq(0, 25, by = 2), las = 2)
mtext("ACF", side = 2, line = 2)

# Plot PACF
pacf(df$x, lag.max = 25, xaxt = "n", main = "", xlab="Lag")
axis(1, at = seq(0, 25, by = 2), las = 2)
mtext("PACF", side = 2, line = 2)
#mtext("Lag", side = 3, line = 1)

par(mfrow = c(1, 1))
dev.off() 

# The pacf plot has two lines outside of the boundaries -> p = 2
# The acf plot has three lines outside of the boundaries -> q = 3
p = 2
q = 3

# Fit an ARMA model (adjust p and q based on ACF/PACF analysis)
arma_model <- Arima(df$x, order = c(p, 0, q))

# View model summary
summary(arma_model)

# Plot residuals
checkresiduals(arma_model)

# Alternatively, manually plot the residuals and ACF of residuals
residuals_arma <- residuals(arma_model)
#plot(residuals_arma, main = "Residuals of ARMA Model")
acf(residuals_arma, main = "ACF of Residuals, ARMA Model (p = 2, q = 3)")

# residuals are within boundaries so the model is good

# If we say p = q = 2 
p = 2
q = 2
arma_model_2_2 <- Arima(df$x, order = c(p, 0, q))
residuals_arma_2_2 <- residuals(arma_model_2_2)
acf(residuals_arma_2_2, main = "ACF of Residuals, ARMA Model (p = 2, q = 2)", xaxt="n")
axis(1, at = seq(0, 40, by = 2), las=2)



# QQ plot to check normality
qqnorm(residuals(arma_model))
qqline(residuals(arma_model), col = "red")

# Shapiro-Wilk test for normality
shapiro.test(residuals(arma_model))



# If we say p = q = 2 
p = 8
q = 8
arma_model_1_1 <- Arima(df$x, order = c(p, 0, q))
residuals_arma_1_1 <- residuals(arma_model_1_1)
acf(residuals_arma_1_1, main = "ACF of Residuals, ARMA Model (p = 2, q = 1)", xaxt="n")
axis(1, at = seq(0, 40, by = 2), las=2)



# QQ plot to check normality
qqnorm(residuals(arma_model))
qqline(residuals(arma_model), col = "red")

# Shapiro-Wilk test for normality
shapiro.test(residuals(arma_model))

# -------------- Calculate LDF


## Make a simple ldf(x,lags)
ldf <- function(x,lags,nBoot=30,plotIt=TRUE,plotFits=FALSE)
{  
  ## Calculate Lag Dependence Functions
  ## by local 2-order polynomial regression with the loess() function,
  ## and leave one out to find the best bandwidth. Finally an approximate
  ## 95% confidence interval is calculated with simple bootstrapping.
  ## Input:
  ## x, is the time series to be analysed
  ## lags, are the values for which the LDF is calculated
  ## nBoot, is the number of bootstrapping samples to make
  ## plotIt, should the ldf be plotted
  ## plotFits, should the smoothed be plotted
  
  ## The result is kept in val
  val <- vector()
  ##
  for(i in 1:length(lags))
  {
    ## Take the k
    k <- lags[i]
    ## print text
    print(paste("Calculating ldf no. ",i," of ",length(lags), sep=""))
    ## Dataframe for modelling: xk is lagged k steps
    D <- data.frame(x=x[-(1:k)],xk=x[-((length(x)-k+1):length(x))])
    ## Leave one out optimization of the bandwidth with loess
    RSSk <- leaveOneOut(D,plotFits)
    ## Calculate the ldf
    RSS <- sum((D$x - mean(D$x))^2)
    val[i] <- (RSS - RSSk) / RSS
  }      
  
  ## Very simple bootstrapping
  iidVal <- vector()
  for(i in 1:nBoot)
  {
    ## Print to entertain the modeller ;-)
    print(paste("Calculating bootstrap no. ",i," of ",nBoot, sep=""))
    ## Bootstrapping to make a confidence band
    xr <- sample(x, min(length(x),100) ,replace=TRUE)
    ## Dataframe for modelling
    DR <- data.frame(x=xr[-1],xk=xr[-length(xr)])
    RSSk <- leaveOneOut(DR)
    ## The ldf is then calculated
    RSS <- sum((DR$x - mean(DR$x))^2)
    (iidVal[i] <- (RSS - RSSk) / RSS)
  }
  
  ## Plot the ldf
  if(plotIt)
  {
    dev.new()
    plot(c(0,lags), c(1,val), type="n", ylim=c(-1,1), ylab="LDF", main="Lag Dependence Functions", xaxt="n", xlab="lag")
    axis(1,c(0,lags))
    abline(0,0,col="gray")
    lines(c(0,lags), c(1,val), type="h")
    ## Draw the approximate 95% confidence interval
    abline(h=quantile(iidVal,0.95), col="blue", lty=2)
  }
  return(val)
}


leaveOneOut <- function(D,plotIt=FALSE)
{
  ## Find the bandwidth giving the best balance between bias and variance
  ## of the model by leave one out.
  
  ## Plot the data
  if(plotIt)
  {
    dev.new()
    par(mfrow=c(1,2))
    plot(D$xk, D$x)
  }
  ## Make the vector of bandwidths which are fitted
  span <- c(seq(0.2, 1, by=0.1),2,4,10)
  ## Matrix for keeping the residuals
  R <- matrix(nrow=nrow(D), ncol=length(span))
  ## Do it for each bandwidth
  for(ii in 1:length(span))
  {
    print(paste("  Fitting for bandwidth",ii,"of",length(span)))
    ## Do the local 2-order polynomial regression one time for each point
    ## leaving the point out while fitting, and then predicting the point.
    for(i in 1:nrow(D))
    {
      R[i,ii] <- D[i,"x"] - predict(loess(x ~ xk, dat=D[-i,], span=span[ii]), D[i,])
    }
  }
  ## Find the best bandwidth
  RSSkAll <- apply(R, 2, function(x){sum(x^2,na.rm=TRUE)})
  ## Take the RRS for the best span value
  spanBest <- span[which.min(RSSkAll)]
  ## Calculate the RSSk
  RSSk <- sum((D$x - predict(loess(x ~ xk, dat=D, span=spanBest), D))^2,na.rm=TRUE)
  ## Plot the fitted function
  if(plotIt)
  {
    DT <- D[order(D$xk),]
    lines(DT$xk, predict(loess(x ~ xk, dat=D, span=spanBest), DT), col=2)
    ## Plot RSSk as a function of the bandwidth
    plot(span, RSSkAll)
  }
  ##
  return(RSSk)
}


ldf(residuals_arma_2_2,1:10)


# Fit a GARCH Model
#install.packages("rugarch")
library(rugarch)

# Specify GARCH model with ARMA terms
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(2, 1), include.mean = TRUE))

# Fit the model
garch_model <- ugarchfit(spec = spec, data = df$x)

# Print model summary
summary(garch_model)

residuals_garch <- residuals(garch_model)
ldf(residuals_garch,1:10)
residuals_garch


# Fit a SETAR model
#install.packages("tsDyn")
library(tsDyn)

ts_data <- ts(df$x)

setar_model <- setar(ts_data, mL = 2, mH = 2, thDelay = 1)

summary(setar_model)

# Plot the fitted model and regime switching
plot(setar_model)


residuals_setar <- residuals(setar_model)
ldf(residuals_setar,1:10)
