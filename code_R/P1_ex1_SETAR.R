rm(list = ls())
gc()
set.seed(123)

## Use rgl
library(rgl)

## Number of samplepoints
n <- 1000
## Uniform distributed x
x <- runif(n,-1,1)
## Errors
r <- rnorm(n)

## Make a time series y with a regime model
y <- rep(NA,n)
y[1] <- r[1]

for(t in 2:n){
    if(y[t-1] <= 0)
      {
        y[t] <- 10 +  0.25*y[t-1] + r[t]
      }
    else
      {
        y[t] <- -10 + 0.25*y[t-1] + r[t]
      }
  }

## Put it into a data.frame, and make x1 and y1 which are lagged one step 
D <- data.frame(y=y[-1], x1=x[-n], y1=y[-n])


# Plot it against the original data
plot(D$y, main="SETAR(2;1;1) Model Forecast", ylab="Forecasted Values", xlab="Time")

