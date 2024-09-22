rm(list=ls())
## Use rgl
library(rgl)

## Number of sample points
n <- 100
## Uniform distributed x
x <- runif(n, -1, 1)
## Errors
r <- rnorm(n)

## Make a time series y with a regime model
y <- rep(NA,n)
y[1] <- r[1]

for(t in 2:n)
  {
    if(x[t-1] < 0)
      {
        y[t] <- 4 + 0.5 * y[t-1] + r[t]
      }
    else
      {
        y[t] <- -4 - 0.5 * y[t-1] + r[t]
      }
  }


## Put it into a data.frame, and make x1 and y1 which are lagged one step 
D <- data.frame(y=y[-1], x1=x[-n], y1=y[-n])

##------------------------------------------------
params_regime1 <- c(1, 0.5)  # Parameters for regime 1
params_regime2 <- c(-1, 0.5)  # Parameters for regime 2

# Simulate the IGAR model
for(t in 2:n) {
  if(runif(1) < 0.5) {  # Stochastic regime switch with equal probability
    y[t] <- params_regime1[1] * y[t-1] + r[t]
  } else {
    y[t] <- params_regime2[1] * y[t-1] + r[t]
  }
}

# Plot the simulated IGAR model
plot.ts(y, main="IGAR Model Simulation", ylab="Values", xlab="Time")


# Hvad er d = 1? Er det laggen i AR term?