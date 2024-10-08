rm(list = ls())
gc()

# Load necessary libraries
library(rgl)
library(tsDyn)

# Define the SETAR(2,1,1) model parameters
set.seed(123)
n <- 3000

r <- rnorm(n)
x <- rep(NA, n)
x[1] <- r[1]

# Simulate the SETAR(2,1,1) model
for(t in 2:n){
    if(x[t-1] <= 0.5)
      {
        x[t] <- 4 +  0.7*x[t-1] + r[t]
      }
    else
      {
        x[t] <- -4 + 0.7*x[t-1] + r[t]
      }
  }


## Parameters for the histogram regression
## Number of intervals 
n.bin <- 50
## The breaks between the intervals 
breaks <- seq(-2,2,len=n.bin+1)
## Initialize
h <- diff(breaks)[1]
lambda <- gamma <- f.hat <- h.hat <- numeric(n.bin)
##----------------------------------------------------------------

##----------------------------------------------------------------
## Cut into intervals conditioned on x_{t-1}
L <- split(x[-1], cut(x[-length(x)],breaks))
## Check if there are at least 5 points in each interval
if(!all(sapply(L,length)>=5)){ print('Stopped: There are less than 5 points in one of the intervals'); break;}
## Calc the hist regressogram, i.e. for each interval
for(i in 1:n.bin)
  {
    x.bin <- L[[i]]
    lambda[i] <- mean(x.bin)
    f.hat[i] <- (n.bin*h)^(-1) * length(x.bin)
    gamma[i] <- sum((x.bin - lambda[i])^2) / length(x.bin)
  }
## Make confidence bands for the cumulated function. Def. (3.10).
## 95% confidence band, c is found in table 3.1
c.alpha <- 1.273
##
Lambda <- cumsum(lambda*h)
for(i in 1:n.bin)
  {
    h.hat[i] <- gamma[i]/f.hat[i];
  }
H.hat <- cumsum(h.hat*h);
##
H.hat.b <- H.hat[n.bin];
Lambda.lower <- Lambda - c.alpha * n.bin^(-0.5) * H.hat.b^(0.5) * (1 + H.hat/H.hat.b);
Lambda.upper <- Lambda + c.alpha * n.bin^(-0.5) * H.hat.b^(0.5) * (1 + H.hat/H.hat.b);
##----------------------------------------------------------------
a = -2
b = 0.5

plot(breaks[-1], Lambda, type = 'l', col = 'blue', lwd = 2, ylim = range(c(6, 9)), main = "Cumulative Conditional Means", xlab = "x", ylab = "Cumulative Mean")
lines(breaks[-1], Lambda.lower, col = 'red', lty = 2)
lines(breaks[-1], Lambda.upper, col = 'red', lty = 2)
grid()

x3 <- 0
x1 <- 0.5
x2 <- 1

før <- 4*(x3-a) + 0.35*(x3^2-a^2)
peak <- 4*(x1-a) + 0.35*(x1^2-a^2)
ned <- -4*x2 + 0.35*(x2^2) + 1.9125
efter <- peak + ned
print(før)
print(peak)
print(efter)
