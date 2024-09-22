rm(list = ls())
gc()

# Load necessary libraries
library(rgl)
library(tsDyn)

# Define the SETAR(2,1,1) model parameters
set.seed(123)
n <- 1000
x <- runif(n, -1, 1)
r <- rnorm(n)
y <- rep(NA, n)
y[1] <- r[1]

# Simulate the SETAR(2,1,1) model
for(t in 2:n){
    if(y[t-1] <= 200)
      {
        y[t] <- 0.2 +  y[t-1] + r[t]
      }
    else
      {
        y[t] <- 10 + 0.95*y[t-1] + r[t]
      }
  }

# Create a data frame
D <- data.frame(y = y[-1], x1 = x[-n], y1 = y[-n])

# Calculate cumulative conditional means
n.bin <- 20
breaks <- seq(min(D$x1), max(D$x1), length.out = n.bin + 1)
h <- diff(breaks)[1]
lambda <- gamma <- f.hat <- h.hat <- numeric(n.bin)

L <- split(D$y, cut(D$x1, breaks))
if (!all(sapply(L, length) >= 5)) {
  stop('Stopped: There are less than 5 points in one of the intervals')
}

for (i in 1:n.bin) {
  x.bin <- L[[i]]
  lambda[i] <- mean(x.bin)
  f.hat[i] <- (n.bin * h)^(-1) * length(x.bin)
  gamma[i] <- sum((x.bin - lambda[i])^2) / length(x.bin)
}

c.alpha <- 1.273
Lambda <- cumsum(lambda * h)
for (i in 1:n.bin) {
  h.hat[i] <- gamma[i] / f.hat[i]
}
H.hat <- cumsum(h.hat * h)
H.hat.b <- H.hat[n.bin]
Lambda.lower <- Lambda - c.alpha * n.bin^(-0.5) * H.hat.b^(0.5) * (1 + H.hat / H.hat.b)
Lambda.upper <- Lambda + c.alpha * n.bin^(-0.5) * H.hat.b^(0.5) * (1 + H.hat / H.hat.b)


# Plot the cumulative conditional means
#plot((x))
plot(breaks[-1], Lambda, type = 'l', col = 'blue', lwd = 2, ylim = range(c(Lambda.lower, Lambda.upper)), main = "Cumulative Conditional Means", xlab = "x", ylab = "Cumulative Mean")
lines(breaks[-1], Lambda.lower, col = 'red', lty = 2)
lines(breaks[-1], Lambda.upper, col = 'red', lty = 2)

print(breaks)