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
for (t in 2:n) {
  if (x[t-1] < 0) {
    y[t] <- 4 + 0.5 * y[t-1] + r[t]
  } else {
    y[t] <- -4 - 0.5 * y[t-1] + r[t]
  }
}

# Create a data frame
D <- data.frame(y = y[-1], x1 = x[-n], y1 = y[-n])

# Fit a local regression model using loess
loess_fit <- loess(y ~ x1, data = D, span = 0.1)  # Adjust span for different bandwidths

# Predict using the loess model
x_grid <- seq(min(D$x1), max(D$x1), length.out = 100)
y_pred <- predict(loess_fit, newdata = data.frame(x1 = x_grid))

# Plot the results
plot(D$x1, D$y, main = "Local Regression Estimate of M(x)", xlab = "x", ylab = "E{X_{t+1} | X_t = x}", pch = 16, col = rgb(0, 0, 0, 0.5))
lines(x_grid, y_pred, col = "blue", lwd = 2)

# Try different bandwidths and compare
loess_fit_1 <- loess(y ~ x1, data = D, span = 0.5)
y_pred_1 <- predict(loess_fit_1, newdata = data.frame(x1 = x_grid))
lines(x_grid, y_pred_1, col = "red", lwd = 2, lty = 2)

loess_fit_2 <- loess(y ~ x1, data = D, span = 0.7)
y_pred_2 <- predict(loess_fit_2, newdata = data.frame(x1 = x_grid))
lines(x_grid, y_pred_2, col = "green", lwd = 2, lty = 3)

legend("topright", legend = c("span=0.1", "span=0.5", "span=0.7"), col = c("blue", "red", "green"), lty = c(1, 2, 3), lwd = 2)