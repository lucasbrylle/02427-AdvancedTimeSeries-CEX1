rm(list = ls())
gc()

# Load necessary libraries
library(rgl)
library(tsDyn)

data = read.csv("DataPart5.csv")
# Fit ARMA model
arma_fit <- arima(data$x, order = c(1, 0, 1))

# Plot the data and ARMA model
plot(data$x, type = 'l', main = "ARMA Model Fit", ylab = "Values", xlab = "Time")
