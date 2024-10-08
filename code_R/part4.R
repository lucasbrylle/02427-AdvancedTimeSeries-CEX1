library(ggplot2)

df = read.csv("~/Documents/Uni/Advanced TSA/comp_ex_1_scrips_2018/DataPart4.csv")
df
summary(df2)

df 
df <- subset(df, Ph >= 0)

W_t = df$W

U_a = df$Ph / (df$Ti- df$Te)



plot(U_a  ~ W_t, xlab="Ua", ylab="Wt", main="Plotting Ua against Wt", ylim=c(130, 250))

plot(W_t ~ U_a, xlab="Ua", ylab="Wt", main="Plotting Ua against Wt", xlim=c(100, 300))

plot(U_a, xlab="t", type="l", ylim=c(0, 450))
lines(W_t, col="red")


# ------------ PLOT DATAFRAME VARIABLES
par(mar = c(5, 4, 4, 4) + 0.1)

plot(df$Ti, xlab = "t", ylab = "Temperature", main="Plotting Dataframe variables: Ti, Te and W", type = "l", col = "red", ylim=c(5, 25))
lines(df$Te, col = "black")
par(new = TRUE)  # Allow a new plot on top of the existing one
plot(df$W, axes = FALSE, ylab = "", xlab = "", type = "l", col = "blue", ylim=c(0, 15))

# Add a y-axis on the right side for W_t
axis(4, at = pretty(range(df$W)), col.axis = "blue", col = "blue")  # Adds a new axis on the right, using 'col.axis' to color tick labels
mtext("Wind speed", side = 4, line = 2.5, col = "blue")  # Label the right y-axis with the desired label and color

# Optionally add a legend to distinguish the lines
legend("topright", legend = c("Ti", "Te", "W"), col = c("red", "black", "blue"), cex = 0.8, lty = 1)

# -------------- PLOT U_a and W_t
par(mar = c(5, 4, 4, 4) + 0.1)

# Plot U_a with its scale on the left y-axis
plot(U_a, xlab = "t", ylab = "Ua(Wt)", main="Plotting Ua(Wt) and Wt", type = "l", col = "black", ylim = c(0, 450))

# Overlay W_t on the same plot using a different scale (0 to 12) on the right y-axis
par(new = TRUE)  # Allow a new plot on top of the existing one

# Plot W_t without the x and y labels, and use a new y-axis range
plot(W_t, type = "l", col = "red", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 13))

# Add a y-axis on the right side for W_t
axis(4, at = pretty(range(W_t)), col="red", col.axis = "red")  # Adds a new axis on the right
mtext("Wt", side = 4, line = 2.5, col="red")  # Label the right y-axis

# Optionally add a legend to distinguish the lines
legend("topright", legend = c("Ua", "Wt"), col = c("black", "red"), cex=0.8, lty = 1)

# --------- FIT SIMPLE MODEL ----------------
summary(W_t)

fit <- loess(U_a ~ W_t)

plot(W_t, U_a, pch = 16, col = "blue", main = "Ua(Wt) ~ Wind Speed", xlab = "Wind Speed (Wt)", ylim=c(50, 280), ylab = "Ua(Wt)")
lines(W_t, predict(fit), col = "red")

# FIT WITH DIFFERENT BANDWIDTHS
bandwidths <- c(0.5, 0.8)

# Define colors for each bandwidth
colors <- rainbow(length(bandwidths) + 1)  # Generate a color palette

plot(W_t, U_a, pch = 16, col = "blue", main = "Ua(Wt) vs Wind Speed", xlab = "Wind Speed (Wt)", ylim=c(50, 280), ylab = "Ua(Wt)")

# Loop through bandwidths and add the fitted curves
for (i in seq_along(bandwidths)) {
  bw <- bandwidths[i]
  fit <- loess(U_a ~ W_t, span = bw)
  lines(W_t, predict(fit), col = colors[i], lwd = 2)
  
}
fit <- loess(U_a ~ W_t, degree=1)
lines(W_t, predict(fit), col = colors[-1], lwd = 2)



# Add a legend to differentiate between bandwidths
legend("topright", legend = paste("Bandwidth:", bandwidths), col = colors, lwd = 2, pch=c(1,2,3))



#install.packages("tidyr")
library(ggplot2)
library(tidyr)

data$U_a <- data$Ph / (data$Ti - data$Te)

# Ensure there are no NaNs or infinities (e.g., due to division by zero)
data <- data[is.finite(data$U_a), ]

loess_fit_1 <- loess(U_a ~ W_t , data = data, span = 0.2)
loess_fit_2 <- loess(U_a ~ W_t, data = data, span = 0.5)
#loess_fit_1 <- loess(U_a ~ W_t, data = data, span = 0.8, parametric="W_t")
#loess_fit_2 <- loess(U_a ~ W_t + data$Ti + data$Te, data = data, span = 0.8, parametric=c("data$Ti", "data$Te"))
loess_fit_3 <- loess(U_a ~ W_t, data = data, span = 0.8)

# Create a sequence of W_t values for prediction
W_t_seq <- seq(min(data$W_t), max(data$W_t), length.out = 100)

# Predict U_a for each model
pred_1 <- predict(loess_fit_1, newdata = data.frame(W_t = W_t_seq))
pred_2 <- predict(loess_fit_2, newdata = data.frame(W_t = W_t_seq))
pred_3 <- predict(loess_fit_3, newdata = data.frame(W_t = W_t_seq))

# Prepare the data for plotting
plot_data <- data.frame(W_t = W_t_seq, 
                        Span_0.2 = pred_1, 
                        Span_0.5 = pred_2, 
                        Span_0.8 = pred_3)

# Reshape the data to a long format
plot_data_long <- gather(plot_data, key = "Span", value = "U_a", -W_t)

# Create the plot
ggplot(data, aes(x = W_t, y = U_a)) +
  geom_point(alpha = 0.5) +
  geom_line(data = plot_data_long, aes(x = W_t, y = U_a, color = Span), size = 1) +
  labs(title = "Estimated U_a as a Function of W_t",
       x = "W_t",
       y = "U_a",
       color = "Bandwidth") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.minor = element_blank() # Remove the grid
        
        ) +
  scale_color_manual(values = c("Span_0.2" = "blue", "Span_0.5" = "red", "Span_0.8" = "green"),
                     labels = c("Span = 0.2", "Span = 0.5", "Span = 0.8")) +
  ylim(125, 250)


