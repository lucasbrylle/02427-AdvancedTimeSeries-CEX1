rm(list = ls())
gc()

# Define the transition matrix for the Markov chain
transition_matrix <- matrix(c(0.6, 0.4,
                              0.05, 0.95), 
                            nrow = 2, byrow = TRUE)

# Function to simulate the Markov chain
simulate_markov_chain <- function(n, transition_matrix) {
  states <- 1:nrow(transition_matrix)
  Jt <- numeric(n)
  Jt[1] <- sample(states, 1)
  for (t in 2:n) {
    Jt[t] <- sample(states, 1, prob = transition_matrix[Jt[t-1], ])
  }
  return(Jt)
}

# Simulate the Markov chain
n <- 1000 
r <- rnorm(n)
Jt <- simulate_markov_chain(n, transition_matrix)

# Simulate the MMAR model
simulate_mmar <- function(n, Jt) {
  y <- numeric(n)
  for (t in 2:n) {
    if (Jt[t] == 1) {
      y[t] <- 4 + 0.7 * y[t-1] + r[t]
    } else {
      y[t] <- -100 + 0.7 * y[t-1] + r[t]
    }
  }
  return(y)
}

# Simulate the MMAR model
y <- simulate_mmar(n, Jt)

# Plot the simulated MMAR model
J_test = Jt[1:100]
plot(y, main="MMAR Model Simulation", ylab="Values", xlab="Time", )