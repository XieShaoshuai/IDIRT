library(ggplot2)

# Prior distribution parameters
alpha_prior_A <- 2
beta_prior_A <- 2
alpha_prior_B <- 2
beta_prior_B <- 2

# Observed data
# Replace these with your actual data
data_A <- peptide$ratio
data_B <- rbinom(100, 1, 0.5)

# Number of trials and number of successes
n_trials_A <- length(data_A)
successes_A <- sum(data_A)
n_trials_B <- length(data_B)
successes_B <- sum(data_B)

# Posterior distribution parameters
alpha_posterior_A <- alpha_prior_A + successes_A
beta_posterior_A <- beta_prior_A + n_trials_A - successes_A
alpha_posterior_B <- alpha_prior_B + successes_B
beta_posterior_B <- beta_prior_B + n_trials_B - successes_B

# Generate samples from the posterior distributions
samples_A <- rbeta(1000, shape1 = alpha_posterior_A, shape2 = beta_posterior_A)
samples_B <- rbeta(1000, shape1 = alpha_posterior_B, shape2 = beta_posterior_B)

x <- seq(0, 1, length.out = 1000)
df <- data.frame(x,samples_A, samples_B)


p <- ggplot(df) +
  geom_density(aes(x = samples_A, colour = 'Posterior A')) +
  geom_density(aes(x = samples_B, colour = 'Posterior B')) +
  labs(x = 'Value', y = 'Density', colour = 'Distribution') +
  theme_minimal()
p


library(bayestestR)
library(rstanarm)

bf <- ttestBF(x = df$samples_B, y = df$samples_B)
d<-p_direction(bf)

pd_to_p(d$pd)



# Calculate the probability that a random sample from B's posterior distribution is greater than a random sample from A's posterior distribution
prob_B_greater_than_A <- mean(samples_A > samples_B)/10000

print(paste('Probability that a random sample from B is greater than a random sample from A:', prob_B_greater_than_A))
