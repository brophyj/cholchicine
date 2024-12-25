# stan model for informative analysis


# file.exists("binom_2_priorCOLCOT.stan")


# prior COLCOT data

alpha1 = as.integer(2379 * 0.071)
beta1 <- 2379 - alpha1
alpha2 = as.integer(2366 * 0.055)
beta2 <- 2366 - alpha2

# Data for the model
data_list <- list(
  n1 = 3534,
  y1 = as.integer(3534 * 0.093),
  n2 = 3528,
  y2 = as.integer(3528 * 0.091),
  alpha1 = alpha1,
  beta1 = beta1,
  alpha2 = alpha2,
  beta2 = beta2
)


# Compile and fit the model
mod2 <- cmdstan_model("binom_2_priorCOLCOT.stan")
fit2 <- mod2$sample(data = data_list, chains = 4, parallel_chains = 4,        
                  refresh = 0, seed = 123)

# Extracting posterior samples
posterior_samples2 <- fit2$draws()

# Plotting
color_scheme_set("blue")
mcmc_trace(posterior_samples2, pars = c("p1", "p2", "rr"), nrow = 3)

# Plot relative risk distributions with no formatting
# mcmc_areas(posterior_samples, pars = "rr", prob = 0.95)

# print summary
fit2$summary()

# Correct extraction of relative risk samples
rr_samples2 <- fit2$draws(variables = "rr")
rr_vector2 <- as.vector(rr_samples2)  # Convert to a simple vector for easier handling

# Calculating probabilities
prob_rr_less_09 <- mean(rr_vector2 < 0.9)
prob_rr_between_09_11 <- mean(rr_vector2 >= 0.9 & rr_vector2 <= 1.1)
prob_rr_greater_11 <- mean(rr_vector2 > 1.1)

# Print the probabilities
cat("Probability RR < 0.9: ", prob_rr_less_09, "\nProbability RR 0.9 to 1.1: ", prob_rr_between_09_11, "\nProbability RR > 1.1: ", prob_rr_greater_11, "\n")

