
# informative prior model

alpha1 = as.integer(2379 * 0.071)
beta1 <- 2379 - alpha1
alpha2 = as.integer(2366 * 0.055)
beta2 <- 2366 - alpha2
# Calculate means for placebo and experimental arms
mean_placebo <- alpha1 / (alpha1 + beta1)
mean_experimental <- alpha2 / (alpha2 + beta2)

# Display the results
mean_placebo
mean_experimental

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


var_placebo <- (alpha1 * beta1) / ((alpha1 + beta1)^2 * (alpha1 + beta1 + 1))  # Beta variance
sd_placebo <- sqrt(var_placebo)  # Beta standard deviation

var_experimental <- (alpha2 * beta2) / ((alpha2 + beta2)^2 * (alpha2 + beta2 + 1))  # Beta variance
sd_experimental <- sqrt(var_experimental)  # Beta standard deviation

mean_rr <- (alpha2 / (alpha2 + beta2)) / (alpha1 / (alpha1 + beta1))

var_log_rr <- (var_experimental / (mean_experimental^2)) + (var_placebo / (mean_placebo^2)) #delta method
sd_log_rr <- sqrt(var_log_rr)


prior_clear <- set_prior(paste0("normal(", log(mean_rr), ", ", sd_log_rr, ")"), class = "Intercept")

dat <- metafor::escalc(measure="RR", ai=e_col, n1i = n_col, ci=e_ctl,  n2i =n_ctl, data=col_df)



clear_model <- brm(
  yi | se(sqrt(vi)) ~ 1,
  data = filter(dat, study == "CLEAR"),  # Only CLEAR data
  prior = prior_clear,
  iter = 20000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99),
  seed = 123
)

draws_clear <- as_draws_df(clear_model)
predicted_clear <- rnorm(
  n = nrow(draws_clear),
  mean = draws_clear$b_Intercept,
  sd = sqrt(draws_clear$sigma^2)
)
predicted_rr_clear <- exp(predicted_clear)
prob_clear_rr_benefit <- mean(predicted_rr_clear <= 0.9)
prob_clear_rr_benefit
