
# hierarchical meta-analytical data prep and plot using brms package
# slightly more informative priors N(0,1) than the default priors
library(tidybayes); library(brms); library(posterior); library(tidyverse)

# data prep
total <- matrix(c(322, (3528-322), 327, (3534-327)), nrow = 2, byrow = TRUE, dimnames = list(c("colchicine", "placebo"), c("Event", "No event")))
epiR::epi.2by2(dat = as.table(total), method = "cohort.count", conf.level = 0.95, units = 100, outcome = "as.columns")

col_df <- data.frame(study=c("CLEAR", "COLCOT"), n_col=c(3528,2366), e_col=c(322,130), n_ctl=c(3534,2379), e_ctl=c(327,169))
dat <- metafor::escalc(measure="RR", ai=e_col, n1i = n_col, ci=e_ctl,  n2i =n_ctl, data=col_df)

# standard (non-bayesian) random effects model
library(metafor)
me.fe <- rma(dat$yi, sei=sqrt(dat$vi), method = "FE") # exponentiation c(exp(me.fe$b), exp(me.fe$ci.lb), exp(me.fe$ci.ub))
me.re <- rma(dat$yi, sei=sqrt(dat$vi), method = "REML") # exponentiation c(exp(me.re$b), exp(me.re$ci.lb), exp(me.re$ci.ub))
plot_Plato <- meta::metabin(col_df$e_col, col_df$n_col, col_df$e_ctl, col_df$n_ctl, sm="RR", method ="I", studlab=dat$Region, prediction=TRUE, comb.random =TRUE)
print(summary(plot_Plato,prediction=TRUE), digits=2)

meta::forest(plot_Plato, studlab = col_df$study)

# Time to go bayesian
# Specify priors for mu and tau (could also ignore and use brms defaults,  set_prior("normal(0,10)", class = "sd"))
prior <- c(
  set_prior("normal(0, 1)", class = "Intercept"),      # Matches Stan's prior for mu
  set_prior("normal(0, 1)", class = "sd")             # Matches Stan's prior for tau
)

# Fit the model with custom priors
brm_out <- brm(
  yi | se(sqrt(vi)) ~ 1 + (1|study), 
  data = dat, 
  iter = 20000, warmup = 2000, cores = 4, refresh = 0,
  control = list(adapt_delta = 0.99), # Improve convergence
  prior = prior,                     # Apply the custom priors
  seed = 123
)

summary(brm_out)


# Extract posterior samples as_draws() is another option
post <- brm_out %>%
  spread_draws(b_Intercept, r_study[study,]) %>%  # variables(brm_out)
  median_qi(condition_mean = b_Intercept + r_study, .width = c(.95)) %>% 
  rename(lower = .lower, rr = condition_mean, upper = .upper) %>% 
  select(study, rr, lower, upper) # Keep relevant columns for study-level output

post[,c(2:4)] <- apply(post[,c(2:4)],2,exp) # change from log(OR) to OR

# Extract posterior draws with posterior package and as_draws_df
draws <- as_draws_df(brm_out)  # variables(draws) to get names

# Mean and 95% CI for the mean effect
mean_intercept <- mean(draws$b_Intercept) # Mean effect (Intercept)
ci_intercept <- quantile(draws$b_Intercept, probs = c(0.025, 0.975)) # Between-study standard deviation
sigma <- draws$sigma                       # Residual standard deviation

# Transform to OR scale
mean_or <- exp(mean_intercept)
ci_or <- exp(ci_intercept)

# Simulate predicted values for the next study (log scale)
predicted_values <- rnorm(n = length(draws$b_Intercept), mean = mean_intercept,
                          sd = sqrt(draws$sd_study__Intercept^2 + sigma^2)
)
predicted_ci <- quantile(predicted_values, probs = c(0.025, 0.975))

# Transform predicted values to OR scale
predicted_mean_or <- exp(mean(predicted_values))
predicted_ci_or <- exp(predicted_ci)


#########

dat <- dat %>% 
  mutate(rr=yi, upper= yi +1.96*sqrt(vi), lower= yi - 1.96*sqrt(vi)) %>% 
  select(study,rr,lower,upper)
dat[,c(2:4)] <- apply(dat[,c(2:4)],2,exp)


post <- rbind(post, dat)

post$lab <- rep(c("Theta", "Y"), each = 2)
post$id <- c(1,2,1,2)

# Create the data frame with desired structure for hierarchical plotting
# Add overall mean and next study
results_df <- tibble(
  study = c("Mean", "Predicted Next Study"),
  lower = c(ci_or[1], predicted_ci_or[1]),
  rr = c(mean_or, predicted_mean_or),
  upper = c(ci_or[2], predicted_ci_or[2])
)
results_df$lab <- c("Mean", "Next")
results_df$id <- c(3,4)

post <- rbind(post, results_df)

ggplot(post, aes(x = forcats::fct_rev(study), y = rr, ymin = lower, ymax = upper, col = lab)) +  
  geom_pointrange(aes(col = lab), position = position_dodge(width = 0.50)) +
  coord_flip() + geom_hline(aes(yintercept = 0.895), lty = 2) +  xlab("") + 
  ylab("")  + theme(legend.position="bottom") + geom_hline(aes(yintercept = 1), lty = 1) +
  scale_colour_discrete(name="", 
                        labels = c("Theta" = bquote("Random effect \n(hierarchical \"shrinking\"):"~exp(theta[J])~" "),
                                   "Y"= bquote("Relative risk \n(observed data):"~exp(Y[J])))) +
  labs(title = "Bayesian forest plot of cholcicine trials",
       subtitle = "Observed and hierarchical individual trial results",
       caption = "Prior tau = normal(0, 1.0)
\nMean prior = normal(0, 1.0) |") +
  theme_bw()

ggsave("output/brms_hier.png", dpi = 600, width = 8)

########## Reproduce same table as from Stan
# Summarize study-level estimates

# Extract and summarize overall mu
mu_values <- brm_out %>% spread_draws(b_Intercept) %>% pull(b_Intercept)
overall_mu <- tibble(study = "Overall (mu)", rr = mean(mu_values, na.rm = TRUE),
                     lower = quantile(mu_values, 0.025, na.rm = TRUE),
                     upper = quantile(mu_values, 0.975, na.rm = TRUE))

# Extract and summarize overall tau
tau_values <- brm_out %>% spread_draws(sd_study__Intercept) %>% pull(sd_study__Intercept)
overall_tau <- tibble(study = "Overall (tau)", rr = mean(tau_values, na.rm = TRUE),
                      lower = quantile(tau_values, 0.025, na.rm = TRUE),
                      upper = quantile(tau_values, 0.975, na.rm = TRUE))

# Combine all summaries
temp <- post[c(1:2),c(1:4)]
temp[,c(2:4)] <- apply(temp[,c(2:4)],2,log) 
overall_post <- bind_rows(temp, overall_mu, overall_tau)
overall_post

# probability of > 10% reduction, log(.9) = -0.11
# mean sd = (.653 + .860)/3.92 = 0.39
pnorm(-.11, -.111, .39)
