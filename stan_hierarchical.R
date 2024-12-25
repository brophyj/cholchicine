library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(cmdstanr, quietly = T)
library(bayesplot,quietly = T)
library(lme4, quietly = T)
library(posterior, quietly = T)
library(brms, quietly = T)

# stan hierarchical model
data <- list(
  n_studies = 2,
  log_or = c(
    log((322 / (3528 - 322)) / (327 / (3534 - 327))),
    log((130 / (2366 - 130)) / (168 / (2379 - 168)))
  ),
  log_or_se = c(
    sqrt(1 / 322 + 1 / (3528 - 322) + 1 / 327 + 1 / (3534 - 327)),
    sqrt(1 / 130 + 1 / (2366 - 130) + 1 / 168 + 1 / (2379 - 168))
  )
)

stan_code <- "
data {
  int<lower=1> n_studies;
  vector[n_studies] log_or;            // Observed log odds ratios
  vector<lower=0>[n_studies] log_or_se; // Standard errors of log odds ratios
}
parameters {
  real mu;                             // Overall log odds ratio
  real<lower=0> tau;                   // Heterogeneity standard deviation
  vector[n_studies] theta;             // Study-specific log odds ratios
}
model {
  mu ~ normal(0, 1);                   // Prior for overall log odds ratio
  tau ~ normal(0, 1);                  // Prior for heterogeneity
  theta ~ normal(mu, tau);             // Hierarchical model
  log_or ~ normal(theta, log_or_se);   // Observed log odds ratios
}
"

model_file <- tempfile(fileext = ".stan")
writeLines(stan_code, model_file)

mod <- cmdstan_model(model_file)

fit <- mod$sample(
  data = data,
  chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000,
  seed = 123
)

fit95 <- fit$summary(NULL, ~quantile(.x, probs = c(0.025, .5, 0.975))) %>% 
  filter(variable %in% c("mu", "theta[1]", "theta[2]")) %>%
  rename(study = variable, RR = `50%`, lower_CI = `2.5%`, upper_CI = `97.5%` ) %>% 
  mutate(lab =c("Mean", "Theta", "Theta"), study = c("Mean", "CLEAR", "COLCOT")) %>% 
  mutate_at(c(2:4), ~ exp(.)) %>%  
  relocate(RR, .after=study)

# useful function posterior <- as_draws_df(fit$draws()) or
# fit$draws(variables = c("theta", "mu", "tau"), format = "df")

# Compute observed OR and 95% CI
observed_data <- data.frame(
  study = c("CLEAR", "COLCOT"),
  RR = exp(data$log_or),
  lower_CI = exp(data$log_or - 1.96 * data$log_or_se),
  upper_CI = exp(data$log_or + 1.96 * data$log_or_se),
  lab = c("Y", "Y")
)

# Extract posterior summaries for shrunk ORs (theta) and overall OR (mu)
# Compute prediction interval
prediction_samples <- fit$draws(variables = c("theta", "mu", "tau"), format = "df") %>%
  mutate(pred_next_study = rnorm(n(), mean = mean(mu), sd = mean(tau)))

prediction_summary <- prediction_samples %>%
  summarise_draws() %>%
  filter(variable == "pred_next_study") %>%
  dplyr::select(c(1,2,4)) %>% 
  mutate(
    RR = exp(mean),
    lower_CI = exp(mean-1.96*sd),
    upper_CI = exp(mean+1.96*sd),
    lab =c("Predicted next study"),
    variable = c("Predicted Next Study")) %>% 
  rename(study = variable) %>% 
  dplyr::select(-c(2,3))

# Combine observed and shrunk data for studies

p.stan <- rbind(fit95, observed_data, prediction_summary )


ggplot(p.stan, aes(x = forcats::fct_rev(study), y = RR, ymin = lower_CI, ymax = upper_CI, col = lab)) +  
  geom_pointrange(aes(col = lab), position = position_dodge(width = 0.50)) +
  coord_flip() + geom_hline(aes(yintercept = 0.887), lty = 2) +  xlab("") + 
  ylab("")  + theme(legend.position="bottom") + geom_hline(aes(yintercept = 1), lty = 1) +
  scale_colour_discrete(name="", 
                        labels = c("Theta" = bquote("Random effect \n(hierarchical \"shrinking\"):"~exp(theta[J])~" "),
                                   "Y"= bquote("Relative risk \n(observed data):"~exp(Y[J])))) +
  labs(title = "Bayesian forest plot of cholcicine trials",
       subtitle = "Observed and hierarchical individual trial results",
       caption = "Prior tau = student_t(3, 0, 0.5)
\nMean prior = normal(0, 1.0) |") +
  theme_bw()

p.stan

# probability of > 10% reduction, RR = 0.9
# mean sd = (1.91-0.404)/3.92 = 0.38
pnorm(0.9, 0.887, .38)

