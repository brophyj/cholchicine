# brms model
pacman::p_load(brms, tidyverse, tidybayes)
data_bin2 <- data.frame(N = c(3535,3528), y = c(329,321), grp2 = as.factor(c("placebo","colchicine"))) 
f = bf(y | trials(N) ~ 0 + grp2)

#get_prior(formula = f,data = data_bin2,family = binomial(link = "identity"))

m <- brm(
  formula = f,
  data = data_bin2,
  family = binomial(link = "identity"),
  prior = c(prior(beta(1, 1), class = b, lb = 0, ub = 1)), # gets rid os a bunch of unhelpful warnings
  chains = 4, warmup = 1000, iter = 2000, seed = 123,
  refresh = 0
)
summary(m)

draws_identity <- m %>% 
  epred_draws(newdata = data_bin2)  # from tidybayes, gives counts

# brms::as_draws_df makes a tibble with proportions
# tidybayes::tidy_draws gives same results

# can also plot model results
# colour scheme
clrs_saguaro <- NatParksPalettes::natparks.pals("Saguaro")
clr_grp20 <- clrs_saguaro[1]
clr_grp21 <- clrs_saguaro[6]
clr_diff <- clrs_saguaro[4]

p1 <- draws_identity %>% 
  ggplot(aes(x = .epred, y = grp2, fill = grp2)) +
  stat_halfeye() +
  scale_fill_manual(values = c(clr_grp20, clr_grp21)) +
  guides(fill = "none") +
  labs(x = "Number of primary outcomes",
       y = NULL) +
  theme_classic() +
  labs(title = "OASIS-9 trial results", subtitle = "Vague non-informative prior")

p2 <- draws_identity %>% 
  ungroup() %>% 
  mutate(grp2 = fct_relevel(grp2, "placebo")) %>% 
  mutate(.epred = .epred /N) %>% 
  compare_levels(.epred, by = "grp2") %>%  # compare_levels() subtracts things using alphabetical order
  ggplot(aes(x = .epred)) +
  stat_halfeye(fill = clr_diff) +
  guides(fill = "none") +
  labs(x = "Difference of outcomes (colchicine - placebo)",
       y = NULL) +
  theme_classic() +
  labs(title = "OASIS-9 trial results", subtitle = "Vague non-informative prior")
