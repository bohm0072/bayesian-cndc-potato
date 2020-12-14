library(tidyverse)

d <- expand_grid(variety = c("A", "B"), biomass = 1:100, draw = 1:5)

d <- d %>% 
  mutate(mu = case_when(
    variety == "A" ~ 4,
    TRUE ~ 2
  )) %>% 
  mutate(beta = case_when(
    variety == "A" ~ -0.1,
    TRUE ~ -0.05
  )) %>% 
  rowwise() %>% 
  mutate(obs = rnorm(1, mean = mu+beta*biomass))



# just plotting the curves
d %>% 
  ggplot(aes(x = biomass, y = obs, color = variety, group = interaction(variety, draw))) +
  geom_line()

# plotting the differences
d %>% 
  select(-mu, -beta) %>% 
  pivot_wider(names_from = variety, values_from = obs) %>% 
  mutate(diff = A - B) %>% 
  ggplot(aes(x = biomass, y = diff, group = draw)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red")
