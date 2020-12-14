library(tidyverse)

d <- expand_grid(variety = c("A", "B"), biomass = 1:100, draw = 1:5)
d <- d %>% 
  mutate(mu = case_when(
    variety == "A" ~ 4,
    TRUE ~ 2
  )) %>% 
  rowwise() %>% 
  mutate(obs = rnorm(1, mean = mu-0.1*biomass))

d %>% 
  print(n=Inf)

# just plotting the curves
d %>% 
  ggplot(aes(x = biomass, y = obs, color = variety, group = interaction(variety, draw))) +
  geom_line()

# plotting the differences
d %>% 
  select(-mu) %>% 
  pivot_wider(names_from = variety, values_from = obs) %>% 
  mutate(diff = A - B) %>% 
  ggplot(aes(x = biomass, y = diff)) +
  geom_line()
