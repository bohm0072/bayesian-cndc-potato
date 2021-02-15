library(tidybayes)
library(tidyverse)

set.seed(25)

tg <- tibble(g1 = rnorm(1000, mean = 5, sd = 1),
             g2 = rnorm(1000, mean = 3, sd = 1))

tg %>% 
  pivot_longer(cols = c(g1,g2)) %>% 
  group_by(name) %>% 
  mutate(median = median(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = value, y = name)) +
  stat_halfeye() +
  geom_vline(aes(xintercept = median), linetype = 2)

tg %>% 
  mutate(diff = g1 - g2) %>% 
  ggplot(aes(x = diff, y = 1)) +
  stat_halfeye()


# now obviously the order of g1 vs g2 is going to matter in terms of calculating the difference- we would likely get a different answer if we shuffled up the orders of the same values in g1 and g2. However, in the context of a Bayesian posterior, the values are inherently linked since they are coming from the same draw. So this is essentially the same exact approach we would take using posterior draw values, and we can see that the direct g1-g2 comparison gives a different answer than comparing the median and credible intervals for each.
