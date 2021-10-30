library(tidyverse)
library(brms)
library(tidybayes)

data_cndc <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd")


m8 <- brm(formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
           Bmax + Si ~ 1 + (1|index),
           alpha1 + alpha2 ~ 1 + (1|group),
           nl = T), 
          data = data_cndc,
prior=c(set_prior("normal(5.2,0.1)", class = "sd", nlpar = "Bmax", group = "index"),
         set_prior("normal(1.2,0.1)", class = "sd", nlpar = "Si", group = "index"),
         set_prior("normal(0.05,0.01)", nlpar = "alpha1", class = "sd", group ="group"),
         set_prior("normal(0.13,0.01)", nlpar = "alpha2", class = "sd", group ="group"),
         set_prior("normal(11.0,0.1)", nlpar = "Bmax", lb = 1),
         set_prior("normal(5.0,0.3)", nlpar = "Si", lb = 0),
         set_prior("normal(5.0,0.1)", nlpar = "alpha1", lb = 0),
         set_prior("normal(0.35,0.01)", nlpar = "alpha2", lb = 0, ub = 1),
         set_prior("student_t(3,1.3,0.1)", class = "sigma")),
iter = 5000,
warmup = 2000,
control = list(adapt_delta = 0.99),
chains = 4, cores = 4)

m8

# write_rds(m8, "brms/models/m8_no_abs_all.rds")

m8 %>%
  spread_draws(b_Si_Intercept, r_index__Si[index,]) %>%
  mutate(index_Si = b_Si_Intercept + r_index__Si) %>%
  ggplot(aes(x = index_Si, y = index)) +
  geom_halfeyeh(alpha = 0.2)

m8 %>%
  spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
  mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`) %>%
  mutate_at(vars(group),as.character) %>% 
  ggplot(aes(x = `group_alpha1`, y = `group`)) +
  stat_halfeye()

m8 %>%
  spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
  mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`) %>%
  mutate_at(vars(group),as.character) %>% 
  ggplot(aes(x = `group_alpha2`, y = `group`)) +
  stat_halfeye()

m8 %>%
  spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
  mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`) %>%
  filter(group %in% c(1,6)) %>% 
  mutate_at(vars(group),as.character) %>% 
  compare_levels(group_alpha2, by = group, comparison = "pairwise") %>% 
  ggplot(aes(x = `group_alpha2`, y = `group`)) +
  stat_halfeye()


m8_prior <- brm(formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                     Bmax + Si ~ 1 + (1|index),
                     alpha1 + alpha2 ~ 1 + (1|group),
                     nl = T), 
          data = m7$data,
          prior=c(set_prior("normal(5.2,0.1)", class = "sd", nlpar = "Bmax", group = "index"),
                  set_prior("normal(1.2,0.1)", class = "sd", nlpar = "Si", group = "index"),
                  set_prior("normal(0.05,0.01)", nlpar = "alpha1", class = "sd", group ="group"),
                  set_prior("normal(0.13,0.01)", nlpar = "alpha2", class = "sd", group ="group"),
                  set_prior("normal(11.0,0.1)", nlpar = "Bmax", lb = 1),
                  set_prior("normal(5.0,0.3)", nlpar = "Si", lb = 0),
                  set_prior("normal(5.0,0.1)", nlpar = "alpha1", lb = 0),
                  set_prior("normal(0.35,0.01)", nlpar = "alpha2", lb = 0, ub = 1),
                  set_prior("student_t(3,1.3,0.1)", class = "sigma")),
          iter = 5000,
          warmup = 2000,
          control = list(adapt_delta = 0.99),
          chains = 4, cores = 4,
          sample_prior = "only")

fmin <- function(x,y){
  pmin(x,y)
}

conditional_effects(m8_prior)

m9_prior <- brm(formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                           Bmax + Si ~ 1 + (1|index),
                           alpha1 + alpha2 ~ 1 + (1|group),
                           nl = T), 
                data = m7$data,
                prior=c(set_prior("normal(5,1)", class = "sd", nlpar = "Bmax", group = "index"),
                        set_prior("normal(1,1)", class = "sd", nlpar = "Si", group = "index"),
                        set_prior("normal(0.05,0.01)", nlpar = "alpha1", class = "sd", group ="group"),
                        set_prior("normal(0.13,0.01)", nlpar = "alpha2", class = "sd", group ="group"),
                        set_prior("normal(10,1)", nlpar = "Bmax", lb = 1),
                        set_prior("normal(5,1)", nlpar = "Si", lb = 0),
                        set_prior("normal(5.0,0.1)", nlpar = "alpha1", lb = 0),
                        set_prior("normal(0.35,0.01)", nlpar = "alpha2", lb = 0, ub = 1),
                        set_prior("student_t(3,1.3,0.1)", class = "sigma")),
                iter = 5000,
                warmup = 2000,
                control = list(adapt_delta = 0.99),
                chains = 4, cores = 4,
                sample_prior = "only")

conditional_effects(m9_prior)

m9 <- brm(formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                           Bmax + Si ~ 1 + (1|index),
                           alpha1 + alpha2 ~ 1 + (1|group),
                           nl = T), 
                data = m7$data,
                prior=c(set_prior("normal(5,1)", class = "sd", nlpar = "Bmax", group = "index"),
                        set_prior("normal(1,1)", class = "sd", nlpar = "Si", group = "index"),
                        set_prior("normal(0.05,0.01)", nlpar = "alpha1", class = "sd", group ="group"),
                        set_prior("normal(0.13,0.01)", nlpar = "alpha2", class = "sd", group ="group"),
                        set_prior("normal(10,1)", nlpar = "Bmax", lb = 1),
                        set_prior("normal(5,1)", nlpar = "Si", lb = 0),
                        set_prior("normal(5.0,0.1)", nlpar = "alpha1", lb = 0),
                        set_prior("normal(0.35,0.01)", nlpar = "alpha2", lb = 0, ub = 1),
                        set_prior("student_t(3,1.3,0.1)", class = "sigma")),
                iter = 5000,
                warmup = 2000,
                control = list(adapt_delta = 0.99),
                chains = 4, cores = 4)
