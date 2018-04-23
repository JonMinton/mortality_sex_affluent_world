rm(list = ls())

library(plotly)
library(tidyverse)
library(readxl)
# Aim of this is to try to fit the siler model by hand 

dta <- read_csv("data/counts.csv")


code_to_country_lookup <- read_excel(
  "support/replication_details.xlsx",
  sheet="code_to_country_lookup"
)

full_country_lookup <- code_to_country_lookup  %>% 
  filter(include_in_full_selection == 1)  %>% 
  .$code


# Original Table 1 : Country, years available, population in (say) 2010 or latest available year



# Derived data ------------------------------------------------------------



dta_selection <- dta %>% 
  filter(sex !="total" & country %in% full_country_lookup) %>% 
  group_by(year, age, sex) %>% 
  summarise(
    population_count = sum(population_count), 
    death_count = sum(death_count)
  ) %>% 
  filter(year >= 1850 )


# Trigonometric quasi-siler function --------------------------------------


# The fundamental idea of the model is to represent mortality change over the life 
# course with three lines over log-mortality.

# The first line is the infantile schedule, with intercept at x = 0 of A_I
# and gradient determined by Theta_I
# The second line is horizontal, with intercept A_R
# The third line is the senescent schedule, with intercept at x = x_max of A_S
# and gradient determined by Theta_S

# The model picks the maximum of the three independent functions at each x value 

# Optim searches over five parameters over the full real number plane 

# These are then mapped to appropriate ranges using either exponential or 
# rescaled logit functions as follows

# A_I = -exp(a_i)
# Theta_I = (pi / 2) * (1 / (1 + exp(-theta_i))) = angler(theta_i)
# A_R = -exp(a_r)
# A_S = -exp(a_s)
# Theta_S = (pi / 2) * (1 / (1 + exp(-theta_s))) = angler(theta_s)

# The infantile schedule is 
# f_i = A_I - x / tan(Theta_I)

# The senescent schedule is 
# f_s = A_S + (x - x_max) / tan(Theta_S)

# And the risk schedule is simply 
# f_r = A_R 

# The combined function is simply 
# Quasi_Siler(x) = max(f_i(x), f_r(x), f_s(x))


# Parameters PAR := {a_i, theta_i, a_r, a_s, theta_s} 
# Parameters are selected via optim by minimising product of RMS of fit to data in 
# year t and sum of squared distances between PAR(t) and PAR(t-1)

# i.e. 
# loss(PAR) := RMS(PAR, data) * SSQ(PAR, PAR_lag)
# Where PAR_lag := {a_i(t-1), theta_i(t-1), a_r(t-1), a_s(t-1), theta_s(t-1)}

# First trial 
#pars <- runif(5, -2, 2)

to_angle <- function(x){
  (pi / 2) * (1 / (1 + exp(-x)))
}

from_angle <- function(y){
  -log(((pi/2) / y) - 1)
}

# Transform from full real to range 0 to -10 (reasonable mort space range)

to_mortspace <- function(x){
  -10 * (1 / (1 + exp(-x)))
}

from_mortspace <- function(y){
  -log((-10 / y) - 1)
}

trans_par <- function(pars){
  out <- c(
    to_mortspace(pars[1]),
    to_angle(pars[2]),
    to_mortspace(pars[3]),
    to_mortspace(pars[4]),
    to_angle(pars[5])
  )
  out 
}


#xfrmed_pars <- trans_par(pars)

quasi_siler <- function(xfrm_pars, max_age = 95){
  x <- 0:max_age
  
  A_I <- xfrm_pars[1]
  Theta_I <- xfrm_pars[2]
  
  A_R <- xfrm_pars[3] 
  
  A_S <- xfrm_pars[4]
  Theta_S <- xfrm_pars[5]
  
  schedule_infant <- A_I - x / tan(Theta_I)
  baseline_risk <- A_R
  schedule_senescent <- A_S + (x - max_age) / tan(Theta_S)
  
  out <- pmax(
    schedule_infant,
    baseline_risk,
    schedule_senescent
  )
  
  out
}

calc_rms <- function(pred_schedule, actual_schedule){
  (actual_schedule - pred_schedule) %>% 
    .^2 %>% 
    mean() %>% 
    .^0.5
}

calc_sq_param_shift <- function(new_pars, old_pars){
  (old_pars - new_pars) %>% 
    .^2 %>% 
    mean()
}

do_siler_rms <- function(pars, observed_schedule){
  xfrmed_pars <- trans_par(pars)
  predicted_schedule <- quasi_siler(xfrmed_pars)
  calc_rms(predicted_schedule, observed_schedule)
}

do_siler_jointloss <- function(pars, observed_schedule, old_pars){
  xfrmed_pars <- trans_par(pars)
  predicted_schedule <- quasi_siler(xfrmed_pars)
  rms <- calc_rms(predicted_schedule, observed_schedule)
  shift_cost <- calc_sq_param_shift(pars, old_pars)
  return(rms * shift_cost)
}

mort_schedules_df <- dta_selection %>% 
  mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
  filter(age <= 95) %>% 
  arrange(sex, year, age) %>% 
  group_by(sex, year) %>% 
  nest() %>% 
  mutate(lmr_schedule = map(data, "lmr")) 


# First year for females 

first_schedule_females <- mort_schedules_df %>% 
  filter(year == min(year)) %>% 
  filter(sex == "female") %>% 
  .[["lmr_schedule"]] %>% .[[1]]

#First year for males

first_schedule_males <- mort_schedules_df %>% 
  filter(year == min(year)) %>% 
  filter(sex == "male") %>% 
  .[["lmr_schedule"]] %>% .[[1]]


# Pull fit from many runs of optim using replicate

calc_starting_fit <- function(par, starting_schedule){
  optim(
    par = par,
    fn = do_siler_rms,
    observed_schedule = starting_schedule
  ) -> tmp
  c(
    value = tmp[["value"]],
    par1 = tmp[["par"]][[1]],
    par2 = tmp[["par"]][[2]],
    par3 = tmp[["par"]][[3]],
    par4 = tmp[["par"]][[4]],
    par5 = tmp[["par"]][[5]]
  )
}

num_initial_replicates <- 10000
num_later_replicates <- 10000


set.seed(5)
many_runs_female <- replicate(
  num_initial_replicates, 
  calc_starting_fit(par = runif(5, -8, 8), starting_schedule = first_schedule_females)
)

set.seed(5)
many_runs_male <- replicate(
  num_initial_replicates, 
  calc_starting_fit(par = runif(5, -8, 8), starting_schedule = first_schedule_males)
)


#install.packages("plotly")

tidy_runs <- many_runs_female %>% 
  t %>% 
  as_tibble %>% 
  mutate(sex = "female") %>% 
  bind_rows(
    many_runs_male %>% 
      t %>% 
      as_tibble %>% 
      mutate(sex = "male")
  ) 

tidy_runs %>%
  mutate(is_male = as.numeric(sex =='male')) %>% 
  plotly::plot_ly(type = 'parcoords',
                  #           color = ~sex, colors = "Set1",
                  #          linetype = ~sex,
                  line = list(color = ~is_male,
                              colors = "Set1"),
                  dimensions = list(
                    list(
                      label = 'Fit', values = ~value),
                    list(range = c(-8, 8),
                         label = 'Infant Intercept', values = ~par1),
                    list(range = c(-8, 8),
                         #                 constraintrange = c(5,6),
                         label = 'Infant Angle', values = ~par2),
                    list(range = c(-8, 8),
                         label = 'Background Intercept', values = ~par3),
                    list(range = c(-8, 8),
                         label = 'Senescent Intercept', values = ~par4),
                    list(range = c(-8, 8),
                         label = 'Senescent Angle', values = ~par5)
                    
                  )
  )

# This plot provides some confidence that the best fitting models are finding reasonably similar 
# estimates for each of the parameters, and also that the range -5 to 5 should not constrain values too severely 
# within this range 
# Maybe use -6 to 6 to be safe? 

# The next stage will feed the estimates from the previous year into the next year's estimates

# It will start with the best estimates for each sex from the first year



# Now to pick the best fits by gender 

tidy_runs %>% 
  group_by(sex) %>% 
  filter(value == min(value))

# # A tibble: 2 x 7
# # Groups:   sex [2]
# value      par1      par2        par3      par4     par5    sex
# <dbl>     <dbl>     <dbl>       <dbl>     <dbl>    <dbl>  <chr>
#   1 0.1834397 -1.223117 0.9318067 -0.07897006 -2.260888 3.036875 female
# 2 0.2025954 -1.294752 0.8536778 -0.04156170 -2.322711 3.075348   male

best_init_pars_female <- tidy_runs %>% 
  filter(sex == "female") %>% 
  filter(value == min(value)) %>% 
  select(par1:par5) %>% 
  gather() %>% 
  deframe()

best_init_pars_male <- tidy_runs %>% 
  filter(sex == "male") %>% 
  filter(value == min(value)) %>% 
  select(par1:par5) %>% 
  gather() %>% 
  deframe()

# n.b. may need to convert to list before passing to optim
# 
# mort_schedules_df_female <- mort_schedules_df %>% filter(sex == "female")
# mort_schedules_df_male <-   mort_schedules_df %>% filter(sex == "male")
# 
# pars_female <- vector("list", length = nrow(mort_schedules_df_female))
# pars_male <- vector("list", length = nrow(mort_schedules_df_male))
# 
# pars_female[[1]] <- best_init_pars_female
# pars_male[[1]] <- best_init_pars_male
# 
# for (i in 2:nrow(mort_schedules_df_female)){
#   lag_par_female <- pars_female[[i-1]]
#   lag_par_male <- pars_male[[i-1]]
#   
#   this_mort_schedule_male   <-  mort_schedules_df_male[["lmr_schedule"]][[i]]
#   this_mort_schedule_female <-  mort_schedules_df_female[["lmr_schedule"]][[i]]
#   
#   set.seed(5)
#   many_runs_male <- replicate(
#     num_later_replicates, 
#     {
#       optim(
#         par = runif(5, -6, 6),
#         do_siler_jointloss,
#         observed_schedule = this_mort_schedule_male,
#         old_pars = lag_par_male
#       ) -> tmp
#       c(
#         value = tmp[["value"]],
#         par1 = tmp[["par"]][[1]],
#         par2 = tmp[["par"]][[2]],
#         par3 = tmp[["par"]][[3]],
#         par4 = tmp[["par"]][[4]],
#         par5 = tmp[["par"]][[5]]
#       )
#     }
#   )
# 
#   set.seed(5)
#   many_runs_female <- replicate(
#     num_later_replicates, 
#     {
#       optim(
#         par = runif(5, -6, 6),
#         do_siler_jointloss,
#         observed_schedule = this_mort_schedule_female,
#         old_pars = lag_par_female
#       ) -> tmp
#       c(
#         value = tmp[["value"]],
#         par1 = tmp[["par"]][[1]],
#         par2 = tmp[["par"]][[2]],
#         par3 = tmp[["par"]][[3]],
#         par4 = tmp[["par"]][[4]],
#         par5 = tmp[["par"]][[5]]
#       )
#     }
#   )
#   
# 
#   tidied_later_runs <- many_runs_female %>% 
#     t %>% 
#     as_tibble %>% 
#     mutate(sex = "female") %>% 
#     bind_rows(
#       many_runs_male %>% 
#         t %>% 
#         as_tibble %>% 
#         mutate(sex = "male")
#     ) 
# 
#   best_this_pars_female <- tidied_later_runs %>% 
#     filter(sex == "female") %>% 
#     filter(value == min(value)) %>% 
#     select(par1:par5) %>% 
#     gather() %>% 
#     deframe()
# 
#   best_this_pars_male <- tidied_later_runs %>% 
#     filter(sex == "male") %>% 
#     filter(value == min(value)) %>% 
#     select(par1:par5) %>% 
#     gather() %>% 
#     deframe()
#   
#   pars_female[[i]] <- best_this_pars_female
#   pars_male[[i]] <- best_this_pars_male
#     
# }
# 
# pars_df <- do.call(bind_rows, pars_female) %>% 
#   mutate(sex = "female") %>% 
#   mutate(year = 1850:2010) %>% 
#   bind_rows(
#     do.call(bind_rows, pars_male) %>% 
#       mutate(sex = "male") %>% 
#       mutate(year = 1850:2010)
#   )
# 
# pars_df %>% 
#   gather(par1:par5, key = "par", value = "value") %>% 
#   ggplot(aes(x = year, y = value, colour = sex)) +
#   geom_line() + facet_wrap(~par, scales = "free_y")


# This approach technically works, but the level of inertia looks far too strong.

# As a first alternative approach I'll try feeding the previous year's values as the 
# starting estimates for the next year 
# As second alternative approach will be to look at ways of systematically 
# weighting the relative importance of the two components before adding their products 

# First alternative - no second weighting, fed forward starting values 


mort_schedules_df_female <- mort_schedules_df %>% filter(sex == "female")
mort_schedules_df_male <-   mort_schedules_df %>% filter(sex == "male")

pars_female <- vector("list", length = nrow(mort_schedules_df_female))
pars_male <- vector("list", length = nrow(mort_schedules_df_male))

pars_female[[1]] <- best_init_pars_female
pars_male[[1]] <- best_init_pars_male

for (i in 2:nrow(mort_schedules_df_female)){
  lag_par_female <- pars_female[[i-1]]
  lag_par_male <- pars_male[[i-1]]
  
  this_mort_schedule_male   <-  mort_schedules_df_male[["lmr_schedule"]][[i]]
  this_mort_schedule_female <-  mort_schedules_df_female[["lmr_schedule"]][[i]]
  
  set.seed(5)
  many_runs_male <- replicate(
    num_later_replicates, 
    {
      optim(
        par = lag_par_male + runif(5, -4, 4),
        do_siler_rms,
        observed_schedule = this_mort_schedule_male
      ) -> tmp
      c(
        value = tmp[["value"]],
        par1 = tmp[["par"]][[1]],
        par2 = tmp[["par"]][[2]],
        par3 = tmp[["par"]][[3]],
        par4 = tmp[["par"]][[4]],
        par5 = tmp[["par"]][[5]]
      )
    }
  )

  
  set.seed(5)
  many_runs_female <- replicate(
    num_later_replicates, 
    {
      optim(
        par = lag_par_female + runif(5, -4, 4),
        do_siler_rms,
        observed_schedule = this_mort_schedule_female
      ) -> tmp
      c(
        value = tmp[["value"]],
        par1 = tmp[["par"]][[1]],
        par2 = tmp[["par"]][[2]],
        par3 = tmp[["par"]][[3]],
        par4 = tmp[["par"]][[4]],
        par5 = tmp[["par"]][[5]]
      )
    }
  )
  
  
  tidied_later_runs <- many_runs_female %>% 
    t %>% 
    as_tibble %>% 
    mutate(sex = "female") %>% 
    bind_rows(
      many_runs_male %>% 
        t %>% 
        as_tibble %>% 
        mutate(sex = "male")
    ) 
  
  best_this_pars_female <- tidied_later_runs %>% 
    filter(sex == "female") %>% 
    filter(value == min(value)) %>% 
    select(par1:par5) %>% 
    gather() %>% 
    deframe()
  
  best_this_pars_male <- tidied_later_runs %>% 
    filter(sex == "male") %>% 
    filter(value == min(value)) %>% 
    select(par1:par5) %>% 
    gather() %>% 
    deframe()
  
  pars_female[[i]] <- best_this_pars_female
  pars_male[[i]] <- best_this_pars_male
  
}

pars_df <- do.call(bind_rows, pars_female) %>% 
  mutate(sex = "female") %>% 
  mutate(year = 1850 + 0:162) %>% 
  bind_rows(
    do.call(bind_rows, pars_male) %>% 
      mutate(sex = "male") %>% 
      mutate(year = 1850 + 0:162)
  )

pars_df %>% 
  gather(par1:par5, key = "par", value = "value") %>% 
  ggplot(aes(x = year, y = value, colour = sex)) +
  geom_line() + facet_wrap(~par, scales = "free_y")


# Now let's transform these pars to get coefficients and angles 

pars_df %>% 
  mutate(
    infant_intercept = to_mortspace(par1),
    infant_angle = to_angle(par2),
    baseline_risk = to_mortspace(par3),
    senescent_intercept = to_mortspace(par4),
    senescent_angle = to_angle(par5)
  ) %>% 
  select(sex:senescent_angle) %>% 
  gather(infant_intercept:senescent_angle, key = "param", value = "value") %>% 
  mutate(param = factor(
    param, 
    levels = c(
      "infant_intercept", "infant_angle", 
      "baseline_risk", 
      "senescent_intercept", "senescent_angle"
      ),
    ordered = T
    )
  ) %>% 
  ggplot(aes(x = year, y = value, colour = sex)) + 
  geom_line() + facet_wrap(~param, scale = "free_y")


# Write out the params

pars_df

pars_df %>% 
  mutate(
    infant_intercept = to_mortspace(par1),
    infant_angle = to_angle(par2),
    baseline_risk = to_mortspace(par3),
    senescent_intercept = to_mortspace(par4),
    senescent_angle = to_angle(par5)
  ) %>% 
  mutate(
    infant_gradient = atan(infant_angle),
    senescent_gradient = atan((pi / 2) - senescent_angle)
  ) %>% 
  select(-infant_angle, -senescent_angle) %>% 
  select(sex:senescent_gradient) %>% 
  gather(infant_intercept:senescent_gradient, key = "param", value = "value") %>% 
  mutate(param = factor(
    param, 
    levels = c(
      "infant_intercept", "infant_gradient", 
      "baseline_risk", 
      "senescent_intercept", "senescent_gradient"
    ),
    ordered = T
  )
  ) %>% 
  ggplot(aes(x = year, y = value, colour = sex)) + 
  geom_line() + facet_wrap(~param, scale = "free_y")

ggsave("figures/siler_best_to_2010.png", dpi = 300, units = "cm", height = 20, width = 20)
write_csv(pars_df, "data/quasi_siler_best_estimates.csv")
# 
# 
# # Let's look at the empirical schedules to get a sense of plausible ranges of values to 
# # consider 
# 
# dta_selection %>% 
#   mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
#   ggplot(aes(x = age, y = lmr, colour = sex)) +
#   geom_point(alpha = 0.05) + 
#   geom_vline(xintercept = 95) + 
#   geom_vline(xintercept = 10)
# 
# # what's the range at age 0, 25, and 95? 
# 
# dta_selection %>% 
#   mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
#   filter(age %in% c(0, 25, 95)) %>% 
#   group_by(age) %>% 
#   summarise(
#     min   = min(lmr),
#     lwr   = quantile(lmr, 0.025),
#     mdlwr = quantile(lmr, 0.25),
#     md    = quantile(lmr, 0.50),
#     mdupr = quantile(lmr, 0.75),
#     upr   = quantile(lmr, 0.975),
#     max   = max(lmr)
#   )
#   
# # Ranges 
# # Infancy 
# # -5.365 to -1.444
# 
# # 25 
# # -7.649 to -2.143
# 
# # Elderly 
# # -1.424 to -0.629
# 
# 
# # Now what about the plausible angles?
# 
# # Ranges to look over are 0 to 10
# # and 40 to 95 
# 
# # First 0 to 10
# 
# dta_selection %>% 
#   mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
#   filter(age <= 10) %>% 
#   group_by(year, sex) %>% 
#   nest() %>% 
#   mutate(mdl = map(data, ~ lm(lmr ~ age, data = .))) %>% 
#   mutate(coef = map(mdl, coefficients)) %>% 
#   mutate(
#     intercept = map_dbl(coef, 1),
#     gradient = map_dbl(coef, 2),
#     angle = atan(1/gradient),
#     angle_deg = (180 / pi) * angle
#     ) %>% 
#   select(year, sex, intercept, gradient, angle, angle_deg) %>% 
#   mutate(type = "infancy") -> angles_infantile
# 
# # Now elderly 
# dta_selection %>% 
#   mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
#   filter(age >= 40 & age <= 95) %>% 
#   group_by(year, sex) %>% 
#   nest() %>% 
#   mutate(mdl = map(data, ~ lm(lmr ~ age, data = .))) %>% 
#   mutate(coef = map(mdl, coefficients)) %>% 
#   mutate(
#     intercept = map_dbl(coef, 1),
#     gradient = map_dbl(coef, 2),
#     angle = atan(1/gradient),
#     angle_deg = (180 / pi) * angle
#   ) %>% 
#   select(year, sex, intercept, gradient, angle, angle_deg) %>% 
#   mutate(type = "elderly") -> angles_elderly
# 
# angles_both <- bind_rows(angles_infantile, angles_elderly)
# # Visualise
# 
# angles_both %>% 
#   select(year, sex, angle_deg, type) %>% 
#   spread(type, angle_deg) %>% 
#   ggplot(aes(x = infancy, y = elderly, colour = sex)) + 
#   geom_point()
# 
# angles_both %>% 
#   select(year, sex, angle_deg, type) %>% 
#   spread(type, angle_deg) %>% 
#   ggplot() +
#   geom_path(aes(x = infancy, y = elderly), data = . %>% filter(sex == "male"), colour = "blue") + 
#   geom_path(aes(x = infancy, y = elderly), data = . %>% filter(sex == "female"), colour = "red")  
#   
# 
# 
# 



# Older material  ---------------------------------------------------------


# siler model 

# The siler model specification is 

# q_x = a1 * exp(-b1*x) + a2 + a3+exp(b3*x)

calc_siler <- function(
  pars = list(
    a1 = 0.01,
    b1 = 3, 
    a2 = 1e-4,
    a3 = 1e-5,
    b3 = 0.1
    ),
  x = 0:95
){
  
  a1 = pars[["a1"]]
  b1 = pars[["b1"]]
  a2 = pars[["a2"]]
  a3 = pars[["a3"]]
  b3 = pars[["b3"]]
  
  pred_q_x = a1*exp(-b1*x) + a2 + a3*exp(b3*x)
  
  return(pred_q_x)
}

calc_rms <- function(
  q_x,
  pred_q_x
){
  
  if (length(q_x) != length(pred_q_x)){
    stop("length of vectors to compare are not equal")
  }
  
  output <- (pred_q_x - q_x) %>% .^2 %>% mean() %>% .^(0.5)
  return(output)
}

# All parameters need to be positive, 
# so search should be over exp(f) rather than f

fit_model_and_get_rms <- function(
  PARS,
  q_x,
  x
){
  pars = lapply(PARS, exp) # convert from -inf to inf space to 0,Inf range 
  pred_q_x <- calc_siler(pars, x)
  calc_rms(q_x, pred_q_x)
}

# First dataset, females, 1950, ages 0:95

dta_selection %>% 
  ungroup() %>% 
  filter(age <= 95) %>% 
  filter(year == 1950) %>% 
  filter(sex == "female") %>% 
  mutate(q_x = death_count / population_count) %>% 
  arrange(age) -> tmp 

tmp$q_x

tmp$age

optim(
  par = list(
    a1 = log(0.01),
    b1 = log(3), 
    a2 = log(1e-4),
    a3 = log(1e-5),
    b3 = log(0.1)
  ),
  fn = fit_model_and_get_rms,
  q_x = tmp$q_x,
  x = tmp$age
) -> opt_out

p2 <- lapply(opt_out$par, exp)

siler_pred <- calc_siler(p2, x = 0:95)


# Now for males in the same year 

dta_selection %>% 
  ungroup() %>% 
  filter(age <= 95) %>% 
  filter(year == 1950) %>% 
  filter(sex == "male") %>% 
  mutate(q_x = death_count / population_count) %>% 
  arrange(age) -> tmp 

tmp$q_x

tmp$age

run_optim <- function(
  q_x, x, 
  init_pars = list(
    a1 = log(0.01),
    b1 = log(3), 
    a2 = log(1e-4),
    a3 = log(1e-5),
    b3 = log(0.1)
    )
  ){
  optim(
    par = init_pars,
    fn = fit_model_and_get_rms,
    q_x = q_x,
    x = x
  ) -> opt_out
  return(opt_out)
}

p2_m <- lapply(opt_out$par, exp)

siler_pred_m <- calc_siler(p2, x = 0:95)

dta_selection %>% 
  ungroup() %>% 
  filter(age <= 95) %>% 
  arrange(year, sex, age) %>% 
  mutate(q_x = death_count / population_count) %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(q_x = map(data, "q_x")) %>% 
  mutate(x = map(data, "age")) %>% 
  mutate(optim_out = map2(q_x, x, run_optim)) %>% 
  select(year, sex, optim_out) %>% 
  mutate(rms = map_dbl(optim_out, "value")) %>% 
  mutate(pars = map(optim_out, function(x) lapply(x[["par"]], exp))) %>% 
  select(-optim_out) -> tmp 

output <- bind_cols(
  tmp %>% select(year, sex, rms),
  pars = bind_rows(tmp$pars)
)

rm(tmp)

output %>% 
  gather(key = "param", value = "value", a1:b3) %>% 
  ggplot(aes(x = year, y = value, colour = sex)) + 
  facet_wrap(~param, scale = "free_y") + 
  geom_line()

# This sort of works but again there seems a lot of instability in places, especially for older 
# periods 

# Trying different configurations of starting values therefore seems important . 
# Let's try to write a wrapper around the earlier code to repeat many times 

rerun_optim <- function(q_x, x, SEED = 5, REPLICATES = 100){
  set.seed(SEED)
  
  randomise_pars_and_run_optim <- function(q_x, x){
    init_vars <- runif(5, -3, 3)
    init_pars = list(
      a1 = init_vars[1],
      b1 = init_vars[2], 
      a2 = init_vars[3],
      a3 = init_vars[4],
      b3 = init_vars[5]
    )
    try(run_optim(q_x, x, init_pars))
  }
  replicate(
    REPLICATES,
    randomise_pars_and_run_optim(q_x, x)
  ) -> output

  return(output)
}


dta_selection %>% 
  ungroup() %>% 
  filter(age <= 95) %>% 
  arrange(year, sex, age) %>% 
  mutate(q_x = death_count / population_count) %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(q_x = map(data, "q_x")) %>% 
  mutate(x = map(data, "age")) %>% 
  mutate(optim_runs = map2(q_x, x, rerun_optim)) -> new_output

new_output$optim_runs[[1]]




# New attempt - run over smoothed data  -----------------------------------

# First attempt: one year lag and lead

dta_selection %>% 
  mutate(lmr = log(death_count / population_count) ) %>% 
  group_by(sex, age) %>% 
  arrange(year) %>% 
  mutate(lag_lmr = lag(lmr, 1, default = NA)) %>% 
  mutate(lead_lmr = lead(lmr, 1, default = NA)) %>% 
  mutate(sm_lmr = case_when(
    is.na(lag_lmr) ~ (lmr + lead_lmr ) / 2,
    is.na(lead_lmr) ~ (lmr + lag_lmr) / 2,
    TRUE ~ (lag_lmr + lmr + lead_lmr) / 3
    )
  ) %T>% sample_n(10) %>% 
  mutate(q_x = exp(sm_lmr)) %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(q_x = map(data, "q_x")) %>% 
  mutate(x = map(data, "age")) %>% 
  mutate(optim_out = map2(q_x, x, run_optim)) %>% 
  select(year, sex, optim_out) %>% 
  mutate(rms = map_dbl(optim_out, "value")) %>% 
  mutate(pars = map(optim_out, function(x) lapply(x[["par"]], exp))) %>% 
  select(-optim_out) -> tmp 

output <- bind_cols(
  tmp %>% select(year, sex, rms),
  pars = bind_rows(tmp$pars)
)

rm(tmp)

output %>% 
  gather(key = "param", value = "value", a1:b3) %>% 
  ggplot(aes(x = year, y = value, colour = sex)) + 
  facet_wrap(~param, scale = "free_y") + 
  geom_line()






