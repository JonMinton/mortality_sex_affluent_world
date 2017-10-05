rm(list = ls())

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
  filter(year >= 1850 & year <=2010)


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



# Let's try to think about this graphically and interactively 
# using the manipulate package. 

# The aim is to fit five parameters such that RMS is minimised

library(manipulate)

estimate_and_plot_log_siler <- function(
  a1, b1, 
  a2, 
  a3, b3

){
  x <- 0:95
  
  prediction <- (a1 - b1 * x) + a2 + (a3 + b3 * x)
  
  df <- data_frame(x = x, pred = prediction)
  
  df %>% 
    ggplot(., aes(x = x, y = pred)) + 
    geom_line() -> p1
  
  print(p1)
}

manipulate(
  estimate_and_plot_log_siler(a1, b1, a2, a3, b3),

  a1 = slider(0, 100),
  b1 = slider(0, 100),
  a2 = slider(0, 20),
  a3 = slider(-100, 100),
  b3 = slider(0, 20)
)




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

# A_I = exp(a_i)
# Theta_I = (pi / 2) * (1 / (1 + exp(-theta_i))) = angler(theta_i)
# A_R = exp(a_r)
# A_S = exp(a_s)
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




