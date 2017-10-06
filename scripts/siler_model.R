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

calc_sq_param_shift <- function(old_pars, new_pars){
  (old_pars - new_pars) %>% 
    .^2 %>% 
    sum()
}

do_siler_rms <- function(pars, observed_schedule){
  xfrmed_pars <- trans_par(pars)
  predicted_schedule <- quasi_siler(xfrmed_pars)
  calc_rms(predicted_schedule, observed_schedule)
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

first_schedule_males <- mort_schedules_df %>% 
  filter(year == min(year)) %>% 
  filter(sex == "male") %>% 
  .[["lmr_schedle"]] %>% .[[1]]

female_start_01 <- optim(
  par = runif(5, -5, 5), 
  fn = do_siler_rms, 
  observed_schedule = first_schedule_females,
  control = list(trace = 1)
  )

female_start_02 <- optim(
  par = runif(5, -5, 5), 
  fn = do_siler_rms, 
  observed_schedule = first_schedule_females,
  control = list(trace = 1)
)

female_start_03 <- optim(
  par = runif(5, -5, 5), 
  fn = do_siler_rms, 
  observed_schedule = first_schedule_females,
  control = list(trace = 1)
)

female_start_04 <- optim(
  par = runif(5, -5, 5), 
  fn = do_siler_rms, 
  observed_schedule = first_schedule_females,
  control = list(trace = 1)
)

female_start_05 <- optim(
  par = runif(5, -5, 5), 
  fn = do_siler_rms, 
  observed_schedule = first_schedule_females,
  control = list(trace = 1)
)

female_fits <- c(
  female_start_01$value,
  female_start_02$value,
  female_start_03$value,
  female_start_04$value,
  female_start_05$value
)

# Pull fit from many runs of optim using replicate

calc_starting_fit <- function(par){
  optim(
    par = par,
    fn = do_siler_rms,
    observed_schedule = first_schedule_females
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

many_runs <- replicate(1000, calc_starting_fit(par = runif(5, -5, 5)))

# reps <- 100 
# i <- 1
# while (i  < reps){
#   this_pars <- runif(5, -5, 5)
#   xfrmed_pars <- trans_par(this_pars)
#   
#   this_schedule <- quasi_siler(xfrmed_pars)
#   
#   plot(this_schedule, type = "l")
#   print(this_pars)
#   browser()
#   
# }



# Let's look at the empirical schedules to get a sense of plausible ranges of values to 
# consider 

dta_selection %>% 
  mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
  ggplot(aes(x = age, y = lmr, colour = sex)) +
  geom_point(alpha = 0.05) + 
  geom_vline(xintercept = 95) + 
  geom_vline(xintercept = 10)

# what's the range at age 0, 25, and 95? 

dta_selection %>% 
  mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
  filter(age %in% c(0, 25, 95)) %>% 
  group_by(age) %>% 
  summarise(
    min   = min(lmr),
    lwr   = quantile(lmr, 0.025),
    mdlwr = quantile(lmr, 0.25),
    md    = quantile(lmr, 0.50),
    mdupr = quantile(lmr, 0.75),
    upr   = quantile(lmr, 0.975),
    max   = max(lmr)
  )
  
# Ranges 
# Infancy 
# -5.365 to -1.444

# 25 
# -7.649 to -2.143

# Elderly 
# -1.424 to -0.629


# Now what about the plausible angles?

# Ranges to look over are 0 to 10
# and 40 to 95 

# First 0 to 10

dta_selection %>% 
  mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
  filter(age <= 10) %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(mdl = map(data, ~ lm(lmr ~ age, data = .))) %>% 
  mutate(coef = map(mdl, coefficients)) %>% 
  mutate(
    intercept = map_dbl(coef, 1),
    gradient = map_dbl(coef, 2),
    angle = atan(1/gradient),
    angle_deg = (180 / pi) * angle
    ) %>% 
  select(year, sex, intercept, gradient, angle, angle_deg) %>% 
  mutate(type = "infancy") -> angles_infantile

# Now elderly 
dta_selection %>% 
  mutate(lmr = ((death_count + 0.5) / (population_count + 0.5)) %>% log) %>% 
  filter(age >= 40 & age <= 95) %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(mdl = map(data, ~ lm(lmr ~ age, data = .))) %>% 
  mutate(coef = map(mdl, coefficients)) %>% 
  mutate(
    intercept = map_dbl(coef, 1),
    gradient = map_dbl(coef, 2),
    angle = atan(1/gradient),
    angle_deg = (180 / pi) * angle
  ) %>% 
  select(year, sex, intercept, gradient, angle, angle_deg) %>% 
  mutate(type = "elderly") -> angles_elderly

angles_both <- bind_rows(angles_infantile, angles_elderly)
# Visualise

angles_both %>% 
  select(year, sex, angle_deg, type) %>% 
  spread(type, angle_deg) %>% 
  ggplot(aes(x = infancy, y = elderly, colour = sex)) + 
  geom_point()

angles_both %>% 
  select(year, sex, angle_deg, type) %>% 
  spread(type, angle_deg) %>% 
  ggplot() +
  geom_path(aes(x = infancy, y = elderly), data = . %>% filter(sex == "male"), colour = "blue") + 
  geom_path(aes(x = infancy, y = elderly), data = . %>% filter(sex == "female"), colour = "red")  
  







# Reasonable range of values 

# Infancy 

# -1.5 to -8.5
# -log(1.5) ~= -0.41
# -log(8.5) ~= -2.14

# Adult 

# -2.5 to -8.0
# -log(2.5) ~= -0.92
# -log(8.0) ~= -2.08

# Elderly 

# -0.1 to -2.0
# -log(0.1) ~=  2.30
# -log(2.0) ~= -0.69


# Attempt with these params 

reps <- 100 
i <- 1
while (i  < reps){
  this_pars <- c(
      a_i = runif(1, -2.14, -0.41),
      theta_i = runif(1, -3, 3),
      a_r = runif(1, -2.08, -0.92),
      a_s = runif(1, -0.69, 2.30),
      theta_s = runif(1, -3, 3)
    )
  
  xfrmed_pars <- trans_par(this_pars)
  
  this_schedule <- quasi_siler(xfrmed_pars)
  
  plot(this_schedule, type = "l")
  print(this_pars)
  print(xfrmed_pars)
  browser()
  
}

