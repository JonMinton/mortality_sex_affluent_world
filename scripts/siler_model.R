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






