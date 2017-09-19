## HPBayes Model application


# The aim now is to learn to apply the 8 paramater Heligman Pollard Mortality model 
# to the HMD.

# Most of the parameters are of substantive interest and interpretability 

# A key way of operationalising the hunch about change from period-based to 
# cohort based trends in convergence then divergence in mortality sex ratios is to 
# see if there's a trend for model fit by period to fall, and model fit by cohort to 
# increase, and further for most of this increase in fit by cohort toe 
# Modelling 


rm(list=ls())

# Prerequisite packages ---------------------------------------------------


pacman::p_load(
  MASS,
  readxl,
  tidyverse,
  grid, ggplot2,
  lattice, 
  latticeExtra,
  HPbayes,
  MortHump
)




# Data --------------------------------------------------------------------


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





# Try example code with HPbayes -------------------------------------------


data(HPprior)

##number of deaths in each age group##
dx <- c(68, 47, 16, 10, 13, 29, 92, 151, 188, 179, 156, 155, 147, 150,
        122, 106, 88, 113, 63, 38, 32, 8)
##number at risk in each age group##
lx <- c(1974, 1906, 1860, 1844, 1834, 1823, 1793, 1700, 1549, 1361,
        1181, 1025, 870, 721, 571, 450, 344, 256, 142, 79, 41, 8)
result <- hp.bm.imis(prior=q0, K=10, nrisk=lx, ndeath=dx)
## End(Not run)

# That was fairly simple, now to try with real data

# 1910, males 

dta_selection %>% 
  ungroup() %>% 
  filter(year == 1910) %>% 
  filter(sex == "male") -> tmp

dx <- as.integer(tmp$death_count + 1)
lx <- as.integer(tmp$population_count + 1)
age <- tmp$age + 0.5

age[1] <- 0.01 # Infant 

result_males <- hp.bm.imis(
  prior = q0,
  K = 10,
  nrisk = lx,
  ndeath = dx,
  age = age
)



# Trying implementation in MortalityLaws package instead ------------------
# https://rdrr.io/github/mpascariu/MortalityLaws/f/inst/doc/MortalityLaws-vignette.pdf


devtools::install_github("mpascariu/MortalityLaws")

require(MortalityLaws)

dta_selection %>% ungroup() %>% filter(year == 1910, sex == "male") %>% arrange(age) -> tmp2
ages <- tmp2$age
deaths = tmp2$death_count
exposure = tmp2$population_count

fit <- MortalityLaw(x = ages,
                    Dx = deaths, # vector with death counts
                    Ex = exposure, # vector containing exposures
                    law = "HP4",
                    how = "poissonL",
                    fit.this.x = 0:90
                    )

# use 
# avaiableLaws()
# to see different curves that can be fit 

# HP4 tended to fit better on the original data used in the 1980 paper
# and to better fit the data this time too 

# The coefficients can be accessed as fit$coefficients 


# Now to vectorise the process 

mod_fit_period <- dta_selection %>% 
  ungroup() %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(hp_mod = map(
    data, 
    function(D) {
      MortalityLaw(
        x = D$age, 
        Dx = D$death_count,
        Ex = D$population_count,
        law = "HP4",
        how = "poissonL",
        fit.this.x = 0:90
      )
    }
    )
  )

params_fit_period <- mod_fit_period %>% 
  mutate(params = map(hp_mod, coefficients)) %>% 
  select(year, sex, params) %>% 
  mutate(params = map(params, enframe)) %>% 
  unnest()

params_fit_period %>%
  filter(name == "G") %>% 
  ggplot(aes(x = year, y = value, group = sex, colour = sex)) + 
  geom_line() 

# Seem to be various problems with fitting this 


mod_fit_period_hp <- dta_selection %>% 
  ungroup() %>% 
  group_by(year, sex) %>% 
  nest() %>% 
  mutate(hp_mod = map(
    data, 
    function(D) {
      MortalityLaw(
        x = D$age, 
        Dx = D$death_count,
        Ex = D$population_count,
        law = "HP",
        how = "poissonL",
        fit.this.x = 0:90
      )
    }
  )
  )

params_fit_period_hp <- mod_fit_period_hp %>% 
  mutate(params = map(hp_mod, coefficients)) %>% 
  select(year, sex, params) %>% 
  mutate(params = map(params, enframe)) %>% 
  unnest()

params_fit_period_hp %>%
  filter(is.finite(value)) %>%
  filter(value < 10E1) %>% 
  ggplot(aes(x = year, y = value, group = sex, colour = sex)) + 
  geom_line() +
  facet_wrap(~name, scale = "free_y") 
