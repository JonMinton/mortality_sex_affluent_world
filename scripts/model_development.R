# Modelling 

# Replication and updating of Rigby & Dorling 2007, Mortality in relation to sex in the affluent world

# Updating using new country selection but ensuring same countries are not included twice


rm(list=ls())

# Prerequisite packages ---------------------------------------------------


pacman::p_load(
  MASS,
  readxl,
  tidyverse,
  grid, ggplot2,
  lattice, 
  latticeExtra
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



dta_selection %>% 
  ungroup() %>%  
  mutate(
    death_count = round(death_count + 1, 0), 
    population_count = round(population_count + 1, 0)
  ) -> dta_for_nb

formulae <- c(
  sex    = "death_count ~ sex + offset(log(population_count))",
  s_a    = "death_count ~ I(age - min(age)) + sex + offset(log(population_count))",
  s_p    = "death_count ~ I(year - min(year)) + sex + offset(log(population_count))",
  s_c    = "death_count ~ I(cohort - min(cohort)) + sex + offset(log(population_count))", 
  s_a_p  = "death_count ~ I(age - min(age)) + I(year - min(year)) + sex + offset(log(population_count))",
  s_A    = "death_count ~ factor(age) + sex + offset(log(population_count))",
  s_P    = "death_count ~ factor(year) + sex + offset(log(population_count))",
  s_C    = "death_count ~ factor(cohort) + sex + offset(log(population_count))"
)

glm.nb(
  death_count ~ sex + offset(log(population_count)),
  data = dta_for_nb
)

glm.nb(
  death_count ~ sex,
  offset = log(population_count),
  data = dta_for_nb
)



# Segmentation analysis ---------------------------------------------------


# Idea - look at gradient and intercept for log mort for males 
# and females and see how this has changed over time 

# do this for the age range 60-95 to start with 

dta_selection %>% 
  ungroup %>% 
  mutate(mr = death_count / population_count) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  filter(age >= 60, age <= 95) %>% 
  group_by(sex, year) %>% 
  nest() %>% 
  mutate(
    mdl = map(data, ~ lm(lmr ~ I(age - min(age)), .))
  ) %>% 
  mutate(
    coef = map(mdl, coef)
  ) %>% 
  mutate(
    intercept = map_dbl(coef, ~ .[[1]]),
    gradient = map_dbl(coef, ~ .[[2]])
  ) -> mdl_old 

mdl_old %>% 
  ggplot(., aes(x = intercept, y = gradient, colour = sex)) +
  geom_point(aes(alpha = year, shape = sex)) + 
  geom_label(aes(label = year), data = . %>% filter(year %in% seq(1850, 2000, by = 25)))

# Let's try a gganimate version of this 

#devtools::install_github("dgrtwo/gganimate")

library(gganimate)


dta_selection %>% 
  ungroup %>% 
  mutate(mr = death_count / population_count) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  filter(age >= 40, age <= 60) %>% 
  group_by(sex, year) %>% 
  nest() %>% 
  mutate(
    mdl = map(data, ~ lm(lmr ~ I(age - min(age)), .))
  ) %>% 
  mutate(
    coef = map(mdl, coef)
  ) %>% 
  mutate(
    intercept = map_dbl(coef, ~ .[[1]]),
    gradient = map_dbl(coef, ~ .[[2]])
  ) -> mdl_middle

p <- ggplot(
  data = mdl_middle, 
  aes(x = intercept, y = gradient, colour = sex)
) + 
  geom_point(aes(frame = year, cumulative = TRUE))

gganimate(p, interval = 0.01, title_frame = T)







