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



dta_ratios <- dta_selection %>% 
  mutate(death_rate = death_count / population_count) %>% 
  select(year, age, sex, death_rate) %>% 
  spread(key=sex, value = death_rate) %>% 
  mutate(sex_ratio = male/ female) %>% 
  select(year, age , sex_ratio) %>% 
  filter(age <= 100)

fn <- function(x){
  smoothed_ratio <- rep(1, 101)
  for (i in 3:98){
    smoothed_ratio[i] <- prod(x$sex_ratio[(i-2):(i+2)]) ^ (1/5)
  }
  smoothed_ratio[-1] <- smoothed_ratio[-1]/smoothed_ratio[-101]
  out <- data.frame(x, smoothed_ratio = smoothed_ratio)
  return(out)
}


dta_ratios_fd <- dta_ratios %>% 
  group_by(year) %>% 
  do(fn(.)) 


# Explore curves by period ------------------------------------------------

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  filter(year %in% seq(1850, 2010, by = 10)) %>% 
  ggplot(., aes(x = age, y = lmr, group = sex, colour = sex)) + 
  geom_point() + facet_wrap(~year)

# Explore excess in lmr by period

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  filter(year %in% seq(1850, 2010, by = 10)) %>% 
  select(year, age, sex, lmr) %>% 
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>% 
  mutate(male_excess = diff_lmr > 0) %>% 
  ggplot(., aes(x = age, y = diff_lmr, colour = male_excess)) + 
  geom_point() + facet_wrap(~year) +
  coord_cartesian(ylim = c(-0.6, 0.6)) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(15, 18, 25, 35, 60), linetype = "dashed")

# Change in diff in lmr with additional year of age - by gender

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  filter(year %in% seq(1850, 2010, by = 10)) %>% 
  group_by(year, sex) %>% 
  arrange(age) %>% 
  mutate(lmr_last = lag(lmr, 1)) %>% 
  mutate(change_lmr = lmr - lmr_last) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(., aes(x = age, y = change_lmr, colour = sex, shape = sex)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(15, 18, 25, 35, 60, 65), linetype = "dashed")


# Bathtubs by cohort 
dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  filter(cohort %in% seq(1800, 1980, by = 10)) %>% 
  ggplot(., aes(x = age, y = lmr, colour = sex, shape = sex)) + 
  geom_point() + 
  facet_wrap(~cohort)

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  filter(cohort %in% seq(1800, 1980, by = 10)) %>%
  select(cohort, age, sex, lmr) %>% 
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>% 
  ggplot(., aes(x = age, y = diff_lmr)) + 
  geom_point() + 
  facet_wrap(~cohort) + 
  geom_hline(yintercept = 0)

# Interest is in cohorts who entered the labour market after 1950

# This means cohorts born from 1930


dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  filter(cohort %in% seq(1930, 1980, by = 2)) %>%
  select(cohort, age, sex, lmr) %>% 
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>% 
  group_by(cohort) %>% 
  mutate(max_diff = max(diff_lmr)) %>%
  mutate(age_max_diff = age[diff_lmr == max_diff]) %>% 
  mutate(age_selection = age >= age_max_diff) %>% 
  ungroup() %>% 
  ggplot(., aes(x = age, y = diff_lmr)) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(~cohort) + 
  geom_hline(yintercept = 0) +
  stat_smooth(data = . %>% filter(age_selection == T), method = "lm", se = F) +
  stat_smooth(data = . %>% filter(age_selection == F), method = "lm", se = F, colour = "red")

# 
# # For each cohort want to extract the following:
# # maximum of three age average rate of change 
# # gradient from maximum age to last observed age 
# 
# calculate_params <- function(DTA){
#   # DTA contains age and diff
#   DTA <- DTA %>% 
#     arrange(age) %>% 
#     mutate(sm_diff_lmr = mean(c(lag(diff_lmr), diff_lmr, lead(diff_lmr)))) 
#   
#   max_diff <- max(DTA$sm_diff_lmr, na.rm = T)
#   age_max_diff <- DTA$age[DTA$sm_diff_lmr == max_diff]
#   
#   max_age <- max(DTA$age, na.rm = T)
#   
#   DTA_SS <- DTA %>% filter(age >= age_max_diff)
#   
#   mdl <- lm(diff_lmr ~ age, data = DTA_SS)
#   fall <- coefficients(mdl)[2]
#   
#   out <- list(age_max_diff, max_diff, fall)
#   out
# }
# 
# 
# dta_selection %>% 
#   ungroup %>% 
#   mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
#   mutate(lmr = log(mr, 10)) %>% 
#   mutate(cohort = year - age) %>% 
#   select(cohort, age, sex, lmr) %>% 
#   spread(sex, lmr) %>% 
#   mutate(diff_lmr = male - female) %>% 
#   filter(cohort >= 1850) %>% 
#   group_by(cohort) %>% 
#   nest() %>% 
#   mutate(params = map(data, safely(calculate_params)))
# 
# 

# Want to look at change in smoothed derivative 

dta_ratios_fd %>% 
  ungroup() %>% 
  filter(year %in% seq(1850, 2010, by = 10)) %>% 
  filter(age >= 5, age <= 95) %>% 
  mutate(divergence = smoothed_ratio > 1) %>% 
  ggplot(., aes(x = age, y = smoothed_ratio, colour = divergence, shape = divergence)) + 
  geom_point() + 
  facet_wrap(~year) + 
  geom_hline(yintercept = 1)


# Now to look at this by cohort 

dta_ratios_fd %>% 
  ungroup() %>% 
  mutate(cohort = year - age) %>% 
  filter(cohort %in% seq(1800, 1980, by = 10)) %>% 
  mutate(divergence = smoothed_ratio > 1) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(., aes(x = age, y = smoothed_ratio, colour = divergence, shape = divergence)) + 
  geom_point() + 
  facet_wrap(~cohort) + 
  geom_hline(yintercept = 1)

dta_ratios_fd %>% 
  ungroup() %>% 
  mutate(cohort = year - age) %>% 
  filter(cohort %in% seq(1930, 1970, by = 2)) %>% 
  mutate(divergence = smoothed_ratio > 1) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(., aes(x = age, y = smoothed_ratio, colour = divergence, shape = divergence)) + 
  geom_point() + 
  facet_wrap(~cohort) + 
  geom_hline(yintercept = 1)


dta_ratios_fd %>% 
  ungroup() %>% 
  mutate(cohort = year - age) %>% 
  mutate(divergence = smoothed_ratio > 1) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(., aes(x = age, y = smoothed_ratio, colour = cohort, shape = divergence)) + 
  geom_jitter(alpha = 0.1) + 
  geom_hline(yintercept = 1) + 
  scale_color_distiller(type = "qual") + 
  coord_cartesian(ylim=c(0.8, 1.2))

dta_ratios_fd %>% 
  ungroup() %>% 
  mutate(divergence = smoothed_ratio > 1) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(., aes(x = age, y = smoothed_ratio, colour = year, shape = divergence)) + 
  geom_jitter(alpha = 0.1) + 
  geom_hline(yintercept = 1) + 
  scale_color_distiller(type = "qual") + 
  coord_cartesian(ylim=c(0.8, 1.2))

dta_ratios_fd %>% 
  ungroup() %>% 
  mutate(divergence = smoothed_ratio > 1) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(., aes(x = age, y = smoothed_ratio, colour = year, shape = divergence)) + 
  geom_jitter(alpha = 0.1) + 
  geom_hline(yintercept = 1) + 
  scale_color_distiller(type = "qual") + 
  coord_cartesian(ylim=c(0.8, 1.2)) 




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







