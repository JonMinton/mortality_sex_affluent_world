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


# Curves by year, fd, 1850 to 1910

avg_smooth_old <- dta_ratios_fd %>% 
  filter(year <= 1910) %>% 
  filter(age >= 5, age <= 95) %>% 
  group_by(age) %>% 
  summarise(
    average_smooth = mean(smoothed_ratio),
    var_smooth = var(smoothed_ratio)
    )


dta_ratios_fd %>% 
  filter(year <= 1910) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(
    ., 
    aes(x = age, y = smoothed_ratio, group = year)     
    ) + 
  geom_line(alpha = 0.05) + 
  geom_hline(yintercept = 1) +
  scale_x_continuous(breaks = seq(5, 95, by = 5)) + 
  scale_y_continuous(limits = c(0.9, 1.2), breaks = seq(0.9, 1.2, by = 0.01)) + 
  ggtitle("Change in sex ratio by age, 1850-1910") + 
  geom_line(
    aes(x = age, y = average_smooth),
    inherit.aes = F,
    data = avg_smooth_old,
    size = 1.2, colour = "blue"
  )

# Now the same over the period 1950-2010

avg_smooth_new <- dta_ratios_fd %>% 
  filter(year >= 1950, year <= 2010) %>% 
  filter(age >= 5, age <= 95) %>% 
  group_by(age) %>% 
  summarise(
    average_smooth = mean(smoothed_ratio),
    var_smooth = var(smoothed_ratio)
  )


dta_ratios_fd %>% 
  filter(year >= 1950, year <= 2010) %>% 
  filter(age >= 5, age <= 95) %>% 
  ggplot(
    ., 
    aes(x = age, y = smoothed_ratio, group = year)     
  ) + 
  geom_line(alpha = 0.05) + 
  geom_hline(yintercept = 1) +
  scale_x_continuous(breaks = seq(5, 95, by = 5)) + 
  scale_y_continuous(limits = c(0.9, 1.2), breaks = seq(0.9, 1.2, by = 0.01)) + 
  ggtitle("Change in sex ratio by age, 1950-2010") + 
  geom_line(
    aes(x = age, y = average_smooth),
    inherit.aes = F,
    data = avg_smooth_new,
    size = 1.2, colour = "red"
  )


# Two averages on same figure

average_two_epochs <- bind_rows(
  avg_smooth_old %>% mutate(epoch = "old"),
  avg_smooth_new %>% mutate(epoch = "new")
) 

ggplot(
  average_two_epochs,
  aes(x = age, y = average_smooth, group = epoch, colour = epoch)
  ) + 
  geom_line(size = 1.2) + 
  geom_hline(yintercept = 1) +
  scale_x_continuous(breaks = seq(5, 95, by = 5)) + 
  scale_y_continuous(limits = c(0.9, 1.2), breaks = seq(0.9, 1.2, by = 0.01)) + 
  ggtitle("Average change in sex ratios over two epochs")  


# Variance of change with age by year within the two epochs

ggplot(
  average_two_epochs,
  aes(x = age, y = var_smooth, group = epoch, colour = epoch)
) + 
  geom_line(size = 1.2) + 
  scale_x_continuous(breaks = seq(5, 95, by = 5)) + 
#  scale_y_continuous(limits = c(0.9, 1.2), breaks = seq(0.9, 1.2, by = 0.01)) + 
  ggtitle("Variance of age-change in sex ratios over two epochs")  


# Explore curves by period ------------------------------------------------

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(decade = cut(
    year, seq(1850, 2010, by = 10), 
    include.lowest = T,
    labels = paste0(seq(1850, 2000, by = 10), "s")
    )
  ) %>% 
  ggplot(., aes(x = age, y = lmr, group = sex, colour = sex)) + 
  geom_point(alpha = 0.1) + 
  facet_wrap(~decade)

# Explore excess in lmr by period

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  select(year, age, sex, lmr) %>% 
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>% 
  mutate(male_excess = diff_lmr > 0) %>% 
  mutate(decade = cut(
    year, seq(1850, 2010, by = 10), 
    include.lowest = T,
    labels = paste0(seq(1850, 2000, by = 10), "s")
  )
  ) %>% 
  ggplot(., aes(x = age, y = diff_lmr, colour = male_excess)) + 
  geom_point(alpha = 0.1) + facet_wrap(~decade) +
  coord_cartesian(ylim = c(-0.6, 0.6)) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(15, 18, 25, 35, 60), linetype = "dashed")



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

# line version of above grouped rather than faceted

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  select(cohort, age, sex, lmr) %>%
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>%
  filter(age <= 90) %>% 
  ggplot(., aes(x = age, y = diff_lmr, group = cohort, alpha = cohort )) + 
  geom_line() + 
  geom_hline(yintercept = 0) + 
  scale_alpha_continuous(limits = c(1800, 1980), range = c(0, 1))


dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  select(cohort, age, sex, lmr) %>%
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>%
  filter(age <= 90) %>% 
  filter(cohort >= 1930) %>% 
  ggplot(., aes(x = age, y = diff_lmr, group = cohort, color = cohort)) + 
  geom_line() + 
  geom_hline(yintercept = 0) + 
#  scale_alpha_continuous(limits = c(1930, 1980), range = c(0, 1)) + 
  scale_color_distiller(palette = "Paired")

# Let's try but faceted by decade

dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  select(cohort, age, sex, lmr) %>%
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>%
  mutate(cohort_decade = 10 * cohort %/% 10) %>% 
  mutate(year_within_cohort = cohort - cohort_decade) %>% 
  filter(age <= 90) %>% 
  filter(cohort >= 1850, cohort <= 1999) %>% 
  ggplot(., aes(x = age, y = diff_lmr, group = cohort, color = factor(year_within_cohort))) + 
  geom_line() + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~cohort_decade) + 
  #  scale_alpha_continuous(limits = c(1930, 1980), range = c(0, 1)) + 
  scale_color_brewer(palette = "Paired", name = "Year\nwithin\ncohort") +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 5)) + 
  scale_y_continuous(breaks = seq(-0.10, 0.60, by = 0.05)) +
  coord_cartesian(ylim = c(-0.1, 0.6)) + 
  labs(title = "Difference in log mortality by cohort", x = "Age in years", y = "Difference in log mortality")

  

# Same sort of thing, but period not cohort 


dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  select(year, age, sex, lmr) %>%
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>%
  mutate(decade = 10 * year %/% 10) %>% 
  mutate(year_within_decade = year - decade) %>% 
  filter(age <= 90) %>% 
  filter(year < 2010) %>% 
  ggplot(., aes(x = age, y = diff_lmr, group = year, color = factor(year_within_decade))) + 
  geom_line() + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~decade) + 
  #  scale_alpha_continuous(limits = c(1930, 1980), range = c(0, 1)) + 
  scale_color_brewer(palette = "Paired", name = "Year\nwithin\ndecade") +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 5)) + 
  scale_y_continuous(breaks = seq(-0.10, 0.60, by = 0.05)) +
  coord_cartesian(ylim = c(-0.1, 0.6)) + 
  labs(title = "Difference in log mortality by period", x = "Age in years", y = "Difference in log mortality")


# Interest in gradient of change between ages 12 and 23 years 

# First by period, then by cohort 

# Note: dots look better here than lines
# By period
dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  select(year,age, sex, lmr) %>% 
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>% 
  group_by(year) %>% 
  arrange(age) %>% 
  mutate(grad_lmr = diff_lmr - lag(diff_lmr)) %>% 
  filter(age >= 3, age <= 50) %>% 
  mutate(decade = 10 * year %/% 10) %>% 
  mutate(year_within_decade = year - decade) %>% 
  filter(year < 2010) %>% 
  ggplot(., aes(x = age, y = grad_lmr, group = year, color = factor(year_within_decade))) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~decade) + 
  scale_color_brewer(palette = "Paired", name = "Year\nwithin\ndecade") +
  scale_x_continuous(limits = c(3, 50), breaks = seq(4, 50, by = 2)) + 
  scale_y_continuous(breaks = seq(-0.05, 0.15, by = 0.01)) +
  coord_cartesian(ylim = c(-0.05, 0.15)) + 
  labs(title = "Change in Difference in log mortality by period", x = "Age in years", y = "Difference in log mortality") 



# By cohort
dta_selection %>% 
  ungroup %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>%
  mutate(cohort = year - age) %>% 
  select(cohort,age, sex, lmr) %>% 
  spread(sex, lmr) %>% 
  mutate(diff_lmr = male - female) %>% 
  group_by(cohort) %>% 
  arrange(age) %>% 
  mutate(grad_lmr = diff_lmr - lag(diff_lmr)) %>% 
  filter(age >= 3, age <= 50) %>% 
  mutate(cohort_decade = 10 * cohort %/% 10) %>% 
  mutate(year_within_cohort = cohort - cohort_decade) %>% 
  filter(cohort < 2000) %>% 
  ggplot(., aes(x = age, y = grad_lmr, group = cohort, color = factor(year_within_cohort))) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~cohort_decade) + 
  scale_color_brewer(palette = "Paired", name = "Year\nwithin\ncohort") +
  scale_x_continuous(limits = c(3, 50), breaks = seq(4, 50, by = 2)) + 
  scale_y_continuous(breaks = seq(-0.05, 0.15, by = 0.01)) +
  coord_cartesian(ylim = c(-0.05, 0.15)) + 
  labs(title = "Change in Difference in log mortality by cohort", x = "Age in years", y = "Difference in log mortality") 



# Piecewise models --------------------------------------------------------


# The purpose of this section is to construct the model segmented by 
# different ages in the life course 
# First this will be by year, then by cohort 
peek <- function(x) print(sample_n(x, 10))

mdls_by_period <- dta_selection %>% 
  ungroup()  %T>% peek %>% 
  mutate(
    age_group = cut(
      age, 
      breaks = c(0, 3, 11, 23, 35, 59, 90, Inf),
      labels = c(
        "0-3",
        "4-11",
        "12-23",
        "24-35",
        "36-59",
        "60-90",
        "older"
      ),
      include.lowest = T
    )
  ) %T>% peek %>%
  filter(age_group != "older") %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %T>% peek %>% 
  group_by(year, age_group, sex) %>% 
  nest() %T>% peek %>% 
  mutate(mdl = map(data, ~ lm(lmr ~ I(age - min(age)), data = .))) %>%
  mutate(
    intercept = map_dbl(mdl, ~ coefficients(.)[[1]]),
    gradient = map_dbl(mdl, ~ coefficients(.)[[2]]),
    fit = map_dbl(mdl, ~ summary(.)[["adj.r.squared"]])
  ) 

# Now to visualise 

ggplot(
  mdls_by_period,
  aes(y = intercept, x = gradient, colour = fit)
  ) + 
  geom_path() + 
  facet_grid(sex ~ age_group, scale = "free_x") + 
  scale_color_distiller(palette = "Paired") + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_text(
    aes(label = year), 
    colour = "black", 
    fontface = "bold",
    check_overlap = T,
    data = . %>% filter(year %in% seq(1850, 2000, by = 25))
    )

# And now to visualise change 

mdls_by_period %>% 
  group_by(year, age_group) %>% 
  mutate(
    dif_intercept = intercept[sex == "male"] - intercept[sex == "female"],
    dif_gradient = gradient[sex == "male"] - gradient[sex == "female"],
    dif_fit = fit[sex == "male"] - fit[sex == "female"]
    ) %>%  
  ggplot(
    aes(y = dif_intercept, x = dif_gradient, colour = dif_fit)
  ) + 
  geom_path() + 
  facet_wrap( ~ age_group, scale = "free", nrow = 1) + 
  scale_color_gradient2(
    low = "red", mid = "grey", high = "blue",
    limits = c(-1, 1)
    ) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_text(
    aes(label = year), 
    colour = "black", 
    fontface = "bold",
    check_overlap = T,
    data = . %>% filter(year %in% seq(1850, 2000, by = 50))
  )


# Now to do the same kind of exercise by cohort 


age_group_selector <- dta_selection %>% 
  ungroup()  %T>% peek %>% 
  mutate(
    age_group = cut(
      age, 
      breaks = c(0, 3, 11, 23, 35, 59, 90, Inf),
      labels = c(
        "0-3",
        "4-11",
        "12-23",
        "24-35",
        "36-59",
        "60-90",
        "older"
      ),
      include.lowest = T
    )
  ) %T>% peek %>%
  filter(age_group != "older") %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %T>% peek %>% 
  separate(
    age_group, into = c("first_age", "last_age"), 
    sep = "-", remove = F
  ) %>%
  mutate(first_age = as.numeric(first_age), 
         last_age = as.numeric(last_age)
  ) %>% 
  mutate(cohort = year - age) %>%
  select(cohort, age, sex, age_group, first_age, last_age) %>% 
  group_by(cohort, age_group) %>% 
  arrange(age) %>% 
  mutate(
    min_found_age = min(age),
    max_found_age = max(age)
  ) %>% 
  mutate(
    age_group_ok = first_age == min_found_age & last_age == max_found_age
  ) %>% 
  select(cohort, age_group, age_group_ok) %>% 
  distinct() %>% ungroup()


dta_for_cohort_models <- dta_selection %>% 
  ungroup()  %T>% peek %>% 
  mutate(
    age_group = cut(
      age, 
      breaks = c(0, 3, 11, 23, 35, 59, 90, Inf),
      labels = c(
        "0-3",
        "4-11",
        "12-23",
        "24-35",
        "36-59",
        "60-90",
        "older"
      ),
      include.lowest = T
    )
  ) %T>% peek %>%
  filter(age_group != "older") %>% 
  mutate(mr = (death_count + 0.5) / (population_count + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  mutate(cohort = year - age) %>% 
  left_join(age_group_selector) %>% 
  filter(age_group_ok) %>% 
  select(cohort, age, age_group, sex, death_count, population_count) %>% 
  arrange(cohort, age, sex) %>% 
  mutate(
    mr = (death_count + 0.5) / (population_count + 0.5),
    lmr = log(mr, 10)
  )

mdls_by_cohort <- dta_for_cohort_models %>% 
  group_by(cohort, age_group, sex) %>% 
  nest() %T>% peek %>% 
  mutate(mdl = map(data, ~ lm(lmr ~ I(age - min(age)), data = .))) %>%
  mutate(
    intercept = map_dbl(mdl, ~ coefficients(.)[[1]]),
    gradient = map_dbl(mdl, ~ coefficients(.)[[2]]),
    fit = map_dbl(mdl, ~ summary(.)[["adj.r.squared"]])
  ) 


ggplot(
  mdls_by_cohort,
  aes(y = intercept, x = gradient, colour = fit)
) + 
  geom_path() + 
  facet_grid(sex ~ age_group, scale = "free") + 
  scale_color_distiller(palette = "Paired") + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_text(
    aes(label = cohort), 
    colour = "black", 
    fontface = "bold",
    check_overlap = T,
    data = . %>% filter(cohort %in% seq(1750, 1975, by = 25))
  )


mdls_by_cohort %>% 
  group_by(cohort, age_group) %>% 
  mutate(
    dif_intercept = intercept[sex == "male"] - intercept[sex == "female"],
    dif_gradient = gradient[sex == "male"] - gradient[sex == "female"],
    dif_fit = fit[sex == "male"] - fit[sex == "female"]
  ) %>%  
  ggplot(
    aes(y = dif_intercept, x = dif_gradient, colour = dif_fit)
  ) + 
  geom_path() + 
  facet_wrap( ~ age_group, scale = "free", nrow = 1) + 
  scale_color_gradient2(
    low = "red", mid = "grey", high = "blue",
    limits = c(-1, 1)
  ) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_text(
    aes(label = cohort), 
    colour = "black", 
    fontface = "bold",
    check_overlap = T,
    data = . %>% filter(cohort %in% seq(1750, 1975, by = 25))
  )

mdls_by_cohort %>% 
  group_by(cohort, age_group) %>% 
  mutate(
    dif_intercept = intercept[sex == "male"] - intercept[sex == "female"],
    dif_gradient = gradient[sex == "male"] - gradient[sex == "female"],
    dif_fit = fit[sex == "male"] - fit[sex == "female"]
  ) %>%  
  ggplot(
    aes(y = dif_intercept, x = dif_gradient, colour = dif_fit)
  ) + 
  geom_path() + 
  facet_wrap( ~ age_group, nrow = 1) + 
  scale_color_gradient2(
    low = "red", mid = "grey", high = "blue",
    limits = c(-1, 1)
  ) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_text(
    aes(label = cohort), 
    colour = "black", 
    fontface = "bold",
    check_overlap = T,
    data = . %>% filter(cohort %in% seq(1750, 1975, by = 25))
  )



# Older material ----------------------------------------------------------




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

mdl_middle %>% 
  ggplot(., aes(x = intercept, y = gradient, colour = sex)) +
  geom_path(aes(alpha = year)) + 
  geom_label(aes(label = year), data = . %>% filter(year %in% seq(1850, 2000, by = 25)))


mdl_both <- bind_rows(
  mdl_middle %>% 
    select(sex, year, intercept, gradient) %>% 
    mutate(age_range = "middle"),

  mdl_old %>% 
    select(sex, year, intercept, gradient) %>% 
    mutate(age_range = "old")
)

mdl_both %>% 
  group_by(year, age_range) %>% 
  mutate(dif_intercept = intercept[sex == "male"] - intercept[sex == "female"]) %>% 
  mutate(dif_gradient = gradient[sex == "male"] - gradient[sex == "female"])  %>% 
  select(year, age_range, dif_intercept, dif_gradient) %>% 
  ungroup() %>% 
  distinct() %>% 
  ggplot(., aes(x = dif_intercept, y = dif_gradient, colour = age_range, group= age_range)) + 
  geom_path(aes(alpha = year)) + 
  geom_label(aes(label = year), data = . %>% filter(year %in% seq(1850, 2000, by = 25))) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)



  


