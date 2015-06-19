# Replication and updating of Rigby & Dorling 2007, Mortality in relation to sex in the affluent world


rm(list=ls())

require(readr)
require(readxl)
require(plyr)
require(stringr)
require(tidyr)
require(dplyr)

require(lattice)
require(latticeExtra)

require(ggplot2)




# Load data ---------------------------------------------------------------

dta <- read_csv("data/counts.csv")

original_countries <- read_excel(
  "support/replication_details.xlsx", 
  sheet = "original_country_selection"
  )

code_to_country_lookup <- read_excel(
  "support/replication_details.xlsx",
  sheet="code_to_country_lookup"
)


code_selection <- code_to_country_lookup %>% 
  filter(in_original_selection == 1) %>% 
  .$code


dta_selection <- dta %>% 
  filter(country %in% code_selection & sex !="Total") %>% 
  group_by(year, age, sex) %>% 
  summarise(
    population_count = sum(population_count), 
    death_count = sum(death_count)
    ) %>% 
  filter(year >= 1850 & year <=2000)


dta_ratios <- dta_selection %>% 
  mutate(death_rate = death_count / population_count) %>% 
  select(year, age, sex, death_rate) %>% 
  spread(key=sex, value = death_rate) %>% 
  mutate(sex_ratio = male/ female) %>% 
  select(year, age , sex_ratio) %>% 
  filter(age <= 100)

contourplot(
    sex_ratio ~ year * age, 
    data=dta_ratios , 
    region=T, 
    par.strip.text=list(cex=1.4, fontface="bold"),
    ylab=list(label="Age in years", cex=1.4),
    xlab=list(label="Year", cex=1.4),
    cex=1.4,
    at=seq(from=0.9, to=3.4, by=0.1),
    col.regions=colorRampPalette(brewer.pal(6, "Reds"))(200),
    main=NULL,
    labels=list(cex=1.2),
    col="blue",
    scales=list(
      x=list(at = seq(from = 1850, to = 2000, by = 10)),
      y=list(at = seq(from = 0, to = 100, by = 10))
    )
  )





