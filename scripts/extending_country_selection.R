# Replication and updating of Rigby & Dorling 2007, Mortality in relation to sex in the affluent world

# Updating using new country selection but ensuring same countries are not included twice


rm(list=ls())

# Prerequisite packages ---------------------------------------------------


require(readr)
require(readxl)
require(plyr)
require(stringr)
require(tidyr)
require(dplyr)

require(lattice)
require(latticeExtra)

require(grid)
require(ggplot2)




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


table_1 <- dta %>% 
  filter(country %in% full_country_lookup) %>%
  filter(sex == "total") %>% 
  group_by(country) %>% 
  summarise(
    min_year = min(year), 
    max_year = max(year), 
    pop_2010 = ifelse(
      max_year >= 2010, 
      sum(population_count[year == 2010]) / 1E6,
      sum(population_count[year == max_year] / 1E6)
    )
  ) %>% 
  arrange(min_year)

write.csv(table_1, "clipboard")
    
# Derived data ------------------------------------------------------------



dta_selection <- dta %>% 
  filter(sex !="total" & country %in% full_country_lookup) %>% 
  group_by(year, age, sex) %>% 
  summarise(
    population_count = sum(population_count), 
    death_count = sum(death_count)
    ) %>% 
  filter(year >= 1850 & year <=2010)





# Original figure 1 -------------------------------------------------------


dta_ratios <- dta_selection %>% 
  mutate(death_rate = death_count / population_count) %>% 
  select(year, age, sex, death_rate) %>% 
  spread(key=sex, value = death_rate) %>% 
  mutate(sex_ratio = male/ female) %>% 
  select(year, age , sex_ratio) %>% 
  filter(age <= 100) %>% 
  mutate(sex_ratio = ifelse(
    sex_ratio > 3.4, 3.399,
    ifelse(
      sex_ratio < 0.9, 0.901,
      sex_ratio
    )
  )
  )


p <- contourplot(
    sex_ratio ~ year * age, 
    data=dta_ratios , 
    region=T, 
    par.strip.text=list(cex=1.4, fontface="bold"),
    ylab=list(label="Age in years", cex=1.4),
    xlab=list(label="Year", cex=1.4),
    cex=1.4,
    at=seq(from=0.9, to=3.4, by=0.1),
    col.regions=colorRampPalette(rev(brewer.pal(6, "Spectral")))(200),
    main=NULL,
    labels=FALSE,
    col="black",
    scales=list(
      x=list(at = seq(from = 1850, to = 2010, by = 10)),
      y=list(at = seq(from = 0, to = 100, by = 10))
    )
  )

png(
  "figures/all_countries/figure1_sex_ratio.png",
  res=300,
  height=20, width=20, units="cm"
)
print(p)
dev.off()

# Original Figure 2 -------------------------------------------------------


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
  do(fn(.)) %>% 
  mutate(smoothed_ratio = ifelse(
    smoothed_ratio > 1.10, 1.099,
    ifelse(
      smoothed_ratio < 0.95, 0.951,
      smoothed_ratio
    )
  )
  )



p <- contourplot(
  smoothed_ratio ~ year * age, 
  data=dta_ratios_fd , 
  region=T, 
  par.strip.text=list(cex=1.4, fontface="bold"),
  ylab=list(label="Age in years", cex=1.4),
  xlab=list(label="Year", cex=1.4),
  cex=1.4,
  at=seq(from=0.95, to=1.10, by=0.01),
  col.regions=colorRampPalette(rev(brewer.pal(6, "Spectral")))(200),
  main=NULL,
  labels=FALSE,
  col="black",
  scales=list(
    x=list(at = seq(from = 1850, to = 2010, by = 10)),
    y=list(at = seq(from = 0, to = 100, by = 10))
  )
)

png(
  "figures/all_countries/figure2_smoothed_first_derivative.png",
  res=300,
  height=20, width=20, units="cm"
)
print(p)
dev.off()




# Original Table 2 --------------------------------------------------------

dta_selection %>% 
  filter(sex !="total" & year %in% 
           c(1860, 1870, 1890, 1900,1910, 1919, 1930, 1940, 
             1950, 1960, 1970, 1980, 1990, 2000, 2010)) %>% 
  group_by(year, sex) %>% 
  summarise(
    population_count = sum(population_count),
    death_count = sum(death_count)          
            ) %>% 
  mutate(
    rate_per_thousand = 1000 * death_count / population_count
    ) %>% 
  select(-population_count, -death_count) %>% 
  spread(key=sex, value = rate_per_thousand) %>% 
  write.table(., file="clipboard")



# Original figure 3 -------------------------------------------------------



dta_selection <- dta %>% 
  filter( sex !="total") %>% 
  group_by(year, age, sex) %>% 
  summarise(
    population_count = sum(population_count), 
    death_count = sum(death_count)
  ) %>% 
  filter(year >= 1841 & year <=2010)


dta_selection %>% 
  filter(age %in% c(10, 30)) %>% 
  mutate(
    mortality_rate = 1000 * death_count / population_count,
    grp = paste(sex, "aged", age)
    ) %>% 
  ggplot(.) + 
  geom_line(aes(x=year, y= mortality_rate, group=grp, colour=grp)) + 
  scale_y_log10(limits=c(0.1, 100.0), breaks=c(0.1, 1.0, 10.0, 100.0)) + 
  theme_minimal() + 
  scale_x_continuous(
    limits = c(1841, 2010), 
    breaks=seq(from = 1841, to = 2010, by = 10)
    ) + 
  labs(y="Mortality/1000/year", x="Year") + 
  scale_colour_discrete(
    guide=guide_legend(title=NULL)
  ) + 
  theme(
    axis.text.x = element_text(angle= 90),
    legend.justification = c(0, 1), 
    legend.position=c(0.7,0.9)
    )


ggsave(filename ="figures/all_countries/figure3_mortper1000peryear.png",
       width=15, height=15, units="cm", dpi = 300
)



# Replicate tables of mortality ratios by birth year ----------------------

dta_selection2 <- dta_selection
dta_selection2 <- dta_selection2 %>% 
  mutate(birth_year = year - age) %>% 
  filter(birth_year %in% seq(from = 1850, to =  2010, by = 5)) %>% 
  filter(age < 100)
  
age_groups <- c(
    "0-4",    "5-9",
    "10-14",  "15-19",
    "20-24",  "25-29",
    "30-34",    "35-39",
    "40-44",    "45-49",
    "50-54",    "55-59",
    "60-64",    "65-69",
    "70-74",    "75-79",
    "80-80",    "85-89",
    "90-94",    "95-99"
)
age_lookup <- data.frame(
  age = 0:99,
  age_group = age_groups[1 + (0:99 %/% 5)]
)

dta_selection2$age_group <- mapvalues(
  dta_selection2$age, 
  from = age_lookup$age, 
  to = as.character(age_lookup$age_group) 
)

dta_selection2$age_group <- factor(
  dta_selection2$age_group,
  levels = age_groups
)

rate_table <- dta_selection2 %>% 
  group_by(birth_year, age_group, sex) %>% 
  summarise(
    population_count = sum(population_count),
    death_count = sum(death_count)
            ) %>% 
  mutate(mortality_rate = death_count / population_count) %>% 
  select(-population_count, -death_count) %>% 
  spread(key=sex, value = mortality_rate) %>% 
  mutate(rate_ratio = male / female) %>% 
  select(-female, -male) %>% 
  spread(key=age_group, value = rate_ratio)

