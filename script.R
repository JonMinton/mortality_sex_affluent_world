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
      x=list(at = seq(from = 1850, to = 2000, by = 10)),
      y=list(at = seq(from = 0, to = 100, by = 10))
    )
  )

png(
  "figures/original/figure1_sex_ratio.png",
  res=300,
  height=20, width=20, units="cm"
)
print(p)
dev.off()

# Original figure 2

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
    x=list(at = seq(from = 1850, to = 2000, by = 10)),
    y=list(at = seq(from = 0, to = 100, by = 10))
  )
)

png(
  "figures/original/figure2_smoothed_first_derivative.png",
  res=300,
  height=20, width=20, units="cm"
)
print(p)
dev.off()



