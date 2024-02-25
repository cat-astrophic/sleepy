# This script validates the mechanism / theory section of the paper

# Loading libraries

library(modelsummary)
library(tidycensus)
library(stargazer)
library(geomander)
library(tinytiger)
library(sandwich)
library(leaflet)
library(ggplot2)
library(xgboost)
library(lmtest)
library(tigris)
library(dplyr)
library(terra)
library(sf)

# Project directory

direc <- 'D:/sleepy/'

# Begin to get data for desired states

all.fips <- counties(state = NULL, year = 2021)

keep.fips <- c('01', '12', '13', '21', '37', '47')

all.fips <- all.fips %>% filter(STATEFP %in% keep.fips)

# Further subset for desired counties

msas <- c('12073', '13215', '47065', '47093')

all.fips$centroid <- st_centroid(all.fips$geometry)

cdf <- all.fips %>% filter(GEOID %in% msas)

d1 <- c()
d2 <- c()
d3 <- c()
d4 <- c()

for (i in 1:nrow(all.fips)) {
  
  print(i)
  
  d1 <- c(d1, st_distance(cdf$centroid[1], all.fips$centroid[i]))
  d2 <- c(d2, st_distance(cdf$centroid[2], all.fips$centroid[i]))
  d3 <- c(d3, st_distance(cdf$centroid[3], all.fips$centroid[i]))
  d4 <- c(d4, st_distance(cdf$centroid[4], all.fips$centroid[i]))
  
}

f1 <- as.integer(d1 < 100000)
f2 <- as.integer(d2 < 100000)
f3 <- as.integer(d3 < 100000)
f4 <- as.integer(d4 < 100000)

all.fips$f1 <- f1
all.fips$f2 <- f2
all.fips$f3 <- f3
all.fips$f4 <- f4

all.fips$d1 <- d1
all.fips$d2 <- d2
all.fips$d3 <- d3
all.fips$d4 <- d4

all.fips$Flag <- as.integer(f1 + f2 + f3 + f4 > 0)

# Remove non-msa cities

all.fips <- all.fips %>% filter(Flag == 1)

# Defining the time zone variable

fancy.pants.tn <- c('47065', '47129', '47013', '47123', '47105', '47089', '47173', '47009',
                    '47107', '47145', '47093', '47029', '47063', '47025', '47143', '47057',
                    '47011', '47121', '47155', '47139', '47001', '47151', '47063', '47067')

fancy.pants.fl <- c('12129', '12037', '12065', '12123', '12079', '12077', '12039', '12073')

fancy.pants.ky <- all.fips[which(all.fips$STATEFP == '21'),]$GEOID
fancy.pants.ga <- all.fips[which(all.fips$STATEFP == '13'),]$GEOID
fancy.pants.nc <- all.fips[which(all.fips$STATEFP == '37'),]$GEOID

fancy.pants <- c(fancy.pants.tn, fancy.pants.ky, fancy.pants.ga, fancy.pants.nc, fancy.pants.fl)

all.fips$TZ <- as.integer(all.fips$GEOID %in% fancy.pants)

# Prepping al.fips for transition to panel data

all.fips <- rbind(all.fips, all.fips, all.fips, all.fips, all.fips, all.fips, all.fips, all.fips, all.fips, all.fips)

# Add house price / rent data to all.fips

values <- c()
years <- c()

for (y in 2015:2021) {
  
  data <- get_acs(geography = 'county', state = c('AL', 'GA', 'FL', 'TN', 'KY', 'NC'), year = y, variables = c('DP04_0089'))
  
  for (i in 1:(nrow(all.fips)/7)) {
    
    tmp <- data %>% filter(GEOID == all.fips$GEOID[i])
    
    values <- c(values, tmp$estimate[1])
    years <- c(years, y)
    
  }
  
}

all.fips$Value <- values
all.fips$Year <- years

# Add population (housing demand) data to all.fips

pops <- c()

for (y in 2015:2021) {
  
  data <- get_acs(geography = 'county', state = c('AL', 'GA', 'FL', 'TN', 'KY', 'NC'), year = y, variables = c('B01003_001'))
  
  for (i in 1:(nrow(all.fips)/7)) {
    
    tmp <- data %>% filter(GEOID == all.fips$GEOID[i])
    
    pops <- c(pops, tmp$estimate[1])
    
  }
  
}

all.fips$Population <- pops

# Add housing units (supply) data to all.fips

units <- c()

for (y in 2015:2021) {
  
  data <- get_acs(geography = 'county', state = c('AL', 'GA', 'FL', 'TN', 'KY', 'NC'), year = y, variables = c('DP04_0001'))
  
  for (i in 1:(nrow(all.fips)/7)) {
    
    tmp <- data %>% filter(GEOID == all.fips$GEOID[i])
    
    units <- c(units, tmp$estimate[1])
    
  }
  
}

all.fips$Units <- units

# Add rent data to all.fips

rent <- c()

for (y in 2015:2021) {
  
  data <- get_acs(geography = 'county', state = c('AL', 'GA', 'FL', 'TN', 'KY', 'NC'), year = y, variables = c('DP04_0134'))
  
  for (i in 1:(nrow(all.fips)/7)) {
    
    tmp <- data %>% filter(GEOID == all.fips$GEOID[i])
    
    rent <- c(rent, tmp$estimate[1])
    
  }
  
}

all.fips$Rent <- rent

# Add income data to all.fips

inc <- c()

for (y in 2015:2021) {
  
  data <- get_acs(geography = 'county', state = c('AL', 'GA', 'FL', 'TN', 'KY', 'NC'), year = y, variables = c('DP03_0062'))
  
  for (i in 1:(nrow(all.fips)/7)) {
    
    tmp <- data %>% filter(GEOID == all.fips$GEOID[i])
    
    inc <- c(inc, tmp$estimate[1])
    
  }
  
}

all.fips$Income <- inc

# Remove counties that are too close to other metro areas

too.close.atl <- c()

atl <- st_sfc(st_point(x = c(-84.3877, 33.7488)))
st_crs(atl) <- 4269

for (i in 1:nrow(all.fips)) {
  
  if (st_distance(atl, all.fips$centroid[i]) < st_distance(cdf$centroid[1], all.fips$centroid[i])) {
    
    if (as.vector(st_distance(atl, all.fips$centroid[i])) < 100000) {
      
      too.close.atl <- c(too.close.atl, all.fips$GEOID[i])
      
    }
    
  }
  
  if (st_distance(atl, all.fips$centroid[i]) < st_distance(cdf$centroid[2], all.fips$centroid[i])) {
    
    if (as.vector(st_distance(atl, all.fips$centroid[i])) < 100000) {
      
      too.close.atl <- c(too.close.atl, all.fips$GEOID[i])
      
    }
    
  }
  
}

too.close.ash <- c()

ash <- st_sfc(st_point(x = c(-82.5515, 35.5951)))
st_crs(ash) <- 4269

for (i in 1:nrow(all.fips)) {
  
  if (st_distance(ash, all.fips$centroid[i]) < st_distance(cdf$centroid[3], all.fips$centroid[i])) {
    
    if (as.vector(st_distance(ash, all.fips$centroid[i])) < 100000) {
      
      too.close.ash <- c(too.close.ash, all.fips$GEOID[i])
      
    }
    
  }
  
}

too.close.nash <- c()

nash <- st_sfc(st_point(x = c(-86.7816, 36.1627)))
st_crs(nash) <- 4269

for (i in 1:nrow(all.fips)) {
  
  if (st_distance(nash, all.fips$centroid[i]) < st_distance(cdf$centroid[1], all.fips$centroid[i])) {
    
    if (as.vector(st_distance(nash, all.fips$centroid[i])) < 100000) {
      
      too.close.nash <- c(too.close.nash, all.fips$GEOID[i])
      
    }
    
  }
  
}

too.close.mon <- c()

mon <- st_sfc(st_point(x = c(-86.3077, 32.3792)))
st_crs(mon) <- 4269

for (i in 1:nrow(all.fips)) {
  
  if (st_distance(mon, all.fips$centroid[i]) < st_distance(cdf$centroid[2], all.fips$centroid[i])) {
    
    if (as.vector(st_distance(mon, all.fips$centroid[i])) < 100000) {
      
      too.close.mon <- c(too.close.mon, all.fips$GEOID[i])
      
    }
    
  }
  
}

all.fips <- all.fips %>% filter(! GEOID %in% too.close.atl) %>% filter(! GEOID %in% too.close.ash) %>% filter(! GEOID %in% too.close.nash) %>% filter(! GEOID %in% too.close.mon)

# Running the regression

rent <- lm(log(Rent) ~ TZ + log(Population) + log(Units) + log(Income) + f1 + f2 + f3 + f4 + factor(STATEFP) + factor(GEOID) + factor(Year), data = all.fips)

rentx <- coeftest(rent, vcov. = vcovCL(rent, type = 'HC1'))

# Placebo tests

set.seed(123456798)

all.fips$TZ2 <- sample(all.fips$TZ)
all.fips$Rent2 <- sample(all.fips$Rent)

rent2 <- lm(log(Rent) ~ TZ2 + log(Population) + log(Units) + log(Income) + f1 + f2 + f3 + f4 + factor(STATEFP) + factor(GEOID) + factor(Year), data = all.fips)

rent2x <- coeftest(rent2, vcov = vcovCL(rent2, type = 'HC1'))

rent3 <- lm(log(Rent2) ~ TZ + log(Population) + log(Units) + log(Income) + f1 + f2 + f3 + f4 + factor(STATEFP) + factor(GEOID) + factor(Year), data = all.fips)

rent3x <- coeftest(rent3, vcov = vcovCL(rent3, type = 'HC1'))

stargazer(rent, rentx, rent2, rent2x, rent3, rent3x, type = 'text', omit = c('GEOID', 'Year', 'STATEFP'), omit.stat = c('f', 'ser'))

# Saving results

write.csv(stargazer(rentx, rent2x, rent3x, omit.stat = c('f', 'ser'), omit = c('GEOID', 'Year', 'STATEFP')), paste0(direc, 'results/mechanism_results.txt'), row.names = FALSE)
write.csv(stargazer(rent, rent2, rent3, omit.stat = c('f', 'ser'), omit = c('GEOID', 'Year', 'STATEFP'), type = 'text'), paste0(direc, 'results/mechanism_results_xxx.txt'), row.names = FALSE)

# Summary statistics

sumdat <- all.fips[,c(34,27,32,33,35,18:21)]

colnames(sumdat) <- c('log(Rent)', 'Eastern Time Zone', 'log(Population)', 'log(Housing Units)', 'log(Income)',
                      'Knoxville, TN Metro Area', 'Chattanooga, TN Metro Area', 'Columbus, GA Metro Area', 'Tallahassee, FL Metro Area')

sumdat$`log(Rent)` <- log(sumdat$`log(Rent)`)
sumdat$`log(Population)` <- log(sumdat$`log(Population)`)
sumdat$`log(Housing Units)` <- log(sumdat$`log(Housing Units)`)
sumdat$`log(Income)` <- log(sumdat$`log(Income)`)

datasummary_skim(sumdat, fmt = '%.3f')

