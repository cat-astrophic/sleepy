# This script runs the econometric analysis and creates figures for my sleepy project

# If you have not updated stargazer, update the filepath and run the accompanying script stargazer_fix.R

# stargazer_fix.R comes directly from :: https://gist.github.com/alexeyknorre/b0780836f4cec04d41a863a683f91b53

# source('C:/Users/User/Documents/stargazer_fix.R')

# Loading libraries

library(modelsummary)
library(stargazer)
library(sandwich)
library(leaflet)
library(ggplot2)
library(lmtest)
library(tigris)
library(dplyr)
library(terra)
library(sf)

# Project directory

direc <- 'D:/sleepy/'

# Reading in the data

rawdf <- read.csv(paste0(direc, 'data/data.csv'))

# Subset to remove observations with missing geospatial data

rawdf <- rawdf %>% filter(is.na(Latitude) == FALSE)

# Clean Sex variable

gender <- c()

for (i in 1:dim(rawdf)[1]) {
  
  if (rawdf$Sex[i] == 'M   ') {
    
    gender <- c(gender, 'M')
    
  } else if (rawdf$Sex[i] == 'F   ') {
    
    gender <- c(gender, 'F')
    
  } else {
    
    gender <- c(gender, rawdf$Sex[i])
    
  }
  
}

rawdf$Gender <- gender
rawdf$Female <- as.integer(rawdf$Gender == 'F')

# Create an observation level ID variable for rawdf

rawdf$ID <- 1:dim(rawdf)[1]

# Creating a runner level ID variable for potential runner level fixed effects

rawdf$Runner_ID <- paste0(rawdf$Name, rawdf$City, rawdf$State)

# Convert df to sf type

rdf <- st_as_sf(rawdf, coords = c('Longitude', 'Latitude'), crs = 4269)

# Getting the Mercer County, ND shapefile

nd <- counties(state = 'ND')
mc <- nd %>% filter(NAME == 'Mercer')

# Determine who is in Mercer County, ND

treated <- st_filter(rdf, mc)
rdf$Treated <- as.integer(rdf$ID %in% treated$ID)
rawdf$Treated <- as.integer(rawdf$ID %in% treated$ID)

# Creating a post-treatment variable

rdf$Post <- as.integer(rdf$Year > 2010)

# Read in a shapefile for US time zones

tz <- st_read(paste0(direc, 'data/Time_Zones/Time_Zones.shp'))
tz <- st_transform(tz, 4269)

# Determine time zones for runners

sf_use_s2(FALSE)
zones <- st_intersects(rdf$geometry, tz$geometry)

poop <- c()

for (i in 1:length(zones)) {
  
  if (is.na(zones[[i]][1]) == FALSE) {
    
    poop <- c(poop, zones[[i]][1])
    
  } else {
    
    poop <- c(poop, 99)
    
  }
  
}

rdf$Zones <- poop

duke <- c()

for (i in 1:length(poop)) {
  
  if (rdf$Treated[i] == 1) {
    
    if (rdf$Post[i] == 1) {
      
      duke <- c(duke, 'Central')
      
    } else {
      
      duke <- c(duke, 'Mountain')
      
    }
    
  } else if (poop[i] == 3) {
    
    duke <- c(duke, 'Central')
    
  } else if (poop[i] == 4) {
    
    duke <- c(duke, 'Mountain')
    
  } else {
    
    duke <- c(duke, 'Other')
    
  }
  
}

rdf$Time_Zones <- duke

# Adjusting the distance variable from miles to meters

rdf$DistanceX <- rdf$Distance * 1609

# Running main regressions

mod1 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf)
mod2 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf[which(rdf$Distance < 180),])
mod3 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf[which(rdf$Distance < 150),])
mod4 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf[which(rdf$Distance < 120),])
mod5 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf[which(rdf$Distance < 90),])
mod6 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf[which(rdf$Distance < 75),])

mod1x <- coeftest(mod1, vcov = vcovCL, cluster = ~City)
mod2x <- coeftest(mod2, vcov = vcovCL, cluster = ~City)
mod3x <- coeftest(mod3, vcov = vcovCL, cluster = ~City)
mod4x <- coeftest(mod4, vcov = vcovCL, cluster = ~City)
mod5x <- coeftest(mod5, vcov = vcovCL, cluster = ~City)
mod6x <- coeftest(mod6, vcov = vcovCL, cluster = ~City)

stargazer(mod1, mod2, mod3, mod4, mod5, mod6, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(mod1x, mod2x, mod3x, mod4x, mod5x, mod6x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Running the regressions with only a central time zone control group

crdf <- rbind(rdf %>% filter(Treated == 1), rdf %>% filter(Time_Zones == 'Central'))
crdf <- crdf[!duplicated(crdf),]

cmod1 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = crdf)
cmod2 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = crdf[which(crdf$Distance < 180),])
cmod3 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = crdf[which(crdf$Distance < 150),])
cmod4 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = crdf[which(crdf$Distance < 120),])
cmod5 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = crdf[which(crdf$Distance < 90),])
cmod6 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = crdf[which(crdf$Distance < 75),])

cmod1x <- coeftest(cmod1, vcov = vcovCL, cluster = ~City)
cmod2x <- coeftest(cmod2, vcov = vcovCL, cluster = ~City)
cmod3x <- coeftest(cmod3, vcov = vcovCL, cluster = ~City)
cmod4x <- coeftest(cmod4, vcov = vcovCL, cluster = ~City)
cmod5x <- coeftest(cmod5, vcov = vcovCL, cluster = ~City)
cmod6x <- coeftest(cmod6, vcov = vcovCL, cluster = ~City)

stargazer(cmod1, cmod2, cmod3, cmod4, cmod5, cmod6, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(cmod1x, cmod2x, cmod3x, cmod4x, cmod5x, cmod6x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Running the regressions with only a mountain time zone control group

mrdf <- rbind(rdf %>% filter(Treated == 1), rdf %>% filter(Time_Zones == 'Mountain'))
mrdf <- mrdf[!duplicated(mrdf),]

mmod1 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = mrdf)
mmod2 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = mrdf[which(mrdf$Distance < 180),])
mmod3 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = mrdf[which(mrdf$Distance < 150),])
mmod4 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = mrdf[which(mrdf$Distance < 120),])

mmod1x <- coeftest(mmod1, vcov = vcovCL(mmod1, type = 'HC0'))
mmod2x <- coeftest(mmod2, vcov = vcovCL(mmod2, type = 'HC0'))
mmod3x <- coeftest(mmod3, vcov = vcovCL(mmod3, type = 'HC0'))
mmod4x <- coeftest(mmod4, vcov = vcovCL(mmod4, type = 'HC0'))

stargazer(mmod1, mmod2, mmod3, mmod4, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(mmod1x, mmod2x, mmod3x, mmod4x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Heterogeneity analysis for gender

rdfmen <- rdf %>% filter(Gender == 'M')
rdfwomen <- rdf %>% filter(Gender == 'F')

menmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + log(Age) + factor(Year) + factor(Runner_ID), data = rdfmen)
fmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + log(Age) + factor(Year) + factor(Runner_ID), data = rdfwomen)

menmodx <- coeftest(menmod, vcov = vcovCL, cluster = ~City)
fmodx <- coeftest(fmod, vcov = vcovCL, cluster = ~City)

stargazer(menmod, fmod, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(menmodx, fmodx, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Running heterogeneity analyses for ability

cutm <- quantile(rdf[which(rdf$Event == 'Marathon'),]$Seconds, c(.5), na.rm = TRUE)
cuth <- quantile(rdf[which(rdf$Event != 'Marathon'),]$Seconds, c(.5), na.rm = TRUE)

rdfm1 <- rdf[which(rdf$Event == 'Marathon'),] %>% filter(Seconds >= cutm)
rdfm2 <- rdf[which(rdf$Event == 'Marathon'),] %>% filter(Seconds < cutm)

rdfh1 <- rdf[which(rdf$Event != 'Marathon'),] %>% filter(Seconds >= cuth)
rdfh2 <- rdf[which(rdf$Event != 'Marathon'),] %>% filter(Seconds < cuth)

q1 <- rbind(rdfm1, rdfh1)
q2 <- rbind(rdfm2, rdfh2)

q1mod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = q1)
q2mod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = q2)

q1modx <- coeftest(q1mod, vcov = vcovCL, cluster = ~City)
q2modx <- coeftest(q2mod, vcov = vcovCL, cluster = ~City)

stargazer(q1mod, q2mod, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(q1modx, q2modx, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Heterogeneity analysis by age

cuta <- quantile(rdf$Age, c(.5), na.rm = TRUE)

a1 <- rdf %>% filter(Age > cuta)
a2 <- rdf %>% filter(Age <= cuta)

a1mod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = a1)
a2mod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = a2)

a1modx <- coeftest(a1mod, vcov = vcovCL, cluster = ~City)
a2modx <- coeftest(a2mod, vcov = vcovCL, cluster = ~City)

stargazer(a1mod, a2mod, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(a1modx, a2modx, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Placebo testing via defining the treatment as the ND Mountain counties

sugar.mountain <- ifelse(rdf$Time_Zones == 'Mountain', 1, 0)
sugar.hell <- ifelse(rdf$State == 'ND', 1, 0)
rdf$Placebo <- sugar.mountain * sugar.hell

pmod1 <- lm(log(Seconds) ~ Placebo*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf)
pmod1x <- coeftest(pmod1, vcov = vcovCL, cluster = ~City)

# Placebo testing via defining the treatment as the westernmost ND Central counties excluding Mercer

sugar.central <- ifelse(rdf$Time_Zones == 'Central', 1, 0)
west <- c()

for (i in 1:nrow(rdf)) {
  
  xxx <- strsplit(rdf$Coordinates[i], ',')[[1]][2]
  xxx <- substr(xxx, 2, nchar(xxx)-1)
  west <- c(west, as.numeric(xxx))
  
}

sugar.west <- ifelse(west < -100.8, 1, 0)

rdf$Placebo2 <- sugar.central * sugar.hell * sugar.west

pmod2 <- lm(log(Seconds) ~ Placebo2*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf)
pmod2x <- coeftest(pmod2, vcov = vcovCL, cluster = ~City)

# Placebo testing via randomizing treatment

set.seed(42069)

vals <- runif(nrow(rdf))
rdf$Placebo3 <- as.integer(vals >= .5)

pmod3 <- lm(log(Seconds) ~ Placebo3*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = rdf)
pmod3x <- coeftest(pmod3, vcov = vcovCL, cluster = ~City)

# Placebo testing via randomizing outcomes

mardf <- rdf %>% filter(Event == 'Marathon')
halfdf <- rdf %>% filter(Event != 'Marathon')

mardf$Thirds <- sample(mardf$Seconds)
halfdf$Thirds <- sample(halfdf$Seconds)

randdf <- rbind(mardf, halfdf)

pmod4 <- lm(log(Thirds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year) + factor(Runner_ID), data = randdf)
pmod4x <- coeftest(pmod4, vcov = vcovCL, cluster = ~City)

# Results from placebo testing

stargazer(pmod1, pmod2, pmod3, pmod4, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(pmod1x, pmod2x, pmod3x, pmod4x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Repeating the entire analysis without runner level fixed effects for comparison

# Running regressions with all observations

fastmod1 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = rdf)
fastmod2 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = rdf[which(rdf$Distance < 180),])
fastmod3 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = rdf[which(rdf$Distance < 150),])
fastmod4 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = rdf[which(rdf$Distance < 120),])
fastmod5 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = rdf[which(rdf$Distance < 90),])
fastmod6 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = rdf[which(rdf$Distance < 75),])

fastmod1x <- coeftest(fastmod1, vcov = vcovCL, cluster = ~City)
fastmod2x <- coeftest(fastmod2, vcov = vcovCL, cluster = ~City)
fastmod3x <- coeftest(fastmod3, vcov = vcovCL, cluster = ~City)
fastmod4x <- coeftest(fastmod4, vcov = vcovCL, cluster = ~City)
fastmod5x <- coeftest(fastmod5, vcov = vcovCL, cluster = ~City)
fastmod6x <- coeftest(fastmod6, vcov = vcovCL, cluster = ~City)

stargazer(fastmod1, fastmod2, fastmod3, fastmod4, fastmod5, fastmod6, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(fastmod1x, fastmod2x, fastmod3x, fastmod4x, fastmod5x, fastmod6x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Running the regressions with only a central time zone control group

crdf <- rbind(rdf %>% filter(Treated == 1), rdf %>% filter(Time_Zones == 'Central'))
crdf <- crdf[!duplicated(crdf),]

cfastmod1 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = crdf)
cfastmod2 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = crdf[which(crdf$Distance < 180),])
cfastmod3 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = crdf[which(crdf$Distance < 150),])
cfastmod4 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = crdf[which(crdf$Distance < 120),])
cfastmod5 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = crdf[which(crdf$Distance < 90),])
cfastmod6 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = crdf[which(crdf$Distance < 75),])

cfastmod1x <- coeftest(cfastmod1, vcov = vcovCL, cluster = ~City)
cfastmod2x <- coeftest(cfastmod2, vcov = vcovCL, cluster = ~City)
cfastmod3x <- coeftest(cfastmod3, vcov = vcovCL, cluster = ~City)
cfastmod4x <- coeftest(cfastmod4, vcov = vcovCL, cluster = ~City)
cfastmod5x <- coeftest(cfastmod5, vcov = vcovCL, cluster = ~City)
cfastmod6x <- coeftest(cfastmod6, vcov = vcovCL, cluster = ~City)

stargazer(cfastmod1, cfastmod2, cfastmod3, cfastmod4, cfastmod5, cfastmod6, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(cfastmod1x, cfastmod2x, cfastmod3x, cfastmod4x, cfastmod5x, cfastmod6x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Running the regressions with only a mountain time zone control group

mrdf <- rbind(rdf %>% filter(Treated == 1), rdf %>% filter(Time_Zones == 'Mountain'))
mrdf <- mrdf[!duplicated(mrdf),]

mfastmod1 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = mrdf)
mfastmod2 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = mrdf[which(mrdf$Distance < 180),])
mfastmod3 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = mrdf[which(mrdf$Distance < 150),])
mfastmod4 <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = mrdf[which(mrdf$Distance < 120),])

mfastmod1x <- coeftest(mfastmod1, vcov = vcovCL, cluster = ~City)
mfastmod2x <- coeftest(mfastmod2, vcov = vcovCL, cluster = ~City)
mfastmod3x <- coeftest(mfastmod3, vcov = vcovCL, cluster = ~City)
mfastmod4x <- coeftest(mfastmod4, vcov = vcovCL, cluster = ~City)

stargazer(mfastmod1, mfastmod2, mfastmod3, mfastmod4, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(mfastmod1x, mfastmod2x, mfastmod3x, mfastmod4x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Heterogeneity analysis for gender

rdfmen <- rdf %>% filter(Gender == 'M')
rdfwomen <- rdf %>% filter(Gender == 'F')

menfastmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + log(Age) + factor(Year), data = rdfmen)
ffastmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + log(Age) + factor(Year), data = rdfwomen)

menfastmodx <- coeftest(menfastmod, vcov = vcovCL, cluster = ~City)
ffastmodx <- coeftest(ffastmod, vcov = vcovCL, cluster = ~City)

stargazer(menfastmod, ffastmod, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(menfastmodx, ffastmodx, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Running heterogeneity analyses for ability

cutm <- quantile(rdf[which(rdf$Event == 'Marathon'),]$Seconds, c(.5), na.rm = TRUE)
cuth <- quantile(rdf[which(rdf$Event != 'Marathon'),]$Seconds, c(.5), na.rm = TRUE)

rdfm1 <- rdf[which(rdf$Event == 'Marathon'),] %>% filter(Seconds >= cutm)
rdfm2 <- rdf[which(rdf$Event == 'Marathon'),] %>% filter(Seconds < cutm)

rdfh1 <- rdf[which(rdf$Event != 'Marathon'),] %>% filter(Seconds >= cuth)
rdfh2 <- rdf[which(rdf$Event != 'Marathon'),] %>% filter(Seconds < cuth)

q1 <- rbind(rdfm1, rdfh1)
q2 <- rbind(rdfm2, rdfh2)

q1fastmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = q1)
q2fastmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = q2)

q1fastmodx <- coeftest(q1fastmod, vcov = vcovCL, cluster = ~City)
q2fastmodx <- coeftest(q2fastmod, vcov = vcovCL, cluster = ~City)

stargazer(q1fastmod, q2fastmod, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(q1fastmodx, q2fastmodx, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Heterogeneity analysis by age

cuta <- quantile(rdf$Age, c(.5), na.rm = TRUE)

a1 <- rdf %>% filter(Age > cuta)
a2 <- rdf %>% filter(Age <= cuta)

a1fastmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = a1)
a2fastmod <- lm(log(Seconds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = a2)

a1fastmodx <- coeftest(a1fastmod, vcov = vcovCL, cluster = ~City)
a2fastmodx <- coeftest(a2fastmod, vcov = vcovCL, cluster = ~City)

stargazer(a1fastmod, a2fastmod, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(a1fastmodx, a2fastmodx, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Placebo testing via defining the treatment as the ND Mountain counties

sugar.mountain <- ifelse(rdf$Time_Zones == 'Mountain', 1, 0)
sugar.hell <- ifelse(rdf$State == 'ND', 1, 0)
rdf$Placebo <- sugar.mountain * sugar.hell

pfastmod1 <- lm(log(Seconds) ~ Placebo*Post + factor(Event) + Female + log(Age) + factor(Year), data = rdf)
pfastmod1x <- coeftest(pfastmod1, vcov = vcovCL, cluster = ~City)

# Placebo testing via defining the treatment as the westernmost ND Central counties excluding Mercer

sugar.central <- ifelse(rdf$Time_Zones == 'Central', 1, 0)
west <- c()

for (i in 1:nrow(rdf)) {
  
  xxx <- strsplit(rdf$Coordinates[i], ',')[[1]][2]
  xxx <- substr(xxx, 2, nchar(xxx)-1)
  west <- c(west, as.numeric(xxx))
  
}

sugar.west <- ifelse(west < -100.8, 1, 0)

rdf$Placebo2 <- sugar.central * sugar.hell * sugar.west

pfastmod2 <- lm(log(Seconds) ~ Placebo2*Post + factor(Event) + Female + log(Age) + factor(Year), data = rdf)
pfastmod2x <- coeftest(pfastmod2, vcov = vcovCL, cluster = ~City)

# Placebo testing via randomizing treatment

set.seed(7734)

vals <- runif(nrow(rdf))
rdf$Placebo3 <- as.integer(vals >= .5)

pfastmod3 <- lm(log(Seconds) ~ Placebo3*Post + factor(Event) + Female + log(Age) + factor(Year), data = rdf)
pfastmod3x <- coeftest(pfastmod3, vcov = vcovCL, cluster = ~City)

# Placebo testing via randomizing outcomes

mardf <- rdf %>% filter(Event == 'Marathon')
halfdf <- rdf %>% filter(Event != 'Marathon')

mardf$Thirds <- sample(mardf$Seconds)
halfdf$Thirds <- sample(halfdf$Seconds)

randdf <- rbind(mardf, halfdf)

pfastmod4 <- lm(log(Thirds) ~ Treated*Post + poly(Distance, 2, raw = TRUE) + factor(Event) + Female + log(Age) + factor(Year), data = randdf)
pfastmod4x <- coeftest(pfastmod4, vcov = vcovCL, cluster = ~City)

# Results from placebo testing

stargazer(pfastmod1, pfastmod2, pfastmod3, pfastmod4, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))
stargazer(pfastmod1x, pfastmod2x, pfastmod3x, pfastmod4x, type = 'text', omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser'))

# Saving results

write.csv(stargazer(fastmod1x, mod1x, fastmod2x, mod2x, fastmod3x, mod3x, fastmod4x, mod4x, fastmod5x, mod5x, fastmod6x, mod6x, omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/ROBUST_combined_main_models.txt'), row.names = FALSE)

write.csv(stargazer(cfastmod1x, cmod1x, cfastmod2x, cmod2x, cfastmod3x, cmod3x, cfastmod4x, cmod4x, cfastmod5x, cmod5x, cfastmod6x, cmod6x, omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/ROBUST_combined_central_models.txt'), row.names = FALSE)

write.csv(stargazer(mfastmod1x, mmod1x, mfastmod2x, mmod2x, mfastmod3x, mmod3x, mfastmod4x, mmod4x, omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/ROBUST_combined_mountain_models.txt'), row.names = FALSE)

write.csv(stargazer(menmodx, fmodx, q1modx, q2modx, a1modx, a2modx, pmod1x, pmod2x, pmod3x, pmod4x, omit = c('Year', 'Runner_ID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/ROBUST_robustness_and_placebo_tests.txt'), row.names = FALSE)

