# This script validates the mechanism / theory section of the paper

# Loading libraries

library(modelsummary)
library(tidycensus)
library(stargazer)
library(geomander)
library(tinytiger)
library(geojsonsf)
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

# Reading in the land cover data

trees <- read.csv(paste0(direc, 'data/county_level_proportions_2001_2011_2021.csv'))

# Defining which counties are border counties and assigning treatment

tz <- geojson_sf(paste0(direc, 'data/NTAD_Time_Zones_873827312334388986.geojson'))

cx <- counties()

fips <- c()

for (i in 1:nrow(trees)) {
  
  if (trees$County[i] < 10000) {
    
    fips <- c(fips, paste0('0', as.character(trees$County[i])))
    
  } else {
    
    fips <- c(fips, as.character(trees$County[i]))
    
  }
  
}

trees$GEOID <- fips

cx <- st_transform(cx, crs = st_crs(tz))

trees <- left_join(trees, cx, by = c('GEOID'))

trees <- st_as_sf(trees)
trees <- st_transform(trees, crs = st_crs(tz))

x <- as.data.frame(cbind(trees$INTPTLON, trees$INTPTLAT))
colnames(x) <- c('LON', 'LAT')
x$LON <- -1 * as.numeric(substr(x$LON, 2, nchar(x$LON)))
x$LAT <- as.numeric(substr(x$LAT, 2, nchar(x$LAT)))
x <- st_as_sf(x, coords = c('LON', 'LAT')) %>% st_set_crs(4326)

sf_use_s2(FALSE)

w <- c() # I was impatient and this was slower than I could handle, so I stopped, restarted, and visualized progress

for (i in 1:(nrow(x)/3)) {
  
  print(i)
  w <- c(w, st_within(x$geometry[i], tz$geometry))
  
}

ww <- c() # I'm not re-running the last loop

for (i in 1:length(w)) {
  
  ww <- c(ww, w[[i]][1])
  
}

trees$TZ_ID <- c(ww, ww, ww)
trees$TZ <- tz$zone[trees$TZ_ID]

# Filter out 2001 and 2011 data

trees <- trees %>% filter(Year > 2011)

# Note :: below is for KY south for the E/C boundary
# Note :: Gulf County, FL omitted because it splits time zones

keep.ec.e <- c(12037, 12077, 12039, 13099, 13061, 13253, 13239, 13259, 13215, 13053, 13145, 13149, 13285, 13233, 13143, 13045, 13115, 13055, 13083, 13087, 13295, 47065, 47143, 47145, 47151, 47129, 21231, 21199, 21147, 21217, 21045, 21123, 21093, 21163)
keep.ec.c <- c(12063, 12013, 1067, 1069, 1005, 1017, 1113, 1081, 1029, 1111, 1019, 1049, 1071, 47007, 47115, 47153, 47035, 47137, 47049, 21207, 21001, 21087, 21085, 21099, 21027)

# Note :: Stanley County, SD and Cherry County, NE omitted because they splits time zones 

keep.cm.c <- c(46021, 46129, 46107, 46065, 46119, 46075, 46121, 46095, 31171, 31117, 31111, 31085, 31087, 20109, 20093, 20203, 20193)
keep.cm.m <- c(46031, 46041, 46071, 46055, 46007, 31005, 31091, 31029, 31057, 31101, 31135, 20199, 20181, 20071, 20075)

ec.data <- trees %>% filter(County %in% c(keep.ec.e, keep.ec.c))
cm.data <- trees %>% filter(County %in% c(keep.cm.c, keep.cm.m))

ec.data$East <- as.integer(ec.data$County %in% keep.ec.e)
cm.data$East <- as.integer(cm.data$County %in% keep.cm.c)

joint <- rbind(ec.data, cm.data)

# Getting ACS data

ren <- c()
pop <- c()
emp <- c()
inc <- c()
hun <- c()
edu <- c()

for (y in 2015:2022) {
  
  covars <- get_acs(geography = 'county', year = y, variables = c('DP04_0134', 'DP03_0062', 'DP03_0009P', 'DP04_0001', 'DP05_0001', 'DP02_0068P'))
  
  for (i in 1:nrow(joint)) {
    
    print(paste0('Checking observation ', i, ' of ', nrow(joint), '.......'))
    
    tmp <- covars[which(covars$GEOID == joint$GEOID[i]),]
    
    tmp.ren <- tmp %>% filter(variable == 'DP04_0134')
    tmp.inc <- tmp %>% filter(variable == 'DP03_0062')
    tmp.emp <- tmp %>% filter(variable == 'DP03_0009P')
    tmp.pop <- tmp %>% filter(variable == 'DP05_0001')
    tmp.hun <- tmp %>% filter(variable == 'DP04_0001')
    tmp.edu <- tmp %>% filter(variable == 'DP02_0068P')
    
    ren <- c(ren, mean(tmp.ren$estimate))
    inc <- c(inc, mean(tmp.inc$estimate))
    pop <- c(pop, mean(tmp.pop$estimate))
    emp <- c(emp, mean(tmp.emp$estimate))
    hun <- c(hun, mean(tmp.hun$estimate))
    edu <- c(edu, mean(tmp.edu$estimate))
    
  }
  
}

joint <- rbind(joint, joint, joint, joint, joint, joint, joint, joint)

joint$Rents <- ren
joint$Population <- pop
joint$Income <- inc
joint$Unemployment <- emp
joint$Housing_Units <- hun
joint$Education <- edu
joint$Year <- c(rep(2015, nrow(joint)/8), rep(2016, nrow(joint)/8), rep(2017, nrow(joint)/8), rep(2018, nrow(joint)/8),
                rep(2019, nrow(joint)/8), rep(2020, nrow(joint)/8), rep(2021, nrow(joint)/8), rep(2022, nrow(joint)/8))

# Defining the time zone boundary

border <- c()

for (i in 1:nrow(joint)) {
  
  if (joint$TZ[i] == 'Eastern') {
    
    border <- c(border, 'Eastern-Central')
    
  }
  
  if (joint$TZ[i] == 'Mountain') {
    
    border <- c(border, 'Central-Mountain')
    
  }
  
  if (joint$TZ[i] == 'Central') {
    
    if (joint$East[i] == 1) {
      
      border <- c(border, 'Central-Mountain')
      
    } else {
      
      border <- c(border, 'Eastern-Central')
      
    }
    
  }
  
}

joint$Border <- border

# Creating an indicator for GA and AL counties

gaal <- c()

for (i in 1:nrow(joint)) {
  
  if (joint$County[i] < 10000) {
    
    gaal <- c(gaal, 1)
    
  } else if (joint$STATEFP[i] == '13') {
    
    gaal <- c(gaal, 1)
    
  } else {
    
    gaal <- c(gaal, 0)
    
  }
  
}

joint$GAAL <- gaal

joint2 <- joint %>% filter(gaal == 0)

# Convert unemployment to [0,1]

joint$Unemployment <- joint$Unemployment / 100
joint2$Unemployment <- joint2$Unemployment / 100

# Running model-verification models

mod1 <- lm(log(Rents) ~ East + log(Housing_Units) + log(Population) + log(Income)
           + Unemployment + factor(Border) + factor(STATEFP) + factor(Year), data = joint)

mod2 <- lm(log(Rents) ~ East + log(Housing_Units) + log(Population) + log(Income)
           + Unemployment + factor(Border) + factor(STATEFP) + factor(Year), data = joint2)

xmod1 <- coeftest(mod1, vcov. = vcovCL(mod1, type = 'HC1'))
xmod2 <- coeftest(mod2, vcov. = vcovCL(mod2, type = 'HC1'))

stargazer(mod1, xmod1, mod2, xmod2, type = 'text', omit = c('STATEFP', 'Year', 'GEOID'), omit.stat = c('f', 'ser'))

# Placebo testing

set.seed(42069)

# Randomizing treatment

vals <- runif(nrow(joint))
joint$Placebo <- as.integer(vals >= .5)

pmod1 <- lm(log(Rents) ~ Placebo + log(Housing_Units) + log(Population) + log(Income)
            + Unemployment + factor(Border) + factor(STATEFP) + factor(Year), data = joint)

xpmod1 <- coeftest(pmod1, vcov. = vcovCL(pmod1, type = 'HC1'))

# Randomizing outcomes

joint$Fake <- sample(joint$Rents)

pmod2 <- lm(log(Fake) ~ East + log(Housing_Units) + log(Population) + log(Income)
            + Unemployment + factor(Border) + factor(STATEFP) + factor(Year), data = joint)

xpmod2 <- coeftest(pmod2, vcov. = vcovCL(pmod2, type = 'HC1'))

# Results

write.csv(stargazer(mod1, mod2, pmod1, pmod2, omit = c('STATEFP'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/rent_placebo.txt'), row.names = FALSE)

write.csv(stargazer(xmod1, xmod2, xpmod1, xpmod2, omit = c('STATEFP'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/rent_placebo_robust.txt'), row.names = FALSE)

stargazer(xmod1, xmod2, xpmod1, xpmod2, type = 'text', omit = c('STATEFP'), omit.stat = c('f', 'ser'))

# Summary statistics of this data

sumdat <- as.data.frame(cbind(joint$East, joint$Rents, joint$Housing_Units, joint$Population, joint$Income, joint$Unemployment))

sumdat$V2 <- log(sumdat$V2)
sumdat$V3 <- log(sumdat$V3)
sumdat$V4 <- log(sumdat$V4)
sumdat$V5 <- log(sumdat$V5)

colnames(sumdat) <- c('East', 'log(Rent)', 'log(Housing Units)', 'log(Population)', 'log(Income)', 'Unemplyoment Rate')

datasummary_skim(sumdat, fmt = '%.3f')

# Creating a figure to motivate this

j1 <- joint %>% filter(STATEFP %in% c('01', '12', '13', '21', '47')) %>% filter(Year == 2022)
j2 <- joint %>% filter(STATEFP %in% c('20', '31', '46')) %>% filter(Year == 2022)

pal1 <- colorNumeric(palette = c('Blue', 'White', 'Red'), domain = j1$Rents)
pal2 <- colorNumeric(palette = c('Blue', 'White', 'Red'), domain = j2$Rents)

leaflet(j1$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 1.0, color = 'black', fillColor = pal1(j1$Rents))
leaflet(j2$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 1.0, color = 'black', fillColor = pal2(j2$Rents))

# Creating a figure to show the sample

j1 <- joint %>% filter(STATEFP %in% c('01', '12', '13', '21', '47')) %>% filter(Year == 2022)
j2 <- joint %>% filter(STATEFP %in% c('20', '31', '46')) %>% filter(Year == 2022)

pal <- colorNumeric(palette = c('red4', 'orange'), domain = j1$East)

leaflet(j1$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 0.666, color = 'black', fillColor = pal(j1$East))
leaflet(j2$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 0.666, color = 'black', fillColor = pal(j2$East))

# Running model-verification models with RUC codes for a reviewer

ruc <- read.csv(paste0(direc, 'data/Ruralurbancontinuumcodes2023.csv'))
ruc <- ruc %>% filter(Attribute == 'RUCC_2023')
ruc <- ruc[,c(1,5)]
colnames(ruc) <- c('County', 'RUC')

joint <- left_join(joint, ruc, by = 'County')
joint2 <- left_join(joint2, ruc, by = 'County')

rucmod1 <- lm(log(Rents) ~ East + factor(RUC) + log(Housing_Units) + log(Population) + log(Income)
              + Unemployment + factor(Border) + factor(STATEFP) + factor(Year), data = joint)

rucmod2 <- lm(log(Rents) ~ East + factor(RUC) + log(Housing_Units) + log(Population) + log(Income)
              + Unemployment + factor(Border) + factor(STATEFP)*factor(Year), data = joint)

xrucmod1 <- coeftest(rucmod1, vcov. = vcovCL(rucmod1, type = 'HC1'))
xrucmod2 <- coeftest(rucmod2, vcov. = vcovCL(rucmod2, type = 'HC1'))

stargazer(rucmod1, xrucmod1, rucmod2, xrucmod2, type = 'text', omit = c('STATEFP', 'Year'), omit.stat = c('f', 'ser'))

