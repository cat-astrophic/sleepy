# This script looks at land cover from the NLCD in an rdd setting around time zone boundaries

# The idea is that there should be a sharp increase in forested land immediately west of a time zone boundary because land is less valuable

# Opposite goes for developed land

# Loading libraries

library(modelsummary)
library(tidycensus)
library(stargazer)
library(geomander)
library(geojsonsf)
library(tinytiger)
library(sandwich)
library(leaflet)
library(ggplot2)
library(xgboost)
library(lmtest)
library(tigris)
library(dplyr)
library(terra)
library(lutz)
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

# Filter out 2001 data

xxx <- trees
trees <- trees %>% filter(Year > 2010)

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

# Getting covariates from the ACS

covars.2021 <- get_acs(geography = 'county', year = 2021, variables = c('DP03_0062', 'DP03_0009P', 'DP03_0019P', 'DP04_0001', 'DP05_0001', 'DP02_0068P', 'DP03_0033', 'DP03_0034', 'DP03_0035', 'DP03_0036', 'DP03_0037', 'DP03_0038', 'DP03_0039', 'DP03_0040', 'DP03_0041', 'DP03_0042', 'DP03_0043', 'DP03_0044', 'DP03_0045'))

covars.2011 <- get_acs(geography = 'county', year = 2011, variables = c('DP03_0062', 'DP03_0009P', 'DP03_0019P', 'DP04_0001', 'DP05_0001', 'DP02_0067P', 'DP03_0033', 'DP03_0034', 'DP03_0035', 'DP03_0036', 'DP03_0037', 'DP03_0038', 'DP03_0039', 'DP03_0040', 'DP03_0041', 'DP03_0042', 'DP03_0043', 'DP03_0044', 'DP03_0045'))

cpi.11.21 <- 1.20

pop <- c()
inc <- c()

for (i in 1:nrow(joint)) {
  
  print(paste0('Checking observation ', i, ' of ', nrow(joint), '.......'))
  
  if (joint$Year[i] == 2021) {
    
    tmp <- covars.2021[which(covars.2021$GEOID == joint$GEOID[i]),]
    
  } else {
    
    tmp <- covars.2021[which(covars.2011$GEOID == joint$GEOID[i]),]
    
  }
  
  tmp.inc <- tmp %>% filter(variable == 'DP03_0062')
  tmp.pop <- tmp %>% filter(variable == 'DP05_0001')
  
  if (joint$Year[i] == 2021) {
    
    inc <- c(inc, mean(tmp.inc$estimate))
    
  } else {
    
    inc <- c(inc, cpi.11.21 * mean(tmp.inc$estimate))
    
  }
  
  pop <- c(pop, mean(tmp.pop$estimate))
  
}

joint$Population <- pop
joint$Income <- inc

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

# Adding in initial values of land cover from 2001

init.dev <- c()
init.agr <- c()
init.for <- c()
init.wet <- c()
init.wat <- c()

trees.01 <- xxx %>% filter(Year == 2001)

for (i in 1:nrow(joint)) {
  
  tmp <- trees.01 %>% filter(GEOID == joint$GEOID[i])
  
  init.dev <- c(init.dev, tmp$Development[1])
  init.agr <- c(init.agr, tmp$Agriculture[1])
  init.for <- c(init.for, tmp$Forests[1])
  init.wet <- c(init.wet, tmp$Wetlands[1])
  init.wat <- c(init.wat, tmp$Water[1])
  
}

joint$Initial_Development <- init.dev
joint$Initial_Agriculture <- init.agr
joint$Initial_Forests <- init.for
joint$Initial_Wetlands <- init.wet
joint$Initial_Water <- init.wat

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

# RDD on E/C and C/M time zone boundaries

mod.all1 <- lm(Development ~ East + factor(Border) + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint)

mod.all2 <- lm(Forests ~ East + factor(Border) + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint)

mod.all3 <- lm(Agriculture ~ East + factor(Border) + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint)

mod.ec1 <- lm(Development ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
              + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint[which(joint$Border == 'Eastern-Central'),])

mod.ec2 <- lm(Forests ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
              + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint[which(joint$Border == 'Eastern-Central'),])

mod.ec3 <- lm(Agriculture ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
              + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint[which(joint$Border == 'Eastern-Central'),])

mod.cm1 <- lm(Development ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
              + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint[which(joint$Border == 'Central-Mountain'),])

mod.cm2 <- lm(Forests ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
              + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint[which(joint$Border == 'Central-Mountain'),])

mod.cm3 <- lm(Agriculture ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
              + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint[which(joint$Border == 'Central-Mountain'),])

# RDD on E/C and C/M time zone boundaries excluding GA and AL counties

xmod.all1 <- lm(Development ~ East + factor(Border) + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
                + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2)

xmod.all2 <- lm(Forests ~ East + factor(Border) + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
                + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2)

xmod.all3 <- lm(Agriculture ~ East + factor(Border) + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
                + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2)

xmod.ec1 <- lm(Development ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2[which(joint2$Border == 'Eastern-Central'),])

xmod.ec2 <- lm(Forests ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2[which(joint2$Border == 'Eastern-Central'),])

xmod.ec3 <- lm(Agriculture ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2[which(joint2$Border == 'Eastern-Central'),])

xmod.cm1 <- lm(Development ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2[which(joint2$Border == 'Central-Mountain'),])

xmod.cm2 <- lm(Forests ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2[which(joint2$Border == 'Central-Mountain'),])

xmod.cm3 <- lm(Agriculture ~ East + factor(Year) + log(Population) + log(Income) + Initial_Development + Initial_Agriculture
               + Initial_Forests + Initial_Water + factor(GEOID) + factor(STATEFP), data = joint2[which(joint2$Border == 'Central-Mountain'),])

# Robust standard errors clustered at the county level

zmod.all1 <- coeftest(mod.all1, vcov. = vcovCL(mod.all1, type = 'HC1'))
zmod.all2 <- coeftest(mod.all2, vcov. = vcovCL(mod.all2, type = 'HC1'))
zmod.all3 <- coeftest(mod.all3, vcov. = vcovCL(mod.all3, type = 'HC1'))

zxmod.all1 <- coeftest(xmod.all1, vcov. = vcovCL(xmod.all1, type = 'HC1'))
zxmod.all2 <- coeftest(xmod.all2, vcov. = vcovCL(xmod.all2, type = 'HC1'))
zxmod.all3 <- coeftest(xmod.all3, vcov. = vcovCL(xmod.all3, type = 'HC1'))

zmod.ec1 <- coeftest(mod.ec1, vcov. = vcovCL(mod.ec1, type = 'HC1'))
zmod.ec2 <- coeftest(mod.ec2, vcov. = vcovCL(mod.ec2, type = 'HC1'))
zmod.ec3 <- coeftest(mod.ec3, vcov. = vcovCL(mod.ec3, type = 'HC1'))

zxmod.ec1 <- coeftest(xmod.ec1, vcov. = vcovCL(xmod.ec1, type = 'HC1'))
zxmod.ec2 <- coeftest(xmod.ec2, vcov. = vcovCL(xmod.ec2, type = 'HC1'))
zxmod.ec3 <- coeftest(xmod.ec3, vcov. = vcovCL(xmod.ec3, type = 'HC1'))

zmod.cm1 <- coeftest(mod.cm1, vcov. = vcovCL(mod.cm1, type = 'HC1'))
zmod.cm2 <- coeftest(mod.cm2, vcov. = vcovCL(mod.cm2, type = 'HC1'))
zmod.cm3 <- coeftest(mod.cm3, vcov. = vcovCL(mod.cm3, type = 'HC1'))

zxmod.cm1 <- coeftest(xmod.cm1, vcov. = vcovCL(xmod.cm1, type = 'HC1'))
zxmod.cm2 <- coeftest(xmod.cm2, vcov. = vcovCL(xmod.cm2, type = 'HC1'))
zxmod.cm3 <- coeftest(xmod.cm3, vcov. = vcovCL(xmod.cm3, type = 'HC1'))

# Main results viz / writing results to file

write.csv(stargazer(mod.all1, mod.all2, mod.all3, mod.ec1, mod.ec2, mod.ec3, mod.cm1, mod.cm2, mod.cm3, omit = c('GEOID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/land_cover_rdd.txt'), row.names = FALSE)

write.csv(stargazer(zmod.all1, zmod.all2, zmod.all3, zmod.ec1, zmod.ec2, zmod.ec3, zmod.cm1, zmod.cm2, zmod.cm3, omit = c('GEOID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/land_cover_rdd_robust.txt'), row.names = FALSE)

write.csv(stargazer(xmod.all1, xmod.all2, xmod.all3, xmod.ec1, xmod.ec2, xmod.ec3, xmod.cm1, xmod.cm2, xmod.cm3, omit = c('GEOID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/land_cover_rdd_no_gaal.txt'), row.names = FALSE)

write.csv(stargazer(zxmod.all1, zxmod.all2, zxmod.all3, zxmod.ec1, zxmod.ec2, zxmod.ec3, zxmod.cm1, zxmod.cm2, zxmod.cm3, omit = c('GEOID'), omit.stat = c('f', 'ser')),
          paste0(direc, 'results/land_cover_rdd_robust_no_gaal.txt'), row.names = FALSE)

stargazer(mod.all1, mod.all2, mod.all3, mod.ec1, mod.ec2, mod.ec3, mod.cm1, mod.cm2, mod.cm3, type = 'text', omit = c('GEOID'), omit.stat = c('f', 'ser'))

stargazer(zmod.all1, zmod.all2, zmod.all3, zmod.ec1, zmod.ec2, zmod.ec3, zmod.cm1, zmod.cm2, zmod.cm3, type = 'text', omit = c('GEOID'), omit.stat = c('f', 'ser'))

stargazer(xmod.all1, xmod.all2, xmod.all3, xmod.ec1, xmod.ec2, xmod.ec3, xmod.cm1, xmod.cm2, xmod.cm3, type = 'text', omit = c('GEOID'), omit.stat = c('f', 'ser'))

stargazer(zxmod.all1, zxmod.all2, zxmod.all3, zxmod.ec1, zxmod.ec2, zxmod.ec3, zxmod.cm1, zxmod.cm2, zxmod.cm3, type = 'text', omit = c('GEOID'), omit.stat = c('f', 'ser'))

# Summary statistics of the data

sumdat <- as.data.frame(cbind(joint$East, joint$Development, joint$Agriculture, joint$Population, joint$Income, joint$Initial_Development, joint$Initial_Agriculture, joint$Initial_Forests, joint$Initial_Water))

sumdat$V4 <- log(sumdat$V4)
sumdat$V5 <- log(sumdat$V5)

colnames(sumdat) <- c('East', 'Developed Land (%)', 'Agricultural Land (%)', 'Log(Population)', 'Log(Income) (2021 USD)', '2001 Developed Land (%)', '2001 Agricultural Land (%)', '2001 Forested Land (%)', '2001 Water (%)')

datasummary_skim(sumdat, fmt = '%.3f')

# Creating leaflets

pal1 <- colorNumeric(palette = c('white', 'red4'), domain = joint$Development)
pal2 <- colorNumeric(palette = c('white', 'green4'), domain = joint$Agriculture)

dev.map <- leaflet(joint$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 1.0, color = 'black', fillColor = pal1(joint$Development))
ag.map <- leaflet(joint$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 1.0, color = 'black', fillColor = pal2(joint$Agriculture))

dev.map
ag.map

# Regional maps for development

EC <- joint %>% filter(STATEFP %in% c('01', '12', '13', '21', '47'))
CM <- joint %>% filter(STATEFP %in% c('20', '31', '46'))

pal3 <- colorNumeric(palette = c('white', 'red4'), domain = EC$Development)
pal4 <- colorNumeric(palette = c('white', 'red4'), domain = CM$Development)

ec.map <- leaflet(EC$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 1.0, color = 'black', fillColor = pal3(EC$Development))
cm.map <- leaflet(CM$geometry) %>% addTiles() %>% addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 1.0, color = 'black', fillColor = pal4(CM$Development))

ec.map
cm.map

