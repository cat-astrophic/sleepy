# This script looks at IRS migration flows and Census wage data to see if time zones affect productivity

# The idea is that people on the good side will have higher wages than those on the wrong side if they commute to the same metro area

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

# Reading in the migration flows data

inflows <- as.data.frame(NULL)
outflows <- as.data.frame(NULL)

yrs <- c('1112', '1213', '1314', '1415', '1516', '1617', '1718', '1819', '1920', '2021')

for (y in yrs) {
  
  tmp.in <- read.csv(paste0(direc, 'data/migration/countyinflow', y, '.csv'))
  tmp.out <- read.csv(paste0(direc, 'data/migration/countyoutflow', y, '.csv'))
  
  colnames(tmp.in) <- c('y2_statefips', 'y2_countyfips', 'y1_statefips', 'y1_countyfips', 'y1_state', 'y1_countyname', 'n1', 'n2', 'agi')
  colnames(tmp.out) <- c('y2_statefips', 'y2_countyfips', 'y1_statefips', 'y1_countyfips', 'y1_state', 'y1_countyname', 'n1', 'n2', 'agi')
  
  tmp.in$Year <- rep(y, nrow(tmp.in))
  tmp.out$Year <- rep(y, nrow(tmp.out))
  
  inflows <- rbind(inflows, tmp.in)
  outflows <- rbind(outflows, tmp.out)
  
}

# Adding FIPS code columns

in.st2 <- ifelse(nchar(inflows$y2_statefips) == 1, '0', '')
in.co2 <- ifelse(nchar(inflows$y2_countyfips) == 2, '0', '')
in.co22 <- ifelse(nchar(inflows$y2_countyfips) == 1, '00', '')

in.st1 <- ifelse(nchar(inflows$y1_statefips) == 1, '0', '')
in.co1 <- ifelse(nchar(inflows$y1_countyfips) == 2, '0', '')
in.co12 <- ifelse(nchar(inflows$y1_countyfips) == 1, '00', '')

out.st2 <- ifelse(nchar(outflows$y2_statefips) == 1, '0', '')
out.co2 <- ifelse(nchar(outflows$y2_countyfips) == 2, '0', '')
out.co22 <- ifelse(nchar(outflows$y2_countyfips) == 1, '00', '')

out.st1 <- ifelse(nchar(outflows$y1_statefips) == 1, '0', '')
out.co1 <- ifelse(nchar(outflows$y1_countyfips) == 2, '0', '')
out.co12 <- ifelse(nchar(outflows$y1_countyfips) == 1, '00', '')

inflows$FIPS_New <- paste0(in.st2, inflows$y2_statefips, in.co2, in.co22, inflows$y2_countyfips)
inflows$FIPS_Old <- paste0(in.st1, inflows$y1_statefips, in.co1, in.co12, inflows$y1_countyfips)

outflows$FIPS_New <- paste0(out.st2, outflows$y2_statefips, out.co2, out.co22, outflows$y2_countyfips)
outflows$FIPS_Old <- paste0(out.st1, outflows$y1_statefips, out.co1, out.co12, outflows$y1_countyfips)

# Subset for desired states

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

all.fips <- all.fips %>% filter(Flag == 1)

all.fips$check <- all.fips$f1 + all.fips$f2 + all.fips$f3 + all.fips$f4

# Creating distance data

dists <- c()

for (i in 1:nrow(all.fips)) {
  
  dists <- c(dists, min(all.fips$d1[i], all.fips$d2[i], all.fips$d3[i], all.fips$d4[i]))
  
}

all.fips$Distance <- dists / 1000

# Creating a full data set for pairwise flows

sources <- c()
sinks <- c()
years <- c()

for (i in 1:length(yrs)) {
  
  tmp1 <- all.fips %>% filter(f1 == 1)
  tmp2 <- all.fips %>% filter(f2 == 1)
  tmp3 <- all.fips %>% filter(f3 == 1)
  tmp4 <- all.fips %>% filter(f4 == 1)
  
  for (j in 1:nrow(tmp1)) {
    
    sources <- c(sources, rep(tmp1$GEOID[j], nrow(tmp1) - 1))
    sinks <- c(sinks, tmp1$GEOID[-j])
    
  }
  
  for (j in 1:nrow(tmp2)) {
    
    sources <- c(sources, rep(tmp2$GEOID[j], nrow(tmp2) - 1))
    sinks <- c(sinks, tmp2$GEOID[-j])
    
  }
  
  for (j in 1:nrow(tmp3)) {
    
    sources <- c(sources, rep(tmp3$GEOID[j], nrow(tmp3) - 1))
    sinks <- c(sinks, tmp3$GEOID[-j])
    
  }
  
  for (j in 1:nrow(tmp4)) {
    
    sources <- c(sources, rep(tmp4$GEOID[j], nrow(tmp4) - 1))
    sinks <- c(sinks, tmp4$GEOID[-j])
    
  }
  
  years <- c(years, rep(yrs[i], length(sources) - length(years)))
  
}

df <- as.data.frame(cbind(years, sources, sinks))
colnames(df) <- c('Year', 'Source', 'Sink')

net.flows1 <- c()
net.flows2 <- c()

for (i in 1:nrow(df)) {
  
  tmp.in <- inflows %>% filter(FIPS_New == df$Sink[i]) %>% filter(FIPS_Old == df$Source[i]) %>% filter(Year == df$Year[i])
  tmp.out <- inflows %>% filter(FIPS_Old == df$Sink[i]) %>% filter(FIPS_New == df$Source[i]) %>% filter(Year == df$Year[i])
  
  net.flows1 <- c(net.flows1, tmp.in$n1[1] - tmp.out$n1[1])
  net.flows2 <- c(net.flows2, tmp.in$n2[1] - tmp.out$n2[1])
  
}

df$Net.Flow.1 <- net.flows1
df$Net.Flow.2 <- net.flows2

df[is.na(df)] <- 0

# Creating a full data set for total net flows

flow.fips <- rep(unique(all.fips$GEOID), length(yrs))
flow.yrs <- c()

for (i in 1:length(yrs)) {
  
  flow.yrs <- c(flow.yrs, rep(yrs[i], length(unique(all.fips$GEOID))))
  
}


df2 <- as.data.frame(cbind(flow.yrs, flow.fips))
colnames(df2) <- c('Year', 'FIPS')

fips.net.flows1 <- c()
fips.net.flows2 <- c()
fips.inflow1 <- c()
fips.inflow2 <- c()

for (i in 1:nrow(df2)) {
  
  tmp.in <- inflows %>% filter(FIPS_New == df2$FIPS[i]) %>% filter(Year == df2$Year[i])
  tmp.out <- inflows %>% filter(FIPS_Old == df2$FIPS[i]) %>% filter(Year == df2$Year[i])
  
  tmp.in <- tmp.in %>% filter(FIPS_Old == '96000')
  tmp.out <- tmp.out %>% filter(FIPS_New != FIPS_Old)
  
  fips.net.flows1 <- c(fips.net.flows1, sum(tmp.in$n1) - sum(tmp.out$n1))
  fips.net.flows2 <- c(fips.net.flows2, sum(tmp.in$n2) - sum(tmp.out$n2))
  
  fips.inflow1 <- c(fips.inflow1, sum(tmp.in$n1))
  fips.inflow2 <- c(fips.inflow2, sum(tmp.in$n2))
  
}

df2$Net.Flow.1 <- fips.net.flows1
df2$Net.Flow.2 <- fips.net.flows2

df2$Inflow.1 <- fips.inflow1
df2$Inflow.2 <- fips.inflow2

# Get population data for FIPS

fips.pops <- c()

for (i in 1:length(unique(all.fips$GEOID))) {
  
  x <- unique(all.fips$GEOID)[i]
  s <- substr(x, 1, 2)
  c <- substr(x, 3, 5)
  
  tmp <- create_tract_table(state = s, county = c, year = 2010)
  fips.pops <- c(fips.pops, sum(tmp$pop))
  
}

pop.df <- as.data.frame(cbind(unique(all.fips$GEOID), fips.pops))
colnames(pop.df) <- c('FIPS', 'POP')

pops1 <- c()

for (i in 1:nrow(df)) {
  
  tmp <- pop.df %>% filter(FIPS == df$Sink[i])
  pops1 <- c(pops1, tmp$POP)
  
}

pops2 <- c()

for (i in 1:nrow(df2)) {
  
  tmp <- pop.df %>% filter(FIPS == df2$FIPS[i])
  pops2 <- c(pops2, tmp$POP)
  
}

df$Population <- as.integer(pops1)
df2$Population <- as.integer(pops2)

# Normalize the Net.Flow variables in df and df2

df$Y1 <- df$Net.Flow.1 / df$Population
df$Y2 <- df$Net.Flow.2 / df$Population

df2$Y1 <- df2$Net.Flow.1 / df2$Population
df2$Y2 <- df2$Net.Flow.2 / df2$Population

df2$In1 <- df2$Inflow.1 / df2$Population
df2$In2 <- df2$Inflow.2 / df2$Population

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

# Adding the time zone variable to df and df2

df.tz <- c()

for (i in 1:nrow(df)) {
  
  df.tz <- c(df.tz, all.fips[which(all.fips$GEOID == df$Sink[i]),]$TZ[1])

}

df2.tz <- c()

for (i in 1:nrow(df2)) {
  
  df2.tz <- c(df2.tz, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$TZ[1])

}

df$TZ <- df.tz
df2$TZ <- df2.tz

# Add the distance variable to df and df2

distf <- c()

for (i in 1:nrow(df)) {
  
  distf <- c(distf, all.fips[which(all.fips$GEOID == df$Sink[i]),]$Distance[1])
  
}

distf2 <- c()

for (i in 1:nrow(df2)) {
  
  distf2 <- c(distf2, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$Distance[1])
  
}

df$Distance <- distf
df2$Distance <- distf2

# Add metro area indicators to df and df2

i1 <- c()
i2 <- c()
i3 <- c()
i4 <- c()

for (i in 1:nrow(df)) {
  
  i1 <- c(i1, all.fips[which(all.fips$GEOID == df$Sink[i]),]$f1[1])
  i2 <- c(i2, all.fips[which(all.fips$GEOID == df$Sink[i]),]$f2[1])
  i3 <- c(i3, all.fips[which(all.fips$GEOID == df$Sink[i]),]$f3[1])
  i4 <- c(i4, all.fips[which(all.fips$GEOID == df$Sink[i]),]$f4[1])
  
}

j1 <- c()
j2 <- c()
j3 <- c()
j4 <- c()

for (i in 1:nrow(df2)) {
  
  j1 <- c(j1, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$f1[1])
  j2 <- c(j2, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$f2[1])
  j3 <- c(j3, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$f3[1])
  j4 <- c(j4, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$f4[1])
  
}

df$M1 <- i1
df$M2 <- i2
df$M3 <- i3
df$M4 <- i4

df2$M1 <- j1
df2$M2 <- j2
df2$M3 <- j3
df2$M4 <- j4

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

df <- df %>% filter(! Sink %in% too.close.atl) %>% filter(! Sink %in% too.close.ash) %>% filter(! Sink %in% too.close.nash) %>% filter(! Sink %in% too.close.mon)

df2 <- df2 %>% filter(! FIPS %in% too.close.atl) %>% filter(! FIPS %in% too.close.ash) %>% filter(! FIPS %in% too.close.nash) %>% filter(! FIPS %in% too.close.mon)

# Creating a State factor for df and df2

s1 <- c()

for (i in 1:nrow(df)) {
  
  s1 <- c(s1, all.fips[which(all.fips$GEOID == df$Sink[i]),]$STATEFP[1])
  
}

s2 <- c()

for (i in 1:nrow(df2)) {
  
  s2 <- c(s2, all.fips[which(all.fips$GEOID == df2$FIPS[i]),]$STATEFP[1])
  
}

df$State <- s1
df2$State <- s2

# Creating a county-to-county distance variable for addressing moving costs

pw.dists <- c()

for (i in 1:nrow(df)) {
  
  a <- all.fips[which(all.fips$GEOID == df$Sink[i]),]$centroid
  b <- all.fips[which(all.fips$GEOID == df$Source[i]),]$centroid
  
  pw.dists <- c(pw.dists, st_distance(a, b))
  
}

df$Travel <- pw.dists / 100000

# See if people prefer moving to / living in one side of the time zone boundary - flows are pairwise across sampled counties only

pw1 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M3 == 1),])
pw2 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1),])
pw3 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M2 == 1),])
pw4 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M4 == 1),])
pw5 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1 | df$M3 == 1),])
pw6 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df)
pw7 <- lm(Y1 ~ TZ + Travel + M3 + M1 + M2 + M4 + factor(Sink) + factor(Source) + factor(Year), data = df)

ppw1 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M3 == 1),])
ppw2 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1),])
ppw3 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M2 == 1),])
ppw4 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M4 == 1),])
ppw5 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1 | df$M3 == 1),])
ppw6 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df)
ppw7 <- lm(Y2 ~ TZ + Travel + M3 + M1 + M2 + M4 + factor(Sink) + factor(Source) + factor(Year), data = df)

pw1x <- coeftest(pw1, vcov. = vcovCL(pw1, type = 'HC1'))
pw2x <- coeftest(pw2, vcov. = vcovCL(pw2, type = 'HC1'))
pw3x <- coeftest(pw3, vcov. = vcovCL(pw3, type = 'HC1'))
pw4x <- coeftest(pw4, vcov. = vcovCL(pw4, type = 'HC1'))
pw5x <- coeftest(pw5, vcov. = vcovCL(pw5, type = 'HC1'))
pw6x <- coeftest(pw6, vcov. = vcovCL(pw6, type = 'HC1'))
pw7x <- coeftest(pw7, vcov. = vcovCL(pw7, type = 'HC1'))

ppw1x <- coeftest(ppw1, vcov. = vcovCL(ppw1, type = 'HC1'))
ppw2x <- coeftest(ppw2, vcov. = vcovCL(ppw2, type = 'HC1'))
ppw3x <- coeftest(ppw3, vcov. = vcovCL(ppw3, type = 'HC1'))
ppw4x <- coeftest(ppw4, vcov. = vcovCL(ppw4, type = 'HC1'))
ppw5x <- coeftest(ppw5, vcov. = vcovCL(ppw5, type = 'HC1'))
ppw6x <- coeftest(ppw6, vcov. = vcovCL(ppw6, type = 'HC1'))
ppw7x <- coeftest(ppw7, vcov. = vcovCL(ppw7, type = 'HC1'))

# Repeat with zeros removed for shits and giggles and shitty giggles

df0 <- df %>% filter(Y1 != 0)

sg1 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M3 == 1),])
sg2 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M1 == 1),])
sg3 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M2 == 1),])
sg4 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M4 == 1),])
sg5 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M1 == 1 | df0$M3 == 1),])
sg6 <- lm(Y1 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0)
sg7 <- lm(Y1 ~ TZ + Travel + M3 + M1 + M2 + M4 + factor(Sink) + factor(Source) + factor(Year), data = df0)

ssg1 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M3 == 1),])
ssg2 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M1 == 1),])
ssg3 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M2 == 1),])
ssg4 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M4 == 1),])
ssg5 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df0[which(df0$M1 == 1 | df0$M3 == 1),])
ssg6 <- lm(Y2 ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + M3 + M1 + M2 + M4 + factor(Year), data = df0)
ssg7 <- lm(Y2 ~ TZ + Travel + M3 + M1 + M2 + M4 + factor(Sink) + factor(Source) + factor(Year), data = df0)

sg1x <- coeftest(sg1, vcov. = vcovCL(sg1, type = 'HC1'))
sg2x <- coeftest(sg2, vcov. = vcovCL(sg2, type = 'HC1'))
sg3x <- coeftest(sg3, vcov. = vcovCL(sg3, type = 'HC1'))
sg4x <- coeftest(sg4, vcov. = vcovCL(sg4, type = 'HC1'))
sg5x <- coeftest(sg5, vcov. = vcovCL(sg5, type = 'HC1'))
sg6x <- coeftest(sg6, vcov. = vcovCL(sg6, type = 'HC1'))
sg7x <- coeftest(sg7, vcov. = vcovCL(sg7, type = 'HC1'))

ssg1x <- coeftest(ssg1, vcov. = vcovCL(ssg1, type = 'HC1'))
ssg2x <- coeftest(ssg2, vcov. = vcovCL(ssg2, type = 'HC1'))
ssg3x <- coeftest(ssg3, vcov. = vcovCL(ssg3, type = 'HC1'))
ssg4x <- coeftest(ssg4, vcov. = vcovCL(ssg4, type = 'HC1'))
ssg5x <- coeftest(ssg5, vcov. = vcovCL(ssg5, type = 'HC1'))
ssg6x <- coeftest(ssg6, vcov. = vcovCL(ssg6, type = 'HC1'))
ssg7x <- coeftest(ssg7, vcov. = vcovCL(ssg7, type = 'HC1'))

# See if people prefer moving to / living in one side of the time zone boundary - total inflows

in1 <- lm(In1 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M3 == 1),])
in2 <- lm(In1 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1),])
in3 <- lm(In1 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M2 == 1),])
in4 <- lm(In1 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M4 == 1),])
in5 <- lm(In1 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1 | df2$M3 == 1),])
in6 <- lm(In1 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2)
in7 <- lm(In1 ~ TZ + M3 + M1 + M2 + M4 + factor(FIPS) + factor(Year), data = df2)

iin1 <- lm(In2 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M3 == 1),])
iin2 <- lm(In2 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1),])
iin3 <- lm(In2 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M2 == 1),])
iin4 <- lm(In2 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M4 == 1),])
iin5 <- lm(In2 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1 | df2$M3 == 1),])
iin6 <- lm(In2 ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2)
iin7 <- lm(In2 ~ TZ + M3 + M1 + M2 + M4 + factor(FIPS) + factor(Year), data = df2)

in1x <- coeftest(in1, vcov. = vcovCL(in1, type = 'HC1'))
in2x <- coeftest(in2, vcov. = vcovCL(in2, type = 'HC1'))
in3x <- coeftest(in3, vcov. = vcovCL(in3, type = 'HC1'))
in4x <- coeftest(in4, vcov. = vcovCL(in4, type = 'HC1'))
in5x <- coeftest(in5, vcov. = vcovCL(in5, type = 'HC1'))
in6x <- coeftest(in6, vcov. = vcovCL(in6, type = 'HC1'))
in7x <- coeftest(in7, vcov. = vcovCL(in7, type = 'HC1'))

iin1x <- coeftest(iin1, vcov. = vcovCL(iin1, type = 'HC1'))
iin2x <- coeftest(iin2, vcov. = vcovCL(iin2, type = 'HC1'))
iin3x <- coeftest(iin3, vcov. = vcovCL(iin3, type = 'HC1'))
iin4x <- coeftest(iin4, vcov. = vcovCL(iin4, type = 'HC1'))
iin5x <- coeftest(iin5, vcov. = vcovCL(iin5, type = 'HC1'))
iin6x <- coeftest(iin6, vcov. = vcovCL(iin6, type = 'HC1'))
iin7x <- coeftest(iin7, vcov. = vcovCL(iin7, type = 'HC1'))

# Okay but what if I log it?

lin1 <- lm(log(In1+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M3 == 1),])
lin2 <- lm(log(In1+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1),])
lin3 <- lm(log(In1+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M2 == 1),])
lin4 <- lm(log(In1+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M4 == 1),])
lin5 <- lm(log(In1+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1 | df2$M3 == 1),])
lin6 <- lm(log(In1+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2)
lin7 <- lm(log(In1+1) ~ TZ + M3 + M1 + M2 + M4 + factor(FIPS) + factor(Year), data = df2)

liin1 <- lm(log(In2+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M3 == 1),])
liin2 <- lm(log(In2+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1),])
liin3 <- lm(log(In2+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M2 == 1),])
liin4 <- lm(log(In2+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M4 == 1),])
liin5 <- lm(log(In2+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2[which(df2$M1 == 1 | df2$M3 == 1),])
liin6 <- lm(log(In2+1) ~ TZ + factor(State) + factor(FIPS) + factor(Year), data = df2)
liin7 <- lm(log(In2+1) ~ TZ + M3 + M1 + M2 + M4 + factor(FIPS) + factor(Year), data = df2)

lin1x <- coeftest(lin1, vcov. = vcovCL(lin1, type = 'HC1'))
lin2x <- coeftest(lin2, vcov. = vcovCL(lin2, type = 'HC1'))
lin3x <- coeftest(lin3, vcov. = vcovCL(lin3, type = 'HC1'))
lin4x <- coeftest(lin4, vcov. = vcovCL(lin4, type = 'HC1'))
lin5x <- coeftest(lin5, vcov. = vcovCL(lin5, type = 'HC1'))
lin6x <- coeftest(lin6, vcov. = vcovCL(lin6, type = 'HC1'))
lin7x <- coeftest(lin7, vcov. = vcovCL(lin7, type = 'HC1'))

liin1x <- coeftest(liin1, vcov. = vcovCL(liin1, type = 'HC1'))
liin2x <- coeftest(liin2, vcov. = vcovCL(liin2, type = 'HC1'))
liin3x <- coeftest(liin3, vcov. = vcovCL(liin3, type = 'HC1'))
liin4x <- coeftest(liin4, vcov. = vcovCL(liin4, type = 'HC1'))
liin5x <- coeftest(liin5, vcov. = vcovCL(liin5, type = 'HC1'))
liin6x <- coeftest(liin6, vcov. = vcovCL(liin6, type = 'HC1'))
liin7x <- coeftest(liin7, vcov. = vcovCL(liin7, type = 'HC1'))

# You know what, let's log these too; I'm [on a] hot [streak]

lpw1 <- lm(log(Y1+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M3 == 1),])
lpw2 <- lm(log(Y1+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1),])
lpw3 <- lm(log(Y1+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M2 == 1),])
lpw4 <- lm(log(Y1+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M4 == 1),])
lpw5 <- lm(log(Y1+1) ~ TZ + Travel + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1 | df$M3 == 1),])
lpw6 <- lm(log(Y1+1) ~ TZ + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df)
lpw7 <- lm(log(Y1+1) ~ TZ + Travel + M3 + M1 + M2 + M4 + factor(Sink) + factor(Source) + factor(Year), data = df)

lppw1 <- lm(log(Y2+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M3 == 1),])
lppw2 <- lm(log(Y2+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1),])
lppw3 <- lm(log(Y2+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M2 == 1),])
lppw4 <- lm(log(Y2+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M4 == 1),])
lppw5 <- lm(log(Y2+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df[which(df$M1 == 1 | df$M3 == 1),])
lppw6 <- lm(log(Y2+1) ~ TZ + Travel + factor(State) + factor(Sink) + factor(Source) + factor(Year), data = df)
lppw7 <- lm(log(Y2+1) ~ TZ + Travel + M3 + M1 + M2 + M4 + factor(Sink) + factor(Source) + factor(Year), data = df)

lpw1x <- coeftest(lpw1, vcov. = vcovCL(lpw1, type = 'HC1'))
lpw2x <- coeftest(lpw2, vcov. = vcovCL(lpw2, type = 'HC1'))
lpw3x <- coeftest(lpw3, vcov. = vcovCL(lpw3, type = 'HC1'))
lpw4x <- coeftest(lpw4, vcov. = vcovCL(lpw4, type = 'HC1'))
lpw5x <- coeftest(lpw5, vcov. = vcovCL(lpw5, type = 'HC1'))
lpw6x <- coeftest(lpw6, vcov. = vcovCL(lpw6, type = 'HC1'))
lpw7x <- coeftest(lpw7, vcov. = vcovCL(lpw7, type = 'HC1'))

lppw1x <- coeftest(lppw1, vcov. = vcovCL(lppw1, type = 'HC1'))
lppw2x <- coeftest(lppw2, vcov. = vcovCL(lppw2, type = 'HC1'))
lppw3x <- coeftest(lppw3, vcov. = vcovCL(lppw3, type = 'HC1'))
lppw4x <- coeftest(lppw4, vcov. = vcovCL(lppw4, type = 'HC1'))
lppw5x <- coeftest(lppw5, vcov. = vcovCL(lppw5, type = 'HC1'))
lppw6x <- coeftest(lppw6, vcov. = vcovCL(lppw6, type = 'HC1'))
lppw7x <- coeftest(lppw7, vcov. = vcovCL(lppw7, type = 'HC1'))

# Viewing results

stargazer(pw1x, pw2x, pw3x, pw4x, pw5x, pw6x, pw7x, type = 'text', omit = c('Sink', 'Source', 'Year', 'State'))
stargazer(ppw1x, ppw2x, ppw3x, ppw4x, ppw5x, ppw6x, ppw7x, type = 'text', omit = c('Sink', 'Source', 'Year', 'State'))

stargazer(sg1x, sg2x, sg3x, sg4x, sg5x, sg6x, sg7x, type = 'text', omit = c('Sink', 'Source', 'Year', 'State'))
stargazer(ssg1x, ssg2x, ssg3x, ssg4x, ssg5x, ssg6x, ssg7x, type = 'text', omit = c('Sink', 'Source', 'Year', 'State'))

stargazer(in1x, in2x, in3x, in4x, in5x, in6x, in7x, type = 'text', omit = c('FIPS', 'Year', 'State'))
stargazer(iin1x, iin2x, iin3x, iin4x, iin5x, iin6x, iin7x, type = 'text', omit = c('FIPS', 'Year', 'State'))

stargazer(lin1x, lin2x, lin3x, lin4x, lin5x, lin6x, lin7x, type = 'text', omit = c('FIPS', 'Year', 'State'))
stargazer(liin1x, liin2x, liin3x, liin4x, liin5x, liin6x, liin7x, type = 'text', omit = c('FIPS', 'Year', 'State'))

stargazer(lpw1x, lpw2x, lpw3x, lpw4x, lpw5x, lpw6x, lpw7x, type = 'text', omit = c('Sink', 'Source', 'Year', 'State'))
stargazer(lppw1x, lppw2x, lppw3x, lppw4x, lppw5x, lppw6x, lppw7x, type = 'text', omit = c('Sink', 'Source', 'Year', 'State'))

# Saving results

write.csv(stargazer(pw1x, pw2x, pw3x, pw4x, pw5x, pw6x, pw7x, omit = c('Sink', 'Source', 'Year', 'State')), paste0(direc, 'results/net_flows_filers.txt'), row.names = FALSE)
write.csv(stargazer(ppw1x, ppw2x, ppw3x, ppw4x, ppw5x, ppw6x, ppw7x, omit = c('Sink', 'Source', 'Year', 'State')), paste0(direc, 'results/net_flows_all.txt'), row.names = FALSE)

write.csv(stargazer(in1x, in2x, in3x, in4x, in5x, in6x, in7x, omit = c('FIPS', 'Year', 'State')), paste0(direc, 'results/inflows_filers.txt'), row.names = FALSE)
write.csv(stargazer(iin1x, iin2x, iin3x, iin4x, iin5x, iin6x, iin7x, omit = c('FIPS', 'Year', 'State')), paste0(direc, 'results/inflows_all.txt'), row.names = FALSE)

# Making a plot where counties are colored by time zone and MSA

mapdat <- all.fips %>% filter(! GEOID %in% too.close.atl) %>% filter(! GEOID %in% too.close.ash) %>% filter(! GEOID %in% too.close.nash) %>% filter(! GEOID %in% too.close.mon)

mapdat$MSA <- ifelse(mapdat$GEOID %in% msas, 1, 0)

color <- c()

for (i in 1:nrow(mapdat)) {
  
  if (mapdat$MSA[i] == 1) {
    
    color <- c(color, 'black')
    
  } else if (mapdat$TZ[i] == 1) {
    
    color <- c(color, 'red')
    
  } else {
    
    color <- c(color, 'orange')
    
  }
  
}

mapdat$Color <- color

map <- leaflet(mapdat$geometry) %>% setView(lng = -85, lat = 33, zoom = 6) %>% addTiles() %>% 
               addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = 0.666, color = 'black', fillColor = mapdat$Color)

map

# Summary statistics of the data

sumdat <- as.data.frame(cbind(df$Net.Flow.1, df$Net.Flow.2, df$Y1, df$Y2, df$TZ, df$M3, df$M1, df$M2, df$M4))
sumdat2 <- as.data.frame(cbind(df2$Inflow.1, df2$Inflow.2, df2$In1, df2$In2, df2$TZ, df2$M3, df2$M1, df2$M2, df2$M4))

colnames(sumdat) <- c('Net Flow (Filers)', 'Net Flow (All)', 'Net Flow per capita (Filers)', 'Net Flow per capita (All)', 'Eastern Time Zone',
                      'Knoxville, TN Metro Area', 'Chattanooga, TN Metro Area', 'Columbus, GA Metro Area', 'Tallahassee, FL Metro Area')

colnames(sumdat2) <- c('Inflow (Filers)', 'Inflow (All)', 'Inflow per capita (Filers)', 'Inflow per capita (All)', 'Eastern Time Zone',
                       'Knoxville, TN Metro Area', 'Chattanooga, TN Metro Area', 'Columbus, GA Metro Area', 'Tallahassee, FL Metro Area')

datasummary_skim(sumdat, fmt = '%.3f')
datasummary_skim(sumdat2, fmt = '%.3f')

# Plotting the pairwise travel distances between county centroids from the pairwise net flows models

df$XTravel <- df$Travel * 100000 / 1000
  
ggplot(data = df, aes(x = XTravel)) +
  geom_histogram(color = 'red4', fill = 'orange', bins = 100) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Histogram of the Distance Between County Centroids within the Same Metro Area') + 
  xlab('Distance (kilometers)') + 
  ylab('')

ggplot(data = df[which(df$Year == '2021'),], aes(x = XTravel)) +
  geom_histogram(color = 'red4', fill = 'orange', bins = 100) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Histogram of the Distance Between County Centroids within the Same Metro Area') + 
  xlab('Distance (kilometers)') + 
  ylab('')

