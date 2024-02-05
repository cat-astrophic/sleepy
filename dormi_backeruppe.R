# This script looks at whether or not growth occurred in the number of businesses in Mercer County post time zone change 

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

# Reading in the CBP data (and subsetting for ND only)

cbp <- as.data.frame(NULL)

yrs <- c('00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
         '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21')

for (y in yrs) {
  
  print(paste0('Reading in data for year 20', y, '.......'))
  
  tmp <- read.csv(paste0(direc, 'data/cbp/cbp', y, 'co.txt'))
  
  colnames(tmp) <- tolower(colnames(tmp))
  tmp <- tmp %>% filter(fipstate == 38)
  tmp$Year <- rep(y, nrow(tmp))
  tmp <- tmp[,which(! colnames(tmp) %in% c('emp_nf', 'qp1_nf', 'ap_nf')),]
  colnames(tmp)[which(colnames(tmp) == 'n.5')] <- 'n1_4'
  
  if (ncol(tmp) == 23) {
    
    tmp <- cbind(tmp[,1:3], rep(NA, nrow(tmp)), tmp[,4:23])
    colnames(tmp)[4] <- 'empflag'
    
  }
  
  if (ncol(tmp) == 21) {
    
    tmp <- cbind(tmp[,1:3], rep(NA, nrow(tmp)), tmp[,4:20], rep(NA, nrow(tmp)), rep(NA, nrow(tmp)), tmp$Year)
    colnames(tmp)[4] <- 'empflag'
    colnames(tmp)[22] <- 'censtate'
    colnames(tmp)[23] <- 'cencty'
    colnames(tmp)[24] <- 'Year'
    
  }
  
  cbp <- rbind(cbp, tmp)
  
}

# Determining highest level naics category entries

cbp$tots <- as.integer(grepl('------', cbp$naics))
cbp$flag <- as.integer(grepl('----', cbp$naics))
cbp$high <- cbp$flag - cbp$tots

cbp <- cbp %>% filter(flag == 1)

# Creating a usable data set

fips <- c()
naics <- c()
year <- c()
jobs <- c()

for (y in unique(cbp$Year)) {
  
  print(paste0('Creating data for 20', y, '.......'))
  
  tmp.y <- cbp %>% filter(Year == y)
  
  for (n in unique(cbp$naics)) {
    
    tmp.n <- tmp.y %>% filter(naics == n)
    
    for (f in unique(cbp$fipscty)) {
      
      tmp <- tmp.n %>% filter(fipscty == f)
      
      fips <- c(fips, f)
      naics <- c(naics, n)
      year <- c(year, y)
      jobs <- c(jobs, tmp$emp[1])
      
    }
    
  }
  
}

df <- as.data.frame(cbind(fips, naics, year, jobs))
colnames(df) <- c('FIPS', 'NAICS', 'Year', 'Jobs')
df <- df %>% filter(FIPS != 999)
df$Year <- as.integer(df$Year)
df$Jobs <- as.integer(df$Jobs)

df[is.na(df)] <- 0

# Adding diff-in-diff variables

df$Post <- as.integer(df$Year > 10)
df$Treated <- as.integer(df$FIPS == 57)

# Creating an HHI variable

hhidf <- df %>% filter(NAICS == '------')

hhi <- c()

for (i in 1:nrow(hhidf)) {
  
  tmp <- df %>% filter(FIPS == hhidf$FIPS[i]) %>% filter(Year == hhidf$Year[i])
  
  shares <- c()
  
  for (j in 1:nrow(tmp)) {shares <- c(shares, (100 * tmp$Jobs[j] / sum(tmp$Jobs))^2)}
  
  hhi <- c(hhi, sum(shares) / 10000)
  
}

hhidf$HHI <- hhi

# Running regressions

mod.hhi <- lm(HHI ~ Treated*Post + factor(Year) + factor(FIPS), data = hhidf)
mod.all <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '------'),])
mod.combo <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = df[which(df$NAICS != '------'),])
mod.11 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '11----'),])
mod.21 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '21----'),])
mod.22 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '22----'),])
mod.23 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '23----'),])
mod.31 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '31----'),])
mod.42 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '42----'),])
mod.44 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '44----'),])
mod.48 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '48----'),])
mod.51 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '51----'),])
mod.52 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '52----'),])
mod.53 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '53----'),])
mod.54 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '54----'),])
mod.55 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '55----'),])
mod.56 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '56----'),])
mod.61 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '61----'),])
mod.62 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '62----'),])
mod.71 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '71----'),])
mod.72 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '72----'),])
mod.81 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '81----'),])
mod.95 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '95----'),])
mod.99 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS), data = df[which(df$NAICS == '99----'),])

mod.hhix <- coeftest(mod.hhi, vcov = vcovCL, cluster = ~FIPS)
mod.allx <- coeftest(mod.all, vcov. = vcovCL, cluster = ~FIPS)
mod.combox <- coeftest(mod.combo, vcov. = vcovCL, cluster = ~FIPS)
mod.11x <- coeftest(mod.11, vcov. = vcovCL, cluster = ~FIPS)
mod.21x <- coeftest(mod.21, vcov. = vcovCL, cluster = ~FIPS)
mod.22x <- coeftest(mod.22, vcov. = vcovCL, cluster = ~FIPS)
mod.23x <- coeftest(mod.23, vcov. = vcovCL, cluster = ~FIPS)
mod.31x <- coeftest(mod.31, vcov. = vcovCL, cluster = ~FIPS)
mod.42x <- coeftest(mod.42, vcov. = vcovCL, cluster = ~FIPS)
mod.44x <- coeftest(mod.44, vcov. = vcovCL, cluster = ~FIPS)
mod.48x <- coeftest(mod.48, vcov. = vcovCL, cluster = ~FIPS)
mod.51x <- coeftest(mod.51, vcov. = vcovCL, cluster = ~FIPS)
mod.52x <- coeftest(mod.52, vcov. = vcovCL, cluster = ~FIPS)
mod.53x <- coeftest(mod.53, vcov. = vcovCL, cluster = ~FIPS)
mod.54x <- coeftest(mod.54, vcov. = vcovCL, cluster = ~FIPS)
mod.55x <- coeftest(mod.55, vcov. = vcovCL, cluster = ~FIPS)
mod.56x <- coeftest(mod.56, vcov. = vcovCL, cluster = ~FIPS)
mod.61x <- coeftest(mod.61, vcov. = vcovCL, cluster = ~FIPS)
mod.62x <- coeftest(mod.62, vcov. = vcovCL, cluster = ~FIPS)
mod.71x <- coeftest(mod.71, vcov. = vcovCL, cluster = ~FIPS)
mod.72x <- coeftest(mod.72, vcov. = vcovCL, cluster = ~FIPS)
mod.81x <- coeftest(mod.81, vcov. = vcovCL, cluster = ~FIPS)
mod.95x <- coeftest(mod.95, vcov. = vcovCL, cluster = ~FIPS)
mod.99x <- coeftest(mod.99, vcov. = vcovCL, cluster = ~FIPS)

# Viewing results

stargazer(mod.hhi, mod.all, mod.combo, mod.11, mod.21, mod.22, mod.23, mod.31, mod.42, mod.44, mod.48, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))
stargazer(mod.51, mod.52, mod.53, mod.54, mod.55, mod.56, mod.61, mod.62, mod.71, mod.72, mod.81, mod.95, mod.99, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))

stargazer(mod.hhix, mod.allx, mod.combox, mod.11x, mod.21x, mod.22x, mod.23x, mod.31x, mod.42x, mod.44x, mod.48x, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))
stargazer(mod.51x, mod.52x, mod.53x, mod.54x, mod.55x, mod.56x, mod.61x, mod.62x, mod.71x, mod.72x, mod.81x, mod.95x, mod.99x, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))


