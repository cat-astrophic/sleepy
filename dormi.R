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

# Adding indicators for time zone status to subset for more nuanced control groups in regressions

mtn <- c()
cent <- c()

mountain.counties <- c(1, 7, 11, 25, 33, 37, 41, 53, 57, 85, 87, 89)

for (i in 1:nrow(df)) {
  
  if (df$FIPS[i] %in% mountain.counties) {
    
    mtn <- c(mtn, 1)
    
  } else {
    
    mtn <- c(mtn, 0)
    
  }
  
  if (! df$FIPS[i] %in% mountain.counties) {
    
    cent <- c(cent, 1)
    
  } else if (df$FIPS[i] == 57) {
    
    cent <- c(cent, 1)
    
  } else {
    
    cent <- c(cent, 0)
    
  }
  
}

df$Mountain <- mtn
df$Central <- cent

# Removing the total jobs category

dfx <- df[which(df$NAICS != '------'),]

# Removing any counties that had zeros

xxx <- dfx %>% group_by(FIPS, Year) %>% summarise(Jobs = sum(Jobs))
snoopy <- which(xxx$Jobs == 0)
worm <- xxx[snoopy,]
droop <- unique(worm$FIPS)
dfx <- dfx %>% filter(! FIPS %in% droop)

# Running regressions

mod1 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfx)
mod2 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfx[which(dfx$Central == 1),])
mod3 <- lm(log(Jobs + 1) ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfx[which(dfx$Mountain == 1),])

mod1x <- coeftest(mod1, vcov = vcovCL, cluster = ~FIPS)
mod2x <- coeftest(mod2, vcov = vcovCL, cluster = ~FIPS)
mod3x <- coeftest(mod3, vcov = vcovCL, cluster = ~FIPS)

# Viewing results

stargazer(mod1, mod2, mod3, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))
stargazer(mod1x, mod2x, mod3x, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))

# Making a data.frame of changes in jobs

dfxx <- dfx %>% filter(Year > 0)

change.jobs <- c()

for (i in 1:nrow(dfxx)) {
  
  tmpa <- dfx %>% filter(FIPS == dfxx$FIPS[i]) %>% filter(NAICS == dfxx$NAICS[i]) %>% filter(Year == dfxx$Year[i] - 1)
  tmpb <- dfx %>% filter(FIPS == dfxx$FIPS[i]) %>% filter(NAICS == dfxx$NAICS[i]) %>% filter(Year == dfxx$Year[i])
  change.jobs <- c(change.jobs, tmpb$Jobs[1] - tmpa$Jobs[1])
  
}

dfxx$Change <- change.jobs

# Running regressions for change in jobs

xmod1 <- lm(Change ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfxx)
xmod2 <- lm(Change ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfxx[which(dfxx$Central == 1),])
xmod3 <- lm(Change ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfxx[which(dfxx$Mountain == 1),])

xmod1x <- coeftest(xmod1, vcov = vcovCL, cluster = ~FIPS)
xmod2x <- coeftest(xmod2, vcov = vcovCL, cluster = ~FIPS)
xmod3x <- coeftest(xmod3, vcov = vcovCL, cluster = ~FIPS)

# Viewing results

stargazer(mod1, mod2, mod3, xmod1, xmod2, xmod3, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))
stargazer(mod1x, mod2x, mod3x, xmod1x, xmod2x, xmod3x, type = 'text', omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year'))

# Saving results

write.csv(stargazer(mod1x, mod2x, mod3x, xmod1x, xmod2x, xmod3x, omit.stat = c('f', 'ser'), omit = c('NAICS', 'FIPS', 'Year')), paste0(direc, 'results/jobs_report.txt'), row.names = FALSE)

# Placebo testing for log(Jobs) :: Randomized Jobs Order

set.seed(42069)

did.coefs <- c()
treat.coefs <- c()
post.coefs <- c()
con.coefs <- c()
ar2 <- c()

for (i in 1:100) {
  
  dfx$Placebo <- sample(dfx$Jobs)
  pmod <- lm(log(Placebo + 1) ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfx)
  did.coefs <- c(did.coefs, pmod$coefficients[length(pmod$coefficients)])
  treat.coefs <- c(treat.coefs, pmod$coefficients[2])
  post.coefs <- c(post.coefs, pmod$coefficients[3])
  con.coefs <- c(con.coefs, pmod$coefficients[1])
  ar2 <- c(ar2, summary(pmod)$adj.r.squared)
  
}

coefs <- c(mean(did.coefs), mean(treat.coefs), mean(post.coefs), mean(con.coefs), mean(ar2))
serrs <- c(sd(did.coefs), sd(treat.coefs), sd(post.coefs), sd(con.coefs), sd(ar2))

placebo.tests <- as.data.frame(cbind(coefs, serrs))

# Placebo testing for Change (in Jobs) :: Randomized Change Order

did.coefsx <- c()
treat.coefsx <- c()
post.coefsx <- c()
con.coefsx <- c()
ar2x <- c()

for (i in 1:100) {
  
  dfxx$Placebo <- sample(dfxx$Change)
  pmod <- lm(Placebo ~ Treated*Post + factor(Year) + factor(FIPS) + factor(NAICS), data = dfxx)
  did.coefsx <- c(did.coefsx, pmod$coefficients[length(pmod$coefficients)])
  treat.coefsx <- c(treat.coefsx, pmod$coefficients[2])
  post.coefsx <- c(post.coefsx, pmod$coefficients[3])
  con.coefsx <- c(con.coefsx, pmod$coefficients[1])
  ar2x <- c(ar2x, summary(pmod)$adj.r.squared)
  
}

coefsx <- c(mean(did.coefsx), mean(treat.coefsx), mean(post.coefsx), mean(con.coefsx), mean(ar2x))
serrsx <- c(sd(did.coefsx), sd(treat.coefsx), sd(post.coefsx), sd(con.coefsx), sd(ar2x))

placebo.testsx <- as.data.frame(cbind(coefsx, serrsx))

# Viewing results from placebo tests

results <- cbind(placebo.tests, placebo.testsx)
results

# Making a nice figure for the paper

nd <- counties(state = 'ND')

plotdf <- dfx %>% filter(Year == 11)
basisdf <- dfx %>% filter(Year == 10)

plotdf <- plotdf %>% group_by(FIPS) %>% summarise(Jobs = sum(Jobs))
basisdf <- basisdf %>% group_by(FIPS) %>% summarise(Jobs = sum(Jobs))

plotdf$Change <- 100 * (plotdf$Jobs - basisdf$Jobs) / basisdf$Jobs

fips <- c()

for (i in 1:nrow(plotdf)) {
  
  f <- as.character(plotdf$FIPS[i])
  
  if (nchar(f) == 3) {f <- paste0('38', f)}
  if (nchar(f) == 2) {f <- paste0('380', f)}
  if (nchar(f) == 1) {f <- paste0('3800', f)}
  
  fips <- c(fips, f)
  
}

plotdf$FIPS <- fips
plotdf[is.infinite(plotdf$Change),]$Change <- 0

plotdf <- st_as_sf(merge(plotdf, nd, by.x = c('FIPS'), by.y = c('GEOID')))

pal <- colorNumeric(palette = 'Reds', domain = plotdf$Change)

jobs_map <- leaflet(plotdf$geometry) %>% setView(lng = -100.33, lat = 47.43, zoom = 6) %>% addTiles() %>% 
  addPolygons(data = plotdf, weight = .2, smoothFactor = 0.5, opacity = 1.0, fillOpacity = .4, color = 'black', fillColor = pal(plotdf$Change))

jobs_map

# Summary statistics of the data

sumdat <- dfx[,4:6]
sumdat$Jobs <- log(sumdat$Jobs + 1)
sumdat$Change <- c(rep(NA, nrow(dfx) - nrow(dfxx)), dfxx$Change)

datasummary_skim(sumdat, fmt = '%.3f')

upper <- quantile(sumdat$Change, c(.995), na.rm = TRUE)
lower <- quantile(sumdat$Change, c(.005), na.rm = TRUE)
sumdat <- sumdat[which(sumdat$Change <= upper & sumdat$Change >= lower),]

ggplot(data = sumdat, aes(x = Change)) +
  geom_histogram(color = 'red4', fill = 'orange', bins = 100) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Histogram of the Change in Jobs at the Industry-County-Year Level') + 
  xlab('Change in Jobs                         ') + 
  ylab('')

