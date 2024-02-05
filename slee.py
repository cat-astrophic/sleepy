# This script performs data prep for the sleepy time project and also creates a figure

# Importing required modules

import pandas as pd
from geopy.geocoders import Nominatim
import geopy.distance
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

# Project directory

direc = 'D:/sleepy/'

# Reading in the marathon data for 2010-2014

mar = pd.DataFrame()

for i in range(2010,2015):
    
    tmp = pd.read_csv(direc + 'data/full/' + str(i) + '.txt', sep = '\t')
    tmp = pd.concat([tmp, pd.Series([i]*len(tmp), name = 'Year')], axis = 1)
    mar = pd.concat([mar, tmp], axis = 0).reset_index(drop = True)

# Reading in the half marathon data for 2010-2014

half = pd.DataFrame()

for i in range(2010,2015):
    
    tmp = pd.read_csv(direc + 'data/half/' + str(i) + '.txt', sep = '\t')
    tmp = pd.concat([tmp, pd.Series([i]*len(tmp), name = 'Year')], axis = 1)
    half = pd.concat([half, tmp], axis = 0).reset_index(drop = True)

# Subset dataframes for relevant data

mar = mar[['Name', 'Sex', 'Age', 'City', 'State', 'Time', 'Year']]
half = half[['Name', 'Sex', 'Age', 'City', 'State', 'Time', 'Year']]

# Location of changes between 'columns' in the txt file that somebody used space to create rather than tabs...

breaks = [0, 6, 15, 20, 44, 63, 66, 70, 74, 82]

# Reading in the marathon data for 2007-2009

mar2 = pd.DataFrame()

for i in range(2007,2010):
    
    text = []
    name = []
    sex = []
    age = []
    city = []
    state = []
    time = []
    year = []
    
    file = open(direc + 'data/full/' + str(i) + '.txt', 'r')
    
    while True:
        
        mess = file.readline()
        text.append(mess)
        
        if not mess:
            
            break
        
    text = text[1:]
    text = [t for t in text if t != '']
    
    for line in text:
        
        name.append(line[20:44])
        sex.append(line[70:74])
        age.append(line[66:70])
        city.append(line[44:63])
        state.append(line[63:66])
        time.append(line[74:82])
        year.append(i)
    
    name = pd.Series(name, name = 'Name')
    sex = pd.Series(sex, name = 'Sex')
    age = pd.Series(age, name = 'Age')
    city = pd.Series(city, name = 'City')
    state = pd.Series(state, name = 'State')
    time = pd.Series(time, name = 'Time')
    year = pd.Series(year, name = 'Year')
    
    tmp = pd.concat([name, sex, age, city, state, time, year], axis = 1)
    mar2 = pd.concat([mar2, tmp], axis = 0).reset_index(drop = True)

# Reading in the marathon data for 2007-2009

half2 = pd.DataFrame()

for i in range(2007,2010):
    
    text = []
    name = []
    sex = []
    age = []
    city = []
    state = []
    time = []
    year = []
    
    file = open(direc + 'data/half/' + str(i) + '.txt', 'r')
    
    while True:
        
        mess = file.readline()
        text.append(mess)
        
        if not mess:
            
            break
        
    text = text[1:]
    text = [t for t in text if t != '']
    
    for line in text:
        
        name.append(line[20:44])
        sex.append(line[70:74])
        age.append(line[66:70])
        city.append(line[44:63])
        state.append(line[63:66])
        time.append(line[74:82])
        year.append(i)
    
    name = pd.Series(name, name = 'Name')
    sex = pd.Series(sex, name = 'Sex')
    age = pd.Series(age, name = 'Age')
    city = pd.Series(city, name = 'City')
    state = pd.Series(state, name = 'State')
    time = pd.Series(time, name = 'Time')
    year = pd.Series(year, name = 'Year')
    
    tmp = pd.concat([name, sex, age, city, state, time, year], axis = 1)
    half2 = pd.concat([half2, tmp], axis = 0).reset_index(drop = True)

# Merge dataframes

mar = pd.concat([mar, mar2], axis = 0).reset_index(drop = True)
half = pd.concat([half, half2], axis = 0).reset_index(drop = True)

mar = pd.concat([mar, pd.Series(['Marathon']*len(mar), name = 'Event')], axis = 1)
half = pd.concat([half, pd.Series(['Half-Marathon']*len(half), name = 'Event')], axis = 1)

df = pd.concat([mar, half], axis = 0).reset_index(drop = True)

# Convert Time column into seconds

def secs_fx(inp):
    
    inp = inp.replace(' ', '')
    s = int(inp[-2:]) + int(inp[2:4])*60 + int(inp[0])*3600
    
    return s

secs = [secs_fx(t) for t in df.Time]

df = pd.concat([df, pd.Series(secs, name = 'Seconds')], axis = 1)

# Set up geocoding feature

geolocator = Nominatim(user_agent = 'myapplication')

# Get runner hometown coordinates

coords = []

for i in range(len(df)):
    
    print('Collecting coordinates for observation ' + str(i+1) + ' of ' + str(len(df)) + '.......')
    
    try:
        
        coords.append(geolocator.geocode(df.City[i] + ', ' + df.State[i])[1])
        
    except:
        
        coords.append(None)

df = pd.concat([df, pd.Series(coords, name = 'Coordinates')], axis = 1)

# Coordinates for the start of the Bismarck Marathon

bm = (46.82009858228105, -100.78278262793603)

# Compute the distance between runners and the race

dists = []

for c in df.Coordinates:
    
    dists.append(geopy.distance.distance(c, bm).mi)

df = pd.concat([df, pd.Series(dists, name = 'Distance')], axis = 1)

# Parse coordiantes because R

lats = [c[0] if c != None else None for c in df.Coordinates]
lons = [c[1] if c != None else None for c in df.Coordinates]

df = pd.concat([df, pd.Series(lats, name = 'Latitude'), pd.Series(lons, name = 'Longitude')], axis = 1)

# Clean city and state names for R

def city_cleaner(dirty):
    
    if len(dirty) < 1:
        
        pass
        
    else:
        
        while (dirty[-1] == ' ') == True:
            
            dirty = dirty[:-1]
            
            if len(dirty) < 1:
                
                break
            
    return dirty

states = [str(s).replace(' ', '') for s in df.State]
cities = [city_cleaner(str(c)) for c in df.City]
locs = [cities[i] + ', ' + states[i] for i in range(len(states))]

df = df[[x for x in df.columns if x != 'City']]
df = df[[x for x in df.columns if x != 'State']]

cities = pd.Series(cities, name = 'City')
states = pd.Series(states, name = 'State')
locs = pd.Series(locs, name = 'Location')

df = pd.concat([df, states, cities, locs], axis = 1)

# Clean runner names for R

names = [city_cleaner(n) for n in df.Name]
df = df[[x for x in df.columns if x != 'Name']]
df = pd.concat([pd.Series(names, name = 'Name'), df], axis = 1)

# Save prepped data to file

df.to_csv(direc + 'data/data.csv', index = False)

# Making the figure for the paper

# Exponential decay function with time zones

def decay(x, y):
    
    d = np.sqrt(x*x + y*y)
        
    if x < -0.5:
        
        decayed = np.exp(-d)
        
    else:
        
        decayed = np.exp(-.5*d)
        
    return decayed

# Generating the data to plot

x = [i / 1000 for i in range(-1000,1001)]*2001

y = []

for i in range(2001):
    
    for j in range(2001):
        
        y.append(x[i])

z = []

for i in range(len(x)):
    
    z.append(decay(x[i], y[i]))

# Making the map

cmap = cm.rainbow
plt.scatter(x, y, c = cmap(z), alpha = 0.5)
plt.show()

