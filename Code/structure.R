
###
###  Set working directory
###

# Your course folder containing the two sub-folders 'Data' and 'Code'.
# setwd("C:/...")

###
###  Load preliminaries
###

source('Code/aux functions.r')

# If these do not work, change to loadPackage() to library() - in this case, you have to install the packages manually before calling library().
loadPackage("R6")
loadPackage("ggplot2")

source('Code/HMDdatClass.r')
source('Code/mortClass.r')
source('Code/fertClass.r')
source('Code/popClass.r')

load('Data/DKfertdat.Rdata')
load('Data/DKpopdat.Rdata')

###
###  OE-data load
### 

HMDobj <- HMDdatClass$new(Otxtfile = 'Data/DNK_Deaths_1x1.txt', Etxtfile = 'Data/DNK_Exposures_1x1.txt')
OEdata <- HMDobj$getOEdata(agelim = c(0,100), timelim = c(1950, 2020))

###
###  Mortality
### 

mort <- mortClass$new(OEdata)
# test <- mort$getMu(); str(test)

###
###  Fertility
### 

fert <- fertClass$new(DKfertdat)
fert$estimate(jumpoffValue = list(SDAgeMother = 4.9)
            , jumpoffSlope = list(MeanAgeMother = 2/25)
            , horizonValue = list(MeanAgeMother = 35, BruttoFer = 2.1)
            , horizon      = list(MeanAgeMother = 2060, BruttoFer = 2060))
# test <- fert$getFer(); str(test)

###
###  Population
### 

mortdat <- mort$getMu(timelim = c(1973, 2250))
fertdat <- fert$getFer(timelim = c(1973, 2250))

pop <- popClass$new("Data/DNK_Pop_1x1.txt")
pop$setMortAndFert(mortdat, fertdat)

# test <- pop$getPop(timelim = c(1973, 2250)); str(test)
