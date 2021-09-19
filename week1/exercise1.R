
# 1.2 ---------------------------------------------------------------------
# Question 1.2: Inspect HMDobj by calling str(HMDobj). Inspect the loaded the data with str(HMDobj$data). 
# Try to overwrite HMDobj$data or anything within it. Hope- fully, you won’t have any luck.

str(HMDobj)

class(HMDobj)

#Check data 
HMDobj$data 

HMDobj$data <- 2 #read only
print(HMDobj)   
methods(print)


# 1.3 ---------------------------------------------------------------------

OEdat <- HMDobj$getOEdata(groups = "Female", agelim = c(100, 110), timelim = c(2000, 2002))
OEdat
str(dat1)
class(dat1)

dat2 <- HMDobj$getSmoothMu(groups = "Female", agelim = c(20,22), timelim = c(1900, 1910))#surf class
str(dat2)
class(dat2)
dat2


# 1.4 ---------------------------------------------------------------------

OEdata <- HMDobj$getOEdata() 
str(OEdata$O)
typeof(OEdata$O)
class(OEdata$O)
is.array(OEdata$O)


calc.OErate <- function(OEdata) {
  # Verify that the input is of class OEdata. Stop if not.
  if(class(OEdata) != "OEdata") stop("not OE class object")
  # Calculate OE-rates
  res <- OEdata$O / OEdata$E
  res[OEdata$E == 0] <- 99
  
  dimnames(res)
  OEdata$agelim
  return(structure(list(groups = OEdata$groups,
                        age = OEdata$agelim,
                        timelim = OEdata$timelim, 
                        mu = res), 
                   class = "surf"))
  # for all groups, ages, and years within input object
  # Output an object of class ’surf’
}
rawOEsurf <- calc.OErate(OEdata)

str(rawOEsurf) # Compare to str(HMDobj$getSmoothMu())
str(rawOEsurf)



#as data frame
surf2plotdat <- function(surf) {
  #surf$mu is array of dim 3
  reshape2::melt(surf$mu, varnames = c("Group", "Age", "Year"), value.name = "Mu")
}
surf2plotdat(surf = rawOEsurf)

library("ggplot2")
library("scales")
theme_set(theme_bw())


plotdat_raw <- surf2plotdat(rawOEsurf)
plotdat_smo <- surf2plotdat(HMDobj$getSmoothMu())
str(plotdat_raw)


plotdat_raw %>% distinct(Year)
subdat_raw <- plotdat_raw[plotdat_raw$Year==2019,]
subdat_smo <- plotdat_smo[plotdat_smo$Year==2019,]

ggplot(mapping = aes(x = Age, y = Mu, color = Group)) +
  geom_point(data = subdat_raw, alpha = 0.6, size = .8) + 
  geom_line(data = subdat_smo) +
  scale_x_continuous(breaks = seq(from = 0, to = 110, by = 10)) +
  scale_y_continuous(trans = log10_trans(), labels = percent_format(accuracy = 0.01), n.breaks = 10) +
  labs(title    = "Historic mortality rates for Denmark in 2019", 
       subtitle = "Smoothed rates superimposed as full lines",
       y        = "Mortality")



# 1.7 ---------------------------------------------------------------------
library(tidyverse)
library(scales)
aidx   <- seq(0, 100, by = 10)
subdat <- subset(plotdat_smo, Age %in% aidx & Year %in% 1950:2020)
ggplot(subdat, aes(x = Year, y = Mu, color = factor(Age))) +
  geom_line(alpha = .8, size = .8) +
  scale_y_continuous(trans = log10_trans(), labels = percent_format(accuracy = 0.01)) +
  labs(title    = "Historic mortality rates for Denmark", 
       subtitle = "Select ages",
       y        = "Mortality") + 
  facet_grid(~Group)



a <- exp(-11)
b <- 0.1

mu <- function(t){
  a*exp(b*t)
}
plot(mu, xlim = range(sample))
rug(sample)

#cum haz
I <- function(t){
  a/b*(exp(b*t)-1)
}

#inv cum haz
Iinv <- function(y){
  log(b/a*y + 1) / b
}
plot(Iinv, xlim = c(0,1))

#sim data
library(survival)
survival::
  x <- rexp(n = 1e5)

sample <- Iinv(y = x)
plot(density(sample));rug(sample(sample, 100))

offset <- 0
dens <- function(x) a*exp(b*(x+offset))*exp(-a/b*(exp(b*(x+offset))-1))
plot(dens, col = 2, add = T, xlim=range(sample))

# 2.3 ---------------------------------------------------------------------
#modal
abline(v = log(b/a)/b, col = 3)
integrate(f = dens, lower = 0, upper = 1000)


# 2.4 ---------------------------------------------------------------------


