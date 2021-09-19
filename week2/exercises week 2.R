#week 2
library(tidyverse)

# exercise 1.1 --------------------------------------------------------------

HMDobj <- HMDdatClass$new(Otxtfile = "Data/USA_Deaths_1x1.txt",
                Etxtfile = "Data/USA_Exposures_1x1.txt")
print(HMDobj)
str(HMDobj)

#get data
OEdata <- HMDobj$getOEdata(timelim = c(1950, 2019))
str(OEdata)

#estimate lee carter model on data
mort <- mortClass$new(OEdata = OEdata)
mort$par

#tip: plot
mort$plotdatEstimates(type = "alpha") %>% class
mort$plotdatEstimates(type = "alpha") %>% filter(Alpha > -40) %>% ggplot(aes(x = Age, y = Alpha, col = Group)) + geom_point()
mort$plotdatEstimates(type = "beta") %>% filter(Beta < 3) %>% ggplot(aes(x = Age, y = Beta, col = Group)) + geom_point()
mort$plotdatEstimates(type = "kappa") %>% filter() %>% ggplot(aes(x = Year, y = Kappa, col = Group)) + geom_point()


# 1.2 ---------------------------------------------------------------------
str(OEdata)
vector(mode = "list", length = 2L)
#get 
surf <- HMDobj$getSmoothMu()
str(surf$mu)

surf$mu["Male",,]
surf$mu["Female",,]

a_x_male <- apply(log(surf$mu["Male",,]), 1, mean)
a_x_female <- rowMeans(log(surf$mu["Female",,]))

# SVD
matrix <- matrix(c(1,2,3,
                   4,3,2,
                   1,1,1,
                   4,5,6), nrow = 4, byrow = T)

svd_obj <- svd(matrix)
svd_obj$u ; diag(svd_obj$d) ; t(svd_obj$v)
svd_obj$u %*% diag(svd_obj$d) %*% t(svd_obj$v)
#
dim(surf$mu["Male",,])
length(a_x_male)
svd_male   <- svd(log(surf$mu["Male",,]) - a_x_male)
svd_female <- svd(log(surf$mu["Female",,]) - a_x_male)
str(svd_male)

### check with mortClass
#alpha
plot(mort$par$Female$alpha, a_x_female);abline(0,1)
plot(mort$par$Male$alpha, a_x_male);abline(0,1)
#beta
mort$par$Female$beta
svd_female$u %>% dim


mort$par$Female$kappa

lc.par2mu(alpha = )


# exercise 2 --------------------------------------------------------------



# Load pop, births and deaths
DNK_pop_raw <- read.csv("Data/DNK_Pop_raw.txt")
str(DNK_pop_raw)

DNK_births_raw <- read.csv("Data/DNK_Births_raw.txt")
str(DNK_births_raw)

DNK_deaths_raw <- read.csv("Data/DNK_Deaths_raw.txt")
str(DNK_deaths_raw)


# Define month endpoints as proportions of the year
b.j <- cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31))/365


#MHD data - load 
HMDdatClass$new
HMDobj <- HMDdatClass$new(Otxtfile = 'Data/DNK_Deaths_1x1.txt', 
                          Etxtfile = 'Data/DNK_Exposures_1x1.txt')

(Etrue <- HMDobj$getOEdata(groups = "female", agelim = 65, timelim = 2020)$E)



# Define the quantities we need for the square in question, cf. Figure 1 in Week sheet 2
P1  <- subset(DNK_pop_raw, Year == 2020 & Age == 65 & Sex == "f")$Population
P2  <- subset(DNK_pop_raw, Year == 2021 & Age == 65 & Sex == "f")$Population
B.U <- subset(DNK_births_raw, Year == 2020-65-1 & Month %in% (1:12))$Births
B.L <- subset(DNK_births_raw, Year == 2020-65   & Month %in% (1:12))$Births
D.L <- subset(DNK_deaths_raw, Year == 2020 & Age == 65 & Sex == "f" & Lexis == "TL")$Deaths
D.U <- subset(DNK_deaths_raw, Year == 2020 & Age == 65 & Sex == "f" & Lexis == "TU")$Deaths

# Define monthly birth fractions, cf. formula (10) in Week sheet 2
b.j   <- cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31))/365
p.L.j <- B.L/sum(B.L)
p.U.j <- B.U/sum(B.U)

# Compute birth means and variances for lower and upper triangle, cf. Q.3.1
bvec1   <- sapply(1:12, function(j) (b.j[j+1]+b.j[j])/2)
bvec2   <- sapply(1:12, function(j) (b.j[j+1]^2+b.j[j]^2+b.j[j+1]*b.j[j])/3)
Bmean.L <- sum(p.L.j*bvec1)
Bmean.U <- sum(p.U.j*bvec1)
Bvar.L  <- sum(p.L.j*bvec2) - Bmean.L^2
Bvar.U  <- sum(p.U.j*bvec2) - Bmean.U^2

# Compute lower and upper exposures, cf. formulas in Q.3.6
E.L <- (1-Bmean.L)*P2 + ((1-Bmean.L)/2 - Bvar.L/(2*(1-Bmean.L)))*D.L
E.U <- Bmean.U*P1 - (Bmean.U/2-Bvar.U/(2*Bmean.U))*D.U

cat("  Our exposure estimate :", E.L+E.U,
    "\nHMD's exposure estimate :", Etrue)
