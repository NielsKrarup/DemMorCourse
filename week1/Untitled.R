#week 1 misc. r-code

#load 
library(tidyverse)


#Fertility data 
DKfertdat
class(DKfertdat)
str(DKfertdat)

#DKPop
DKpopdat
str(DKpopdat)

DKpopdat %>% distinct(Year)
DKpopdat %>% filter(Year == max(Year) & Age == min(Age))

DKpopdat %>% distinct(Group)
DKpopdat %>% filter(Group == "Total")
DKpopdat %>% group_by(Year, Group) %>% filter(N != 0) %>% filter(Age == max(Age))

DKpop

which.min(c(1:10,-2))
sample <- sample(1:6, 4, T)
sample
which.min(which(sample == 5))


paste(data.frame(x = "A", y = 9))
paste(df)
