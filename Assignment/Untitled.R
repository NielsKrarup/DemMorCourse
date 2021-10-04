tmp <- 
  surf_obj$mu[,,dimnames(surf_obj$mu)$Year %in% as.character(1836:2020)] -
  surf_obj$mu[,,dimnames(surf_obj$mu)$Year %in% as.character(1835:2019)]

tmp[,,5]

#sum array 

apply(X = OEdata$O, MARGIN = c(1,3), FUN = sum)

arr_long <- reshape2::melt(surf_obj$mu) %>% arrange(Year, Group, Age)
arr_long %>% mutate(age_mod = Age %% 10)

arr_long %>% filter(Year == 1839)

0:100 %% 10

ggplot(arr_long %>% mutate(age_mod = Age %% 10) %>% filter(Year < 1900),
       aes(x = Year, y = value, col = Age)) + 
  geom_line(aes(group = Age, col = factor(age_mod)), alpha = 0.4) +
  facet_wrap(vars(Group), scales = "free") + 
  #scale_y_continuous(trans = "log10") +
  guides(col="none")

ggplot(arr_long %>% filter(Age ==0)%>% mutate(Age_mod = Age %% 10), aes(x = Year, y = value)) + 
  geom_line(aes(group = factor(Group):factor(Age, col = factor(Age_mod)), alpha = 0.4) +
  facet_wrap(vars(Age_grp), scales = "free") + 
  #facet_wrap(vars(Group), scales = "free") + 
  
  #scale_y_continuous(trans = "log10") +
  guides(col="none")

  
  surf_obj$mu[,1,1:2]
  