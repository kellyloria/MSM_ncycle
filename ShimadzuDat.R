## Load packages
se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

date <-"20220602" 

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","xts","dygraphs",
         "nrlmetab","cowplot"), require, character.only=T) #"LakeMetabolizer"

dat <- read.csv("/Users/kellyloria/Downloads/shimadzuDat_DOC2.csv", sep=",")
datDOC <- dat %>% 
  group_by(site) %>%
  summarise(
    DOC_mgL_mean = mean(DOC_mgL),
    DOC_mgL_sd = sd(DOC_mgL),
    n = n(),
    DOC_mgL_se = DOC_mgL_sd / sqrt(n)) %>%
  select(site,DOC_mgL_mean, DOC_mgL_sd, DOC_mgL_se)

    
dat <- read.csv("/Users/kellyloria/Downloads/shimadzuDat_TN2.csv", sep=",")
str(dat)
names(dat)
datTN <- dat %>% 
  group_by(site) %>%
  summarise(
    TN_ugL_mean = mean(TN_ugL, na.rm=T),
    TN_ugL_sd = sd(TN_ugL,  na.rm=T),
    n = n(),
    TN_ugL_se = TN_ugL_sd / sqrt(n))%>%
  select(site,TN_ugL_mean, TN_ugL_sd, TN_ugL_se)

dat <- read.csv("/Users/kellyloria/Downloads/shimadzuDat_TC.csv", sep=",")
str(dat)
names(dat)
datTC <- dat %>%
  group_by(site) %>%
  summarise(
    TC_mgL_mean = mean(TC_mgL, na.rm=T),
    TC_mgL_sd = sd(TC_mgL,  na.rm=T),
    n = n(),
    TC_mgL_se = TC_mgL_sd / sqrt(n)) %>%
  select(site,TC_mgL_mean, TC_mgL_sd, TC_mgL_se)

shimDat<- datDOC %>% 
  #full_join(datTC) %>%
  full_join(datTN)

save.dat.fxn <- function(dat){
  write.csv(dat, paste("/Users/kellyloria/Documents/UNR/shimadzu_out/Shimadzu_",date,"_.csv",sep =""))
}

save.dat.fxn(shimDat)





str(dat)
names(dat)

##
dat_q <- dat %>% 
  #select(Type, Anal., Sample, Name,Time,Conc.,Result) %>% 
  subset(Type=="Unknown")
  #rename(datetime = X.1, do.obs = Dissolved.Oxygen) %>% 
  #mutate(datetime = mdy_hm(datetime))  
str(dat_q)
  
  #mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S")) 
