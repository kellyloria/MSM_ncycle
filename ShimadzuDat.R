## Load packages
se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

date <-"20230501" 

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","xts","dygraphs",
         "nrlmetab","cowplot"), require, character.only=T) #"LakeMetabolizer"

dat <- read.csv("/Users/kellyloria/Downloads/DOC0501.csv", sep=",")

str(dat)
names(dat)

dat<- subset(dat, flag=="A")


datDOC <- dat %>% 
  group_by(site) %>%
  dplyr::summarise(
    DOC_mgL_mean = mean(DOC_mgL),
    DOC_mgL_sd = sd(DOC_mgL),
    n = n(),
    DOC_mgL_se = DOC_mgL_sd / sqrt(n)) %>%
  dplyr::select(site,DOC_mgL_mean, DOC_mgL_sd, DOC_mgL_se)

    
dat <- read.csv("/Users/kellyloria/Downloads/TN0501.csv", sep=",")

str(dat)
names(dat)

dat<- subset(dat, flag=="A")

datTN <- dat %>% 
  group_by(site) %>%
  dplyr::summarise(
    TN_ugL_mean = mean(TN_ugL, na.rm=T), # normally ugL
    TN_ugL_sd = sd(TN_ugL,  na.rm=T), # normally ugL
    n = n(),
    TN_ugL_se = TN_ugL_sd / sqrt(n))%>%
  dplyr::select(site,TN_ugL_mean, TN_ugL_sd, TN_ugL_se)

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

shimDat <- datDOC %>% 
  #full_join(datTC) %>%
  full_join(datTN)

save.dat.fxn <- function(dat){
  write.csv(dat, paste("/Users/kellyloria/Documents/UNR/shimadzu_out/Shimadzu_",date,"_.csv",sep =""))
}

save.dat.fxn(shimDat)




TNcheck <- ggplot(shimDat, aes(x = site, y = TN_ugL_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin =  TN_ugL_mean- TN_ugL_se, ymax = TN_ugL_mean + TN_ugL_se), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab(expression(paste("TN (", mu, "L)")))

ggsave("ShimTNQAC.png", plot = TNcheck, width = 7, height = 5, units = "in")


DOCcheck <- ggplot(shimDat, aes(x = site, y = DOC_mgL_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin =  DOC_mgL_mean- DOC_mgL_se, ymax = DOC_mgL_mean + DOC_mgL_se), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("DOC mgL")

ggsave("ShimDOCQAC.png", plot = DOCcheck, width = 7, height = 5, units = "in")




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
