## Load packages
se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

date <-"20240108_DOC_MA" 

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","xts","dygraphs",
         "nrlmetab","cowplot"), require, character.only=T) #"LakeMetabolizer"

dat <- read.csv("/Users/kellyloria/Documents/UNR/shimadzu_out/raw/2024_01_08_in.csv", sep=",")

str(dat)
names(dat)

datDOC <- dat %>% 
  group_by(Site) %>%
  #filter(flag=="N")%>%
  dplyr::summarise(
    DOC_mgL_mean = mean(Values),
    DOC_mgL_sd = sd(Values),
    n = n(),
    DOC_mgL_se = DOC_mgL_sd / sqrt(n)) %>%
  dplyr::select(Site, DOC_mgL_mean, DOC_mgL_sd, DOC_mgL_se) %>%
  mutate(DOC_flag = ifelse(DOC_mgL_se > 0.25 * DOC_mgL_mean, ">25% variance", "good"))

    
dat <- read.csv("/Users/kellyloria/Downloads/TN_0801.csv", sep=",")

str(dat)
names(dat)

datTN <- dat %>% 
  group_by(Site) %>%
  #filter(flag=="N")%>%
  dplyr::summarise(
    TN_ugL_mean = mean(Values, na.rm=T), # normally ugL
    TN_ugL_sd = sd(Values,  na.rm=T), # normally ugL
    n = n(),
    TN_ugL_se = TN_ugL_sd / sqrt(n))%>%
  dplyr::select(Site,TN_ugL_mean, TN_ugL_sd, TN_ugL_se)%>%
  mutate(TN_flag = ifelse(TN_ugL_se > 0.25 * TN_ugL_mean, ">25% variance", "good"))

# 
# dat <- read.csv("/Users/kellyloria/Downloads/shimadzuDat_TC.csv", sep=",")
# str(dat)
# names(dat)
# datTC <- dat %>%
#   group_by(site) %>%
#   summarise(
#     TC_mgL_mean = mean(TC_mgL, na.rm=T),
#     TC_mgL_sd = sd(TC_mgL,  na.rm=T),
#     n = n(),
#     TC_mgL_se = TC_mgL_sd / sqrt(n)) %>%
#   select(site,TC_mgL_mean, TC_mgL_sd, TC_mgL_se)

shimDat <- datDOC %>% 
  #full_join(datTC) %>%
  full_join(datTN)

save.dat.fxn <- function(dat){
  write.csv(dat, paste("/Users/kellyloria/Documents/UNR/shimadzu_out/Shimadzu_",date,"_.csv",sep =""))
}

save.dat.fxn(shimDat)




TNcheck <- ggplot(shimDat, aes(x = Site, y = TN_ugL_mean, color=TN_flag)) +
  geom_point() +
  geom_errorbar(aes(ymin =  TN_ugL_mean- TN_ugL_se, ymax = TN_ugL_mean + TN_ugL_se), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab(expression(paste("TN (", mu, "L)")))

#ggsave("ShimTNQAC_v1.png", plot = TNcheck, width = 9, height = 5, units = "in")


DOCcheck <- ggplot(shimDat, aes(x = Site, y = DOC_mgL_mean, color=DOC_flag)) +
  geom_point() +
  geom_errorbar(aes(ymin =  DOC_mgL_mean- DOC_mgL_se, ymax = DOC_mgL_mean + DOC_mgL_se), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("DOC mgL")

# ggsave("ShimDOCQAC_v1.png", plot = DOCcheck, width = 9, height = 5, units = "in")




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
