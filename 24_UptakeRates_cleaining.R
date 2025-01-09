
library(ggplot2)
library(tidyverse)
library(dplyr)

## think of filtering out the negatice sample slo//pes 

GBL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output_v2.rds")
names(GBL_nh3_datq)
GBL_nh3_datq$up_method <- "NH4"
str(GBL_nh3_datq)

GBL_no3_plot <- ggplot(GBL_nh3_datq, aes(x = TMR_NH3)) +
  geom_point(aes(y = Uadd), col = "red") +
  geom_line(aes(y = Uadd), col = "red") +
  theme_bw() + facet_wrap(~ date) 
GBL_no3_plot

GBL_no3_plot <- ggplot(GBL_nh3_datq%>%filter(Uadd>0), aes(x = nh3_N_add_int)) +
  geom_point(aes(y = Uadd), col = "red") +
  geom_line(aes(y = Uadd), col = "red") +
  theme_bw() + facet_wrap(~ date) 
GBL_no3_plot

GBL_no3_plot <- ggplot(GBL_nh3_datq%>%filter(sw>0), aes(x = TMR_NH3)) +
  geom_point(aes(y = sw), col = "red") +
  geom_line(aes(y = sw), col = "red") +
  theme_bw() + facet_wrap(~ date) 
GBL_no3_plot


GBL_nh3_day <- GBL_nh3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NH3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))

GBU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_output_v2.rds")
names(GBU_nh3_datq)
str(GBU_nh3_datq)
GBU_nh3_datq$up_method <- "NH4"

GBU_nh3_day <- GBU_nh3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NH3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))

BWL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output_v2.rds")
names(BWL_nh3_datq)
str(BWL_nh3_datq)
BWL_nh3_datq$up_method <- "NH4"


BWL_nh3_day <- BWL_nh3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NH3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))


BWU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output_v2.rds")
names(BWU_nh3_datq)
str(BWU_nh3_datq)
BWU_nh3_datq$up_method <- "NH4"

BWU_nh3_day <- BWU_nh3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NH3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))



## example analysis GBL
covariat_dat_GBL <- covariat_datq%>%
  filter(Site=="GBL")%>%
  left_join(GBL_nh3_day, by=c("date"))

GBU_nh3_day1 <- GBU_nh3_day %>%
  filter(!date==as.Date("2022-06-23"))

GBU_nh3_day1 <- GBU_nh3_day1 %>%
  filter(!date==as.Date("2022-07-22"))

covariat_dat_GBU <- covariat_datq%>%
  filter(Site=="GBU")%>%
  left_join(GBU_nh3_day1, by=c("date"))

BWL_nh3_day1 <- BWL_nh3_day %>%
  filter(!date==as.Date("2022-05-26"))

covariat_dat_BWL <- covariat_datq%>%
  filter(Site=="BWL")%>%
  left_join(BWL_nh3_day1, by=c("date"))

BWU_nh3_day1 <- BWU_nh3_day %>%
  filter(!date==as.Date("2023-07-18"))

BWU_nh3_day1 <- BWU_nh3_day1 %>%
  filter(!date==as.Date("2023-08-10"))

covariat_dat_BWU <- covariat_datq%>%
  filter(Site=="BWU")%>%
  left_join(BWU_nh3_day1, by=c("date"))

#####
####


GBL_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_output_v2.rds")
names(GBL_NO3_datq)
GBL_NO3_datq$up_method <- "NO3"
str(GBL_NO3_datq)

GBL_NO3_day <- GBL_NO3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NO3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))

GBU_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_output_v2.rds")
names(GBU_NO3_datq)
str(GBU_NO3_datq)
GBU_NO3_datq$up_method <- "NO3"

GBU_NO3_day <- GBU_NO3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NO3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))

BWL_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_output_v2.rds")
names(BWL_NO3_datq)
str(BWL_NO3_datq)
BWL_NO3_datq$up_method <- "NO3"

BWL_NO3_day <- BWL_NO3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NO3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))


BWU_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_output_v2.rds")
names(BWU_NO3_datq)
str(BWU_NO3_datq)
BWU_NO3_datq$up_method <- "NO3"

BWU_NO3_day <- BWU_NO3_datq%>%
  filter(sw>0 & Uadd>0)%>%
  group_by(date, up_method, ln_injectate_ratio) %>%
  summarise(
    sw_m= mean(sw, na.rm=T),
    sw_sd= sd(sw, na.rm=T),
    Uadd_m=mean(Uadd, na.rm=T),
    Uadd_sd=sd(Uadd, na.rm=T),
    Uadd_int_m=mean(Uadd_int, na.rm=T),
    Uadd_int_sd=sd(Uadd_int, na.rm=T),
    Vf_add_int_m=mean(Vf_add_int, na.rm=T),
    Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
    TMR_Cl = max(TMR_Cl, na.rm = T),
    TMR_N = max(TMR_NO3, na.rm = T)) %>%
  mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
         percent_recovery = c(
           (((sample_recovery*-1)- (ln_injectate_ratio*-1))
            /(ln_injectate_ratio*-1))*100))




## example analysis GBL
GBL_NO3_day1 <- GBL_NO3_day %>%
  filter(!date==as.Date("2021-06-23"))

GBL_NO3_day1 <- GBL_NO3_day1 %>%
  filter(!date==as.Date("2023-06-01"))

covariat_dat_GBL_n3 <- covariat_datq%>%
  filter(Site=="GBL")%>%
  left_join(GBL_NO3_day1, by=c("date"))

GBU_NO3_day1 <- GBU_NO3_day %>%
  filter(!date==as.Date("2021-06-23"))

GBU_NO3_day1 <- GBU_NO3_day1 %>%
  filter(!date==as.Date("2022-04-07"))

GBU_NO3_day1 <- GBU_NO3_day1 %>%
  filter(!date==as.Date("2022-10-03"))

covariat_dat_GBU_n3 <- covariat_datq%>%
  filter(Site=="GBU")%>%
  left_join(GBU_NO3_day1, by=c("date"))


BWL_NO3_day1 <- BWL_NO3_day %>%
  filter(!date==as.Date("2022-08-24"))

covariat_dat_BWL_n3 <- covariat_datq%>%
  filter(Site=="BWL")%>%
  left_join(BWL_NO3_day1, by=c("date"))

BWU_NO3_day1 <- BWU_NO3_day %>%
  filter(!date==as.Date("2023-07-18"))

covariat_dat_BWU_n3 <- covariat_datq%>%
  filter(Site=="BWU")%>%
  left_join(BWU_NO3_day, by=c("date"))
