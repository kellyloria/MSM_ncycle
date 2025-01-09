# Lets start calculating uptake rates
## 24_GBL_NO3_BTC_output.rds

library(drc)
## remember: maximum response d = V max - maximum response
## remember: maximum response e = km - half saturation constant
library(ggplot2)
library(tidyverse)
library(dplyr)

###
site_colors <- c(
  "BWL" = "#3283a8",
  "BWU" = "#3258a8",
  "GBL" = "#a67d17",
  "BWL" = "#a65d17"
)
##

# GBL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output_v2.rds")
# names(GBL_nh3_datq)
# str(GBL_nh3_datq)
# GBL_nh3_datq$site <- "GBL"


# BWL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output_v2.rds")
# names(BWL_nh3_datq)
# str(BWL_nh3_datq)
# BWL_nh3_datq$site <- "BWL"

BWL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output_v2.rds")
names(BWL_nh3_datq)
str(BWL_nh3_datq)
BWL_nh3_datq$site <- "BWL"

BWU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output.rds")
names(BWU_nh3_datq)
str(BWU_nh3_datq)

############
dates<- unique(BWL_nh3_datq$date)
BWL_nh3_datq$Uadd_int1 <- ifelse(is.na(BWL_nh3_datq$Uadd_int) | is.nan(BWL_nh3_datq$Uadd_int), 0.000001, BWL_nh3_datq$Uadd_int)


### 1
BWL_nh3_datq_220526 <- BWL_nh3_datq%>%
  filter(date==as.Date("2022-05-26"))

BWL_nh3_datq_220526<-BWL_nh3_datq_220526[c(-10,-11,-12),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_220526, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_220526$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_1<- ggplot(BWL_nh3_datq_220526, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000),lty=3) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_220526$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_220526$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "grey25") +
  facet_grid(.~date)


sw <- mean(na.omit(BWL_nh3_datq_220526$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_220526$sw))

v <- mean(na.omit(BWL_nh3_datq_220526$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_220526$Vf_add_int))



### 2
BWL_nh3_datq_220824 <- BWL_nh3_datq%>%
  filter(date==as.Date("2022-08-24"))

BWL_nh3_datq_220824<-BWL_nh3_datq_220824[c(-21, -22,-23),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_220824, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_220824$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_2<- ggplot(BWL_nh3_datq_220824, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_220824$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_220824$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "grey30")+
  facet_grid(.~date)

sw <- mean(na.omit(BWL_nh3_datq_220824$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_220824$sw))

v <- mean(na.omit(BWL_nh3_datq_220824$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_220824$Vf_add_int))



###3 
BWL_nh3_datq_221012 <- BWL_nh3_datq%>%
  filter(date==as.Date("2022-10-12"))

#BWL_nh3_datq_221012<-BWL_nh3_datq_221012[c(-16),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_221012, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_221012$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_3<- ggplot(BWL_nh3_datq_221012, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_221012$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_221012$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


Uadd_plot_sw<- ggplot(BWL_nh3_datq_221012, aes(x = TMR_NH3*1000, y = sw*1000, color = site)) +
  geom_point(size = 2, shape = 17)

#
BWL_nh3_datq_221012<-BWL_nh3_datq_221012[c(-1,-2),]


sw <- mean(na.omit(BWL_nh3_datq_221012$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_221012$sw))

v <- mean(na.omit(BWL_nh3_datq_221012$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_221012$Vf_add_int))



### 4
BWL_nh3_datq_221121 <- BWL_nh3_datq%>%
  filter(date==as.Date("2022-11-21"))

BWL_nh3_datq_221121<-BWL_nh3_datq_221121[c(-19,-20, -21),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_221121, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_221121$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_4<- ggplot(BWL_nh3_datq_221121, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  annotate("text", x = max(BWL_nh3_datq_221121$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_221121$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)

BWL_nh3_datq_221121a<-BWL_nh3_datq_221121[c(-1,-2,-3),]


Uadd_plot_sw<- ggplot(BWL_nh3_datq_221121a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()

sw <- mean(na.omit(BWL_nh3_datq_221121a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_221121a$sw))

v <- mean(na.omit(BWL_nh3_datq_221121a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_221121a$Vf_add_int))




### 5
BWL_nh3_datq_221219 <- BWL_nh3_datq%>%
  filter(date==as.Date("2022-12-19"))

BWL_nh3_datq_221219<-BWL_nh3_datq_221219[c(-17,-18, -19),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_221219, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_221219$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_5<- ggplot(BWL_nh3_datq_221219, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_221219$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_221219$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)




BWL_nh3_datq_221219a<-BWL_nh3_datq_221219[c(-1,-2,-3,-4),]
Uadd_plot_sw<- ggplot(BWL_nh3_datq_221219a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()

sw <- mean(na.omit(BWL_nh3_datq_221219a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_221219a$sw))
v <- mean(na.omit(BWL_nh3_datq_221219a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_221219a$Vf_add_int))




### 6
BWL_nh3_datq_230215 <- BWL_nh3_datq%>%
  filter(date==as.Date("2023-02-15"))

BWL_nh3_datq_230215<-BWL_nh3_datq_230215[c(-13,-14, -15, -16, -17,-18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_230215, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_230215$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_6<- ggplot(BWL_nh3_datq_230215, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_230215$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_230215$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_nh3_datq_230215a<-BWL_nh3_datq_230215[c(-3,-4, -5),]
Uadd_plot_sw<- ggplot(BWL_nh3_datq_230215a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()

sw <- mean(na.omit(BWL_nh3_datq_230215a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_230215a$sw))
v <- mean(na.omit(BWL_nh3_datq_230215a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_230215a$Vf_add_int))


### 7
BWL_nh3_datq_230405 <- BWL_nh3_datq%>%
  filter(date==as.Date("2023-04-05"))

BWL_nh3_datq_230405<-BWL_nh3_datq_230405[c(-15,-16,-17,-18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_230405, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_230405$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_7<- ggplot(BWL_nh3_datq_230405, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_230405$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_230405$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_nh3_datq_230405a<-BWL_nh3_datq_230405[c(-1,-2,-3,-4, -5),]
Uadd_plot_sw<- ggplot(BWL_nh3_datq_230405a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()

sw <- mean(na.omit(BWL_nh3_datq_230405a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_230405a$sw))
v <- mean(na.omit(BWL_nh3_datq_230405a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_230405a$Vf_add_int))


### 8
BWL_nh3_datq_23078 <- BWL_nh3_datq%>%
  filter(date==as.Date("2023-07-18"))

BWL_nh3_datq_23078<-BWL_nh3_datq_23078[c(-11,-12,-13,-14),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_23078, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_23078$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_8<- ggplot(BWL_nh3_datq_23078, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_23078$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_23078$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_nh3_datq_23078a<-BWL_nh3_datq_23078[c(-1,-2),]
Uadd_plot_sw<- ggplot(BWL_nh3_datq_23078a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_nh3_datq_23078a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_23078a$sw))
v <- mean(na.omit(BWL_nh3_datq_23078a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_23078a$Vf_add_int))




### 9
BWL_nh3_datq_230810 <- BWL_nh3_datq%>%
  filter(date==as.Date("2023-08-10"))

BWL_nh3_datq_230810<-BWL_nh3_datq_230810[c(-13,-14,-15,-16),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_230810, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_230810$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_9<- ggplot(BWL_nh3_datq_230810, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_230810$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_230810$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_nh3_datq_230810a<-BWL_nh3_datq_230810[c(-1,-2,-3),]
Uadd_plot_sw<- ggplot(BWL_nh3_datq_230810a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_nh3_datq_230810a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_230810a$sw))
v <- mean(na.omit(BWL_nh3_datq_230810a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_230810a$Vf_add_int))




### 10
BWL_nh3_datq_230925 <- BWL_nh3_datq%>%
  filter(date==as.Date("2023-09-25"))

BWL_nh3_datq_230925<-BWL_nh3_datq_230925[c(-17,-18,-19),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWL_nh3_datq_230925, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWL_nh3_datq_230925$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_10<- ggplot(BWL_nh3_datq_230925, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_nh3_datq_230925$TMR_NH3 * 1000), 
           y = min(BWL_nh3_datq_230925$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



BWL_nh3_datq_230925a<-BWL_nh3_datq_230925[c(-1,-2,-3,-4,-5,-6),]
Uadd_plot_sw<- ggplot(BWL_nh3_datq_230925a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_nh3_datq_230925a$sw))
sw_sd <- sd(na.omit(BWL_nh3_datq_230925a$sw))
v <- mean(na.omit(BWL_nh3_datq_230925a$Vf_add_int))
v_sd <- sd(na.omit(BWL_nh3_datq_230925a$Vf_add_int))




library(ggpubr)

GBL_nh4_grid <- ggarrange(Uadd_plot_1,
                          Uadd_plot_2,
                          Uadd_plot_3,
                          Uadd_plot_4,
                          Uadd_plot_5,
                          Uadd_plot_6,
                          Uadd_plot_7,
                          Uadd_plot_8,
                          Uadd_plot_9,
                          Uadd_plot_10,
                      ncol = 3, nrow = 4,
                      common.legend = TRUE, 
                      legend = "bottom")


ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/BWL_nh4_mmfits_grid.png", plot = GBL_nh4_grid, width = 8.5, height = 7, units = "in")



