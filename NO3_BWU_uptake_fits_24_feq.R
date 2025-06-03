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
  "BWU" = "#a65d17"
)

BWU_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_output_v2.rds")
names(BWU_no3_datq)
str(BWU_no3_datq)
BWU_no3_datq$site <- "BWU"

############
dates<- unique(BWU_no3_datq$date)
BWU_no3_datq$Uadd_int1 <- ifelse(is.na(BWU_no3_datq$Uadd_int) | is.nan(BWU_no3_datq$Uadd_int), 0.000001, BWU_no3_datq$Uadd_int)


### 1
BWU_no3_datq_230718 <- BWU_no3_datq%>%
  filter(date==as.Date("2023-07-18"))

BWU_no3_datq_230718<-BWU_no3_datq_230718[c(-23,-24,-25),]

# model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWU_no3_datq_230718, fct = MM.2())
# summary(model.drm1)
# 
# mm2 <- data.frame(TMR_NO3 = seq(0, max(BWU_no3_datq_230718$TMR_NO3), length.out = 100))
# mm2$Uadd <- predict(model.drm1, newdata = mm2)
# mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_1<- ggplot(BWU_no3_datq_230718, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  #geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  facet_grid(.~date)


### 2 
BWU_no3_datq_230810 <- BWU_no3_datq%>%
  filter(date==as.Date("2023-08-10"))

BWU_no3_datq_230810<-BWU_no3_datq_230810[c(-9,-10,-11,-12,-13,-14),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWU_no3_datq_230810, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWU_no3_datq_230810$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_2<- ggplot(BWU_no3_datq_230810, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWU_no3_datq_230810$TMR_NO3 * 1000), 
           y = min(BWU_no3_datq_230810$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



BWL_nh3_datq_230925a<-BWU_no3_datq_230810[c(-1,-2,-3,-4,-5,-6),]
Uadd_plot_sw<- ggplot(BWU_no3_datq_230810, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWU_no3_datq_230810$sw))
sw_sd <- sd(na.omit(BWU_no3_datq_230810$sw))
v <- mean(na.omit(BWU_no3_datq_230810$Vf_add_int))
v_sd <- sd(na.omit(BWU_no3_datq_230810$Vf_add_int))


### 5
BWU_no3_datq_230925 <- BWU_no3_datq%>%
  filter(date==as.Date("2023-09-25"))

BWU_no3_datq_230925<-BWU_no3_datq_230925[c(-10,-11,-12,-13,-14,-15),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWU_no3_datq_230925, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWU_no3_datq_230925$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_3<- ggplot(BWU_no3_datq_230925, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWU_no3_datq_230925$TMR_NO3 * 1000), 
           y = min(BWU_no3_datq_230925$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



BWU_no3_datq_230925a<-BWU_no3_datq_230925[c(-1,-2),]
Uadd_plot_sw<- ggplot(BWU_no3_datq_230925, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWU_no3_datq_230925a$sw))
sw_sd <- sd(na.omit(BWU_no3_datq_230925a$sw))
v <- mean(na.omit(BWU_no3_datq_230925a$Vf_add_int))
v_sd <- sd(na.omit(BWU_no3_datq_230925a$Vf_add_int))





library(ggpubr)

GBL_nh4_grid <- ggarrange(Uadd_plot_1,
                          Uadd_plot_2,
                          Uadd_plot_3,
                      ncol = 3, nrow = 1,
                      common.legend = TRUE, 
                      legend = "bottom")


ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/BWU_no3_mmfits_grid_25.png", plot = GBL_nh4_grid, width = 8.5, height = 2.5, units = "in")



