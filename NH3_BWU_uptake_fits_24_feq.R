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

BWU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output_v2.rds")
names(BWU_nh3_datq)
str(BWU_nh3_datq)
BWU_nh3_datq$site <- "BWU"

############
dates<- unique(BWU_nh3_datq$date)
BWU_nh3_datq$Uadd_int1 <- ifelse(is.na(BWU_nh3_datq$Uadd_int) | is.nan(BWU_nh3_datq$Uadd_int), 0.000001, BWU_nh3_datq$Uadd_int)


### 1
BWU_nh3_datq_220824 <- BWU_nh3_datq%>%
  filter(date==as.Date("2022-08-24"))

BWU_nh3_datq_220824<-BWU_nh3_datq_220824[c(-22,-23,-24,-25),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWU_nh3_datq_220824, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWU_nh3_datq_220824$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_1<- ggplot(BWU_nh3_datq_220824, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWU_nh3_datq_220824$TMR_NH3 * 1000), 
           y = min(BWU_nh3_datq_220824$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)




BWU_nh3_datq_220824a<-BWU_nh3_datq_220824[c(-1,-2,-3,-4,-5,-6, -7,-8),]
Uadd_plot_sw<- ggplot(BWU_nh3_datq_220824a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWU_nh3_datq_220824a$sw))
sw_sd <- sd(na.omit(BWU_nh3_datq_220824a$sw))
v <- mean(na.omit(BWU_nh3_datq_220824a$Vf_add_int))
v_sd <- sd(na.omit(BWU_nh3_datq_220824a$Vf_add_int))


###3 
BWU_nh3_datq_221012 <- BWU_nh3_datq%>%
  filter(date==as.Date("2022-10-12"))

BWU_nh3_datq_221012<-BWU_nh3_datq_221012[c(-1,-2,-3,-14,-15),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWU_nh3_datq_221012, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWU_nh3_datq_221012$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_2<- ggplot(BWU_nh3_datq_221012, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWU_nh3_datq_221012$TMR_NH3 * 1000), 
           y = min(BWU_nh3_datq_221012$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



BWU_nh3_datq_221012a<-BWU_nh3_datq_221012[c(-1,-2,-3,-4),]
Uadd_plot_sw<- ggplot(BWU_nh3_datq_221012a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWU_nh3_datq_221012a$sw))
sw_sd <- sd(na.omit(BWU_nh3_datq_221012a$sw))
v <- mean(na.omit(BWU_nh3_datq_221012a$Vf_add_int))
v_sd <- sd(na.omit(BWU_nh3_datq_221012a$Vf_add_int))

### 3
BWU_nh3_datq_23078 <- BWU_nh3_datq%>%
  filter(date==as.Date("2023-07-18"))

BWU_nh3_datq_23078<-BWU_nh3_datq_23078[c(-10,-11,-12,-13,-14, -15),]
# 
# model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWU_nh3_datq_23078, fct = MM.2())
# summary(model.drm1)
# 
# mm2 <- data.frame(TMR_NH3 = seq(0, max(BWU_nh3_datq_23078$TMR_NH3), length.out = 100))
# mm2$Uadd <- predict(model.drm1, newdata = mm2)
# mm2$site <- "BWU"
# 
# ## summaries to add to plot:
# summary_drm <- summary(model.drm1)
# est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
# std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
# est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
# std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])
# 
# label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_3<- ggplot(BWU_nh3_datq_23078, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
 # geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  facet_grid(.~date)



### 4
BWU_nh3_datq_230810 <- BWU_nh3_datq%>%
  filter(date==as.Date("2023-08-10"))

BWU_nh3_datq_230810<-BWU_nh3_datq_230810[c(-7,-8,-9,-10,-11,-12,-13,-14,-15,-16, -17,-18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWU_nh3_datq_230810, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWU_nh3_datq_230810$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_4<- ggplot(BWU_nh3_datq_230810, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWU_nh3_datq_230810$TMR_NH3 * 1000), 
           y = min(BWU_nh3_datq_230810$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)




BWU_nh3_datq_230810a<-BWU_nh3_datq_230810[c(-1),]
Uadd_plot_sw<- ggplot(BWU_nh3_datq_230810a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWU_nh3_datq_230810a$sw))
sw_sd <- sd(na.omit(BWU_nh3_datq_230810a$sw))
v <- mean(na.omit(BWU_nh3_datq_230810a$Vf_add_int))
v_sd <- sd(na.omit(BWU_nh3_datq_230810a$Vf_add_int))

### 5
BWU_nh3_datq_230925 <- BWU_nh3_datq%>%
  filter(date==as.Date("2023-09-25"))

BWU_nh3_datq_230925<-BWU_nh3_datq_230925[c(-1,-14,-15,-16,-17,-18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = BWU_nh3_datq_230925, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(BWU_nh3_datq_230925$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_5<- ggplot(BWU_nh3_datq_230925, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWU_nh3_datq_230925$TMR_NH3 * 1000), 
           y = min(BWU_nh3_datq_230925$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)




library(ggpubr)

GBL_nh4_grid <- ggarrange(Uadd_plot_1,
                          Uadd_plot_2,
                          Uadd_plot_3,
                          Uadd_plot_4,
                          Uadd_plot_5,
                      ncol = 3, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")


ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/BWU_nh4_mmfits_grid.png", plot = GBL_nh4_grid, width = 8, height = 4.2, units = "in")



