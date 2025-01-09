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
  "GBU" = "#a65d17"
)
##


GBU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_output_v2.rds")
names(GBU_nh3_datq)
str(GBU_nh3_datq)
GBU_nh3_datq$site <- "GBU"


############
unique(GBU_nh3_datq$date)
GBU_nh3_datq$Uadd_int1 <- ifelse(is.na(GBU_nh3_datq$Uadd_int) | is.nan(GBU_nh3_datq$Uadd_int), 0.000001, GBU_nh3_datq$Uadd_int)


### 1
GBU_nh3_datq_220623 <- GBU_nh3_datq%>%
  filter(date==as.Date("2022-06-23"))

GBU_nh3_datq_220623<-GBU_nh3_datq_220623[c(-11,-12),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBU_nh3_datq_220623, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBU_nh3_datq_220623$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_1<- ggplot(GBU_nh3_datq_220623, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBU_nh3_datq_220623$TMR_NH3 * 1000), 
           y = min(GBU_nh3_datq_220623$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black") +
  facet_grid(.~date)


GBU_nh3_datq_220623a<-GBU_nh3_datq_220623[c(11,12,13,14,15),]
Uadd_plot_sw<- ggplot(GBU_nh3_datq_220623a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBU_nh3_datq_220623a$sw))
sw_sd <- sd(na.omit(GBU_nh3_datq_220623a$sw))
v <- mean(na.omit(GBU_nh3_datq_220623a$Vf_add_int))
v_sd <- sd(na.omit(GBU_nh3_datq_220623a$Vf_add_int))



### 2
GBU_nh3_datq_220722 <- GBU_nh3_datq%>%
  filter(date==as.Date("2022-07-22"))

GBU_nh3_datq_220722<-GBU_nh3_datq_220722[c(-17,-18,-19, -20,-21, -22),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBU_nh3_datq_220722, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBU_nh3_datq_220722$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_2<- ggplot(GBU_nh3_datq_220722, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000), lty=3) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBU_nh3_datq_220722$TMR_NH3 * 1000), 
           y = min(GBU_nh3_datq_220722$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "grey30")+
  facet_grid(.~date)


GBU_nh3_datq_220722a<-GBU_nh3_datq_220722[c(-1,-2,-14,-15,-16),]
Uadd_plot_sw<- ggplot(GBU_nh3_datq_220722a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBU_nh3_datq_220722a$sw))
sw_sd <- sd(na.omit(GBU_nh3_datq_220722a$sw))
v <- mean(na.omit(GBU_nh3_datq_220722a$Vf_add_int))
v_sd <- sd(na.omit(GBU_nh3_datq_220722a$Vf_add_int))





###3 
GBU_nh3_datq_221003 <- GBU_nh3_datq%>%
  filter(date==as.Date("2022-10-03"))

GBU_nh3_datq_221003<-GBU_nh3_datq_221003[c(-15,-16,-17,-18,-19, -20),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBU_nh3_datq_221003, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBU_nh3_datq_221003$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_3<- ggplot(GBU_nh3_datq_221003, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBU_nh3_datq_221003$TMR_NH3 * 1000), 
           y = min(GBU_nh3_datq_221003$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)

GBU_nh3_datq_221003a<-GBU_nh3_datq_221003[c(-1,-2),]
Uadd_plot_sw<- ggplot(GBU_nh3_datq_221003a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBU_nh3_datq_221003a$sw))
sw_sd <- sd(na.omit(GBU_nh3_datq_221003a$sw))
v <- mean(na.omit(GBU_nh3_datq_221003a$Vf_add_int))
v_sd <- sd(na.omit(GBU_nh3_datq_221003a$Vf_add_int))



### 4
GBU_nh3_datq_230615 <- GBU_nh3_datq%>%
  filter(date==as.Date("2023-06-15"))

GBU_nh3_datq_230615<-GBU_nh3_datq_230615[c(-12,-13, -14),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBU_nh3_datq_230615, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBU_nh3_datq_230615$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_4<- ggplot(GBU_nh3_datq_230615, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  annotate("text", x = max(GBU_nh3_datq_230615$TMR_NH3 * 1000), 
           y = min(GBU_nh3_datq_230615$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



GBU_nh3_datq_230615a<-GBU_nh3_datq_230615[c(-1,-2,-3),]
Uadd_plot_sw<- ggplot(GBU_nh3_datq_230615a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBU_nh3_datq_230615a$sw))
sw_sd <- sd(na.omit(GBU_nh3_datq_230615a$sw))
v <- mean(na.omit(GBU_nh3_datq_230615a$Vf_add_int))
v_sd <- sd(na.omit(GBU_nh3_datq_230615a$Vf_add_int))



### 5
GBU_nh3_datq_230710 <- GBU_nh3_datq%>%
  filter(date==as.Date("2023-07-10"))

GBU_nh3_datq_230710<-GBU_nh3_datq_230710[c(-11,-12,-13, -14, -15, -16),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBU_nh3_datq_230710, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBU_nh3_datq_230710$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_5<- ggplot(GBU_nh3_datq_230710, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBU_nh3_datq_230710$TMR_NH3 * 1000), 
           y = min(GBU_nh3_datq_230710$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


GBU_nh3_datq_230615a<-GBU_nh3_datq_230710[c(-1,-2,-3),]
Uadd_plot_sw<- ggplot(GBU_nh3_datq_230710, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBU_nh3_datq_230615a$sw))
sw_sd <- sd(na.omit(GBU_nh3_datq_230615a$sw))
v <- mean(na.omit(GBU_nh3_datq_230615a$Vf_add_int))
v_sd <- sd(na.omit(GBU_nh3_datq_230615a$Vf_add_int))




### 9
GBU_nh3_datq_230808 <- GBU_nh3_datq%>%
  filter(date==as.Date("2023-08-08"))

GBU_nh3_datq_230808<-GBU_nh3_datq_230808[c(-11,-12,-13,-14, -15, -16),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBU_nh3_datq_230808, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBU_nh3_datq_230808$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBU"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_6<- ggplot(GBU_nh3_datq_230808, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBU_nh3_datq_230808$TMR_NH3 * 1000), 
           y = min(GBU_nh3_datq_230808$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



GBU_nh3_datq_230615a<-GBU_nh3_datq_230808[c(-1,-2,-3),]
Uadd_plot_sw<- ggplot(GBU_nh3_datq_230808, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBU_nh3_datq_230615a$sw))
sw_sd <- sd(na.omit(GBU_nh3_datq_230615a$sw))
v <- mean(na.omit(GBU_nh3_datq_230615a$Vf_add_int))
v_sd <- sd(na.omit(GBU_nh3_datq_230615a$Vf_add_int))



library(ggpubr)

GBL_nh4_grid <- ggarrange(Uadd_plot_1,
                          Uadd_plot_2,
                          Uadd_plot_3,
                          Uadd_plot_4,
                          Uadd_plot_5,
                          Uadd_plot_6,
                      ncol = 3, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")


ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/GBU_nh4_mmfits_grid.png", plot = GBL_nh4_grid, width = 8, height = 4, units = "in")



