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

GBL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output_v2.rds")
names(GBL_nh3_datq)
str(GBL_nh3_datq)
GBL_nh3_datq$site <- "GBL"

############
unique(GBL_nh3_datq$date)
GBL_nh3_datq$Uadd_int1 <- ifelse(is.na(GBL_nh3_datq$Uadd_int) | is.nan(GBL_nh3_datq$Uadd_int), 0.000001, GBL_nh3_datq$Uadd_int)


### 2
GBL_nh3_datq_220623 <- GBL_nh3_datq%>%
  filter(date==as.Date("2022-06-23"))

GBL_nh3_datq_220623<-GBL_nh3_datq_220623[c(-8,-9,-10,-11,-12),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_220623, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_220623$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_2<- ggplot(GBL_nh3_datq_220623, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBL_nh3_datq_220623$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_220623$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black") +
  facet_grid(.~date)


GBL_nh3_datq_220623a<-GBL_nh3_datq_220623[c(-1, -2),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_220623a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_220623a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_220623a$sw))
v <- mean(na.omit(GBL_nh3_datq_220623a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_220623a$Vf_add_int))



### 1
GBL_nh3_datq_220407 <- GBL_nh3_datq%>%
  filter(date==as.Date("2022-04-07"))

GBL_nh3_datq_220407<-GBL_nh3_datq_220407[c(-12,-13,-14,-15, -16,-17, -18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_220407, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_220407$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_1<- ggplot(GBL_nh3_datq_220407, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  annotate("text", x = max(GBL_nh3_datq_220407$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_220407$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black") +
  facet_grid(.~date)


GBL_nh3_datq_220407a<-GBL_nh3_datq_220407[c(-1, -2,-3),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_220407a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_220407a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_220407a$sw))
v <- mean(na.omit(GBL_nh3_datq_220407a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_220407a$Vf_add_int))



### 3
GBL_nh3_datq_220103 <- GBL_nh3_datq%>%
  filter(date==as.Date("2022-10-03"))

GBL_nh3_datq_220103<-GBL_nh3_datq_220103[c(-16,-17,-18,-19,-20,-21, -22,-23, -24),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_220103, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_220103$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_3<- ggplot(GBL_nh3_datq_220103, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBL_nh3_datq_220103$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_220103$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


GBL_nh3_datq_220103a<-GBL_nh3_datq_220103[c(-1, -2,-3),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_220103a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_220103a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_220103a$sw))
v <- mean(na.omit(GBL_nh3_datq_220103a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_220103a$Vf_add_int))





### 4
GBL_nh3_datq_221104 <- GBL_nh3_datq%>%
  filter(date==as.Date("2022-11-04"))

GBL_nh3_datq_221104<-GBL_nh3_datq_221104[c(-15,-16,-17,-18, -19,-20, -21),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_221104, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_221104$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_4<- ggplot(GBL_nh3_datq_221104, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBL_nh3_datq_221104$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_221104$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



GBL_nh3_datq_221104a<-GBL_nh3_datq_221104[c(-1, -2,-3,-4,-5),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_221104a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_221104a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_221104a$sw))
v <- mean(na.omit(GBL_nh3_datq_221104a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_221104a$Vf_add_int))

### 5
GBL_nh3_datq_221212 <- GBL_nh3_datq%>%
  filter(date==as.Date("2022-12-12"))

GBL_nh3_datq_221212<-GBL_nh3_datq_221212[c(-12,-13,-14,-15, -16,-17, -18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_221212, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_221212$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_5<- ggplot(GBL_nh3_datq_221212, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  annotate("text", x = max(GBL_nh3_datq_221212$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_221212$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


GBL_nh3_datq_221212a<-GBL_nh3_datq_221212[c(-1, -2,-3),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_221212a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_221212a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_221212a$sw))
v <- mean(na.omit(GBL_nh3_datq_221212a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_221212a$Vf_add_int))


### 6
GBL_nh3_datq_230327 <- GBL_nh3_datq%>%
  filter(date==as.Date("2023-03-27"))

GBL_nh3_datq_230327<-GBL_nh3_datq_230327[c(-7,-8,-9,-10,-11,-12,-13,-14 -15,-16,-17,-18,-19),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_230327, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_230327$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_6<- ggplot(GBL_nh3_datq_230327, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  #geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  # annotate("text", x = max(GBL_nh3_datq_230327$TMR_NH3 * 1000), 
  #          y = min(GBL_nh3_datq_230327$Uadd_int1 * 1000), 
  #          label = label_text, 
  #          hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


GBL_nh3_datq_230327a<-GBL_nh3_datq_230327[c(-1, -2,-3,-4),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_230327a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_230327a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_230327a$sw))
v <- mean(na.omit(GBL_nh3_datq_230327a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_230327a$Vf_add_int))




### 7
GBL_nh3_datq_230615 <- GBL_nh3_datq%>%
  filter(date==as.Date("2023-06-15"))

GBL_nh3_datq_230615<-GBL_nh3_datq_230615[c(-7,-8,-9,-10,-11,-12,-13, -14),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_230615, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_230615$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_7<- ggplot(GBL_nh3_datq_230615, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  annotate("text", x = max(GBL_nh3_datq_230615$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_230615$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)




GBL_nh3_datq_230615a<-GBL_nh3_datq_230615[c(-1, -2),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_230615a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_230615a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_230615a$sw))
v <- mean(na.omit(GBL_nh3_datq_230615a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_230615a$Vf_add_int))




### 8
GBL_nh3_datq_230710 <- GBL_nh3_datq%>%
  filter(date==as.Date("2023-07-10"))

GBL_nh3_datq_230710<-GBL_nh3_datq_230710[c(-10,-11,-12,-13, -14, -15),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_230710, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_230710$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_8<- ggplot(GBL_nh3_datq_230710, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBL_nh3_datq_230710$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_230710$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


GBL_nh3_datq_230710a<-GBL_nh3_datq_230710[c(-1, -2, -3),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_230710a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_230710a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_230710a$sw))
v <- mean(na.omit(GBL_nh3_datq_230710a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_230710a$Vf_add_int))




### 9
GBL_nh3_datq_230808 <- GBL_nh3_datq%>%
  filter(date==as.Date("2023-08-08"))

GBL_nh3_datq_230808<-GBL_nh3_datq_230808[c(-11,-12,-13,-14, -15, -16, -17, -18, -19),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NH3, data = GBL_nh3_datq_230808, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NH3 = seq(0, max(GBL_nh3_datq_230808$TMR_NH3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "GBL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.3f ± %.3f", est_d, std_error_d)

Uadd_plot_9<- ggplot(GBL_nh3_datq_230808, aes(x = TMR_NH3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NH3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 17) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NH[4]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(GBL_nh3_datq_230808$TMR_NH3 * 1000), 
           y = min(GBL_nh3_datq_230808$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


GBL_nh3_datq_230808a<-GBL_nh3_datq_230808[c(-1, -2, -3),]
Uadd_plot_sw<- ggplot(GBL_nh3_datq_230808a, aes(x = TMR_NH3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(GBL_nh3_datq_230808a$sw))
sw_sd <- sd(na.omit(GBL_nh3_datq_230808a$sw))
v <- mean(na.omit(GBL_nh3_datq_230808a$Vf_add_int))
v_sd <- sd(na.omit(GBL_nh3_datq_230808a$Vf_add_int))



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
                      ncol = 3, nrow = 3,
                      common.legend = TRUE, 
                      legend = "bottom")


ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/GBL_nh4_mmfits_grid_25.png", plot = GBL_nh4_grid, width = 8.5, height = 6.5, units = "in")
