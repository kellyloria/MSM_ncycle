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


# BWL_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output_v2.rds")
# names(BWL_no3_datq)
# str(BWL_no3_datq)
# BWL_no3_datq$site <- "BWL"

BWL_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_output_v2.rds")
names(BWL_no3_datq)
str(BWL_no3_datq)
BWL_no3_datq$site <- "BWL"

BWU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output.rds")
names(BWU_nh3_datq)
str(BWU_nh3_datq)

############
dates<- unique(BWL_no3_datq$date)
BWL_no3_datq$Uadd_int1 <- ifelse(is.na(BWL_no3_datq$Uadd_int) | is.nan(BWL_no3_datq$Uadd_int), 0.000001, BWL_no3_datq$Uadd_int)

### 1
BWL_no3_datq_210728 <- BWL_no3_datq%>%
  filter(date==as.Date("2021-07-28"))

BWL_no3_datq_210728<-BWL_no3_datq_210728[c(-11,-12,-13),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_210728, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_210728$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_1<- ggplot(BWL_no3_datq_210728, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_210728$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_210728$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black") +
  facet_grid(.~date)


BWL_no3_datq_210728a<-BWL_no3_datq_210728[c(-1,-2,-3,-4),]
Uadd_plot_sw<- ggplot(BWL_no3_datq_210728a, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_210728a$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_210728a$sw))
v <- mean(na.omit(BWL_no3_datq_210728a$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_210728a$Vf_add_int))


### 2
BWL_no3_datq_220526 <- BWL_no3_datq%>%
  filter(date==as.Date("2022-05-26"))

BWL_no3_datq_220526<-BWL_no3_datq_220526[c(-8,-9,-10,-11,-12),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_220526, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_220526$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_2<- ggplot(BWL_no3_datq_220526, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  #geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000),lty=3) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  # annotate("text", x = max(BWL_no3_datq_220526$TMR_NO3 * 1000), 
  #          y = min(BWL_no3_datq_220526$Uadd_int1 * 1000), 
  #          label = label_text, 
  #          hjust = 1, vjust = 0, size = 4, color = "grey25") +
  facet_grid(.~date)


### 3
BWL_no3_datq_220824 <- BWL_no3_datq%>%
  filter(date==as.Date("2022-08-24"))

BWL_no3_datq_220824<-BWL_no3_datq_220824[c(-11, -12,-13),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_220824, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_220824$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_3<- ggplot(BWL_no3_datq_220824, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_220824$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_220824$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_no3_datq_220824a<-BWL_no3_datq_220824[c(-1),]
Uadd_plot_sw<- ggplot(BWL_no3_datq_220824a, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_220824a$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_220824a$sw))
v <- mean(na.omit(BWL_no3_datq_220824a$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_220824a$Vf_add_int))




### 4 
BWL_no3_datq_221012 <- BWL_no3_datq%>%
  filter(date==as.Date("2022-10-12"))

BWL_no3_datq_221012<-BWL_no3_datq_221012[c(-18,-19,-20),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_221012, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_221012$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_4<- ggplot(BWL_no3_datq_221012, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_221012$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_221012$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_no3_datq_221012a<-BWL_no3_datq_221012[c(-1),]
Uadd_plot_sw<- ggplot(BWL_no3_datq_221012, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_221012$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_221012$sw))
v <- mean(na.omit(BWL_no3_datq_221012$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_221012$Vf_add_int))


### 5
BWL_no3_datq_221121 <- BWL_no3_datq%>%
  filter(date==as.Date("2022-11-21"))

BWL_no3_datq_221121<-BWL_no3_datq_221121[c(-15, -16, -17, -18, -19),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_221121, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_221121$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_5<- ggplot(BWL_no3_datq_221121, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  annotate("text", x = max(BWL_no3_datq_221121$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_221121$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



BWL_no3_datq_221121a<-BWL_no3_datq_221121[c(-1, -2),]
Uadd_plot_sw<- ggplot(BWL_no3_datq_221121a, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_221121a$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_221121a$sw))
v <- mean(na.omit(BWL_no3_datq_221121a$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_221121a$Vf_add_int))




### 5
BWL_no3_datq_221219 <- BWL_no3_datq%>%
  filter(date==as.Date("2022-12-19"))

BWL_no3_datq_221219<-BWL_no3_datq_221219[c(-1,-15,-16,-17,-18, -19,-20),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_221219, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_221219$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_6 <- ggplot(BWL_no3_datq_221219, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_221219$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_221219$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



BWL_no3_datq_221121a<-BWL_no3_datq_221121[c(-1, -2),]
Uadd_plot_sw<- ggplot(BWL_no3_datq_221219, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_221219$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_221219$sw))
v <- mean(na.omit(BWL_no3_datq_221219$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_221219$Vf_add_int))

### 7
BWL_no3_datq_230215 <- BWL_no3_datq%>%
  filter(date==as.Date("2023-02-15"))

BWL_no3_datq_230215<-BWL_no3_datq_230215[c(-13,-14, -15, -16, -17,-18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_230215, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_230215$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_7 <- ggplot(BWL_no3_datq_230215, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_230215$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_230215$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)

BWL_no3_datq_221121a<-BWL_no3_datq_221121[c(-1, -2),]
Uadd_plot_sw<- ggplot(BWL_no3_datq_230215, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_230215$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_230215$sw))
v <- mean(na.omit(BWL_no3_datq_230215$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_230215$Vf_add_int))



### 8
BWL_no3_datq_230405 <- BWL_no3_datq%>%
  filter(date==as.Date("2023-04-05"))

BWL_no3_datq_230405<-BWL_no3_datq_230405[c(-11,-12, -13),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_230405, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_230405$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_8<- ggplot(BWL_no3_datq_230405, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_230405$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_230405$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)

BWL_no3_datq_230405a<-BWL_no3_datq_230405[c(-1, -2, -3),]

Uadd_plot_sw<- ggplot(BWL_no3_datq_230405a, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_230405a$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_230405a$sw))
v <- mean(na.omit(BWL_no3_datq_230405a$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_230405a$Vf_add_int))



### 8
BWL_no3_datq_23078 <- BWL_no3_datq%>%
  filter(date==as.Date("2023-07-18"))

BWL_no3_datq_23078<-BWL_no3_datq_23078[c(-11,-12,-13),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_23078, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_23078$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_9<- ggplot(BWL_no3_datq_23078, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_23078$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_23078$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)



Uadd_plot_sw<- ggplot(BWL_no3_datq_23078, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_23078$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_23078$sw))
v <- mean(na.omit(BWL_no3_datq_23078$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_23078$Vf_add_int))



### 10
BWL_no3_datq_230810 <- BWL_no3_datq%>%
  filter(date==as.Date("2023-08-10"))

BWL_no3_datq_230810<-BWL_no3_datq_230810[c(-14,-15,-16,-17),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_230810, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_230810$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_10<- ggplot(BWL_no3_datq_230810, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_230810$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_230810$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_no3_datq_230810a<-BWL_no3_datq_230810[c(-1),]

Uadd_plot_sw<- ggplot(BWL_no3_datq_230810a, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_230810a$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_230810a$sw))
v <- mean(na.omit(BWL_no3_datq_230810a$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_230810a$Vf_add_int))


### 10
BWL_no3_datq_230925 <- BWL_no3_datq%>%
  filter(date==as.Date("2023-09-25"))

BWL_no3_datq_230925<-BWL_no3_datq_230925[c(-15,-16,-17,-18),]

model.drm1 <- drm (Uadd_int1 ~ TMR_NO3, data = BWL_no3_datq_230925, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(TMR_NO3 = seq(0, max(BWL_no3_datq_230925$TMR_NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)
mm2$site <- "BWL"

## summaries to add to plot:
summary_drm <- summary(model.drm1)
est_d <- c(summary_drm$coefficients["d:(Intercept)", "Estimate"]*1000)
std_error_d <- c(summary_drm$coefficients["d:(Intercept)", "Std. Error"]*1000)
est_e <- c(summary_drm$coefficients["e:(Intercept)", "Estimate"])
std_error_e <- c(summary_drm$coefficients["e:(Intercept)", "Std. Error"])

label_text <- sprintf("Uadd: %.4f ± %.4f", est_d, std_error_d)

Uadd_plot_11<- ggplot(BWL_no3_datq_230925, aes(x = TMR_NO3*1000, y = Uadd_int1*1000, color = site)) +
  geom_line(data = mm2, aes(x = TMR_NO3*1000, y = Uadd*1000)) +
  geom_point(size = 2, shape = 19) +  # Color will come from site mapping
  scale_color_manual(values = site_colors) +
  labs(y=expression(U[add]~(μg~L^-1~s^-1)), x= expression(TMR~NO[3]~(μg~s^-1~L^-3))) + 
  theme_bw() +
  annotate("text", x = max(BWL_no3_datq_230925$TMR_NO3 * 1000), 
           y = min(BWL_no3_datq_230925$Uadd_int1 * 1000), 
           label = label_text, 
           hjust = 1, vjust = 0, size = 4, color = "black")+
  facet_grid(.~date)


BWL_no3_datq_230925a<-BWL_no3_datq_230925[c(-1, -2, -3, -4),]

Uadd_plot_sw<- ggplot(BWL_no3_datq_230925a, aes(x = TMR_NO3*1000, y = sw, color = site)) +
  geom_point()
Uadd_plot_sw

sw <- mean(na.omit(BWL_no3_datq_230925a$sw))
sw_sd <- sd(na.omit(BWL_no3_datq_230925a$sw))
v <- mean(na.omit(BWL_no3_datq_230925a$Vf_add_int))
v_sd <- sd(na.omit(BWL_no3_datq_230925a$Vf_add_int))





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
                          Uadd_plot_11,
                      ncol = 3, nrow = 4,
                      common.legend = TRUE, 
                      legend = "bottom")


ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/BWL_no3_mmfits_grid.png", plot = GBL_nh4_grid, width = 8.5, height = 7, units = "in")



