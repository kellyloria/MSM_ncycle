#=========== Preliminaries
rm(list=ls())
# load packages
library(tidyverse)
library(lubridate)
library(plotly)
library(devtools)
library(nrlmetab)
library(zoo)
library(suncalc)
library(readxl)
library(patchwork)
library(gridExtra)
library(dplyr)
library(reshape)
library(scales)

# /Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/BTC_out/plot_out


## N uptake figures

rawdat = list.files(paste("./BTC_out/plot_output/plot_out",sep=""), full.names = T)

rawdat_A = list.files(paste("./BTC_out/Archived output/",sep=""), full.names = T)


GBL_06 <- read.csv("./BTC_out/plot_output/GBL_BTC_20210603.csv")
GBL_06$Site<- "GBU"

## BWL NH4 ##
BW_NH4_08 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL220824.csv")
BW_NH4_08$Site<-"BWL"

BW_NH4_11 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL221121.csv")
BW_NH4_11$Site<-"BWL"  

BW_NH4_05 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL220526.csv")
BW_NH4_05$Site<-"BWL" 

BW_NH4_10 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL221012.csv")
BW_NH4_10$Site<-"BWL" 

BW_NH4_0723 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL230718.csv")
BW_NH4_0723$Site<-"BWL" 

BW_NH4_0423 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL230405.csv")
BW_NH4_0423$Site<-"BWL" 

BW_NH4_0223 <- read.csv("./BTC_out/plot_output/BWL_NH4_BTC_BWL230215.csv")
BW_NH4_0223$Site<-"BWL" 

BWL_NH4<- rbind(BW_NH4_05, BW_NH4_08, BW_NH4_10, BW_NH4_11, BW_NH4_0723, BW_NH4_0423, BW_NH4_0223)

## 
BW_NO3_10 <- read.csv("./BTC_out/plot_output/BWL_NO3_BTC_BWL221012.csv")
BW_NO3_10$Site<-"BWL"

BW_NO3_05 <- read.csv("./BTC_out/plot_output/BWL_NO3_BTC_BWL220526.csv")
BW_NO3_05$Site<-"BWL"

BW_NO3_11 <- read.csv("./BTC_out/plot_output/BWL_NO3_BTC_BWL221121.csv")
BW_NO3_11$Site<-"BWL"

BW_NO3_07 <- read.csv("./BTC_out/plot_output/BW_BTC_20210728.csv")
BW_NO3_07$Site<-"BWL"

BW_NO3_07 <- read.csv("./BTC_out/plot_output/BW_BTC_20210728.csv")
BW_NO3_07$Site<-"BWL"

BW_NO3_0723 <- read.csv("./BTC_out/plot_output/BWL_NO3_BTC_BWL230718.csv")
BW_NO3_0723$Site<-"BWL" 

BW_NO3_0423 <- read.csv("./BTC_out/plot_output/BWL_NO3_BTC_BWL230405.csv")
BW_NO3_0423$Site<-"BWL" 

BW_NO3_0223 <- read.csv("./BTC_out/plot_output/BWL_NO3_BTC_BWL230215.csv")
BW_NO3_0223$Site<-"BWL" 


BWL_NO3 <- rbind(BW_NO3_10,BW_NO3_05,BW_NO3_11,BW_NO3_07, BW_NO3_0723, BW_NO3_0423, BW_NO3_0223)





# INC_08 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221012.csv")
# INC_08[3,10]<-NA
# INC_08$Site<-"Incline"
# GBU_07 <- read.csv("./BTC_out/GB_BTC_20210722.csv")
# GBU_07$Site<-"GBL"
# GEN_07 <- read.csv("./BTC_out/GEN_BTC_20210728.csv")
# GEN_07$Site<-"General"


## plot in 4 loop
tbl <-
  list.files(path = "./BTC_out/plot_output/",
             pattern = "*.csv", 
             full.names = T) %>% 
  map_df(~read_csv(., col_types = cols(.default = "c"))) %>%
  mutate(
    Datetime= as.POSIXct(as.character(datetime), format="%m/%d/%y %H:%M"),
    date= as.Date(Datetime))
#   
# # Define the mapping of old values to new values
# value_mapping <- c("BWL" = "BWL", "BWL_NH4" = "BWL", "BWL_NO3" = "BWL",
#                    "BWU" = "BWU", "GB" = "GBL", "GBL" = "GBL",
#                    "GBL_NH4" = "GBL", "GBL_NO3" = "GBL", "GBU" = "GBU",
#                    "GBU_NH4" = "GBU")
# 
# # Replace values in the "Site" column using case_when
# tbl <- tbl %>% 
#   mutate(Site = case_when(
#     Site %in% names(value_mapping) ~ value_mapping[Site],
#     TRUE ~ Site  # Keep unchanged if not in value_mapping
#   ))

# Print the modified tbl structure
str(tbl)



tbl$NH4 <- as.numeric(tbl$NH4)
tbl$Cl <- as.numeric(tbl$Cl)
tbl$kw <- as.numeric(tbl$kw)
tbl$sw <- as.numeric(tbl$sw)
tbl$Uadd <- as.numeric(tbl$Uadd)
tbl$NO3 <- as.numeric(tbl$NO3)

unique(tbl$Site)

tbld<- as.data.frame(tbl)
names(tbld)
str(tbl)


tbl_ave <- tbld %>%
  dplyr::group_by(Site, date, sp) %>% 
  dplyr::summarise(
          NH4_m = mean(NH4, na.rm = T),
          NH4_sd = sd(NH4, na.rm = T),
          Cl_m = mean(Cl, na.rm = T),
          Cl_sd = sd(Cl, na.rm = T),
          kw_m= mean(kw, na.rm = T),
          kw_sd= sd(kw, na.rm = T),
          sw_m=mean(sw, na.rm = T),
          sw_sd=sd(sw, na.rm = T),
          Uadd_m=mean(Uadd, na.rm = T),
          Uadd_sd=sd(Uadd, na.rm = T),
          NO3_m= mean(NO3, na.rm = T),
          NO3_sd= sd(NO3, na.rm = T))
    



tbl_ave1 <- tbl_ave %>%
  mutate(label = case_when(
    Site %in% c("BWU", 
                "BWL"
    ) ~ "west",
    Site %in% c("GBL", 
                "GBU"
    ) ~ "east",
    TRUE ~ "unknown"  # Assign "unknown" for unmatched site IDs (optional)
  ))


# RS_dat



Uaddplotlog <- ggplot(tbl_ave1, aes(x=date, y=log(Uadd_m +1), colour=Site, shape=sp)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line() +
  geom_point(size=3) +
  geom_pointrange(aes(ymin=log(Uadd_m+1)-log(Uadd_sd+1), ymax=log(Uadd_m+1)+log(Uadd_sd+1))) +
  labs(y= expression(paste("log(Uadd) (mg m^-2 min^-1)")),
       x= expression(paste("Date")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    #"#2a7191", 
    "#77bad9",
    "#DAA520", 
    #"#DAA520",
    #"#d1bb84",
    "#d1bb84"),0.5)) +
  theme_bw() +
  facet_grid(label~.)


Uaddplot <- ggplot(tbl_ave, aes(x=date, y=(Uadd_m), colour=Site, shape=sp)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line() +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=(Uadd_m)-(Uadd_sd), ymax=(Uadd_m)+(Uadd_sd))) +
  labs(y= expression(paste("Uadd (mg m^-2 min^-1)")),
       x= expression(paste("Date")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520",
    "#d1bb84"),0.7)) +
  theme_bw() +
  facet_grid(sp~.)

# ggsave(plot = Uaddplot, filename = paste("./figures/Uaddplot.jpeg",sep=""),width=5.5,height=4.5,dpi=300)

swplot <- ggplot(tbl_ave1, aes(x=date, y= sw_m, colour=Site, shape=sp)) +
  ylim(0,260) +
  geom_line() +
  geom_point(size = 3) +
   geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd)) +
  labs(y = expression(paste("Sw (m)")),
       x= expression(paste("Date")))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520",
    "#d1bb84"),0.7)) +
  theme_bw() + 
  facet_grid(label~.)


Uaddplot_b <- ggplot(tbl_ave1, aes(x=date, y= Uadd_m, colour=Site, shape=sp)) +
  #ylim(0,300) +
  geom_line() +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=Uadd_m-Uadd_sd, ymax=Uadd_m+Uadd_sd)) +
  labs(y= expression(paste("Uadd (mg m^-2 min^-1)")),
       x= expression(paste("Date")))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520",
    "#d1bb84"),0.7)) +
  theme_bw() + 
  facet_grid(label~.)

dim(tbl_ave1)


BTC_grid <- ggarrange(Uaddplot_b,
                            swplot,
                            ncol = 2, nrow = 1,
                            common.legend = TRUE, 
                      labels = c("A", "B"),
                            widths = c(1,1, 1),
                            legend = "bottom")

# ggsave(plot = BTC_grid, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/BTC_grid_grid.png",sep=""),width=6.5,height=3.25,dpi=300)




##
# Daily data
RAFDM <- read.csv("/Users/kellyloria/Downloads/RSAFMD.csv")
RAFDM_W <- RAFDM %>%
  mutate(date = as.Date((date), format ="%Y-%m-%d")) %>%
  subset(Site== "BWU" |Site== "BWL") 
RAFDM_W$label<- "west"

RAFDM_E <- RAFDM %>%
  mutate(date = as.Date((date), format ="%Y-%m-%d")) %>%
  subset(Site== "GBU"|Site== "GBL") 
RAFDM_E$label <- "east"

RAFDM_all <- rbind(RAFDM_W,RAFDM_E)
names(RAFDM_all)
RAFDM_all$AFDM..mg.cm..2. <- as.numeric(RAFDM_all$AFDM..mg.cm..2.)
summary(RAFDM_all)

RAFDM_all <- RAFDM_all%>%
  dplyr::select(Site, label, date, AFDM..mg.cm..2.)


## 
# GPP 

## BWL 
metabdat <- list.files(paste("/Users/kellyloria/Documents/UNR/MSMmetab/23_output/",sep=""), full.names = T) %>%
  lapply(read_csv) %>%
  bind_rows()
names(metabdat)
summary(metabdat)

hist(metabdat$GPP_mean)
hist(metabdat$GPP_daily_mean)

metdat <- metabdat %>%
  mutate(date= as.Date(date, format="%m/%d/%y")) %>%
  dplyr::select(Site, shore, date, 
                GPP_mean, GPP_sd,
                ER_mean, ER_sd)
summary(metdat)

### 

newdat <- left_join(tbl_ave1, RAFDM_all, by = c("date", "label", "Site"))
newdat <- left_join(newdat, metdat, by = c("date", "Site"))

newdat_KNO3 <- newdat %>%
  filter(sp == "KNO3")

newdat_KNO3BL<- newdat_KNO3 %>%
  filter(Site == "BWL")
summary(newdat_KNO3BL)



newdat_KNO3GBL<- newdat_KNO3 %>%
  filter(Site == "GBL")
summary(newdat_KNO3GBL)


newdat_NH4Cl <- newdat %>%
  filter(sp == "NH4Cl")

newdat_NH4ClBL<- newdat_NH4Cl %>%
  filter(Site == "BWL")
summary(newdat_NH4ClBL)

newdat_NH4ClGB<- newdat_NH4Cl %>%
  filter(Site == "GBL")
summary(newdat_NH4ClGB)


AFDM_sw <- ggplot(newdat, aes(x=AFDM..mg.cm..2., y= sw_m, colour=Site, shape=sp)) +
  ylim(0,360) +
 # geom_line() +
  geom_smooth(method="lm", se=F)+
  geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd)) +
  geom_point(size = 3) +
  #geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd)) +
  labs(y = expression(paste("Sw (m)")),
       x= expression(paste("Algal AFDM (mg/cm2) ")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() +
  facet_grid(sp~.)


GPP_sw <- ggplot(newdat, aes(x=GPP_mean, y= sw_m)) +
  ylim(0,250) + geom_smooth(method="lm", se=F, color="grey50")+
  geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd,colour=Site)) +
  geom_errorbarh(aes(xmin=GPP_mean-GPP_sd, xmax=GPP_mean+GPP_sd, colour=Site)) +
  geom_point(size = 3, aes(colour=Site, shape=sp)) +
   labs(y = expression(paste("Sw (m)")),
        x= expression(paste("GPP")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() 

ER_sw <- ggplot(newdat, aes(x=ER_mean, y= sw_m, colour=Site, shape=sp)) +
  ylim(0,250) +
  geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd)) +
  geom_errorbarh(aes(xmin=ER_mean-ER_sd, xmax=ER_mean+ER_sd)) +
  geom_point(size = 3) +  
  labs(y = expression(paste("Sw (m)")),
       x= expression(paste("ER")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() 


AFDM_udd <- ggplot(newdat, aes(x=AFDM..mg.cm..2., y= Uadd_m, colour=Site, shape=sp)) +
  #ylim(0,20) + xlim(0,1500) +
  # geom_line() +
  geom_smooth(method="lm", se=F)+
  geom_pointrange(aes(ymin=Uadd_m-Uadd_sd, ymax=Uadd_m+Uadd_sd)) +
  geom_point(size = 3) +
  #geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd)) +
  # labs(y = expression(paste("Sw (m)")),
  #      x= expression(paste("Algal AFDM (mg/cm2) ")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() +
  facet_grid(sp~.)

newdat$AFDM_udd_ymin <- pmax(newdat$Uadd_m - newdat$Uadd_sd, 0)


GPP_udd <- ggplot(newdat, aes(x=GPP_mean, y= Uadd_m)) +
  geom_pointrange(aes(ymin=AFDM_udd_ymin, ymax=Uadd_m+Uadd_sd, colour=Site)) +
  geom_errorbarh(aes(xmin=GPP_mean-GPP_sd, xmax=GPP_mean+GPP_sd, colour=Site)) +
  geom_point(size = 3, aes(colour=Site, shape=sp)) +
  geom_smooth(method="lm", se=F, color="grey50")+
  labs(y = expression(paste("Uadd (mg m^-2 min^-1)")),
       x= expression(paste("GPP")))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() 

ER_udd <- ggplot(newdat, aes(x=ER_mean, y= Uadd_m, colour=Site, shape=sp)) +
  geom_pointrange(aes(ymin=AFDM_udd_ymin, ymax=Uadd_m+Uadd_sd)) +
  geom_errorbarh(aes(xmin=ER_mean-ER_sd, xmax=ER_mean+ER_sd)) +
  geom_point(size = 3) +
  labs(y = expression(paste("Uadd (mg m^-2 min^-1)")),
       x= expression(paste("ER")))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() 

# ggsave(plot = AFDM_sw, filename = paste("./figures/udd_AFDM.jpeg",sep=""),width=5.5,height=4.5,dpi=300)


## total grid:
NAupGrid_grid <- ggarrange(
  GPP_udd,
  ER_udd,
  GPP_sw,
  ER_sw, 
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = TRUE, 
  legend = "bottom")

# ggsave(plot = NAupGrid_grid, filename = paste("./figures/Metab_BTC.png",sep=""),width=5,height=5.25,dpi=300)


gpp_up<- lmer(Uadd_m~ scale(GPP_mean) + (1|sp), data=newdat)
summary(gpp_up)

gpp_up<- lm(Uadd_m~ scale(ER_mean), data=newdat)
summary(gpp_up)

gpp_sw<- lm(sw_m~ scale(ER_mean), data=newdat)
summary(gpp_sw)

gpp_sw<- lm(sw_m~ scale(GPP_mean), data=newdat)
summary(gpp_sw)



AFDM_swNH <- ggplot(newdat_NH4Cl, aes(x=AFDMmgcm2, y= sw_m, colour=Site, shape=sp)) +
  #ylim(0,20) + xlim(0,1500) +
  # geom_line() +
  geom_smooth(method="lm", se=F)+
  geom_point(size = 3) +
  #geom_pointrange(aes(ymin=sw_m-sw_sd, ymax=sw_m+sw_sd)) +
  labs(y = expression(paste("Sw (m)")),
       x= expression(paste("Algal AFDM (mg/cm2) ")))+
  #scale_shape_manual(values = c(17, 15, 17, 17, 15, 17, 15))+
  scale_colour_manual(values = alpha(c(
    "#3283a8", 
    "#77bad9",
    "#DAA520", 
    "#d1bb84"), 0.5)) +
  theme_bw() 

library(ggpubr)
## total grid:
all_grid <- ggarrange(swplot, 
                      Uaddplotlog, 
                      ncol = 1, nrow = 2,
                      labels = c("A", "B"),
                      hjust = c(-5, -5),
                      vjust = c(18, 18),
                      widths = c(1,1, 0.8),
                      common.legend = TRUE, 
                      legend = "right")


# ggsave(plot = all_grid, filename = paste("./figures/Uaddswplot2.jpeg",sep=""),width=5.5,height=6,dpi=300)








BWL_NH4 <- tbl %>%
  subset(Site=="BWL_NH4")

BWL_NO3 <- tbl %>%
  subset(Site=="BWL_NO3")

BWU_NH4 <- tbl %>%
  subset(Site=="BWU_NH4")

GBL_NH4 <- tbl %>%
  subset(Site=="GBL_NH4")

GBL_NO3 <- tbl %>%
  subset(Site=="GBL_NO3")

GBU_NH4 <- tbl %>%
  subset(Site=="GBU_NH4")

GBU_NO3 <- tbl %>%
  subset(Site=="GBU_NO3")


# 
# plot_dat <- data.frame(Site = NA, date=as.Date(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
# for (i in 2:nrow(datq)) {
#   temp_dat <- datq[c(i-1,i),]
#   slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
#   kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
#   datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
#   NH4<- datq$Nh4_C[i]
#   Cl<- datq$Cl_mgL[i]
#   temp_out <- data.frame(Site = "BWL_NH4", 
#                          stamps = paste(i, i-1, sep = "-"), 
#                          slope_sample = slope_sample, 
#                          kw=kw, 
#                          datetime=datetime,
#                          NH4=NH4,
#                          Cl=Cl)
#   out <- rbind(out, temp_out)
# }
# 
# ## Cadd geometric mean of background concetrations 
# out <- out[c(-1,-2, -3),]
# out$sw <- -1/(out$kw)
# out$Uadd <- Q*Cadd/out$sw*w




Uadd_plotBW <- ggplot(BWL_NH410, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  #ggtitle("Blackwood L Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))



#### Nutrient uptake rate for GB
library(ggplot2)
library(ggformula) # optional
library(MASS)
library(drc)
library(gridExtra)


rawdat = list.files(paste("./BTC_out/",sep=""), full.names = T)

# plot BWL individually #
BWL_NH410 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221012.csv")

model.drm1 <- drm (Uadd ~ NH4, data = BWL_NH410, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_NH410$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

## BW 

Uadd_plotBW <- ggplot(BWL_NH410, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  #ggtitle("Blackwood L Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))


# # plot BWL individually #
BWL_NH408 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL220824.csv")

model.drm1 <- drm (Uadd ~ NH4, data = BWL_NH408, fct = MM.2())
summary(model.drm1)

mm3 <- data.frame(NO3 = seq(0, max(BWL_NH408$NH4), length.out = 100))
mm3$Uadd <- predict(model.drm1, newdata = mm3)

# # plot BWL individually #
BWL_NH408 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL220824.csv")

model.drm1 <- drm (Uadd ~ NH4, data = BWL_NH408, fct = MM.2())
summary(model.drm1)

mm3 <- data.frame(NO3 = seq(0, max(BWL_NH408$NH4), length.out = 100))
mm3$Uadd <- predict(model.drm1, newdata = mm3)

summary(model.drm1)
(model.drm1)


# # plot BWL individually #
BWL_NH411 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221121.csv")

model.drm1 <- drm(Uadd ~ NH4, data = BWL_NH408, fct = MM.2())
summary(model.drm1)

mm3 <- data.frame(NO3 = seq(0, max(BWL_NH408$NH4), length.out = 100))
mm3$Uadd <- predict(model.drm1, newdata = mm3)

summary(model.drm1)
(model.drm1)








## BW 

Uadd_plotBW2 <- ggplot() +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(data= BWL_NH410, aes(x=NH4, y=Uadd), size = 3, shape= 15, col = "#3283a8") +
  geom_line(data = mm3, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(data= BWL_NH408, aes(x=NH4, y=Uadd), size = 3, shape= 17, col = "#3283a8") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  #ggtitle("Blackwood L Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))





# gen
model.drm1 <- drm (Uadd ~ NO3, data = GEN_07, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GEN_07$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGen <- ggplot(GEN_07, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#226b5b") +
  geom_point(size = 3, shape= 17, col = "#226b5b") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg g m^-2 min^-1)"))) +
  ggtitle("General Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))



## GB ##
model.drm1 <- drm (Uadd ~ NO3, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)



Uadd_plotGB <- ggplot(GBL_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))


## GBU
model.drm1 <- drm (Uadd ~ NO3, data = GBU_07, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBU_07$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)


Uadd_plotGBU <- ggplot(GBU_07, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a84e32") +
  geom_point(size = 3, shape= 17, col = "#a69117") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))




## inc
model.drm1 <- drm (Uadd ~ NO3, data = INC_08, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(INC_08$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)


Uadd_plotINC <- ggplot(INC_08, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a84e32") +
  geom_point(size = 3, shape= 17, col = "#a84e32") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Incline Creek") +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))


library(ggpubr)

## total grid:
all_grid <- ggarrange(Uadd_plotBW, 
                      Uadd_plotGB, 
                      Uadd_plotINC,
                      Uadd_plotGen, 
                      Uadd_plotGBU, 
                      ncol = 3, nrow = 2)
                      #labels = c("A", "","", "B","","", "C", "","", "D", "", ""),
                      #hjust = c(-5, -5, -5,-5),
                      #widths = c(1,1, 0.8),
                      #common.legend = TRUE, 
                      #legend = "bottom")


# ggsave(plot = all_grid, filename = paste("./figures/NO3_Allgrid.png",sep=""),width=9.5,height=5.25,dpi=300)

Alluptake <- rbind(GBL_06, GBU_07, INC_08, GEN_07, BW_07)
Alluptake$date<- as.Date(Alluptake$datetime)

Alluptake_a <- Alluptake %>%
  group_by(Site, date) %>%
  summarise(
    sw_a = mean(sw),
    sw_sd = sd(sw),
    Uadd_a = mean(Uadd, na=T),
    Uadd_sd = sd(Uadd, na=T))



tsac_tb <- read.csv("/Users/kellyloria/Documents/UNR/Ncycle/TASCC\ table.csv")
tsac_tb <- tsac_tb %>%
  mutate(date = as.Date((Date), format ="%m/%d/%y"))

  
GBL_d <- read.csv("/Users/kellyloria/Desktop/LakeTahoeNS/plotDat/GBL_daily.csv")
GBL_d <- GBL_d %>%
  mutate(date = as.Date((date)))
str(GBL_d)

GBL_d$Site <- "GBL" 
BWL_d <- read.csv("/Users/kellyloria/Desktop/LakeTahoeNS/plotDat/BWL_daily.csv")
BWL_d<- BWL_d %>%
  mutate(date = as.Date((date)))
BWL_d$Site <- "BWL" 

metabout<- rbind(GBL_d,BWL_d)



metabplot <- tsac_tb %>%
  dplyr::full_join(metabout,  by = "date")




AFDM <- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/CleanDat/22AFDM.csv")
AFDM <- AFDM %>%
  mutate(date = as.Date((date), format ="%Y-%m-%d"))

AFDM[7,3]<-"2021-07-28"

Alluptake_a<- Alluptake_a%>%
  full_join(tsac_tb)





SW_supply <- ggplot(metabplot, aes(x=mean.Udd..mg.m2.min., y= GPP_mean, colour=Site.x)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=sw_a-sw_sd, ymax=sw_a+sw_sd))+
  labs(x = expression(paste("GPP")),
       y= "Sw (m)") +
  # scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
  #                                "#226b5b", "#a84e32")) +
  theme_classic() 



SW_AFDM <- ggplot(Alluptake_a, aes(x=RS_AFDM, y= sw_a, colour=Site)) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=sw_a-sw_sd, ymax=sw_a+sw_sd))+
  labs(x = expression(paste("Rock scrape AFDM")),
       y= "Sw (m)")+
  scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
                                 "#226b5b", "#a84e32")) +
  theme_classic() 




SW_supply <- ggplot(Alluptake_a, aes(x=NO3.supply..g.day., y= sw_a, colour=Site)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=sw_a-sw_sd, ymax=sw_a+sw_sd))+
  labs(x = expression(paste("Nitrate supply (g/day)")),
       y= "Sw (m)")+ theme_classic()
  # scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
  #                                "#226b5b", "#a84e32")) +
  theme_classic() 



SW_AFDM <- ggplot(Alluptake_a, aes(x=RS_AFDM, y= sw_a, colour=Site)) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=sw_a-sw_sd, ymax=sw_a+sw_sd))+
  labs(x = expression(paste("Rock scrape AFDM")),
       y= "Sw (m)")+
  scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
                                 "#226b5b", "#a84e32")) +
  theme_classic() 


Uadd_supply <- ggplot(Alluptake_a, aes(x=NO3.supply..g.day., y= Uadd_a, colour=Site)) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=Uadd_a-Uadd_sd, ymax=Uadd_a+Uadd_sd))+
  labs(x = expression(paste("Nitrate supply (g/day)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)")))+
  scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
                                 "#226b5b", "#a84e32")) +
  theme_classic() 

Uadd_AFDM <- ggplot(Alluptake_a, aes(x=RS_AFDM, y= Uadd_a, colour=Site)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=Uadd_a-Uadd_sd, ymax=Uadd_a+Uadd_sd))+
  labs(x = expression(paste("Rock scrape AFDM")),
       y= expression(paste("Uadd (mg m^-2 min^-1)")))+
  scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
                                 "#226b5b", "#a84e32")) +
  theme_classic() 


## total grid:
all_grid2 <- ggarrange(SW_supply, 
                      SW_AFDM, 
                      Uadd_supply,
                      Uadd_AFDM, 
                      ncol = 2, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")


# ggsave(plot = all_grid2, filename = paste("./figures/Summary_Allgrid2.png",sep=""),width=6,height=6,dpi=300)

supp_AFDM <- ggplot(Alluptake_a, aes(x=NO3.supply..g.day., y= RS_AFDM, colour=Site)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  #geom_pointrange(aes(ymin=Uadd_a-Uadd_sd, ymax=Uadd_a+Uadd_sd))+
  labs(y = expression(paste("Rock scrape AFDM")),
       x= expression(paste("Nitrate supply (g/day)")))+
  scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
                                 "#226b5b", "#a84e32")) +
  theme_classic() 


#ggsave(plot = supp_AFDM, filename = paste("./figures/NO3supp.png",sep=""),width=4.25,height=3,dpi=300)

