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

## N uptake figures

rawdat = list.files(paste("./BTC_out",sep=""), full.names = T)
GBL_06 <- read.csv("./BTC_out/GBL_BTC_20210603.csv")
GBL_06$Site<- "GBU"
BW_07 <- read.csv("./BTC_out/BW_BTC_20210728.csv")
BW_07$Site<-"BWL"
INC_08 <- read.csv("./BTC_out/INC_BTC_20210805.csv")
INC_08[3,10]<-NA
INC_08$Site<-"Incline"
GBU_07 <- read.csv("./BTC_out/GB_BTC_20210722.csv")
GBU_07$Site<-"GBL"
GEN_07 <- read.csv("./BTC_out/GEN_BTC_20210728.csv")
GEN_07$Site<-"General"


#### Nutrient uptake rate for GB
library(ggplot2)
library(ggformula) # optional
library(MASS)
library(drc)
library(gridExtra)


model.drm1 <- drm (Uadd ~ NO3, data = BW_07, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BW_07$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotBW <- ggplot(BW_07, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Blackwood L Creek") +
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



tsac_tb <- read.csv("/Users/kellyloria/Documents/UNR/Ncycle/TASCCtable_plot.csv")
tsac_tb <- tsac_tb %>%
  mutate(Date = as.Date((Date), format ="%m/%d/%y"))

AFDM <- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/CleanDat/22AFDM.csv")
AFDM <- AFDM %>%
  mutate(date = as.Date((date), format ="%Y-%m-%d"))

AFDM[7,3]<-"2021-07-28"

Alluptake_a<- Alluptake_a%>%
  full_join(tsac_tb)




SW_supply <- ggplot(Alluptake_a, aes(x=NO3.supply..g.day., y= sw_a, colour=Site)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_smooth(method=lm, se=F, color="gray") +
  geom_point(size = 3) +
  geom_pointrange(aes(ymin=sw_a-sw_sd, ymax=sw_a+sw_sd))+
  labs(x = expression(paste("Nitrate supply (g/day)")),
       y= "Sw (m)")+
  scale_colour_manual(values = c("#3283a8", "#a67d17", "#a69117",
                                 "#226b5b", "#a84e32")) +
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

