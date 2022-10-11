library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)
library(lme4)
library(lmerTest)
library(gridExtra)
library(cowplot)

## ---------------------------

# Read in the nitrogen uptake assay data:
BW_NA <- read.csv("/Users/kellyloria/Documents/UNR/Ncycle/NA21_dat/NA21_BW20210728.csv")
BW_NA$datetime <- as.POSIXct(as.character(BW_NA$datetime), format="%Y-%m-%d %H:%M:%S") ## modify the format to match your data
# 
summary(BW_NA)

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


## SPC data from HOBO
BW_Hobo <-read.csv("/Users/kellyloria/Documents/UNR/Ncycle/NA21_dat/GeneralandBW_20210727_20775520_4.csv", skip=1)
summary(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                        "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                        "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
range(BW_Hobo$DateTime)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
BW_Hobo <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 15:30:00") & DateTime <= as.POSIXct("2021-07-28 17:01:50"))

##
## Reach morphology estimates:
##

## Equation:
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 15:31:00") & DateTime <= as.POSIXct("2021-07-28 15:39:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo$DateTime), BW_Hobo$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-28 15:40:00") #Lolomai 
peak_time <- BW_Hobo[which.max(BW_Hobo$SpCond),]$DateTime 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds
v <- 250/time_diff_sec
v

#############################
## Estimate mean depth (z) ##
#############################
## effective depth (z) can be estimated using the following equation
##
## z = Q/(w*v)
## where z is effective depth (m)
## Q is discharge (m^3/sec)
## w is average width (m)
## v is velocity (m/sec)
## Enter average width measurement in m
w <- mean(5.49, 5.12, 3.98, 4.79, 6.1, 
          4.32, 4.45, 3.78, 6.02, 5.33)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## double check flow from USGS
library(dataRetrieval)
siteNo <- "10336660"
pCode <- "00060" #cubic ft per second
start.date <- "2021-07-27"
end.date <- "2021-07-29"

BWUSGS <- readNWISuv(siteNumbers = siteNo,
                     parameterCd = pCode,
                     startDate = start.date,
                     endDate = end.date)

BWUSGS <- renameNWISColumns(BWUSGS)
names(BWUSGS)

BWUSGS$Flow_InstLsec<- c(BWUSGS$Flow_Inst*28.3168466)




##
## Left join hobo conductivity estimates with samples:
##

BW_NA <- left_join(BW_NA, BW_Hobo[c("DateTime", "SpCond", "TempC")],
                   by= c("datetime"="DateTime"))

summary(BW_NA)

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

##
## Calculations Notes:
##
# 1. Correct for background concentrations (_C):
#     "Results" in mg N/L
BW_NA$NO3_C <- (BW_NA$Results-0.005) 

BW_NA$SpCond_C <- c(BW_NA$SpCond - 71)

#No Cl samples so Cl approx.
BW_NA$Cl_mgL <- ((0.05/0.105)*BW_NA$SpCond_C)

qplot(Cl_mgL, NO3_C, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 250g KNO3 in 10L carboy
NO3mgL <- 250 * (1000) * (62/101) *(1/10)
# Carboy concentrations 2000g NaCl in 10L carboy
NaClmgL <- 2000 * (1000) * (1/10)

NO3mgL/NaClmgL

#mass recovery= 
BW_NA$NtoNaCl <-  BW_NA$NO3_C/BW_NA$Cl_mgL

qplot(datetime, NtoNaCl, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

BW_NA$massR <- (NO3mgL/NaClmgL)- BW_NA$NtoNaCl
BW_NA$massRPer <- (1-((NO3mgL/NaClmgL)- BW_NA$NtoNaCl)/(NO3mgL/NaClmgL)) * 100


library(dplyr)
BW_NA$massRPerSum<- c(0.000000,  
                      2.377923, 
                      6.016597,
                      13.89078,
                      25.98728,
                      35.24685,
                      45.33801,
                      54.95329,
                      63.86841,
                      77.49025,
                      92.76113,
                      109.2268,
                      130.4606)
                      
                      
qplot(datetime, massRPer, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).

injectC1 <- c(log(NO3mgL/NaClmgL), 0)
streamdist <- c(0, 250)

lm(injectC1~streamdist)
qplot(streamdist, injectC1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


injectC2 <- c(log(NO3mgL/NaClmgL), log(0.001824644))
lm(injectC2~streamdist)
qplot(streamdist, injectC2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC3 <- c(log(NO3mgL/NaClmgL), log(0.002792052))
lm(injectC3~streamdist)
qplot(streamdist, injectC3, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


injectC4 <- c(log(NO3mgL/NaClmgL), log(0.006042070))
lm(injectC4~streamdist)
qplot(streamdist, injectC4, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC5 <- c(log(NO3mgL/NaClmgL), log(0.009281968))
lm(injectC5~streamdist)
qplot(streamdist, injectC5, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC6 <- c(log(NO3mgL/NaClmgL), log(0.007105118))
lm(injectC6~streamdist)
qplot(streamdist, injectC6, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC7 <- c(log(NO3mgL/NaClmgL), log(0.007743216))
lm(injectC7~streamdist)
qplot(streamdist, injectC7, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC8 <- c(log(NO3mgL/NaClmgL), log(0.007378061))
lm(injectC8~streamdist)
qplot(streamdist, injectC8, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC9 <- c(log(NO3mgL/NaClmgL), log(0.006840811))
lm(injectC9~streamdist)
qplot(streamdist, injectC9, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC10 <- c(log(NO3mgL/NaClmgL), log(0.010452402))
lm(injectC10~streamdist)
qplot(streamdist, injectC10, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC11 <- c(log(NO3mgL/NaClmgL), log(0.011717756))
lm(injectC11~streamdist)
qplot(streamdist, injectC11, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC12 <- c(log(NO3mgL/NaClmgL), log(0.012634535))
lm(injectC12~streamdist)
qplot(streamdist, injectC12, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC13 <- c(log(NO3mgL/NaClmgL), log(0.016293296))
lm(injectC13~streamdist)
qplot(streamdist, injectC13, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

Kw = c(0.01027, 
       -0.01496, 
       -0.01325, 
       -0.01017, 
       -0.008449, 
       -0.009518,
       -0.009174, 
       -0.009367, 
       -0.00967, 
       -0.007974,
       -0.007517,
       -0.007216,
       -0.006198)

BW_NA$kw<-Kw
BW_NA$sw<- -1/(BW_NA$kw)

qplot(NO3_C, sw, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# back ground corrected N concentrations: Cadd 
Cadd <- 0.05

# Uadd: added areal uptake rate
BW_NA$Uadd <- Q*Cadd/BW_NA$sw*w

# Vf: the uptake velocity (L/T) or uptake efficiency
BW_NA$Vf <- BW_NA$Uadd/Cadd
# injection time 1 hour to 99% recovery

##
## Plotting
##
library(ggplot2)
library(ggformula) # optional
library(MASS)
library(drc)
library(gridExtra)

BW_NA$site <- "Blackwood"
BW_plot <- BW_NA[,c("datetime",
                    "SpCond_C",
                    "Cl_mgL",
                    "NO3_C",
                    "massR",
                    "kw",
                    "sw",
                    "Uadd",
                    "Vf",
                    "site")]


BW_plot1 <- subset(BW_plot, datetime >= as.POSIXct("2021-07-28 16:00:00"))


plot_grid(
  ggplot(BW_plot1, aes(NO3_C, Uadd)) + geom_point(),
  ggplot(BW_plot1, aes(NO3_C, sw)) + geom_point(),
  ggplot(BW_NA, aes(NO3_C, NtoNaCl)) + geom_point(),
  ncol=1, align="hv")



#### M-M curve fit -- Error here.
model.drm1 <- drm (Uadd ~ NO3_C, data = BW_plot1, fct = MM.2())
summary(model.drm1)

mm1<- data.frame(Uadd = seq(0, max(BW_plot$Uadd), length.out = 100))
mm2 <- data.frame(NO3_C = seq(0, max(BW_plot$NO3_C), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm1)

summary(model.drm1)
(model.drm1)


Uadd_plotBW <- ggplot(BW_plot, aes(x=NO3_C, y=Uadd)) + 
  #ylim(0,20) + xlim(0,1500) + 
  geom_line(data = mm2, aes(x = NO3_C, y = Uadd), colour = "black") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate ugL")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))

max(BW_plot$sw)
max(na.omit(BW_plot$Vfaddint))
max(BW_plot$Uadd)






















# OLD CODE 
# 
# concentrations <- glm(NO3_C~Cl_mgL, data=BW_NA)
# summary(concentrations)
# BW_NA$Cgrab = (summary(concentrations)$coefficients)[2,1] * BW_NA$SpCond_C + (summary(concentrations)$coefficients)[1,1]
# 
# ###
# ## TMR ?
# BW_NA$TMR_N = (BW_NA$NO3_C)*Q #*((as.numeric(GB_NA$datetime)*60)- 97619148000)
# BW_NA$TMR_Cl = (BW_NA$Cl_ugL)*Q #*((as.numeric(GB_NA$datetime)*60)- 97619148000)
# 
# slug<- log(((250*1000)/(2000*1000)))#*Q)
# 
# BW_NA$TMR_rat <- log(BW_NA$TMR_N/BW_NA$TMR_Cl) 
# BW_NA$dist <- 250
# 
# ## Remove blank and head of reach .. so first 2 rows
# BW_NA<- BW_NA[-c(1,1),]
# 
# #new df
# distance <- c(0, 
#               BW_NA$dist)
# 
# TMR_con <- c(slug,
#              BW_NA$TMR_rat)
# 
# sample_ID <- c("AAA",
#                BW_NA$Sample_ID)
# 
# newdf <- data_frame(distance, TMR_con, sample_ID)
# 
# kw_glm <- glm(TMR_con~scale(distance)*sample_ID, data=newdf)
# summary(kw_glm)
# 
# (summary(kw_glm)$coefficients)[3:14,]
# BW_NA$kw = c(-1.0160,
#              -1.6892,
#              -0.9172,
#              -0.6877,
#              -0.7551,
#              -0.6691,
#              -0.7109,
#              -3.3743,
#              -0.3691,
#              -0.2548,
#              -0.1504,
#              -0.1404)
# 
# 
# # sw: Added nutrient uptake length calculated with the slug BTC-integrated approach (m) 
# BW_NA$sw <- -1/BW_NA$kw 
# 
# # Uadd-int: Added nutrient areal uptake rate calculated with the slug BTC-integrated approach (ug m^-2 sec^-1)
# BW_NA$Uadd = ((Q * (BW_NA$NO3_C))/ (BW_NA$sw * (w)))
# 
# # Vf-add-int: Added nutrient uptake velocity calculated with the dynamic TASCC approach (m sec^-1)
# BW_NA$Vfaddint = BW_NA$Uadd/ (BW_NA$NO3_C)
# 
# 
# library(ggplot2)
# library(ggformula) # optional
# library(MASS)
# library(drc)
# library(gridExtra)
# 
# 
# BW_NA$site <- "Blackwood"
# BW_plot <- BW_NA[,c("datetime",
#                     "SpCond_C",
#                     "Cl_ugL",
#                     "NO3_C",
#                     "TMR_rat",
#                     "kw",
#                     "sw",
#                     "Uadd",
#                     "Vfaddint",
#                     "site")]
# 
# 
# #### Nutrient uptake rate for GB
# model.drm1 <- drm (Uadd ~ NO3_C, data = BW_plot, fct = MM.2())
# summary(model.drm1)
# 
# mm2 <- data.frame(NO3_C = seq(0, max(BW_plot$NO3_C), length.out = 100))
# mm2$Uadd <- predict(model.drm1, newdata = mm2)
# 
# summary(model.drm1)
# (model.drm1)
# 
# 
# Uadd_plotBW <- ggplot(BW_plot, aes(x=NO3_C, y=Uadd)) + 
#   #ylim(0,20) + xlim(0,1500) + 
#   geom_line(data = mm2, aes(x = NO3_C, y = Uadd), colour = "black") +
#   geom_point(size = 3, shape= 17, col = "#a67d17") +
#   labs(x = expression(paste("Nitrate ugL")),
#        y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
#   theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))
# 
# max(BW_plot$sw)
# max(na.omit(BW_plot$Vfaddint))
# max(BW_plot$Uadd)
