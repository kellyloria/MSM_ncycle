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
library(drc)

## ---------------------------

# Read in the nitrogen uptake assay data:
BW_NA <- read.csv("./NA21_dat/NA21_BW20210728.csv")
BW_NA$datetime <- as.POSIXct(as.character(BW_NA$datetime), format="%Y-%m-%d %H:%M:%S") ## modify the format to match your data
# 
summary(BW_NA)

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


## SPC data from HOBO
BW_Hobo <-read.csv("./NA21_dat/GeneralandBW_20210727_20775520_4.csv", skip=1)
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
end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

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
reachL = 250

BW_NA$NO3_C <- (BW_NA$Results-0) 

BW_NA$SpCond_C <- c(BW_NA$SpCond - 71)

#No Cl samples so Cl approx.
BW_NA$Cl_mgL <- ((0.05/0.105)*BW_NA$SpCond_C)

qplot(Cl_mgL, NO3_C, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 250g KNO3 in 10L carboy
NO3mgL <- 250 * (1000) * (62/101) * (1/10)
# Carboy concentrations 2000g NaCl in 10L carboy
NaClmgL <- 2000 * (1000) * (35.45/58.44) *(1/10)
carboy <- NO3mgL/NaClmgL


#mass recovery= 
#BW_NA$MR_No3 <- (NO3mgL-(BW_NA$NO3_C/NO3mgL)) *100
BW_NA$NtoCllog <-  log(BW_NA$NO3_C/BW_NA$Cl_mgL)

qplot(datetime, NtoCllog, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
BW_NA$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NO3=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(BW_NA)) {
  temp_dat <- GB_NA[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((BW_NA$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NO3<- BW_NA$NO3_C[i]
  Cl<- BW_NA$Cl_mgL[i]
  temp_out <- data.frame(Site = strsplit(BW_NA$Sample_ID, split = "_")[[1]][1], 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

Cadd<- 0.05
out <- out[-1,]
out$sw<- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w
out$Uadd <- Q*Cadd/out$sw*w

BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 16:00:00") & DateTime <= as.POSIXct("2021-07-28 17:01:00"))

BW_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(BW_Hobo1, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# ggsave(plot = BW_uptake, filename = paste("./figures/BW220728_v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BW_BTC_20210728.csv", row.names = TRUE)


# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)

#### M-M curve fit -- Error here.
library(dr4pl)
model.drm1 <- drc::drm (Uadd ~ NO3, data = out, fct = MM.3())
summary(model.drm1)

mm1<- data.frame(Uadd = seq(0, max(out$Uadd), length.out = 100))
mm2 <- data.frame(NO3_C = seq(0, max(out$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm1)

summary(model.drm1)
(model.drm1)


Uadd_plotBW <- ggplot(out, aes(x=NO3, y=Uadd)) + 
  #ylim(0,20) + xlim(0,1500) + 
  geom_line(data = mm2, aes(x = NO3_C, y = Uadd), colour = "black") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate ugL")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))

max(out$sw)
max(out$Uadd)



























### old method below
# plan to run in loop: 
lm(injectC1~streamdist)
qplot(streamdist, injectC1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


injectC2 <- c(log(carboy), log(0.001824644))
lm(injectC2~streamdist)
qplot(streamdist, injectC2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC3 <- c(log(carboy), log(0.002792052))
lm(injectC3~streamdist)
qplot(streamdist, injectC3, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


injectC4 <- c(log(carboy), log(0.006042070))
lm(injectC4~streamdist)
qplot(streamdist, injectC4, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC5 <- c(log(carboy), log(0.009281968))
lm(injectC5~streamdist)
qplot(streamdist, injectC5, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC6 <- c(log(carboy), log(0.007105118))
lm(injectC6~streamdist)
qplot(streamdist, injectC6, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC7 <- c(log(carboy), log(0.007743216))
lm(injectC7~streamdist)
qplot(streamdist, injectC7, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC8 <- c(log(carboy), log(0.007378061))
lm(injectC8~streamdist)
qplot(streamdist, injectC8, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC9 <- c(log(carboy), log(0.006840811))
lm(injectC9~streamdist)
qplot(streamdist, injectC9, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC10 <- c(log(carboy), log(0.010452402))
lm(injectC10~streamdist)
qplot(streamdist, injectC10, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC11 <- c(log(carboy), log(0.011717756))
lm(injectC11~streamdist)
qplot(streamdist, injectC11, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC12 <- c(log(carboy), log(0.012634535))
lm(injectC12~streamdist)
qplot(streamdist, injectC12, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

injectC13 <- c(log(carboy), log(0.016293296))
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

BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 16:00:00") & DateTime <= as.POSIXct("2021-07-28 17:01:00"))

BW_uptake<- plot_grid(
  ggplot(BW_plot1, aes(NO3_C, Uadd)) + geom_point(), 
  ggplot(BW_plot1, aes(NO3_C, sw)) + geom_point(),
  ggplot(BW_NA, aes(datetime, NtoNaCl)) + geom_point(),
  ggplot(BW_Hobo1, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

ggsave(plot = BW_uptake, filename = paste("./figures/BW220728.png",sep=""),width=4,height=7,dpi=300)


# N supply:
L= 250
N_supp <-(86400*Q*(Cadd*0.001))/(w*L)

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






###
###
###

# Check the salt volume... 
# stream widths...
# reach length 







####################
## BWL 2022-05-26 ##
BW_Hobo <-read.csv("./NA22_dat/BWL_20220526/20775523_BWL20220526BOR.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
                      "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-05-26 10:15:00") & DateTime <= as.POSIXct("2022-05-26 11:20:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range-- the second salt slug for NO3:
BW_Hobo2 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-05-26 12:07:00") & DateTime <= as.POSIXct("2022-05-26 12:27:10"))

qplot(DateTime, Cond, data = BW_Hobo2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
##
## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo2, DateTime >= as.POSIXct("2022-05-26 12:08:00") & DateTime <= as.POSIXct("2022-05-26 12:13:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*1763 # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo2$DateTime), BW_Hobo2$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2022-05-26 12:16:25") #Lolomai 
peak_time <- BW_Hobo2[which.max(BW_Hobo2$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- reachL/time_diff_sec
v
#############################
## Estimate mean depth (z) ##
#############################
## Enter average width measurement in m
w <- mean(c(6.0, 5.2,10.7, 8, 7, 9, 12.2, 10.7,7.1,6.1,5,5.3))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z





####################
## BWL 2022-08-24 ##
BW_Hobo <-read.csv("./NA22_dat/BWL_20220824/BWLNH4BOR_20775520_20220824.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-08-24 10:05:00") & DateTime <= as.POSIXct("2022-08-24 11:36:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-08-24 10:05:00") & DateTime <= as.POSIXct("2022-08-24 10:17:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*1500) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-08-24 10:18:05") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(4,5.5,4.7,6.8,5.9,8.1,6.3,6.1,5.7,7.7,9,9.3,8.2,7.8,7))
## Calculate effective depth
z <- c(Q/1000)/(w*v)
z

### NO3 pulse
BW_Hobo <-read.csv("./NA22_dat/BWL_20220824/BWL_NO3BOR_20775523_20220824andBWU\ NH4HOR.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%m/%d/%y %H:%M:%S") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo2 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-08-24 11:56:00") & DateTime <= as.POSIXct("2022-08-24 14:00:00"))

qplot(DateTime, Cond, data = BW_Hobo2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-08-24 11:56:30") & DateTime <= as.POSIXct("2022-08-24 12:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000 # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo$DateTime), BW_Hobo$SpCond, bg_SpCond, SpCond_mass)

#inj_time <- as.POSIXct("2021-07-28 15:40:00") #Lolomai 
peak_time <- BW_Hobo[which.max(BW_Hobo$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c() #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- mean()
## Calculate effective depth
z <- (Q/1000)/(w*v)
z



####################
## BWU 2022-08-24 ##
####################

BW_Hobo <-read.csv("./NA22_dat/BWU_20220824/BWL_NO3BOR_20775523_20220824BWUNH4HOR.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))



BW_Hobo <-read.csv("./NA22_dat/BWU_20220824/BWL_NO3BOR_20775523_20220824BWUNH4HOR.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))






# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-08-24 10:05:00") & DateTime <= as.POSIXct("2022-08-24 11:36:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-08-24 10:05:00") & DateTime <= as.POSIXct("2022-08-24 10:17:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*1500) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-08-24 10:18:05") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(4,5.5,4.7,6.8,5.9,8.1,6.3,6.1,5.7,7.7,9,9.3,8.2,7.8,7))
## Calculate effective depth
z <- c(Q/1000)/(w*v)
z





####################
## BWL 2022-10-12 ##
BW_Hobo <-read.csv("./NA22_dat/BWL_20221012/20775520_13.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-12 10:49:00") & DateTime <= as.POSIXct("2022-10-12 13:30:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


BW_Hobo2 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-12 10:49:00") & DateTime <= as.POSIXct("2022-10-12 13:40:00"))

qplot(DateTime, Cond, data = BW_Hobo2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:

## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-12 10:49:00") & DateTime <= as.POSIXct("2022-10-12 11:08:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*600)
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-12 11:02:20") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(175) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(6.2,8.5,7,6.4,10,9.7,8.3,4.8,
            5.4,5,6,7.1,6.2,4.6,7))
## Calculate effective depth
z <- (Q/1000)/(w*v)
z



# Adjust the time range:
BW_Hobo2 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-12 13:31:00") & DateTime <= as.POSIXct("2022-10-12 14:42:00"))

qplot(DateTime, Cond, data = BW_Hobo2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:

## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo2, DateTime >= as.POSIXct("2022-10-12 13:31:00") & DateTime <= as.POSIXct("2022-10-12 13:36:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700)
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo2$DateTime), BW_Hobo2$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-12  13:16:20") #Lolomai 
peak_time <- BW_Hobo2[which.max(BW_Hobo2$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(150) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(6.2,8.5,7,6.4,10,9.7,8.3,4.8,
            5.4,5,6,7.1,6.2,4.6,7))
## Calculate effective depth
z <- (Q/1000)/(w*v)
z






## BWL 2022-10-12 ## NO3 pulse
BW_Hobo <-read.csv("./NA22_dat/BWL_20221012/BWLBOR20775520_14.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%m/%d/%y %H:%M:%S") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-12 03:32:00") & DateTime <= as.POSIXct("2022-10-12 05:17:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-10-12 03:32:00") & DateTime <= as.POSIXct("2022-10-12 03:45:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*600 # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-12 02:40:00") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(75) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- mean()
## Calculate effective depth
z <- (Q/1000)/(w*v)
z



## BWU 2022-10-12 ##
BW_Hobo <-read.csv("./NA22_dat/BWL_20221012/BWLBOR20775520_14.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-12 15:20:00") & DateTime <= as.POSIXct("2022-10-12 17:17:20"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-10-12 15:20:00") & DateTime <= as.POSIXct("2022-10-12 15:31:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*600) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-12 15:17:00") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*3600 
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(75) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(0.7,1.4,1.2,2.2,6,8.2,1.8,2.8, 1.2,1.3,3.5,1.1,1.4,2,3.8))
## Calculate effective depth
z <- c(Q/1000)/(w*v)
z

