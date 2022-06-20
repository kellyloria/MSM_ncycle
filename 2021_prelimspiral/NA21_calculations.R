library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(rstan)
library(unitted)
library(zoo)
library(lubridate)
library(dataRetrieval)

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR/SummerResearch2021')){
  inputDir<- '/Users/kellyloria/Documents/UNR/SummerResearch2021/'
  outputDir<- '/Users/kellyloria/Documents/UNR/SummerResearch2021/DO_downloads/' 
}
## ---------------------------



#http://www.mcglynnlab.com/uploads/1/0/6/4/10645747/covino_et_al_land_use_-_land_cover_and_scale_influences_on_in-stream_nitrogen_uptake_kinetics_jgr-b_2012.pdf

# Read in the nitrogen uptake assay data:
GB_NA <- read.csv(paste0(inputDir,("/NA21/NA21_GB20210722.csv")))
GB_NA$datetime <- as.POSIXct(as.character(GB_NA$datetime), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# 

qplot(datetime, Results, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


## SPC data from HOBO
GB_Hobo <-read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/GlenbrookNA20210722_20775520_4.csv", skip=1)
summary(GB_Hobo)


# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
GB_Hobo <- GB_Hobo[,c("Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
              "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(GB_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
GB_Hobo$DateTime <- as.POSIXct(as.character(GB_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# 
range(GB_Hobo$DateTime)

GB_Hobo$SpCond <- GB_Hobo$Cond/(1-(25-GB_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = GB_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
GB_Hobo <- subset(GB_Hobo, DateTime >= as.POSIXct('2021-07-22 13:30:00') & DateTime <= as.POSIXct('2021-07-22 15:58:00'))
#

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
sub_bg <- subset(GB_Hobo, DateTime >= as.POSIXct("2021-07-22 13:30:00") & DateTime <= as.POSIXct("2021-07-22 13:35:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(GB_Hobo$DateTime), GB_Hobo$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-22 13:40:00") #Lolomai 
peak_time <- GB_Hobo[which.max(GB_Hobo$SpCond),]$DateTime 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds
v <- 160/time_diff_sec
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
w <- mean(1.1, 2.0, 1.12, 1.2, 0.72, 0.43, 0.21, 0.61, 0.52, 0.55)
## Calculate effective depth in m
z <- (Q/1000)/(w*v)
z



##
## Left join hobo conductivity estimates with samples:
##

GB_NA <- left_join(GB_Hobo, GB_NA[c("datetime", "Results", "YSI_spc")],
                   by= c("DateTime"="datetime"))

summary(GB_NA)

qplot(DateTime, Results, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

##
## Calculations Notes:
##

# Nutrient uptake length (Sw) can be estimated by measuring the longitudinal 
#   decline of biologically active tracer (nutrient) relative to conservative 
#   tracer with distance downstream during the plateau portion of a constant-rate 
#   tracer addition experiment (Stream Solute Workshop 1990).

# longitudinal uptake rate of added nutrient (kw-add-plat) (in L^-1)
# 

# The plateau approach uptake length of added nutrient (Sw-add-plat) (in L):
#       Sw-add-plat = –1/kw-add-plat

# plateau approach added nutrient areal uptake rates (Uadd-plat):
#      Uadd-plat = Q * [NO3-Nadd-plat]/Sw-add-plat * w

# Uptake velocities (Vf-add-plat):
#      Vf-add-plat = Uadd-plat/[NO3-Nadd-plat]

# Uadd-plat is the plateau approach areal uptake rate of added nutrient (M * L–2 * T–1)
# Amount of dissloved N in amount of water and time 

# Q is stream discharge (L^3 * T^–1)

# [NO3-Nadd-plat] is the geometric mean of background corrected NO3-N concentrations of longitudinal grab samples collected across the stream reach during constant-rate plateau conditions (M*L^–3)



# TMR =Q∫TC(t)dt

# TMR is the tracer mass recovery (M)
# TC is the time- integrated tracer concentrations (M*T*L^–3) of background corrected Cl and NO3-N
# Sw-add-int was calculated by plotting the natural log of the injectate NO3-N:Cl ratio and the BTC-integrated NO3-N:Cl [i.e., TMR(NO3-N):TMR(Cl)] ratio against stream distance
# The slope of the line derived from these data is the BTC-integrated longitudinal uptake rate of added nutrient (kw- add-int), and Sw-add-int is the negative inverse of kw-add-int.


# 1. Correct for background concentrations (_C):
GB_NA$NO3_C <- GB_NA$Results-0.02
GB_NA$NO3_C <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0.001) # Na's produced in TMR calculations if 0

GB_NA$SpCond_C <- GB_NA$SpCond- bg_SpCond
GB_NA$SpCond_C <-replace(GB_NA$SpCond_C, GB_NA$SpCond_C<0, 0.001)

#No Cl samples so Cl approx.
GB_NA$Cl_mgL <- (0.05/0.105)*GB_NA$SpCond_C

# Possible over correction for molarity:
# GB_NA$Cl_mgL <- (0.05/0.105)*GB_NA$SpCond_C
# GB_NA$Cl_M <- GB_NA$Cl_mgL * 35.453 
# GB_NA$NO3_M <- GB_NA$NO3_C * 101.1032 

GB_NA$NO3_Clratio <- GB_NA$NO3_C/GB_NA$SpCond_C
GB_NA$NO3_Clratio <-replace(GB_NA$NO3_Clratio, GB_NA$NO3_Clratio<0, 0.001)



# mg per L injectate:
NO3_ClratioInjTMR <- ((250*101.1032)/(2000*35.453))/ (18.95)
                      
plot(GB_NA$DateTime, GB_NA$NO3_C)
plot(GB_NA$DateTime, GB_NA$Cl_mgL)



###
# Restrict time to remove above head of reach samples from kw:
GB_NAom<- na.omit(GB_NA)
GB_NAom <- subset(GB_NAom, DateTime >= as.POSIXct("2021-07-22 14:06:00"))

plot(GB_NAom$DateTime, GB_NAom$NO3_Clratio)

# TMR:
GB_NAom$TMR = (Q1*(GB_NAom$NO3_Clratio))
TMR1<- NO3_ClratioInj*Q1
#GB_NA$TMR = (v)*((GB_NA$NO3_Clratio))


# kw and sw:
# new df
distance <- c(0, 
              160
              )

# Different TMR methods
# concentrations <- c(1.77095, 
#                     #0.17473858 
#                     #0.0001 
#                     #-0.01787905  
#                     #0.02777912 
#                     #0.03785915  
#                     #0.04208335  
#                     #0.04947852  
#                     #0.07153914
#                     #0.08165646  
#                     #0.08647761  
#                     #0.02602138  
#                     #0.11906101  
#                     #0.21277225  
#                     #0.16200644 
#                     #0.14758292  
#                     #0.17142364
#                     #0.20604545  
#                     #0.19862176  
#                     #0.18160820  
#                     0.21183411
#                    )
# 
# concentrations2 <- c(0.035172588, 
#                      0.0004261029
#                      #0.0001  
#                      #0.005591573  
#                      #0.007620551  
#                      #0.008470828  
#                      #0.009959378  
#                      #0.014399893 
#                      #0.016436376  
#                      #0.017406810  
#                      #0.005237763  
#                      #0.023965422  
#                      #0.042828268  
#                      #0.032609775  
#                      #0.029706510  
#                      #0.034505335
#                      #0.041474252  
#                      #0.039979961  
#                      #0.036555354  
#                      #0.042639434 
# )
# 
# 
# TMR3 <- c(0.0005703506, 
#                      #0.0004261029
#                      #1.600000e-05  
#                      #1.150986e-05  
#                      #3.838739e-06 
#           #8.946517e-06 
#           #1.219288e-05 
#           #1.355333e-05 
#           #1.593500e-05 
#           #2.303983e-05
#           #2.629820e-05 
#           #2.785090e-05 
#           #8.380421e-06 
#           #3.834468e-05 
#           #6.852523e-05 
#           #5.217564e-05 
#           #4.753042e-05 
#           #5.520854e-05
#           #6.635880e-05 
#           #6.396794e-05 
#           #5.848857e-05 
#           6.822309e-05
# )


TMR4 <- c(3.013415e-05, 
          #4.036065e-07 
          #1.346098e-07 
          #3.137199e-06 
          #4.275574e-06 
          #4.752629e-06 
          #5.587792e-06 
          #8.079181e-06 
          #9.221767e-06
          #9.766237e-06 
          2.938691e-06 
          #1.344600e-05 
          #2.402916e-05 
          #1.829599e-05 
          #1.666709e-05 
          #1.935951e-05 
          #2.326948e-05
          #2.243109e-05 
          #2.050969e-05 
          #2.392321e-05

          )

newdf <- data_frame(distance, TMR4)

# bp <- ggplot(newdf, aes(x=distance, y=log(TMR3))) + 
#   geom_point(aes(color=ID)) 

Kwslope.mod <- lm(log(TMR4)~scale(distance), data=newdf)
summary(Kwslope.mod)


# kw: Added nutrient longitudinal uptake rate calculated with the slug BTC-integrated approach (m^-1)
# kw <- c(
#   -2.9369,
#   -2.9371,
#   -2.9384,
#        -2.938,
#        -2.719,
#        -2.644,
#        -2.530,
#        -2.269,
#        -2.1756,
#        -2.1350,
#        -2.984,
#        -1.9089,
#        -1.498,
#        -1.6911,
#        -1.7571,
#        -1.6512,
#        -1.5211,
#        -1.5471,
#        -1.6104,
#        -1.5015
#        )
# 
# kw2 <- c(
#         -2.527,
#         -2.760,
#         -3.536,
#         -2.938,
#         -2.719,
#         -2.644,
#         -2.530,
#         -2.269,
#         -2.176,
#         -2.135,
#         -2.984,
#         -1.909,
#         -1.498,
#         -1.691,
#         -1.757,
#         -1.651,
#         -1.521,
#         -1.547,
#         -1.610,
#         -1.502)

kw4 <- c( # TMR calculated with mg/L
  -3.05,
  -3.826,
  -1.60,
  -1.381,
  -1.381,
  -1.192,
  -0.9308,
  -0.8373,
  -0.7967,
  -1.646,
  -0.5706,
  -0.1601,
  -0.3528,
  -0.4188,
  -0.3129,
  -0.1828,
  -0.2087,
  -0.2721,
  -0.1632 
         )

GB_NAom$kw <- kw4

# sw: Added nutrient uptake length calculated with the slug BTC-integrated approach (m)       
sw <- -1/(kw4) # mg to ug correction?
GB_NAom$sw <- sw
plot(GB_NAom$NO3_C, GB_NAom$sw)

# Uadd-int: Added nutrient areal uptake rate calculated with the slug BTC-integrated approach (ug m^-2 sec^-1)
GB_NAom$Uadd = ((Q * (GB_NAom$NO3_C*1000))/ (GB_NAom$sw * (w)))

# Vf-add-int: Added nutrient uptake velocity calculated with the dynamic TASCC approach (m sec^-1)
GB_NAom$Vfaddint = GB_NAom$Uadd/ (GB_NAom$NO3_C*1000)
plot(GB_NAom$NO3_C, GB_NAom$Uadd)
plot(GB_NAom$NO3_C, GB_NAom$Vfaddint)

##---------------------------
##---------------------------
## 
## Blackwood creek 
##
##---------------------------
##---------------------------

#http://www.mcglynnlab.com/uploads/1/0/6/4/10645747/covino_et_al_land_use_-_land_cover_and_scale_influences_on_in-stream_nitrogen_uptake_kinetics_jgr-b_2012.pdf

# Read in the nitrogen uptake assay data:
BW_NA <- read.csv(paste0(inputDir,("/NA21/NA21_BW20210728.csv")))
BW_NA$datetime <- as.POSIXct(as.character(BW_NA$datetime), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# 
summary(BW_NA)

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))




## SPC data from HOBO
BW_Hobo1 <-read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/NA21/GeneralandBW_20210727_20775520_4.csv", skip=1)
summary(BW_Hobo1)


# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo1 <- BW_Hobo1[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo1) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo1$DateTime <- as.POSIXct(as.character(BW_Hobo1$DateTime), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# 
range(BW_Hobo1$DateTime)

summary(BW_Hobo1)


qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

BW_Hobo <- subset(BW_Hobo1, DateTime >= as.POSIXct("2021-07-28 15:30:00") & DateTime <= as.POSIXct("2021-07-28 17:01:50"))

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# BW_Hobo$DateTime <- as.POSIXct(round_date(
#   as.POSIXct(BW_Hobo$DateTime, format="%Y-%m-%d %H:%M:%OS"), hour, unit="5 seconds"))



#BW_NA$SpCond <- BW_NA$YSI_spc/(1-(25-21.7)*0.021/100)

qplot(DateTime, SpCond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:

## Equation:
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 15:30:00") & DateTime <= as.POSIXct("2021-07-28 15:40:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)


## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo$DateTime), BW_Hobo$SpCond, bg_SpCond, SpCond_mass)
## Summary


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-28 15:40:00") #Lolomai 
peak_time <- BW_Hobo[which.max(BW_Hobo$SpCond),]$DateTime 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds 
# 250 m upstream
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

##
##
##
## YOU ARE HERE !!!!
#
#
#
#
#
#


## left join cond data to NA data
BW_NA1 <- left_join(BW_NA, BW_Hobo[c("DateTime", "TempC", "SpCond")],
                   by= c("datetime"="DateTime"))

summary(BW_NA1)

# 1. Correct for background concentrations (_C):
BW_NA1$NO3_C <- BW_NA1$Results-0.005
BW_NA1$SpCond_C <- BW_NA1$SpCond- 71
BW_NA1$Cl_mgL <- (0.05/0.105)*BW_NA1$SpCond_C
#BW_NA$Cl_M <- BW_NA$Cl_mgL * 35.453 
#BW_NA$NO3_M <- BW_NA$NO3_C * 101.1032 

BW_NA1$NO3_Clratio <- BW_NA1$NO3_C/BW_NA1$Cl_mgL
NO3_ClratioInj <- (((250*101.1032)/(2000*35.453)))/ (3.78541*18.05)

plot(BW_NA1$datetime, BW_NA1$NO3_C)
plot(BW_NA1$datetime, BW_NA1$Cl_mgL)

plot(BW_NA1$datetime, BW_NA1$NO3_Clratio)

###
# Restrict time to remove above reach samples:
#BW_NA <- subset(BW_NA, datetime >= as.POSIXct("2021-07-28 15:32:00"))

# TMR:
BW_NA1$TMR = Q*((BW_NA1$NO3_Clratio))
TMR1 <- NO3_ClratioInj*Q

# kw and sw:
# new df
distance <- c(0, 
              500)

concentrations <- c(0.2693697,
                    #0.00000000 
                    #0.28262868 
                    #0.14415854 
                    #0.31196269 
                    #0.39245627 
                    #0.36684974 
                    #0.39979589 
                    #0.38345048 
                    #0.02672888
                    #0.53967592 
                    #0.60500837 
                    #0.67159835 
                    #0.78061719
)


newdf <- data_frame(distance, concentrations)

# bp <- ggplot(newdf, aes(x=distance, y=log(concentrations), group=ID)) + 
#   geom_point(aes(color=ID)) 

Kwslope.mod <- lm(log(concentrations)~scale(distance), data=newdf)
summary(Kwslope.mod)


# kw: Added nutrient longitudinal uptake rate calculated with the slug BTC-integrated approach (m^-1)


# kw1<- c(
#   -2.395,
#   -3.187,
#   -1.661,
#   -1.353,
#   -1.273,
#   -1.244,
#   -1.252,
#   -3.129,
#   -1.023,
#   -0.9728,
#   -0.9594,
#   -1.06)
  

kw2 <- c(
  0.3719,
  0.03398,
  -0.4421,
  0.1038,
  0.2661,
  
  
  
  
  
  
  
  -1.755,
  -0.229,
  0.07909,
  0.1591,
  0.1882,
  0.1804,
  -1.696,
  0.4094,
  0.4596,
  0.473,
  0.3719
  )


# sw: Added nutrient uptake length calculated with the slug BTC-integrated approach (m)       
sw <- -1/(kw2)
BW_NA$sw <- sw
plot(BW_NA$NO3_C, BW_NA$sw)

lm(sw~NO3_M, data=BW_NA)

# Uadd-int: Added nutrient areal uptake rate calculated with the slug BTC-integrated approach (mg m^-2 sec^-1)
BW_NA$Uadd = (Q * (BW_NA$NO3_C*1000))/ (BW_NA$sw * (w))

# Vf-add-int: Added nutrient uptake velocity calculated with the dynamic TASCC approach (m sec^-1)
BW_NA$Vfaddint = BW_NA$Uadd/ BW_NA$NO3_C


plot(BW_NA$NO3_C, BW_NA$Uadd)
plot(BW_NA$NO3_C, BW_NA$Vfaddint)


##
## Visualizations:
##


library(ggplot2)
library(ggformula) # optional
library(MASS)
library(drc)
library(gridExtra)


# goal of making uptake plots of uptake lengeth SW, upake rate U add, and uptake velocity V add

# New df:
GB_NAom$site <- "Glenbrook"
GB_plot <- GB_NAom[,c("DateTime",
                      "SpCond_C",
                      "NO3_C",
                      "NO3_Clratio",
                      "TMR",
                      "kw",
                      "sw",
                      "Uadd",
                      "Vfaddint",
                      "site")]



BW_NA$site <- "Blackwood"
BW_plot <- BW_NA[,c("datetime",
                      "SpCond_C",
                      "NO3_C",
                      "NO3_Clratio",
                      "TMR",
                      "kw",
                      "sw",
                      "Uadd",
                      "Vfaddint",
                      "site")]
                      
colnames(BW_plot) <- c("DateTime",
                       "Cl_M",
                       "NO3_C",
                       "NO3_Clratio",
                       "TMR",
                       "kw",
                       "sw",
                       "Uadd",
                       "Vfaddint",
                       "site")



# # Trying dif line fits
# vCurve <- nls(Vfaddint ~ a*NO3_C / ( b+NO3_C), 
#               data = BW_plot, list(a = 1, b = 1))
# 
# summary(vCurve)
# 
# vlinear.model <-lm(Vfaddint ~ NO3_C, BW_plot)
# vlog.model <-lm(log(Vfaddint) ~ scale(NO3_C), BW_plot)
# vexp.model <-lm(Vfaddint ~ exp(NO3_C), BW_plot)
# 
# vlog.model.df <- data.frame(x = BW_plot$NO3_C,
#                             y = exp(fitted(vlog.model)))
# 
# vadd_plot <- ggplot(BW_plot, aes(x=NO3_C, y=Vfaddint)) +
#   geom_line(data = vlog.model.df, aes(x, y, color = "Log Model"), size = 1, linetype = 2, color="black")  +
#   geom_point(size = 3, col = "#60a6bf") +
#   labs(x = expression(paste("Nitrate concentration (mol)")),
#        y= expression(paste("Vadd (mm min^-1)"))) +
#   theme_classic() +
#   annotate("text", x = c(19,17,17), y = c(8,7.65,7.25), label = c("Vmax = 8.221 mm min^-1", "r^2 = 0.394", "p = 0.017"))


##
## Nutrient uptake rate

uCurve2 <- nls(Uadd ~ a*NO3_C / ( b+NO3_C), 
               data = BW_plot, list(a = 1, b = 1))

summary(uCurve2)

ulinear.model <-lm(Uadd ~ NO3_M, BW_plot)
ulog.model <-lm(log(Uadd) ~ scale(NO3_M), BW_plot)
uexp.model <-lm(Uadd ~ exp(NO3_M), BW_plot)

ulog.model.df <- data.frame(x = BW_plot$NO3_M,
                           y = exp(fitted(ulog.model)))


## mm curve:
model.drm <- drm (Uadd ~ NO3_C, data = BW_plot, fct = MM.2())
summary(model.drm)

mml <- data.frame(NO3_C= seq(0, max(BW_plot$NO3_C), length.out = 100))
mml$Uadd <- predict(model.drm, newdata = mml)

# 
# ggplot(BW_plot, aes(x = NO3_C, y = Uadd)) +
#   theme_bw() +
#   ggtitle("Michaelis-Menten kinetics") +
#   geom_point(alpha = 0.5) +
#   geom_line(data = mml, aes(x = NO3_C, y = Uadd), colour = "red")
# 


Uadd_plotBW <- ggplot(BW_plot, aes(x=NO3_C, y=Uadd)) +
  geom_line(data = mml, aes(x = NO3_C, y = Uadd), linetype = 2, colour = "black")+
  #ylim(0,5) + xlim(0,250) + 
  geom_point(size = 3, col = "#60a6bf") +
  labs(x = expression(paste("Nitrate ugL ")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() + 
  annotate("text", x = c(65,55,55), y = c(5, 4.5, 4), label = c("Umax = 6081.72", "Km = 50.16", "p > 0.05 "))

max(BW_plot$sw)
max(BW_plot$Vfaddint)
max(BW_plot$Uadd)



#### Nutrient uptake rate for GB
uCurve1 <- nls(Uadd ~ a*NO3_M / ( b+NO3_M), 
              data = GB_plot, list(a = 1, b = 1))

summary(uCurve1)

residuals(uCurve1)
plot(residuals(uCurve1) ~ predict(uCurve1))
shapiro.test(residuals(uCurve1))
confint(uCurve1)


ulinear.model <-lm(Uadd ~ NO3_C, GB_plot)
ulog.model <-lm(log(Uadd+1) ~ (NO3_C), GB_plot)
uexp.model <-lm(Uadd ~ exp(NO3_C), GB_plot)

ulog.model.df <- data.frame(x = GB_plot$NO3_C,
                            y = exp(fitted(ulog.model)))


model.drm1 <- drm (Uadd ~ NO3_C, data = GB_plot, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3_C = seq(0, max(GB_plot$NO3_C), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)



Uadd_plotGB <- ggplot(GB_plot, aes(x=NO3_C, y=Uadd)) + 
  #ylim(0,20) + xlim(0,1500) + 
  geom_line(data = mm2, aes(x = NO3_C, y = Uadd), colour = "black") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate ugL")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() + 
  annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))

max(GB_plot$sw)
max(GB_plot$Vfaddint)
max(GB_plot$Uadd)


# 
# ModelDat <- grid.arrange(Uadd_plotBW, Uadd_plotGB, ncol=2)
#                          #ncol = 2, heights=c(1,1))
# 
# #ggsave(paste0(outputDir,"/MSM_Nuptake.jpeg"), ModelDat, scale = 0.8, width =20, height = 10, units = c("cm"), dpi = 500)
# 

