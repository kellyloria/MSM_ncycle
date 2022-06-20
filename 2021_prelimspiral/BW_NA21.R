library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)
library(lme4)
library(lmerTest)

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR/SummerResearch2021')){
  inputDir<- '/Users/kellyloria/Documents/UNR/SummerResearch2021/'
  outputDir<- '/Users/kellyloria/Documents/UNR/SummerResearch2021/DO_downloads/' 
}
## ---------------------------

# Read in the nitrogen uptake assay data:
BW_NA <- read.csv(paste0(inputDir,("/NA21/NA21_BW20210728.csv")))
BW_NA$datetime <- as.POSIXct(as.character(BW_NA$datetime), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# 

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


## SPC data from HOBO
BW_Hobo <-read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/NA21/GeneralandBW_20210727_20775520_4.csv", skip=1)
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


##
## Left join hobo conductivity estimates with samples:
##

BW_NA <- left_join(BW_NA, BW_Hobo[c("DateTime", "SpCond")],
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
BW_NA$NO3_C <- (BW_NA$Results-0.005) *1000

BW_NA$SpCond_C <- c(BW_NA$SpCond - 71)

#No Cl samples so Cl approx.
BW_NA$Cl_ugL <- ((0.05/0.105)*BW_NA$SpCond_C)*1000


qplot(Cl_ugL, NO3_C, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

concentrations <- glm(NO3_C~Cl_mgL, data=BW_NA)
summary(concentrations)
BW_NA$Cgrab = (summary(concentrations)$coefficients)[2,1] * BW_NA$SpCond_C + (summary(concentrations)$coefficients)[1,1]

###
## TMR ?
BW_NA$TMR_N = (BW_NA$NO3_C)*Q #*((as.numeric(GB_NA$datetime)*60)- 97619148000)
BW_NA$TMR_Cl = (BW_NA$Cl_ugL)*Q #*((as.numeric(GB_NA$datetime)*60)- 97619148000)

slug<- log(((250*1000)/(2000*1000)))#*Q)

BW_NA$TMR_rat <- log(BW_NA$TMR_N/BW_NA$TMR_Cl) 
BW_NA$dist <- 250

## Remove blank and head of reach .. so first 2 rows
BW_NA<- BW_NA[-c(1,1),]

#new df
distance <- c(0, 
              BW_NA$dist)

TMR_con <- c(slug,
             BW_NA$TMR_rat)

sample_ID <- c("AAA",
               BW_NA$Sample_ID)

newdf <- data_frame(distance, TMR_con, sample_ID)

kw_glm <- glm(TMR_con~scale(distance)*sample_ID, data=newdf)
summary(kw_glm)

(summary(kw_glm)$coefficients)[3:14,]
BW_NA$kw = c(-1.0160,
             -1.6892,
             -0.9172,
             -0.6877,
             -0.7551,
             -0.6691,
             -0.7109,
             -3.3743,
             -0.3691,
             -0.2548,
             -0.1504,
             -0.1404)


# sw: Added nutrient uptake length calculated with the slug BTC-integrated approach (m) 
BW_NA$sw <- -1/BW_NA$kw 

# Uadd-int: Added nutrient areal uptake rate calculated with the slug BTC-integrated approach (ug m^-2 sec^-1)
BW_NA$Uadd = ((Q * (BW_NA$NO3_C))/ (BW_NA$sw * (w)))

# Vf-add-int: Added nutrient uptake velocity calculated with the dynamic TASCC approach (m sec^-1)
BW_NA$Vfaddint = BW_NA$Uadd/ (BW_NA$NO3_C)


library(ggplot2)
library(ggformula) # optional
library(MASS)
library(drc)
library(gridExtra)


BW_NA$site <- "Blackwood"
BW_plot <- BW_NA[,c("datetime",
                    "SpCond_C",
                    "Cl_ugL",
                    "NO3_C",
                    "TMR_rat",
                    "kw",
                    "sw",
                    "Uadd",
                    "Vfaddint",
                    "site")]


#### Nutrient uptake rate for GB
model.drm1 <- drm (Uadd ~ NO3_C, data = BW_plot, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3_C = seq(0, max(BW_plot$NO3_C), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

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
