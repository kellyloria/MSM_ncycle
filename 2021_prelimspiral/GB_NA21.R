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
GB_Hobo$DateTime <- as.POSIXct(as.character(GB_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
range(GB_Hobo$DateTime)

# 
GB_Hobo$SpCond <- GB_Hobo$Cond/(1-(25-GB_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = GB_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
GB_Hobo <- subset(GB_Hobo, DateTime >= as.POSIXct('2021-07-22 13:30:00') & DateTime <= as.POSIXct('2021-07-22 15:58:00'))

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

GB_NA <- left_join(GB_NA, GB_Hobo[c("DateTime", "SpCond")],
                   by= c("datetime"="DateTime"))

summary(GB_NA)

qplot(datetime, Results, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

##
## Calculations Notes:
##
# 1. Correct for background concentrations (_C):
GB_NA$NO3_C <- (GB_NA$Results-0.021) *1000
GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0


GB_NA[2,7] = 0.05
GB_NA$SpCond_C <- c(GB_NA$SpCond - 413.2548)
GB_NA$SpCond_C <-replace(GB_NA$SpCond_C, GB_NA$SpCond_C<0, 0)

#No Cl samples so Cl approx.
GB_NA$Cl_mgL <- ((0.05/0.105)*GB_NA$SpCond_C)*1000


qplot(Cl_mgL, NO3_CC, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

concentrations <- glm(NO3_CC~Cl_mgL, data=GB_NA)
summary(concentrations)
GB_NA$Cgrab = (summary(concentrations)$coefficients)[2,1] * GB_NA$SpCond_C + (summary(concentrations)$coefficients)[1,1]

###
## TMR ?
GB_NA$TMR_N = (GB_NA$NO3_CC)*Q #*((as.numeric(GB_NA$datetime)*60)- 97619148000)
GB_NA$TMR_Cl = (GB_NA$Cl_mgL)*Q #*((as.numeric(GB_NA$datetime)*60)- 97619148000)

slug<- log(((250*1000)/(2000*1000)))#*Q)

GB_NA$TMR_rat <- log(GB_NA$TMR_N/GB_NA$TMR_Cl) 
GB_NA$dist <- 160

## Remove blank and head of reach .. so first 2 rows

GB_NA<- GB_NA[-c(1,2),]

#new df
distance <- c(0, 
              GB_NA$dist)

TMR_con <- c(slug,
             GB_NA$TMR_rat)

sample_ID <- c("AAA",
               GB_NA$Sample_ID)

newdf <- data_frame(distance, TMR_con, sample_ID)

newdf$TMR_con <-replace(newdf$TMR_con, TMR_con==-Inf, -7)
newdf$TMR_con <-replace_na(newdf$TMR_con, -7)

kw_glm <- glm(TMR_con~scale(distance)*sample_ID, data=newdf)
summary(kw_glm)

(summary(kw_glm)$coefficients)[3:24,]
GB_NA$kw = c(-2.86367561,
             -2.86367561,
             -2.02328686,
             -1.72692921,
             -1.63516153,
             -1.48794594,
             -1.12417899,
             -0.99589352,
             -0.94538373,
             -0.94538373,
             -0.62773154,
             -0.04435973,
             -0.31571081,
             -0.40766076,
             -0.25373352,
             -0.06376129,
             -0.09365,
             -0.17296,
             -0.16296)


# sw: Added nutrient uptake length calculated with the slug BTC-integrated approach (m) 
GB_NA$sw <- -1/GB_NA$kw 

# Uadd-int: Added nutrient areal uptake rate calculated with the slug BTC-integrated approach (ug m^-2 sec^-1)
GB_NA$Uadd = ((Q * (GB_NA$NO3_CC))/ (GB_NA$sw * (w)))

# Vf-add-int: Added nutrient uptake velocity calculated with the dynamic TASCC approach (m sec^-1)
GB_NA$Vfaddint = GB_NA$Uadd/ (GB_NA$NO3_CC)
             

library(ggplot2)
library(ggformula) # optional
library(MASS)
library(drc)
library(gridExtra)


GB_NA$site <- "Glenbrook"
GB_plot <- GB_NA[,c("datetime",
                      "SpCond_C",
                      "Cl_mgL",
                      "NO3_CC",
                      "TMR",
                      "kw",
                      "sw",
                      "Uadd",
                      "Vfaddint",
                      "site")]


#### Nutrient uptake rate for GB
model.drm1 <- drm (Uadd ~ NO3_CC, data = GB_plot, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3_CC = seq(0, max(GB_plot$NO3_CC), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)


Uadd_plotGB <- ggplot(GB_plot, aes(x=NO3_CC, y=Uadd)) + 
  #ylim(0,20) + xlim(0,1500) + 
  geom_line(data = mm2, aes(x = NO3_CC, y = Uadd), colour = "black") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate ugL")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))

max(GB_plot$sw)
max(na.omit(GB_plot$Vfaddint))
max(GB_plot$Uadd)
