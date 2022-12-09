library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)
library(lme4)
library(lmerTest)
library(cowplot)

## ---------------------------
## GBU 2021-07-22
# Read in the nitrogen uptake assay data:
GB_NA <- read.csv("./NA21_dat/NA21_GB20210722.csv")
GB_NA$datetime <- as.POSIXct(as.character(GB_NA$datetime), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# 

qplot(datetime, Results, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


## SPC data from HOBO
GB_Hobo <-read.csv("./NA21_dat/GlenbrookNA20210722_20775520_4.csv", skip=1)
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

qplot(DateTime, Cond, data = GB_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

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
end_time <-as.POSIXct("2021-07-22 15:58:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- 160
v <- reachL/time_diff_sec
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

GB_NA <- subset(GB_NA, datetime >= as.POSIXct('2021-07-22 13:40:00'))


qplot(datetime, Results, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))



##
## Calculations Notes:
##
# 1. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
GB_NA$NO3_C <- (GB_NA$Results) 

#GB_NA[2,7] = 0.05
GB_NA$SpCond_C <- c(GB_NA$SpCond - 406)
#GB_NA$SpCond_C <-replace(GB_NA$SpCond_C, GB_NA$SpCond_C<0, 0)

#No Cl samples so Cl approx.
GB_NA$Cl_mgL <- ((0.05/0.105)*GB_NA$SpCond_C)

qplot(Cl_mgL, NO3_C, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 250g KNO3 in 10L carboy
NO3mgL <- 250 * (1000) * (62/101) *(1/10)
# Carboy concentrations 2000g NaCl in 10L carboy
NaClmgL <- 2000 * (1000) * (35.45/58.44) * (1/10)
carboy <- NO3mgL/NaClmgL

#mass recovery= 
GB_NA$NtoNaCl <-  GB_NA$NO3_C/GB_NA$Cl_mgL
GB_NA$NtoNaCllog <-  log(GB_NA$NO3_C/GB_NA$Cl_mgL)

qplot(datetime, NtoNaCl, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

GB_NA$massR <- (carboy)- GB_NA$NtoNaCl
GB_NA$massRPer <- (1-((carboy)- GB_NA$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
GB_NA$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NO3=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(GB_NA)) {
  temp_dat <- GB_NA[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((GB_NA$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NO3<- GB_NA$NO3_C[i]
  Cl<- GB_NA$Cl_mgL[i]
  temp_out <- data.frame(Site = strsplit(GB_NA$Sample_ID, split = "_")[[1]][1], 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

Cadd <- 0.021

out <- out[-1,]
out$sw<- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(GB_Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# ggsave(plot = BW_uptake, filename = paste("./figures/BW220728_v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GB_BTC_20210722.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)

#### M-M curve fit -- Error here.
library(dr4pl)
model.drm1 <- drc::drm (Uadd ~ NO3, data = out, fct = MM.3())
summary(model.drm1)

mm1<- data.frame(Uadd = seq(0, max(out$Uadd), length.out = 100))
mm1 <- data.frame(NO3_C = seq(0, max(out$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm1)



Uadd_plotBW <- ggplot(out, aes(x=NO3, y=Uadd)) + 
  #ylim(0,20) + xlim(0,1500) + 
  geom_line(data = mm2[-1,], aes(x = NO3_C, y = Uadd), colour = "black") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate ugL")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))

max(out$sw)
max(out$Uadd)


## ---------------------------
## GBU 2021-06-03

# Glenbrook 2021-06-23 lower
# Read in the nitrogen uptake assay data:
GB_NA <- read.csv("./NA21_dat/samples/GBL20210603sample.csv")
GB_NA$datetime <- as.POSIXct(paste(GB_NA$date, GB_NA$time), format = "%m/%d/%y %H:%M:%S")


qplot(datetime, NO3, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

# 
# ## SPC data from HOBO
GB_Hobo <-read.csv("./NA21_dat/GlenbrookNutAssay1_20775509.csv")
summary(GB_Hobo)
# 
# 
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat)
GB_Hobo <- GB_Hobo[,c("Date.Time",
                      "Low.Range",
                      "Temp")]
# 
colnames(GB_Hobo) <- c("DateTime","Cond","TempC")
# # Convert DateTime
GB_Hobo$DateTime <- as.POSIXct(as.character(GB_Hobo$DateTime), format="%m/%d/%y %H:%M")
range(GB_Hobo$DateTime)
GB_Hobo<-GB_Hobo[!duplicated(GB_Hobo$DateTime), ]
# 
GB_Hobo$SpCond <- GB_Hobo$Cond/(1-(25-GB_Hobo$TempC)*0.021/100)

qplot(DateTime, SpCond, data = GB_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:


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
sub_bg <- subset(GB_Hobo, DateTime <= as.POSIXct("2021-06-03 12:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*250
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(GB_Hobo$DateTime), GB_Hobo$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-06-03 12:10:00") #Lolomai 
peak_time <- GB_Hobo[which.max(GB_Hobo$SpCond),]$DateTime 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- 50
v <- reachL/time_diff_sec
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
w <- mean(c(1.1, 2.0, 1.12, 1.2, 1.72, 2.43, 1.21, 1.61, 0.52, 1.55))
## Calculate effective depth in m
#z <- (Q/1000)/(w*v)
z <- 0.3125



##
## Left join hobo conductivity estimates with samples:
##

GB_NA <- left_join(GB_NA, GB_Hobo[c("DateTime", "SpCond")],
                   by= c("datetime"="DateTime"))

summary(GB_NA)

#GB_NA <- subset(GB_NA, datetime >= as.POSIXct('2021-07-22 13:40:00'))


qplot(datetime, NO3, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))



##
## Calculations Notes:
##
# 1. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
GB_NA$NO3_C <- c(GB_NA$NO3) 

#GB_NA[2,7] = 0.05
GB_NA$SpCond_C <- c(GB_NA$SpCond - 326)
#GB_NA$SpCond_C <-replace(GB_NA$SpCond_C, GB_NA$SpCond_C<0, 0)

#No Cl samples so Cl approx.
GB_NA$Cl_mgL <- ((0.05/0.105)*GB_NA$SpCond_C)

qplot(Cl_mgL, NO3_C, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 250g KNO3 in 10L carboy
NO3mgL <- 50 * (1000) * (62/101) *(1/10)
# Carboy concentrations 2000g NaCl in 10L carboy
NaClmgL <- 250 * (1000) * (35.45/58.44) * (1/10)
carboy <- NO3mgL/NaClmgL

#mass recovery= 
GB_NA$NtoNaCl <-  GB_NA$NO3_C/GB_NA$Cl_mgL
GB_NA$NtoNaCllog <-  log(GB_NA$NO3_C/GB_NA$Cl_mgL)

qplot(datetime, NtoNaCl, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

GB_NA$massR <- (carboy)- GB_NA$NtoNaCl
GB_NA$massRPer <- (1-((carboy)- GB_NA$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
GB_NA$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NO3=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(GB_NA)) {
  temp_dat <- GB_NA[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((GB_NA$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NO3<- GB_NA$NO3_C[i]
  Cl<- GB_NA$Cl_mgL[i]
  temp_out <- data.frame(Site = "GB", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

Cadd <- 0.003

out <- out[-1,]
out$sw<- -1/(out$kw)
out$Uadd <- (Q*0.25)*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(GB_Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# write.csv(x = out, file = "./BTC_out/GBL_BTC_20210603.csv", row.names = TRUE)

N_supp <-(86400*(Q*0.25)*(Cadd*0.001))/(w*reachL)



## ---------------------------
## GBU 2021-06-03

# Glenbrook 2021-06-23 lower
# Read in the nitrogen uptake assay data:
GB_NA <- read.csv("./NA21_dat/samples/GBL20210603sample.csv")
GB_NA$datetime <- as.POSIXct(paste(GB_NA$date, GB_NA$time), format = "%m/%d/%y %H:%M:%S")


qplot(datetime, NO3, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

# 
# ## SPC data from HOBO
GB_Hobo <-read.csv("./NA21_dat/GlenbrookNutAssay1_20775509.csv")
summary(GB_Hobo)
# 
# 
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat)
GB_Hobo <- GB_Hobo[,c("Date.Time",
                      "Low.Range",
                      "Temp")]
# 
colnames(GB_Hobo) <- c("DateTime","Cond","TempC")
# # Convert DateTime
GB_Hobo$DateTime <- as.POSIXct(as.character(GB_Hobo$DateTime), format="%m/%d/%y %H:%M")
range(GB_Hobo$DateTime)
GB_Hobo<-GB_Hobo[!duplicated(GB_Hobo$DateTime), ]
# 
GB_Hobo$SpCond <- GB_Hobo$Cond/(1-(25-GB_Hobo$TempC)*0.021/100)

qplot(DateTime, SpCond, data = GB_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:


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
sub_bg <- subset(GB_Hobo, DateTime <= as.POSIXct("2021-06-03 12:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*250
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(GB_Hobo$DateTime), GB_Hobo$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-06-03 12:10:00") #Lolomai 
peak_time <- GB_Hobo[which.max(GB_Hobo$SpCond),]$DateTime 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- 50
v <- reachL/time_diff_sec
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
w <- mean(c(1.1, 2.0, 1.12, 1.2, 1.72, 2.43, 1.21, 1.61, 0.52, 1.55))
## Calculate effective depth in m
#z <- (Q/1000)/(w*v)
z <- 0.3125



##
## Left join hobo conductivity estimates with samples:
##

GB_NA <- left_join(GB_NA, GB_Hobo[c("DateTime", "SpCond")],
                   by= c("datetime"="DateTime"))

summary(GB_NA)

#GB_NA <- subset(GB_NA, datetime >= as.POSIXct('2021-07-22 13:40:00'))


qplot(datetime, NO3, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))



##
## Calculations Notes:
##
# 1. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
GB_NA$NO3_C <- c(GB_NA$NO3) 

#GB_NA[2,7] = 0.05
GB_NA$SpCond_C <- c(GB_NA$SpCond - 326)
#GB_NA$SpCond_C <-replace(GB_NA$SpCond_C, GB_NA$SpCond_C<0, 0)

#No Cl samples so Cl approx.
GB_NA$Cl_mgL <- ((0.05/0.105)*GB_NA$SpCond_C)

qplot(Cl_mgL, NO3_C, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 250g KNO3 in 10L carboy
NO3mgL <- 50 * (1000) * (62/101) *(1/10)
# Carboy concentrations 2000g NaCl in 10L carboy
NaClmgL <- 250 * (1000) * (35.45/58.44) * (1/10)
carboy <- NO3mgL/NaClmgL

#mass recovery= 
GB_NA$NtoNaCl <-  GB_NA$NO3_C/GB_NA$Cl_mgL
GB_NA$NtoNaCllog <-  log(GB_NA$NO3_C/GB_NA$Cl_mgL)

qplot(datetime, NtoNaCl, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

GB_NA$massR <- (carboy)- GB_NA$NtoNaCl
GB_NA$massRPer <- (1-((carboy)- GB_NA$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
GB_NA$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NO3=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(GB_NA)) {
  temp_dat <- GB_NA[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((GB_NA$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NO3<- GB_NA$NO3_C[i]
  Cl<- GB_NA$Cl_mgL[i]
  temp_out <- data.frame(Site = "GBU", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

Cadd <- 0.003

out <- out[-1,]
out$sw<- -1/(out$kw)
out$Uadd <- (Q*0.25)*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(GB_Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# write.csv(x = out, file = "./BTC_out/GBL_BTC_20210603.csv", row.names = TRUE)

N_supp <-(86400*(Q*0.25)*(Cadd*0.001))/(w*reachL)






###########
################################################################################################
###########






####################
## GBL 2022-04-07 ##
####################
BW_Hobo <-read.csv("./NA22_dat/GBL_20220407/Copy\ of\ 20775523_GBLBOR20220407_NO3.csv", skip=1)
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
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-04-07 10:00:00") & DateTime <= as.POSIXct("2022-04-07 11:00:00"))

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
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-04-07 10:00:00") & DateTime <= as.POSIXct("2022-04-07 10:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*1141 # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2021-04-07 10:16:20") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2021-04-07 11:00:00")
time_diff_sec <- as.numeric(peak_time - inj_time)
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(70) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- mean()
## Calculate effective depth
z <- (Q/1000)/(w*v)
z






#####################
## GBU 2022-04-07 ##
#####################
BW_Hobo <-read.csv("./NA22_dat/GBU_20220407/20775520_GBUBOR20220407_NH4.csv", skip=1)
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

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-04-07 13:50:00") & DateTime <= as.POSIXct("2022-04-07 14:22:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }

sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-04-07 13:50:00") & DateTime <= as.POSIXct("2022-04-07 13:58:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*712.2  # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-04-07 13:59:20") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-04-07 14:22:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
wft <- mean()
w <- wft*0.3048
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

#################### 
## GBL 2022-06-23 ##

BW_Hobo <-read.csv("./NA22_dat/GBL_20220623/GBLBOR20775520_220623NH4.csv", skip=1)
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
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-06-23 10:40:00") & DateTime <= as.POSIXct("2022-06-23 11:30:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-06-23 10:40:00") & DateTime <= as.POSIXct("2022-06-23 10:46:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*800 # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

#inj_time <- as.POSIXct("2021-07-28 15:40:00") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-06-23 10:48:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
wft <- mean(c(4.0833333, 8.75, 6.583333,4.833333, 9,
          3.666667, 5.333333, 6, 9, 8, 4.75, 5.75,
          5.583333))

w<- wft*0.3048
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## ---------------------------
#################### 
## GBU 2022-06-23 ##

BW_Hobo <-read.csv("./NA22_dat/GBU_20220623/GBuplowermixBORNA_220623.csv", skip=1)
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
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-06-23 12:45:00") & DateTime <= as.POSIXct("2022-06-23 13:30:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-04-07 10:00:00") & DateTime <= as.POSIXct("2022-04-07 10:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*800 # NOT sure
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
w <- mean(6.666667, 7.166667, 10.5, 3.5, 4.66667, 6, 4.333, 7, 5.583333, 6.0833333, 4.75, 4.666667,3.4166, 4)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


#################### 
## GBL 2022-10-03 ##

BW_Hobo <-read.csv("./NA22_dat/GBL_20221003/BORGBLandGBU20775520_12.csv", skip=1)
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
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-03 10:40:00") & DateTime <= as.POSIXct("2022-10-03 11:42:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-10-03 10:40:00") & DateTime <= as.POSIXct("2022-10-03 10:50:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-03 10:50:10") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-10-03 11:42:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(1.3,1.75,1,2.1,1.5,1.5,1.25,
          1,1.4,2,1.7,1.8,1.5,1.6,1.3))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z











# Adjust the time NO3 range:
BW_Hobo2 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-03 12:45:00") & DateTime <= as.POSIXct("2022-10-03 13:50:00"))

qplot(DateTime, Cond, data = BW_Hobo2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo2, DateTime >= as.POSIXct("2022-10-03 12:45:00") & DateTime <= as.POSIXct("2022-10-03 12:55:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*800) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo2$DateTime), BW_Hobo2$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-03 12:53:20") #Lolomai 
peak_time <- BW_Hobo2[which.max(BW_Hobo2$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean()
## Calculate effective depth
z <- c(Q/1000)/(w*v)
z


#################### 
## GBU 2022-10-03 ##

BW_Hobo <-read.csv("./NA22_dat/GBL_20221003/BORGBLandGBU20775520_12.csv", skip=1)
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

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

# Adjust the time range:
BW_Hobo1 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-03 14:24:00") & DateTime <= as.POSIXct("2022-10-03 15:15:00"))

qplot(DateTime, Cond, data = BW_Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

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
sub_bg <- subset(BW_Hobo1, DateTime >= as.POSIXct("2022-10-03 14:24:00") & DateTime <= as.POSIXct("2022-10-03 14:31:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo1$DateTime), BW_Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-03 14:31:15") #Lolomai 
peak_time <- BW_Hobo1[which.max(BW_Hobo1$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(1.4,1.2,0.8,0.7,1.6,
          0.5,0.8,1.3,1,1.1,1.4,
          1.3,1.8,1.2,1.2))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


# Adjust the time range:
BW_Hobo2 <- subset(BW_Hobo, DateTime >= as.POSIXct("2022-10-03 15:20:00") & DateTime <= as.POSIXct("2022-10-03 16:32:00"))

qplot(DateTime, Cond, data = BW_Hobo2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


## Reach morphology estimates:
##

## (1) Determine the background conductivity
sub_bg <- subset(BW_Hobo2, DateTime >= as.POSIXct("2022-10-03 15:15:00") & DateTime <= as.POSIXct("2022-10-03 15:38:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*800) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo2$DateTime), BW_Hobo2$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-03  15:52:05") #Lolomai 
peak_time <- BW_Hobo2[which.max(BW_Hobo2$SpCond),]$DateTime 
#end_time <-as.POSIXct("2021-07-28 17:01:50")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(1.4,1.2,0.8,0.7,1.6,
            0.5,0.8,1.3,1,1.1,1.4,
            1.3,1.8,1.2,1.2))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z
