library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)
library(lme4)
library(lmerTest)
library(cowplot)
library(drc)
library(dr4pl)

###
## Reach morphology estimates:
##
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }


##########################
# BWU 2022-08-24
# 


##
##########################
### BWL 2022-08-24 ##

Hobo <-read.csv("./NA22_dat/BWL_20220824/BWLNH4BOR_20775520_20220824.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%OSZ") 
range(Hobo$DateTime)
str(Hobo)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-08-24 10:05:00") & DateTime <= as.POSIXct("2022-08-24 11:36:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-08-24 10:05:00") & DateTime <= as.POSIXct("2022-08-24 10:15:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*1500) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-08-24 10:18:05") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-08-24 11:35:30")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(4, 5.5,4.7,6.8,5.9,
            8.1,6.3,6.1,5.7,7.7,
            9,9.3,8.2,7.8,7))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWL_20220824/BWL20220824_NH4v2.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:23),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(24:26),c(4)])

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)


datq[6,11]= 83.62996
datq[8,11]= 127.8487
datq[9,11]= 131.4467
datq[10,11]= 140.8758
datq[11,11]= 142.2619



datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
Nh4mgL <- 300 * (1000) * (18.04/53.491) *(1/10)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 1500 * (1000) * (35.45/58.44) * (1/10)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
datq$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(datq)) {
  temp_dat <- datq[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- datq$Nh4_C[i]
  Cl<- datq$Cl_mgL[i]
  temp_out <- data.frame(Site = "BWL_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2,-3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL220824v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL220824.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$NO3_mgNL))
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)


##########################
### BWL 2022-10-12 ###

Hobo <-read.csv("./NA22_dat/BWL_20221012/20775520_13.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%OSZ") 
range(Hobo$DateTime)
str(Hobo)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-10-12 10:55:00") & DateTime <= as.POSIXct("2022-10-12 13:00:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-10-12 10:49:00") & DateTime <= as.POSIXct("2022-10-12 11:08:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*600)
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-12 11:02:20") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
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

## NH4 sample data ## 
dat <- read.csv("./NA22_dat//BWL_20221012/BWL_20221012_NH4.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:15),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(16:17),c(4)])

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)


datq[2,9]= 56.72769
datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 12L carboy
Nh4mgL <- 300 * (1000) * (18.04/53.491) *(1/12)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 700 * (1000) * (35.45/58.44) * (1/12)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
datq$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(datq)) {
  temp_dat <- datq[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- datq$Nh4_C[i]
  Cl<- datq$Cl_mgL[i]
  temp_out <- data.frame(Site = "BWL_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2, -3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

getwd()

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL221012.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL221012.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))


##########################
### BWU 2022-08-24 ##

BW_Hobo <-read.csv("./NA22_dat/BWU_20220824/BWL_NO3BOR_20775523_20220824BWUNH4HOR.csv", skip=1)
summary(BW_Hobo)
names(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%m/%d/%y %H:%M:%S") #format="%Y-%m-%dT%H:%M:%SZ") 
range(BW_Hobo$DateTime)
str(BW_Hobo)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))



Hobo <-read.csv("./NA22_dat/BWU_20220824/BWL_NO3BOR_20775523_20220824BWUNH4HOR.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-08-24 14:15:00") & DateTime <= as.POSIXct("2022-08-24 14:40:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWU_20220824/BWU20220824_NH4v2.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(dat, datetime <= as.POSIXct("2022-08-24 13:53:00")) #Lolomai
bg_SpCond <- mean(na.omit(sub_bg$ysi_SPC))
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*1500) # NOT sure
## Calculate Q
## Units = L/sec
dat_fl <-na.omit((dat[,c(8,10)]))
Q <- Qint(as.numeric(dat_fl$datetime), na.omit(dat_fl$ysi_SPC), bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-08-24 13:41:10") #Lolomai 
peak_time <- dat_fl[which.max(dat_fl$ysi_SPC),]$datetime 
end_time <-as.POSIXct("2022-08-24 15:15:30")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(4, 2.6,2.8,4.9,4.5,
            3.3,2.4,2.5,3.4,3.2,
            2.1,4.9,3,4.6,5.1))
## Calculate effective depth
z <- c(Q/1000)/(w*v)
z

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:25),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(26),c(4)])

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq$SpCond_C <- c(datq$ysi_SPC  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 12L carboy
Nh4mgL <- 300 * (1000) * (18.04/53.491) *(1/10)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 1500 * (1000) * (35.45/58.44) * (1/10)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
datq$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(datq)) {
  temp_dat <- datq[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- datq$Nh4_C[i]
  Cl<- datq$Cl_mgL[i]
  temp_out <- data.frame(Site = "BWU_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2,-3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(dat_fl, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWU220824v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWU_NH4_BTC_BWU220824.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(na.omit(Hobo$TempC))
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt <- mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)

####################
## BWL 2022-10-12 ##
# 
# Hobo <-read.csv("./NA22_dat/BWL_20221012/20775523_17HOR.csv", skip=1)
# summary(Hobo)
# names(Hobo)
# 
# # modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
# Hobo <- Hobo[,c("Date.Time..GMT.07.00",
#                 "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
#                 "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523.")]
# 
# colnames(Hobo) <- c("DateTime","Cond","TempC")
# # Convert DateTime
# Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%m/%d/%y %H:%M:%S") 
# range(Hobo$DateTime)
# str(Hobo)
# 
# # 
# Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)
# 
qplot(DateTime, TempC, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
# 
# # Adjust the time range:
# Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-08-24 14:15:00") & DateTime <= as.POSIXct("2022-08-24 14:40:00"))
# 
# qplot(DateTime, Cond, data = Hobo, geom="point") +
#   theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWL_20221012/BWL_20221012_NH4v2.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(dat, datetime <= as.POSIXct("2022-10-12 11:20:00")) #Lolomai
bg_SpCond <- mean(na.omit(sub_bg$ysi_SPC))
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*600) # 
## Calculate Q
## Units = L/sec
dat_fl <-na.omit((dat[,c(8,10)]))
Q <- Qint(as.numeric(dat_fl$datetime), na.omit(dat_fl$ysi_SPC), bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-12 10:59:50") #Lolomai 
peak_time <- dat_fl[which.max(dat_fl$ysi_SPC),]$datetime 
end_time <-as.POSIXct("2022-10-12 12:05:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(175) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(6.2,8.5,7,6.4,10,9.7,
            8.3,4.8,5.4,5,6,7.1,
            6.2,4.6,7))
## Calculate effective depth
z <- c(Q/1000)/(w*v)
z

### BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:15),]
## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(16:19),c(4)])

# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq$SpCond_C <- c(datq$ysi_SPC  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 12L carboy
Nh4mgL <- 400 * (1000) * (18.04/53.491) *(1/12)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 600 * (1000) * (35.45/58.44) * (1/12)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
datq$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(datq)) {
  temp_dat <- datq[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- datq$Nh4_C[i]
  Cl<- datq$Cl_mgL[i]
  temp_out <- data.frame(Site = "BWL_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2,-3,-4),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(dat_fl, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL221012.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL221012.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
#mean(na.omit(Hobo$TempC))
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt <- mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)

####################
## BWL 2022-11-21 ##

Hobo <-read.csv("./NA22_dat/BWL_20221121/BWL_BOR_20221121_20775520_19.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%OSZ") 
range(Hobo$DateTime)
str(Hobo)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-11-21 10:00:00") & DateTime <= as.POSIXct("2022-11-21 11:30:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-11-21 10:00:00") & DateTime <= as.POSIXct("2022-11-21 10:15:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*750) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-11-21 10:16:20") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-11-21 11:13:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(150) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(5.1,6.6,5.9,5,5.2,
            3,6.4,8.9,8.2,7.1,
            5,7,6.1,9.4,8.3))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWL_20221121/BWL20221121_NH4.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:21),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(22:23),c(4)])

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - 0
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)


datq[2,11]= 41.41997
datq[5,11]= 42.25523
datq[11,11]= 77.73627

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
Nh4mgL <- 375 * (1000) * (18.04/53.491) *(1/7.5)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 750 * (1000) * (35.45/58.44) * (1/7.5)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
datq$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(datq)) {
  temp_dat <- datq[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- datq$Nh4_C[i]
  Cl<- datq$Cl_mgL[i]
  temp_out <- data.frame(Site = "BWL_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2,-3,-4,-5),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*0.0001/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL221121.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL221121.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(0*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$NO3_mgNL))
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)


##
####################
## BWL 2022-11-21 ##

Hobo <-read.csv("./NA22_dat/BWL_20221121/BWL_BOR_20221121_20775520_19.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%OSZ") 
range(Hobo$DateTime)
str(Hobo)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-11-21 11:50:00") & DateTime <= as.POSIXct("2022-11-21 13:35:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-11-21 11:50:00") & DateTime <= as.POSIXct("2022-11-21 12:15:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*750) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-11-21 12:42:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-11-21 13:35:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(150) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(5.1,6.6,5.9,5,5.2,
            3,6.4,8.9,8.2,7.1,
            5,7,6.1,9.4,8.3))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWL_20221121/BWL20221121_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:19),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1:2),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq[14,11]= 61.32323


datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 375 * (1000) * (62/101) *(1/7.5)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 750 * (1000) * (35.45/58.44) * (1/7.5)
carboy <- NmgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$NO3_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$NO3_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
datq$carboy <- log(carboy)

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NO3=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(datq)) {
  temp_dat <- datq[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((datq$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NO3<- datq$NO3_C[i]
  Cl<- datq$Cl_mgL[i]
  temp_out <- data.frame(Site = "BWL_NO3", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*0.0001/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL_NO3_221121.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NO3_BTC_BWL221121.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$Nh4_mgNL))
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)

