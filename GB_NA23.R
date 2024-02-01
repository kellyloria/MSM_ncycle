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

## Reach morphology estimates:
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }

##########################
### GBL 2022-03-27 ##
Hobo <-read.csv("./NA22_dat/GBL_20230327/20775520_21.csv", skip=1)
summary(Hobo)
names(Hobo)

Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]
# # modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
# Hobo <- Hobo[,c("Date.Time..GMT.07.00",
#                 "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
#                 "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523.")]
colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)
# Adjust the time range:
Hobo1 <- subset(Hobo, DateTime >= as.POSIXct("2023-03-27 12:00:00") & DateTime <= as.POSIXct("2023-03-27 12:40:30"))

qplot(DateTime, Cond, data = Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo1, DateTime >= as.POSIXct("2023-03-27 12:00:00") & DateTime <= as.POSIXct("2023-03-27 12:05:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*345 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo1$DateTime), Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-03-27 12:15:10") #Lolomai 
peak_time <- Hobo1[which.max(Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-03-27 12:40:10")
time_diff_sec <- as.numeric(peak_time - inj_time) * 60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- reachL/time_diff_sec
v
## Enter average width measurement in m
w <- mean(c(2.8, 3.6, 3.6, 4.5, 2.9, 4.3, 1.9, 2.3, 2.1, 1.5, 1.7, 2.1, 2.0, 1.8, 2.4))
## Calculate effective depth
z <- (Q/1000)/(w*v)
z
## NH4 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20230327/GBL_20230327_NH4.csv")
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
datq <- dat[c(1:18),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- 0
# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations N in 6 L carboy
Nh4mgL <- 110 * (1000) * (18.04/53.491) *(1/6)
# Carboy concentrations NaCl in 6 L carboy
NaClmgL <- 345 * (1000) * (35.45/58.44) * (1/6)
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
  temp_out <- data.frame(Site = "GBL_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2, -3, -5),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(datq, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# write.csv(x = out, file = "./BTC_out/GBL_BTC_NH4_20230327.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo1$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<-mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)


###### NO3 ############
# Adjust the time range:
Hobo1 <- subset(Hobo, DateTime >= as.POSIXct("2023-03-27 13:15:00") & DateTime <= as.POSIXct("2023-03-27 14:00:30"))

qplot(DateTime, Cond, data = Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo1, DateTime >= as.POSIXct("2023-03-27 13:15:00") & DateTime <= as.POSIXct("2023-03-27 13:20:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*365 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo1$DateTime), Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-03-27 13:36:00") #Lolomai 
peak_time <- Hobo1[which.max(Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-03-27 13:55:00")
time_diff_sec <- as.numeric(peak_time - inj_time) * 60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- reachL/time_diff_sec
v
## Enter average width measurement in m
w <- mean(c(2.8, 3.6, 3.6, 4.5, 2.9, 4.3, 1.9, 2.3, 2.1, 1.5, 1.7, 2.1, 2.0, 1.8, 2.4))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## Read in sample data:
dat <- read.csv("./NA22_dat//GBL_20230327/GBL_20230327_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:18),]

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(25,26),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

# No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations N in 6 carboy
NmgL <- 105 * (1000) * (62/101) *(1/6)
# Carboy concentrations NaCl in 6L carboy
NaClmgL <- 365 * (1000) * (35.45/58.44) * (1/6)
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
  temp_out <- data.frame(Site = "GBL_NO3", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w

GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL_NO3_220407.png",sep=""),width=4,height=7,dpi=300)
# write.csv(x = out, file = "./BTC_out/GBL_NO3_BTC_230327.csv", row.names = TRUE)

# estimate N supply:
Nalt <- 0.001
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo1$TempC)
mean(na.omit(datq$DOC_mgL))
Nalt<- mean(na.omit(datq$Nh4_mgNL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)

############################
######## 2022-06-15 ########
### GBL NH4 2022-06-15 ##
Hobo <-read.csv("./NA22_dat/GBL_20230615/20775523_25BOR.csv", skip=1)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat)
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
                "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523.")]
colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)
# Adjust the time range:
Hobo1 <- subset(Hobo, DateTime >= as.POSIXct("2023-06-15 11:15:00") & DateTime <= as.POSIXct("2023-06-15 11:35:00"))

qplot(DateTime, Cond, data = Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo1, DateTime >= as.POSIXct("2023-06-15 11:15:00") & DateTime <= as.POSIXct("2023-06-15 11:20:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*700 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo1$DateTime), Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-06-15 11:22:00") #Lolomai 
peak_time <- Hobo1[which.max(Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-06-15 11:30:00")
time_diff_sec <- as.numeric(peak_time - inj_time) * 60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(75) #
v <- reachL/time_diff_sec
v
## Enter average width measurement in m
w <- mean(c(2.91, 3.42, 3.57, 4.01, 4.32, 3.70, 2.15, 2.22, 1.79, 1.94, 2.34, 2.12, 1.99, 2.58))
## Calculate effective depth
z <- (Q/1000)/(w*v)
z
## NH4 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20230615/GBL_20230615_NH4.csv")
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
Cadd <- dat[c(25),c(4)]
# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq[4,11]= c(201.4118)
datq[5,11]= c(208.9371)
datq[6,11]= c(199.9065)
datq[7,11]= c(171.8061)


datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations  N in 6 L carboy
Nh4mgL <- 200 * (1000) * (18.04/53.491) *(1/8)
# Carboy concentrations  NaCl in 6L carboy
NaClmgL <- 700 * (1000) * (35.45/58.44) * (1/8)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

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
  temp_out <- data.frame(Site = "GBL_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 

out <- out[c(-1,-2, -14, -15),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(datq, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# write.csv(x = out, file = "./BTC_out/GBL_BTC_NH4_20230615.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo1$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<-mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)


###### NO3 ############
# Adjust the time range:
Hobo1 <- subset(Hobo, DateTime >= as.POSIXct("2023-06-15 12:05:00") & DateTime <= as.POSIXct("2023-06-15 12:30:30"))

qplot(DateTime, Cond, data = Hobo1, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo1, DateTime >= as.POSIXct("2023-06-15 12:05:00") & DateTime <= as.POSIXct("2023-06-15 12:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*700 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo1$DateTime), Hobo1$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-06-15 12:15:00") #Lolomai 
peak_time <- Hobo1[which.max(Hobo1$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-06-15 12:30:00")
time_diff_sec <- as.numeric(peak_time - inj_time) * 60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(75) #
v <- reachL/time_diff_sec
v
## Enter average width measurement in m
w <- mean(c(2.91, 3.42, 3.57, 4.01, 4.32, 3.70, 2.15, 2.22, 1.79, 1.94, 2.34, 2.12, 1.99, 2.58))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## Read in sample data:
dat <- read.csv("./NA22_dat/GBL_20230615/GBL_20230615_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:15),]

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(25),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

# No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations N in 6 carboy
NmgL <- 125 * (1000) * (62/101) *(1/8)
# Carboy concentrations NaCl in 6L carboy
NaClmgL <- 700 * (1000) * (35.45/58.44) * (1/8)
carboy <- NmgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$NO3_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$NO3_C/datq$Cl_mgL)

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
  temp_out <- data.frame(Site = "GBL_NO3", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w
out <- out[c(-1, -3, -15),]

GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo1$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<-mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)

