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
### GBL 2022-04-07 ##

Hobo <-read.csv("./NA22_dat/GBL_20220407/Copy\ of\ 20775523_GBLBOR20220407_NO3.csv", skip=1)
summary(Hobo)
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

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-04-07 10:00:00") & DateTime <= as.POSIXct("2022-04-07 11:00:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-04-07 10:00:00") & DateTime <= as.POSIXct("2022-04-07 10:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*1141 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2021-04-07 10:16:20") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2021-04-07 10:40:10")
time_diff_sec <- as.numeric(peak_time - inj_time)
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(70) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- c(2.05)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


## NH4 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20220407/GBL20220407_NH4v2.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:18),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(19:20),c(4)])

# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq$SpCond_C <- c(datq$ysi_SPC  - c(428))
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations N in 12L carboy
Nh4mgL <- 189.2 * (1000) * (18.04/53.491) *(1/9.46353)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 704.6 * (1000) * (35.45/58.44) * (1/9.46353)
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
out <- out[c(-1,-2),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(datq, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GB221104v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_BTC_NH4_20220407.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<-mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)

######################
## NO3 sample data ## 
dat <- read.csv("./NA22_dat/GBU_20220407/GBL20220407_NO3.csv")
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
datq <- dat[c(1:17),]

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(18:19),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

# No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 308.0 * (1000) * (62/101) *(1/9.5)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 1141.0 * (1000) * (35.45/58.44) * (1/9.5)
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
out <- out[c(-1,-2),]
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

# write.csv(x = out, file = "./BTC_out/GBL_NO3_BTC_220407.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)




####################
## GBU 2022-04-07 ##
Hobo <-read.csv("./NA22_dat/GBU_20220407/20775520_GBUBOR20220407_NO3.csv", skip=1)
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

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-04-07 12:40:00") & DateTime <= as.POSIXct("2022-04-07 13:15:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-04-07 12:40:00") & DateTime <= as.POSIXct("2022-04-07 12:45:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*1060.8 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2021-04-07 12:42:30") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime #
end_time <-as.POSIXct("2021-04-07 13:03:00")
time_diff_sec <- as.numeric(peak_time - inj_time)
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- c(1.85)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


## NO3 sample data ## 
dat <- read.csv("./NA22_dat/GBU_20220407/GBU20220407_NO3.csv")
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
datq <- dat[c(1:19),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1:2),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 300.6 * (1000) * (62/101) *(1/9.5)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 1060.8 * (1000) * (35.45/58.44) * (1/9.5)
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
  temp_out <- data.frame(Site = "GBU_NO3", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2, -3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBU_NO3_220407.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBU_NO3_BTC_220407.csv", row.names = TRUE)

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


########################
## GBU NO3 2022-06-23 ##

Hobo <-read.csv("./NA22_dat/GBU_20220623/GBuplowermixBORNA_220623.csv", skip=1)
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

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-06-23 13:30:00") & DateTime <= as.POSIXct("2022-06-23 14:13:00"))
qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-06-23 13:30:00") & DateTime <= as.POSIXct("2022-06-23 13:40:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*800 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-06-23 13:54:45") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime #
end_time <-as.POSIXct("2022-06-23 14:13:00")
time_diff_sec <- as.numeric(peak_time - inj_time) *60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
wft <- mean(c(6.666667, 7.166667, 10.5, 3.5, 4.666667,
       4.583333, 6, 4.3333, 7, 5.583333, 6.0833, 4.75,
       4.66667, 3.41667, 4))
w <- 0.3048*wft
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


## NO3 sample data ## 
dat <- read.csv("./NA22_dat/GBL20230623/GBU20220623_NO3.csv")
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
## Cadd geometric mean of background concentrations 
Cadd <- c(0.01)
datq <- dat

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)
datq[2,6] <-0.010

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 250 * (1000) * (62/101) *(1/7.6)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 800 * (1000) * (35.45/58.44) * (1/7.6)
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
out <- out[c(-1,-2, -3, -4,-5),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w

GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL_NO3_220623.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_NO3_BTC_220623.csv", row.names = TRUE)

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
###
###
###
###
###

## GBU 2022-06-23 NH4

Hobo <-read.csv("./NA22_dat/GBU_20220623/GBuplowermixBORNA_220623.csv", skip=1)
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

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-06-23 12:45:00") & DateTime <= as.POSIXct("2022-06-23 13:30:00"))
qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-06-23 12:45:00") & DateTime <= as.POSIXct("2022-06-23 12:48:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*800 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-06-23 12:50:30") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime #
end_time <-as.POSIXct("2022-06-23 13:12:00")
time_diff_sec <- as.numeric(peak_time - inj_time) *60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
wft <- mean(c(6.666667, 7.166667, 10.5, 3.5, 4.666667,
              4.583333, 6, 4.3333, 7, 5.583333, 6.0833, 4.75,
              4.66667, 3.41667, 4))
w <- 0.3048*wft
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


## NH4 sample data ## 
dat <- read.csv("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/NA22_dat/GBU_20220623/GBU20220623_NH4v2.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%m/%d/%y %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

# 1. select the sample selection for: GBL_NH4
## Cadd geometric mean of background concentrations 
Cadd <- mean(0.0028,0.0157)

# GBL BTC ###
# 1. select the sample selection for: GBL_NH4
#datq <- dat[c(-2, -4, -6, -10),]
datq <- dat

# leftjoin 

## Cadd geometric mean of background concentrations 

# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations N in 12L carboy
Nh4mgL <- 250 * (1000) * (18.04/53.491) *(1/9.46353)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 800 * (1000) * (35.45/58.44) * (1/9.46353)
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
  temp_out <- data.frame(Site = "GBU_NH4", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NH4=NH4,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2, -3, -4),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(datq, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GB221104v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBU_BTC_NH4_20220623.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(0.021*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<-mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)


###
##
##


## NO3 sample data ## 
dat <- read.csv("./NA22_dat/GBU_20220623/GBU20220623_NO3.csv")
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
## Cadd geometric mean of background concentrations 
Cadd <- mean(na.omit(dat[c(1:3),c(6)]))
datq <- dat

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 250 * (1000) * (62/101) *(1/7.6)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 800 * (1000) * (35.45/58.44) * (1/7.6)
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
  temp_out <- data.frame(Site = "GBU_NO3", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

## Cadd geometric mean of background concetrations 
out <- out[c(-1,-2, -3, -4,-5),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w

GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBU_NO3_220623.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBU_NO3_BTC_220623.csv", row.names = TRUE)

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
##
##

#################### 
## GBL NH4 2022-06-23 ##

Hobo <-read.csv("./NA22_dat/GBL_20220623/GBLBOR20775520_220623NH4.csv", skip=1)
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

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-06-23 10:30:00") & DateTime <= as.POSIXct("2022-06-23 11:35:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-06-23 10:35:00") & DateTime <= as.POSIXct("2022-06-23 10:40:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*800) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-06-23 10:45:20") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-06-23 11:18:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
wft <- mean(c(4.167, 3.75, 3.583, 4.83, 3.667, 5.33, 5, 4.75, 5.75, 5.583))
w <- wft * 0.3048
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


### NH4 samples ###
GBL_NH4 <- read.csv("./NA22_dat/GBL20230623/GBL20220623_NH4.csv")
GBL_NH4$datetime <- as.POSIXct(paste(GBL_NH4$date, GBL_NH4$time), format = "%Y-%m-%d %H:%M:%S")
str(GBL_NH4)

GBL_NH4 <- left_join(GBL_NH4, Hobo[c("DateTime", "SpCond")],
                     by= c("datetime"="DateTime"))

summary(GBL_NH4)

qplot(datetime, Nh4_mgNL, data = GBL_NH4, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
GBL_NH4v2 <- GBL_NH4[c(1:11),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- c(0.015) # mean(GBL_NH4[c(1,11),c(4)])

# 2. Correct for background concentrations (_C):
GBL_NH4v2$Nh4_C <- (GBL_NH4v2$Nh4_mgNL) - Cadd
GBL_NH4v2$Nh4_C <-replace(GBL_NH4v2$Nh4_C, GBL_NH4v2$Nh4_C<0, 0)

GBL_NH4v2$SpCond_C <- c(GBL_NH4v2$SpCond  - bg_SpCond)
GBL_NH4v2$SpCond_C <-replace(GBL_NH4v2$SpCond_C, GBL_NH4v2$SpCond_C<0, 0)

## I'm here ###
#No Cl samples so Cl approx.
GBL_NH4v2$Cl_mgL <- ((0.05/0.105)*GBL_NH4v2$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = GBL_NH4v2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 6L carboy
Nh4mgL <- 250 * (1000) * (18.04/53.491) *(1/6)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 800 * (1000) * (35.45/58.44) * (1/6)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
GBL_NH4v2$NtoNaCl <-  GBL_NH4v2$Nh4_C/GBL_NH4v2$Cl_mgL
GBL_NH4v2$NtoNaCllog <-  log(GBL_NH4v2$Nh4_C/GBL_NH4v2$Cl_mgL)

qplot(datetime, NtoNaCllog, data = GBL_NH4v2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

GBL_NH4v2$massR <- (carboy)- GBL_NH4v2$NtoNaCl
GBL_NH4v2$massRPer <- (1-((carboy)- GBL_NH4v2$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
GBL_NH4v2$carboy <- log(carboy)

dat<- GBL_NH4v2

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(dat)) {
  temp_dat <- dat[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((dat$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- dat$Nh4_C[i]
  Cl<- dat$Cl_mgL[i]
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
out <- out[c(-1,-2,-5),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL220103v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_BTC_NH4_20220623.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(dat$PO4_ugL))
mean(na.omit(dat$DOC_mgL))
N_alt<-mean(na.omit(dat$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)









##
##
#################### 
## GBL NH4 2022-10-03 ##

Hobo <-read.csv("./NA22_dat/GBL_20221003/BORGBLandGBU20775520_12.csv", skip=1)
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

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-10-03 10:40:00") & DateTime <= as.POSIXct("2022-10-03 11:43:30"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-10-03 10:40:00") & DateTime <= as.POSIXct("2022-10-03 10:45:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-03 10:50:10") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-10-03 11:43:00")
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


### NH4 samples ###
GBL_NH4 <- read.csv("./NA22_dat/GBL_20221003/GBL20221003_NH4.csv")
GBL_NH4$datetime <- as.POSIXct(paste(GBL_NH4$date, GBL_NH4$time), format = "%Y-%m-%d %H:%M:%S")
str(GBL_NH4)

GBL_NH4 <- left_join(GBL_NH4, Hobo[c("DateTime", "SpCond")],
                     by= c("datetime"="DateTime"))

summary(GBL_NH4)

qplot(datetime, Nh4_mgNL, data = GBL_NH4, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
GBL_NH4v2 <- GBL_NH4[c(1:23),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(GBL_NH4[c(25:27),c(4)])

# 2. Correct for background concentrations (_C):
GBL_NH4v2$Nh4_C <- (GBL_NH4v2$Nh4_mgNL) - Cadd
GBL_NH4v2$Nh4_C <-replace(GBL_NH4v2$Nh4_C, GBL_NH4v2$Nh4_C<0, 0)

GBL_NH4v2$SpCond_C <- c(GBL_NH4v2$SpCond  - bg_SpCond)
GBL_NH4v2$SpCond_C <-replace(GBL_NH4v2$SpCond_C, GBL_NH4v2$SpCond_C<0, 0)

## I'm here ###
#No Cl samples so Cl approx.
GBL_NH4v2$Cl_mgL <- ((0.05/0.105)*GBL_NH4v2$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = GBL_NH4v2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 6L carboy
Nh4mgL <- 400 * (1000) * (18.04/53.491) *(1/6)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 700 * (1000) * (35.45/58.44) * (1/6)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
GBL_NH4v2$NtoNaCl <-  GBL_NH4v2$Nh4_C/GBL_NH4v2$Cl_mgL
GBL_NH4v2$NtoNaCllog <-  log(GBL_NH4v2$Nh4_C/GBL_NH4v2$Cl_mgL)

qplot(datetime, NtoNaCllog, data = GBL_NH4v2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

GBL_NH4v2$massR <- (carboy)- GBL_NH4v2$NtoNaCl
GBL_NH4v2$massRPer <- (1-((carboy)- GBL_NH4v2$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
GBL_NH4v2$carboy <- log(carboy)

dat<- GBL_NH4v2

## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NH4=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(dat)) {
  temp_dat <- dat[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((dat$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NH4<- dat$Nh4_C[i]
  Cl<- dat$Cl_mgL[i]
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
out <- out[c(-1, -2, -3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL220103v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_BTC_NH4_20221003.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(dat$PO4_ugL))
mean(na.omit(dat$DOC_mgL))
N_alt<-0.006
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)

## 2022-10-03 NO3 ###
## NO3 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20221003/GBL20221003_NO3.csv")
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
datq <- dat[c(1:21),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(22:23),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq[12,11]= 820.127
datq[21,11]= 388.6032

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 400 * (1000) * (62/101) *(1/6)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 800 * (1000) * (35.45/58.44) * (1/6)
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
datq<-datq[c(-7),]

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
out <- out[c(-1,-3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL_NO3_221003.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_NO3_BTC_221003.csv", row.names = TRUE)

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

#### M-M curve fit -- Error here.
library(dr4pl)
library(drc)

model.drm1 <- drc::drm (Uadd ~ NH4, data = plot_out, fct = MM.3())
summary(model.drm1)

mm2 <- data.frame(NH4_C = seq(0, max(plot_out$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


Uadd_plotBW <- ggplot(plot_out, aes(x=NH4, y=Uadd)) + 
  #ylim(0,20) + xlim(0,1500) + 
  geom_line(data = mm2[-1,], aes(x = NH4_C, y = Uadd), colour = "black") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate ugL")),
       y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
  theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))

#########################
## GBL NH4 2022-10-03 ##
#########################

Hobo <-read.csv("./NA22_dat/GBL_20221212/GBL_221212_20775520_BOR.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-12-12 12:58:00") & DateTime <= as.POSIXct("2022-12-12 14:15:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-12-12 12:58:00") & DateTime <= as.POSIXct("2022-12-12 13:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*500) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-12-12 13:27:15") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-12-12 14:13:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(50) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(1.8, 2.2, 1.3, 1.2, 1.7, 2, 1.8, 2, 2))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## 2022-10-03 NO3 ###
## NO3 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20221212/GBL20221212_NO3.csv")
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
datq <- dat[c(1:19),]

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(20:21),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 150 * (1000) * (62/101) *(1/6)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 500 * (1000) * (35.45/58.44) * (1/6)
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
out <- out[c(-1,-2,-3,-19),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL_NO3_221212.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_NO3_BTC_221212.csv", row.names = TRUE)

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




##########################
### GBU 2022-10-03 NH4 ###

Hobo <-read.csv("./NA22_dat/GBL_20221003/BORGBLandGBU20775520_12.csv", skip=1)
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

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-10-03 14:24:00") & DateTime <= as.POSIXct("2022-10-03 15:15:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-10-03 14:24:00") & DateTime <= as.POSIXct("2022-10-03 14:31:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-10-03 14:32:15") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-10-03 15:10:10")
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

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/GBU_20221003/GBU20221003_NH4v2.csv")
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
datq <- dat[c(1:19),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(20:21),c(4)])

# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)


datq[19,11]= 418.951
datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 6L carboy
Nh4mgL <- 200 * (1000) * (18.04/53.491) *(1/6)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 700 * (1000) * (35.45/58.44) * (1/6)
carboy <- Nh4mgL/NaClmgL

# mass recovery = 
datq$NtoNaCl <-  datq$Nh4_C/datq$Cl_mgL
datq$NtoNaCllog <-  log(datq$Nh4_C/datq$Cl_mgL)

qplot(datetime, NtoNaCllog, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

datq$massR <- (carboy)- datq$NtoNaCl
datq$massRPer <- (1-((carboy)- datq$NtoNaCl)/(carboy)) * 100

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
  temp_out <- data.frame(Site = "GBU_NH4", 
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


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# ggsave(plot = GB_uptake, filename = paste("./figures/GBU220103v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBU_BTC_20221003v2.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<- mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)


plot_out<- out[c(2:21),]

##########################
### GBL 2022-11-04 ##

Hobo <-read.csv("./NA22_dat/GBL_20221104/20775523_19_BOR.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
                "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523.")]

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
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-11-04 11:10:00") & DateTime <= as.POSIXct("2022-11-04 12:35:50"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-11-04 11:10:00") & DateTime <= as.POSIXct("2022-11-04 11:30:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-11-04 11:26:40") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-11-04 12:34:30")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(75) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(1.42,1.22,0.82,0.72,1.62,
            0.52,0.82,1.32,1,1.1,1.4,
            1.3,1.8,1.2,1.2))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20221104/GBL20221104_NH4v2.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:20),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(21:22),c(4)])

# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)


datq[17,11]= 349.1365
datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((18.04/53.491)*datq$SpCond_C)

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
out <- out[c(-1,-2,-6),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# ggsave(plot = GB_uptake, filename = paste("./figures/GBU221104.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_BTC_NH4_20221104.csv", row.names = TRUE)

# estimate N supply:
Cadd <- 0.019
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))

N_alt<- mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)


##########################
## GBL 2022-12-21 ##
Hobo <-read.csv("./NA22_dat/GBL_20221212/GBL_221212_20775520_BOR.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-12-12 11:50:00") & DateTime <= as.POSIXct("2022-12-12 13:00:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-12-12 11:50:00") & DateTime <= as.POSIXct("2022-12-12 11:57:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*500 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-12-12 11:54:30") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-12-12 12:40:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(70) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- mean(c(1.8, 2.2, 1.3, 1.2, 1.7, 2, 1.8, 2, 2))
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


## NH4 sample data ## 
dat <- read.csv("./NA22_dat/GBL_20221212/GBL20221212_NH4.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

### GBL BTC ###
# 1. select the sample selection for: GBL_NH4
datq <- dat[c(1:18),]

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(19:20),c(4)])

# 2. Correct for background concentrations (_C):
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq[1,11]= 318.4727
datq[2,11]= 318.634
datq[8,11]=376.859

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations N in 12L carboy
Nh4mgL <- 150 * (1000) * (18.04/53.491) *(1/9.46353)
# Carboy concentrations 700g NaCl in 6L carboy
NaClmgL <- 500 * (1000) * (35.45/58.44) * (1/9.46353)
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
out <- out[c(-1,-2,-3,-4),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(datq, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBL221212.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_BTC_221212.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
N_alt<-mean(na.omit(datq$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)



###### 2022-12-12 #######
## NO3 sample data ## 
Hobo <-read.csv("./NA22_dat/GBL_20221212/GBL_221212_20775520_BOR.csv", skip=1)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(as.character(Hobo$DateTime), format="%Y-%m-%dT%H:%M:%SZ") 
# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-12-12 13:00:00") & DateTime <= as.POSIXct("2022-12-12 14:15:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-12-12 13:00:00") & DateTime <= as.POSIXct("2022-12-12 13:10:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*500 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-12-12 13:25:30") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-12-12 14:15:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(70) #
v <- reachL/time_diff_sec
v

## Enter average width measurement in m
w <- mean(c(1.8, 2.2, 1.3, 1.2, 1.7, 2, 1.8, 2, 2))
## Calculate effective depth
z <- (Q/1000)/(w*v)
z




dat <- read.csv("./NA22_dat/GBL_20221212/GBL20221212_NO3.csv")
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
datq <- dat[c(1:19),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1:2),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C <0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NmgL <- 150 * (1000) * (62/101) *(1/9.5)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <- 500 * (1000) * (35.45/58.44) * (1/9.5)
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
out <- out[c(-1,-2, -3, -5,-19),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w

GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBU_NO3_220407.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBL_NO3_BTC_221212.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- 0.039
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)



