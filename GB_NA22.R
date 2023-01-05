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
#end_time <-as.POSIXct("2021-04-07 11:00:00")
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
out <- out[c(-1,-2, -3),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(datq, aes(datetime, ysi_SPC)) + geom_point(),
  ncol=1, align="hv")
GB_uptake

# ggsave(plot = GB_uptake, filename = paste("./figures/GBU221104v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/GBU_BTC_20221104_MW.csv", row.names = TRUE)

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
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-10-03 10:40:00") & DateTime <= as.POSIXct("2022-10-03 11:42:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
##
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-10-03 10:40:00") & DateTime <= as.POSIXct("2022-10-03 10:50:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100*700) # NOT sure
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

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


### NH4 samples ###
GBL_NH4 <- read.csv("./NA22_dat/GBL_20221003/GBL20221003_NH4v2.csv")
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
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
GBL_NH4v2$Nh4_C <- (GBL_NH4v2$Nh4_mgNL) - Cadd
GBL_NH4v2$Nh4_C <-replace(GBL_NH4v2$Nh4_C, GBL_NH4v2$Nh4_C<0, 0)


#GB_NA[2,7] = 0.05
GBL_NH4v2$SpCond_C <- c(GBL_NH4v2$SpCond  - bg_SpCond)
GBL_NH4v2$SpCond_C <-replace(GBL_NH4v2$SpCond_C, GBL_NH4v2$SpCond_C<0, 0)

## I'm here ###
#No Cl samples so Cl approx.
GBL_NH4v2$Cl_mgL <- ((0.05/0.105)*GBL_NH4v2$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = GBL_NH4v2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g KNO3 in 6L carboy
Nh4mgL <- 300 * (1000) * (18.04/53.491) *(1/6)
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
out <- out[c(-1,-2,-3),]
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

# write.csv(x = out, file = "./BTC_out/GBL_BTC_20221003v2.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
mean(out$sw)
mean(out$Uadd)
mean(Hobo$TempC)
mean(na.omit(dat$PO4_ugL))
mean(na.omit(dat$DOC_mgL))
N_alt<-mean(na.omit(dat$NO3_mgNL))
N_supp <-(86400*Q*(N_alt*0.001))/(w*reachL)

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
reachL <- c(50) #
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

summary(dat)

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

# write.csv(x = out, file = "./BTC_out/GBU_BTC_20221104.csv", row.names = TRUE)

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







