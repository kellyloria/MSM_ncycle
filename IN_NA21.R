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
## GBU 2021-06-03
# Read in the nitrogen uptake assay data:
GN_NA <- read.csv("./NA21_dat/samples/IncL20210805sample.csv")
GN_NA$datetime <- as.POSIXct(paste(GN_NA$date, GN_NA$time), format = "%m/%d/%y %H:%M:%S")

  
  
qplot(datetime, NO3, data = GN_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))





## SPC data from HOBO
BW_Hobo <-read.csv("./NA21_dat/InclineNassay20775515_0.csv", skip=1)
summary(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775515..SEN.S.N..20775515.",
                      "Temp...C..LGR.S.N..20775515..SEN.S.N..20775515.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
range(BW_Hobo$DateTime)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
BW_Hobo <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-08-05 8:45:00") & DateTime <= as.POSIXct("2021-08-05 10:02:00"))

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
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-08-05 8:45:00") & DateTime <= as.POSIXct("2021-08-05 9:00:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo$DateTime), BW_Hobo$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-08-05 09:10:00") #Lolomai 
peak_time <- BW_Hobo[which.max(BW_Hobo$SpCond),]$DateTime 

time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
v <- c(100/time_diff_sec)
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
w <- c(1.905)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z



##
## Left join hobo conductivity estimates with samples:
##

GN_NA <- left_join(GN_NA, BW_Hobo[c("DateTime", "SpCond")],
                   by= c("datetime"="DateTime"))

summary(GN_NA)

#GB_NA <- subset(GB_NA, datetime >= as.POSIXct('2021-07-22 13:40:00'))


qplot(datetime, NO3, data = GN_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

GN_NA <- subset(GN_NA,  datetime>= as.POSIXct("2021-08-05 9:11:00") & datetime <= as.POSIXct("2021-08-05 09:45:00"))




##
## Calculations Notes:
##
# 1. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
GN_NA$NO3_C <- (GN_NA$NO3) 

#GB_NA[2,7] = 0.05
GN_NA$SpCond_C <- c(GN_NA$SpCond - 90.00)
#GB_NA$SpCond_C <-replace(GB_NA$SpCond_C, GB_NA$SpCond_C<0, 0)

#No Cl samples so Cl approx.
GN_NA$Cl_mgL <- ((0.05/0.105)*GN_NA$SpCond_C)


qplot(Cl_mgL, NO3_C, data = GN_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 250g KNO3 in 10L carboy
NO3mgL <- 246 * (1000) * (62/101) *(1/10)
# Carboy concentrations 2000g NaCl in 10L carboy
NaClmgL <- 2000 * (1000) * (35.45/58.44) * (1/10)
carboy <- NO3mgL/NaClmgL

#mass recovery= 
GN_NA$NtoNaCl <-  GN_NA$NO3_C/GN_NA$Cl_mgL
GN_NA$NtoNaCllog <-  log(GN_NA$NO3_C/GN_NA$Cl_mgL)

qplot(datetime, NtoNaCl, data = GN_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

GB_NA$massR <- (carboy)- GB_NA$NtoNaCl
GB_NA$massRPer <- (1-((carboy)- GB_NA$NtoNaCl)/(carboy)) * 100


# The added longitudinal uptake rate(kw-dyn) was calculated by plotting the logged N:Cl of the injectate and each grab sample against stream distance 
# and then calculating the slope between each pair of points (injectate sample and each grab sample).
GN_NA$carboy <- log(carboy)

reachL<- c(100)
## way of iterating slope change between the row values
out <- data.frame(Site = NA, datetime=as.POSIXct(NA), NO3=NA, Cl= NA, stamps = NA, slope_sample=NA, kw = NA)
for (i in 2:nrow(GN_NA)) {
  temp_dat <- GN_NA[c(i-1,i),]
  slope_sample <- (temp_dat$NtoNaCllog[2]-temp_dat$NtoNaCllog[1])/(as.numeric(temp_dat$datetime[2] - temp_dat$datetime[1]))
  kw <- (temp_dat$carboy[2]-temp_dat$NtoNaCllog[1])/(as.numeric(0-reachL))
  datetime<- as.POSIXct((GN_NA$datetime[i]), format="%Y-%m-%d %H:%M:%S") 
  NO3<- GN_NA$NO3_C[i]
  Cl<- GN_NA$Cl_mgL[i]
  temp_out <- data.frame(Site = "INC", 
                         stamps = paste(i, i-1, sep = "-"), 
                         slope_sample = slope_sample, 
                         kw=kw, 
                         datetime=datetime,
                         NO3=NO3,
                         Cl=Cl)
  out <- rbind(out, temp_out)
}

Cadd <- 0.011

out <- out[c(-1),]
out$sw<- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


GB_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(BW_Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")

# ggsave(plot = BW_uptake, filename = paste("./figures/BW220728_v2.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/INC_BTC_20210805.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
mean(BW_Hobo$TempC)
mean(out$sw)
mean(0.04770792, 0.04588845, 0.04707783, 0.04771765, 0.04715911, 0.04636185,
     0.04478110, 0.04794743, 0.04778461, 0.04659000, 0.04579577, 0.04861030)
