
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
# BWL 2022-05-26
# BWL20220526_NH4 "./NA22_dat/BWL_20220526/BWL20220526_NH4.csv"

Hobo <-read.csv("./NA22_dat/BWL_20230215/20775520_21_NO3uptake.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520." )]

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
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-02-15 11:30:00") & DateTime <= as.POSIXct("2023-02-15 12:30:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-02-15 11:30:00") & DateTime <= as.POSIXct("2023-02-15 11:40:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 800) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-02-15 11:55:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-02-15 12:15:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(7.5, 8.2, 10.6, 8.02, 11.7, 3.8, 
            7.5, 10.5, 7.5, 8.3, 3.4, 
            3.1, 2.5, 3.5, 2.7))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWL_20230215/BWL_20230215_NH4.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:18),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1,2,3),c(4)])

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq[2,12]= c(0.003)
datq[3,12]= c(0.004)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
Nh4mgL <- 217.7 * (1000) * (18.04/53.491) *(1/8)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  800 * (1000) * (35.45/58.44) * (1/8)
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
out <- out[c(-1,-8, -11, -12, -13, -14, -15),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL230215.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
Nalt <- mean(na.omit(datq$NO3_mgNL))
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)
#######################

### NO3 #################
Hobo <-read.csv("./NA22_dat/BWL_20220526/20775523_BWL20220526BOR.csv", skip=1)
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
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2022-05-26 10:10:00") & DateTime <= as.POSIXct("2022-05-26 11:50:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2022-05-26 10:10:00") & DateTime <= as.POSIXct("2022-05-26 10:20:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 1740) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2022-05-26 10:41:10") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2022-05-26 10:51:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(8, 10.2,7.2,9,8.7,
            10.5,9.5,8.3,9.2,8.1,
            7,9.5,10.1,8.5,9))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

##  sample data ## 
dat <- read.csv("./NA22_dat/BWL_20220526/BWL20220526_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:11),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1, 13:15),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C<0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NO3mgL <- 452.1 * (1000) * (62/101) *(1/8)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  1763.3 * (1000) * (35.45/58.44) * (1/8)
carboy <- NO3mgL/NaClmgL

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
  NO3<- datq$NO3_mgNL[i]
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
out <- out[c(-1,-2,-11),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL_NO3_220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NO3_BTC_BWL220526.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$Nh4_mgNL))
Nalt <- 0.004
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)

######
### NO3 #################
Hobo <-read.csv("./NA22_dat/BWL_20230215/20775520_21_NO3uptake.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.08.00",
                "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520." )]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
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
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-02-15 12:40:00") & DateTime <= as.POSIXct("2023-02-15 13:50:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))




## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-02-15 12:40:00") & DateTime <= as.POSIXct("2023-02-15 12:50:00")) #Lolomai
bg_SpCond <- mean(na.omit(sub_bg$SpCond))
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 800) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-02-15 13:19:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-02-15 13:41:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(7.5, 8.2, 10.6, 8.02, 
             11.7, 3.8, 7.5, 10.5, 
             7.5, 8.3, 3.4, 3.1, 
             2.5, 3.5, 2.7))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

##  sample data ## 
dat <- read.csv("./NA22_dat/BWL_20230215/BWL_20230215_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:18),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1,19),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C<0, 0)

datq[16,11]= c(30.64700)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NO3mgL <- 215 * (1000) * (62/101) *(1/7.8)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  800 * (1000) * (35.45/58.44) * (1/7.8)
carboy <- NO3mgL/NaClmgL

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
  NO3<- datq$NO3_mgNL[i]
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
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL_NO3_220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NO3_BTC_BWL23021.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$Nh4_mgNL))
Nalt <- 0.004
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)

###########################
##########################
# BWL 2022-04-06

Hobo <-read.csv("./NA22_dat/BWL_20230405/20775523_24.csv", skip=1)
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

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-04-05 11:25:00") & DateTime <= as.POSIXct("2023-04-05 12:35:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-04-05 11:25:00") & DateTime <= as.POSIXct("2023-04-05 11:30:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 1001) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-04-05 11:45:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-04-05 12:22:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(6.8, 7.09, 8.24))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BWL_20230405/BWL_20230405_NH4.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:19),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- c(0.001)

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq[2,12]= c(0.003)
datq[3,12]= c(0.004)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
Nh4mgL <- 215 * (1000) * (18.04/53.491) *(1/8)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  1001 * (1000) * (35.45/58.44) * (1/8)
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
out <- out[c(-1,-2,-3, -4, -5, -14,-15,-16),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL230405.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
Nalt <- mean(na.omit(datq$NO3_mgNL))
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)
#######################


### NO3 #################
Hobo <-read.csv("./NA22_dat/BWL_20230405/20775523_24.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
                "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523.")]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime


colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
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
Hoboq <- subset(Hobo, DateTime >= as.POSIXct("2023-04-05 13:00:00") & DateTime <= as.POSIXct("2023-04-05 13:40:00"))

qplot(DateTime, Cond, data = Hoboq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-04-05 13:00:00") & DateTime <= as.POSIXct("2023-04-05 13:01:00")) #Lolomai
bg_SpCond <- mean(na.omit(sub_bg$SpCond))
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 1064) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hoboq$DateTime), Hoboq$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-04-05 13:18:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-04-05 13:45:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(90) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(7.377))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

##  sample data ## 
dat <- read.csv("./NA22_dat/BWL_20230405/BWL_20230405_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:14),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- c(0.003)

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C<0, 0)

datq[16,11]= c(30.64700)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NO3mgL <- 229 * (1000) * (62/101) *(1/7.8)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  1064 * (1000) * (35.45/58.44) * (1/7.8)
carboy <- NO3mgL/NaClmgL

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
  NO3<- datq$NO3_mgNL[i]
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
out <- out[c(-1,-2,-3,-4,-14),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL_NO3_220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NO3_BTC_BWL230405.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$Nh4_mgNL))
Nalt <- 0.004
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)

###################
##########################
# BWL 2022-04-06

Hobo <-read.csv("./NA22_dat/BW_20230718/20775523_BOR.csv", skip=1)
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

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 10:00:00") & DateTime <= as.POSIXct("2023-07-18 10:50"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 10:00:00") & DateTime <= as.POSIXct("2023-07-18 10:15:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 1000) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-07-18 10:20:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-07-18 10:44:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(5.56, 6.01, 10.62, 14.23, 10.97, 9.84, 10.34, 11.12, 15.87, 12.00, 7.78, 10.15, 9.98, 9.23, 8.42))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

## NH4 sample data ## 
dat <- read.csv("./NA22_dat/BW_20230718/BWL_20230718_NH4.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

dat[5,11]= c(26.98148)
dat[13,4]= c(0.00565)
# 26.98148

qplot(datetime, Nh4_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:14),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(0.00540, 0.00110)

# 2. Correct for background concentrations (_C):
#GB_NA$NO3_C <- (GB_NA$Results-0.021) 
#GB_NA$NO3_CC <-replace(GB_NA$NO3_C, GB_NA$NO3_C<0, 0) # Na's produced in TMR calculations if 0
datq$Nh4_C <- (datq$Nh4_mgNL) - Cadd
datq$Nh4_C <-replace(datq$Nh4_C, datq$Nh4_C<0, 0)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, Nh4_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
Nh4mgL <- 250 * (1000) * (18.04/53.491) *(1/9)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  1000 * (1000) * (35.45/58.44) * (1/9)
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
out <- out[c(-1,-2,-3,-13, -14),]
out$sw <- -1/(out$kw)
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NH4, sw)) + geom_point(),
  ggplot(out, aes(NH4, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NH4/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NH4_BTC_BWL230718.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
mean(na.omit(datq$DOC_mgL))
Nalt <- mean(na.omit(datq$NO3_mgNL))
Nalt
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)
#######################
## NO3 ##
# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 10:45:00") & DateTime <= as.POSIXct("2023-07-18 11:16:30"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 10:45:00") & DateTime <= as.POSIXct("2023-07-18 10:50:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 1000) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-07-18 10:54:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-07-18 11:17:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(5.56, 6.01, 10.62, 14.23, 10.97, 9.84, 10.34, 11.12, 15.87, 12.00, 7.78, 10.15, 9.98, 9.23, 8.42))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z

##  sample data ## 
dat <- read.csv("./NA22_dat/BW_20230718/BWL_20230718_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

dat <- left_join(dat, Hobo[c("DateTime", "SpCond")],
                 by= c("datetime"="DateTime"))

summary(dat)

qplot(datetime, NO3_mgNL, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("1 min"))

# 1. select the sample selection for:
datq <- dat[c(1:14),]
# leftjoin 

## Cadd geometric mean of background concentrations 
Cadd <- mean(dat[c(1,14),c(6)])

# 2. Correct for background concentrations (_C):
datq$NO3_C <- (datq$NO3_mgNL) - Cadd
datq$NO3_C <-replace(datq$NO3_C, datq$NO3_C<0, 0)

#datq[16,11]= c(30.64700)

datq$SpCond_C <- c(datq$SpCond  - bg_SpCond)
datq$SpCond_C <-replace(datq$SpCond_C, datq$SpCond_C<0, 0)

#No Cl samples so Cl approx.
datq$Cl_mgL <- ((0.05/0.105)*datq$SpCond_C)

qplot(Cl_mgL, NO3_C, data = datq, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Carboy concentrations 300g in 10 carboy
NO3mgL <- 285 * (1000) * (62/101) *(1/9)
# Carboy concentrations 1500 NaCl in 6L carboy
NaClmgL <-  1000 * (1000) * (35.45/58.44) * (1/9)
carboy <- NO3mgL/NaClmgL

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
  NO3<- datq$NO3_mgNL[i]
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
out$Uadd <- Q*Cadd/out$sw*w


BW_uptake<- plot_grid(
  ggplot(out, aes(NO3, sw)) + geom_point(),
  ggplot(out, aes(NO3, Uadd)) + geom_point(), 
  ggplot(out, aes(datetime, log(NO3/Cl))) + geom_point(),
  ggplot(Hobo, aes(DateTime, SpCond)) + geom_point(),
  ncol=1, align="hv")
BW_uptake

# ggsave(plot = BW_uptake, filename = paste("./figures/BWL_NO3_220525.png",sep=""),width=4,height=7,dpi=300)

# write.csv(x = out, file = "./BTC_out/BWL_NO3_BTC_BWL230718.csv", row.names = TRUE)

# estimate N supply:
N_supp <-(86400*Q*(Cadd*0.001))/(w*reachL)
N_supp
mean(na.omit(out$sw))
mean(na.omit(out$Uadd))
mean(Hobo$TempC)
mean(na.omit(datq$PO4_ugL))
Nalt <- mean(na.omit(datq$Nh4_mgNL))
Nalt <- 0.005
mean(na.omit(datq$DOC_mgL))

N_supp_alt <-(86400*Q*(Nalt*0.001))/(w*reachL)

###########################