
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
# BWL20230718_NH4 "./NA22_dat/BWL_20220526/BWL20220526_NH4.csv"

Hobo <-read.csv("./NA22_dat/BWU_20230718/20775523_BOR.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
                "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523." )]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(Hobo$DateTime, format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 11:50:10") & DateTime <= as.POSIXct("2023-07-18 12:25:10"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 11:51:00") & DateTime <= as.POSIXct("2023-07-18 12:00:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 700) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-07-18 12:03:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-02-15 12:25:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(6.14, 4.18, 6.45, 4.43, 4.12, 4.78, 7.34, 9.81, 6.43, 10.67, 7.23, 4.88, 5.29, 6.01, 6.26))
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
### 2023-07-18 NO3 ###

# BWL20230718_NH4 "./NA22_dat/BWL_20220526/BWL20220526_NH4.csv"

Hobo <-read.csv("./NA22_dat/BWU_20230718/20775523_BOR.csv", skip=1)
summary(Hobo)
names(Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Hobo <- Hobo[,c("Date.Time..GMT.07.00",
                "Full.Range..μS.cm..LGR.S.N..20775523..SEN.S.N..20775523.",
                "Temp...C..LGR.S.N..20775523..SEN.S.N..20775523." )]

colnames(Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Hobo$DateTime <- as.POSIXct(Hobo$DateTime, format="%Y-%m-%dT%H:%M:%SZ") 
range(Hobo$DateTime)
str(Hobo)

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# 
Hobo$SpCond <- Hobo$Cond/(1-(25-Hobo$TempC)*0.021/100)

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 12:30:00") & DateTime <= as.POSIXct("2023-07-18 13:08:00"))

qplot(DateTime, Cond, data = Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

## Reach morphology estimates:
## (1) Determine the background conductivity
sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 12:30:00") & DateTime <= as.POSIXct("2023-07-18 12:40:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 700) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Hobo$DateTime), Hobo$SpCond, bg_SpCond, SpCond_mass)

inj_time <- as.POSIXct("2023-07-18 12:46:00") #Lolomai 
peak_time <- Hobo[which.max(Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2023-07-18 13:07:30")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(100) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(6.14, 4.18, 6.45, 4.43, 4.12, 4.78, 7.34, 9.81, 6.43, 10.67, 7.23, 4.88, 5.29, 6.01, 6.26))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


#######################
### 2023-08-10 NO3 ###

# BWL20230718_NH4 "./NA22_dat/BWL_20220526/BWL20220526_NH4.csv"

dat <- read.csv("./NA22_dat/BWU_20230810/BWU_20230810_NO3.csv")
dat$datetime <- as.POSIXct(paste(dat$date, dat$time), format = "%Y-%m-%d %H:%M:%S")
str(dat)

qplot(datetime, ysi_SPC, data = dat, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
Hobo <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 12:30:00") & DateTime <= as.POSIXct("2023-07-18 13:08:00"))

datq <- dat[c(1:16),]

## Reach morphology estimates:
## (1) Determine the background conductivity
#sub_bg <- subset(Hobo, DateTime >= as.POSIXct("2023-07-18 12:30:00") & DateTime <= as.POSIXct("2023-07-18 12:40:00")) #Lolomai
bg_SpCond <- c(33.5)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- c(2100* 1000) 
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(datq$datetime), datq$ysi_SPC, bg_SpCond, SpCond_mass)
Q
inj_time <- as.POSIXct("2023-08-10 13:43:40") #Lolomai 
peak_time <- datq[which.max(datq$ysi_SPC),]$datetime 
end_time <-as.POSIXct("2023-08-10 14:12:00")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL <- c(150) #
v <- c(reachL/time_diff_sec)
v

## Enter average width measurement in m
w <- mean(c(5.02, 6.53, 3.55, 5.13, 4.22, 6.83, 5.61, 4.70, 4.35, 4.20, 7.89, 5.39, 4.95, 4.20, 5.13))
w
## Calculate effective depth
z <- (Q/1000)/(w*v)
z
