# Dillution gaging notes:
InclineSS<- read.csv("./NA21_dat/DillutionEst/InclineNassay20775515_0.csv", skip=1)
names(InclineSS)

Inc_Hobo# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Inc_Hobo <- InclineSS[,c("Date.Time..GMT.07.00",
                         "Full.Range..μS.cm..LGR.S.N..20775515..SEN.S.N..20775515.",
                         "Temp...C..LGR.S.N..20775515..SEN.S.N..20775515.")]

colnames(Inc_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Inc_Hobo$DateTime <- as.POSIXct(as.character(Inc_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
range(Inc_Hobo$DateTime)

# 
Inc_Hobo$SpCond <- Inc_Hobo$Cond/(1-(25-Inc_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = Inc_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
Inc_Hobo <- subset(Inc_Hobo, DateTime >= as.POSIXct('2021-08-05 08:31:00') & DateTime <= as.POSIXct('2021-08-05 10:0:00 '))

qplot(DateTime, Cond, data = Inc_Hobo, geom="point") +
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
sub_bg <- subset(Inc_Hobo, DateTime >= as.POSIXct('2021-08-05 08:31:00') & DateTime <= as.POSIXct('2021-08-05 09:00:00')) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Inc_Hobo$DateTime), Inc_Hobo$SpCond, bg_SpCond, SpCond_mass)
Q

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-08-05 09:10:00") #Lolomai 
peak_time <- Inc_Hobo[which.max(Inc_Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2021-03-26 13:31:51")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL<-c(100)
v <- c(reachL/time_diff_sec)
v
## Enter average width measurement in m
w <- 1.905
w
## Calculate effective depth in m
z <- c(Q/1000)/(w*v)
z


# Dillution gaging notes:
GBSS<- read.csv("./NA21_dat/DillutionEst/GlenbrookNutAssay1_20775509.csv", header=T)
names(GBSS)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
GB_Hobo <- GBSS[,c("Date.Time",
                         "Low.Range",
                         "Temp")]

colnames(GB_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
GB_Hobo$DateTime <- as.POSIXct(as.character(GB_Hobo$DateTime), format="%m/%d/%y %H:%M") 
range(GB_Hobo$DateTime)

# 
GB_Hobo$SpCond <- GB_Hobo$Cond/(1-(25-GB_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = GB_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))


## Reach morphology estimates:
##

## Equation:
## (1) Determine the background conductivity
sub_bg <- subset(GB_Hobo, DateTime >= as.POSIXct('2021-06-03 12:00:00') & DateTime <= as.POSIXct('2021-06-03 12:10:00')) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Inc_Hobo$DateTime), Inc_Hobo$SpCond, bg_SpCond, SpCond_mass)
Q

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-03-26 13:25:00") #Lolomai 
peak_time <- Inc_Hobo[which.max(Inc_Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2021-03-26 13:31:51")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL<-100
v <- reachL/time_diff_sec
v
## Enter average width measurement in m
w <- 1.905
w
## Calculate effective depth in m
z <- (Q/1000)/(w*v)
z



# Dillution gaging notes:
InclineSS<- read.csv("./NA21_dat/DillutionEst/InclineNassay20775515_0.csv", skip=1)
names(InclineSS)

Inc_Hobo# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
Inc_Hobo <- InclineSS[,c("Date.Time..GMT.07.00",
                         "Full.Range..μS.cm..LGR.S.N..20775515..SEN.S.N..20775515.",
                         "Temp...C..LGR.S.N..20775515..SEN.S.N..20775515.")]

colnames(Inc_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
Inc_Hobo$DateTime <- as.POSIXct(as.character(Inc_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
range(Inc_Hobo$DateTime)

# 
Inc_Hobo$SpCond <- Inc_Hobo$Cond/(1-(25-Inc_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = Inc_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
Inc_Hobo <- subset(Inc_Hobo, DateTime >= as.POSIXct('2021-08-05 08:31:00') & DateTime <= as.POSIXct('2021-08-05 10:0:00 '))

qplot(DateTime, Cond, data = Inc_Hobo, geom="point") +
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
sub_bg <- subset(Inc_Hobo, DateTime >= as.POSIXct('2021-08-05 08:31:00') & DateTime <= as.POSIXct('2021-08-05 09:00:00')) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(Inc_Hobo$DateTime), Inc_Hobo$SpCond, bg_SpCond, SpCond_mass)
Q

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-03-26 13:25:00") #Lolomai 
peak_time <- Inc_Hobo[which.max(Inc_Hobo$SpCond),]$DateTime 
end_time <-as.POSIXct("2021-03-26 13:31:51")
time_diff_sec <- as.numeric(peak_time - inj_time)*60
time_tota_sec <- (as.numeric(end_time - inj_time)) * 3600 # minutes

## Velocity = distance in meters/time in seconds
reachL<-100
v <- reachL/time_diff_sec
v
## Enter average width measurement in m
w <- 1.905
w
## Calculate effective depth in m
z <- (Q/1000)/(w*v)
z
