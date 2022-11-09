##
## Nitrogen Uptake Assay
## ** NEED to double check morph and salt measurements before adding as a column to the dataset. 
##


## ---------------------------
# Blackwood 2021-06-09 

#######################
## Read in HOBO Data ##
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/BW_NAssay_20775520_1.csv", skip=1, header=T)

dat$Site <- "Blackwood"
dat$Seiral <- 20775520

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
              "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-06-09 12:30:00") & timestamp <= as.POSIXct("2021-06-09 13:36:00")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))

#write.csv(dat1, (paste0(outputDir,"/MSM21NA_Blackwood_20210609.csv")))


################
## Estimate Q ##
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-06-09 12:35:00") & timestamp <= as.POSIXct("2021-06-09 12:44:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*500 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-06-09 12:50:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 365.76/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
w <- (7.057292)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z
## Example Summary

####################################
## Read in sample collection data ##
####################################
dat1$time <- (format(as.POSIXlt(dat1$timestamp),format = '%T'))


sample <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/Blackwood20210609_SampleCollection.csv", skip=1, header=T)
sample <- sample[c(1:23),]
sample$date <- as.Date(sample$date, format="%m/%d/%y")


#sample$time <- as.POSIXct(as.character(sample$time), format="%H:%M:%S") ## modify the format to match your data
sample$timestamp <- as.POSIXct(paste(sample$date, sample$time), format="%Y-%m-%d %H:%M:%S")


summary(sample)

sample1 <- left_join(dat1, sample[c("time", "sample.no", "YSI.SPC")],
                       by = c("time" = "time"))
summary(sample1)

sample2 <- left_join(sample, dat1[c("timestamp", "Full.Range", "TempC", "SpCond")],
                     by = c("timestamp" = "timestamp"))
summary(sample2)


qplot(timestamp, SpCond, data = sample2, geom="point", ylab ="Flow", color = (Keep)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

qplot(timestamp, SpCond, data = sample1, geom="point", ylab ="Flow", color = factor(sample.no)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))









## ---------------------------
# Blackwood 2021-06-30 

#######################
## Read in HOBO Data ##
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/20775520_BlackwoodNAssay_20210630.csv", skip=1, header=T)

dat$Site <- "Blackwood"
dat$Seiral <- 20775520

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
              "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-06-30 12:44:00") & timestamp <= as.POSIXct("2021-06-30 13:32:30")) ## modify to match your data
ggplot(dat1, aes(timestamp, Full.Range))+geom_point()
ggplot(dat1, aes(timestamp, TempC))+geom_point()
range(na.omit(dat1$timestamp))

#write.csv(dat1, (paste0(outputDir,"/MSM21NA_Blackwood_20210630.csv")))


################
## Estimate Q ##
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-06-30 12:20:00") & timestamp <= as.POSIXct("2021-06-30 12:30:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*240 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-06-30 12:30:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 365.76/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
w <- 20
## Calculate effective depth
z <- (Q/4000)/(w*v)
z
## Example Summary

####################################
## Read in sample collection data ##
####################################
# *** NOT A CLEAR SIGNAL so discard all samples 




###
# Blackwood 2021-07-28 upper

#######################
## Read in HOBO Data ##
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/GeneralNAssay_UP20775520_0.csv", skip=1, header=T)
dat$Site <- "General"
dat$Seiral <- 20775520

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
              "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 10:51:00")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))

#write.csv(dat1, (paste0(outputDir,"/MSM21NA_GeneralUP_20210715.csv")))


################
## Estimate Q ##
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 09:55:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*2000 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-15 10:00:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 50/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
wft <- mean(20.5,10.5,16,9,9.5,18,16.5,10.5)
w<- 4.418729143
## Calculate effective depth
z <- (Q/1000)/(w*v)
z
###


























## ---------------------------
# General 2021-07-15 upper

#######################
## Read in HOBO Data ##
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/GeneralNAssay_UP20775520_0.csv", skip=1, header=T)
dat$Site <- "General"
dat$Seiral <- 20775520

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
              "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 10:51:00")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))

#write.csv(dat1, (paste0(outputDir,"/MSM21NA_GeneralUP_20210715.csv")))


################
## Estimate Q ##
################
## Equation
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 09:55:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*2000 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-15 10:00:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 50/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
w <- mean(20.5,10.5,16,9,9.5,18,16.5,10.5)
## Calculate effective depth
z <- (Q/50)/(w*v)
z





## ---------------------------
# General 2021-07-15 lower

#######################
## Read in HOBO Data ##
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/GeneralNAssay_LOW20775520_0.csv", skip=1, header=T)
dat$Site <- "General"
dat$Seiral <- 20775515

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775515..SEN.S.N..20775515.",
              "Temp...C..LGR.S.N..20775515..SEN.S.N..20775515.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%m/%d/%y %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:55:00") & timestamp <= as.POSIXct("2021-07-15 10:50:00")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))

#write.csv(dat1, (paste0(outputDir,"/MSM21NA_GeneralLow_20210715.csv")))


################
## Estimate Q ##
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 09:55:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*2000 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-15 10:00:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 100/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
wft <- mean(19.5,9.5,15,8,8.5,17,15.5,9.5) 

w<-3.74364

## Calculate effective depth
z <- (Q/1000)/(w*v)
z


####################################
## Read in sample collection data ##
####################################
sample <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/General_20210715_samplecollection.csv", skip=1, header=T)
sample <- sample[c(1:18),]
sample$date <- as.Date(sample$date, format="%m/%d/%y")

sample$timestamp <- as.POSIXct(paste(sample$date, sample$time), format="%Y-%m-%d %H:%M:%S")
summary(sample)

sample2 <- left_join(sample, dat1[c("timestamp", "Full.Range", "TempC", "SpCond")],
                     by = c("timestamp" = "timestamp"))
summary(sample2)

qplot(timestamp, SpCond, data = sample2, geom="point", ylab ="Flow", color = (Keep)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))



## ---------------------------
# Glenbrook 2021-06-03 lower

#######################
## Read in HOBO Data ##
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/Glenbrook20210603_Nuptake.csv",skip=1, header=T)
dat$Site <- "Glenbrook"
dat$Seiral <- 20775509

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00",
              "Full.Range..μS.cm..LGR.S.N..20775509..SEN.S.N..20775509.",
              "Temp...C..LGR.S.N..20775509..SEN.S.N..20775509.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-06-03 12:00:00") & timestamp <= as.POSIXct("2021-06-03 12:59:59")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat1$timestamp))

#write.csv(dat1, (paste0(outputDir,"/MSM21NA_Glenbrook_20210603.csv")))


################
## Estimate Q ## ** CHECK SALT AMOUNTS
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-06-03 12:00:00") & timestamp <= as.POSIXct("2021-06-03 12:05:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*200 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-06-03 12:10:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 365.76/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
w <- mean(20.5,10.5,16,9,9.5,18,16.5,10.5)
## Calculate effective depth
z <- (Q/375)/(w*v)
z


####################################
## Read in sample collection data ##
####################################
sample <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/Glenbrook20210603_samplecollection.csv", skip=1, header=T)
#sample <- sample[c(1:18),]
sample$date <- as.Date(sample$date, format="%m/%d/%y")

sample$timestamp <- as.POSIXct(paste(sample$date, sample$time), format="%Y-%m-%d %H:%M:%S")
summary(sample)

sample2 <- left_join(sample, dat1[c("timestamp", "Full.Range", "TempC", "SpCond")],
                     by = c("timestamp" = "timestamp"))
summary(sample2)

qplot(timestamp, SpCond, data = sample2, geom="point", ylab ="Flow", color = (Keep)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

######



## ---------------------------
# Glenbrook 2021-06-23 lower

#######################
## Read in HOBO Data ## # Unfinished stream analysis decided to keep samples 
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/Glenbrook20210623_NuptakeREAL.csv", skip=1, header=T)
dat$Site <- "Glenbrook"
dat$Seiral <- 20775509

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00" ,
              "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
              "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-06-23 10:08:01") & timestamp <= as.POSIXct("2021-06-23 12:40:00")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))


dat2 <- subset(dat, timestamp >= as.POSIXct("2021-06-23 11:40:00") & timestamp <= as.POSIXct("2021-06-23 12:40:00")) ## modify to match your data
ggplot(dat2, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))


dat3 <- subset(dat, timestamp >= as.POSIXct("2021-06-23 12:55:00") & timestamp <= as.POSIXct("2021-06-23 14:25:00")) ## modify to match your data
ggplot(dat3, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat1$timestamp))



#write.csv(dat1, (paste0(outputDir,"/MSM21NA_GeneralLow_20210715.csv")))


################
## Estimate Q ##
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 09:55:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*2000 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-15 10:00:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 365.76/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
w <- mean(20.5,10.5,16,9,9.5,18,16.5,10.5)
## Calculate effective depth
z <- (Q/100)/(w*v)
z


####################################
## Read in sample collection data ##
####################################
sample <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/General_20210715_samplecollection.csv", skip=1, header=T)
sample <- sample[c(1:18),]
sample$date <- as.Date(sample$date, format="%m/%d/%y")

#sample$time <- as.POSIXct(as.character(sample$time), format="%H:%M:%S") ## modify the format to match your data
sample$timestamp <- as.POSIXct(paste(sample$date, sample$time), format="%Y-%m-%d %H:%M:%S")
summary(sample)

sample2 <- left_join(sample, dat1[c("timestamp", "Full.Range", "TempC", "SpCond")],
                     by = c("timestamp" = "timestamp"))
summary(sample2)

qplot(timestamp, SpCond, data = sample2, geom="point", ylab ="Flow", color = (Keep)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

#
#
#



## ---------------------------
# Incline 2021-06-23 lower

#######################
## Read in HOBO Data ## # Unfinished stream analysis decided to keep samples 
#######################

# Import
dat <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/Icline/N_assay/20210707/20775515_InclineNassay_0.csv", skip=1, header=T)
dat$Site <- "Incline"
dat$Seiral <- 20775515

# Subset
# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
dat <- dat[,c("Site",
              "Seiral",
              "Date.Time..GMT.07.00" ,
              "Full.Range..μS.cm..LGR.S.N..20775515..SEN.S.N..20775515.",
              "Temp...C..LGR.S.N..20775515..SEN.S.N..20775515.")]

colnames(dat) <- c("Site","Seiral","timestamp","Full.Range","TempC")
# Convert DateTime
dat$timestamp <- as.POSIXct(as.character(dat$timestamp), format="%y/%m/%d %H:%M:%S") ## modify the format to match your data
# Convert electrical conductance to specific conductance if not already automated
dat$SpCond <- dat$Full.Range/(1-(25-dat$TempC)*0.021/100)
# Visualize
ggplot(dat, aes(timestamp, SpCond))+geom_point()
# Subset based on plot to the slug window
# Pull out slug window:
dat1 <- subset(dat, timestamp >= as.POSIXct("2021-07-07 10:10:00") & timestamp <= as.POSIXct("2021-07-07 10:57:00")) ## modify to match your data
ggplot(dat1, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))


dat2 <- subset(dat, timestamp >= as.POSIXct("2021-07-07 11:40:00") & timestamp <= as.POSIXct("2021-07-07 12:33:00")) ## modify to match your data
ggplot(dat2, aes(timestamp, SpCond))+geom_point()
range(na.omit(dat$timestamp))


#write.csv(dat2, (paste0(outputDir,"/MSM21NA_InclineUP_20210707.csv")))


################
## Estimate Q ##
################
## Equation
Qint<-function(time, cond, bkg, condmass){
  condcorr<-cond-bkg
  ##below routine integrates
  ydiff<- condcorr[-1]+ condcorr[-length(condcorr)]
  condint<-sum(diff(time)*ydiff/2)
  Q<-condmass/condint
  Q }
## (1) Determine the background conductivity
## Select area before or after the salt wave which is constant
## for at least 30 minutes and take the average
sub_bg <- subset(dat, timestamp >= as.POSIXct("2021-07-15 09:40:00") & timestamp <= as.POSIXct("2021-07-15 09:55:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
## 1 g salt in 1 L of water gives cond=2100 uS / cm
# Oak @ Lolomai: 500 g 
SpCond_mass <- 2100*2000 # 500 g NaCl
## Calculate Q!
## Units = L/sec
Q <- Qint(as.numeric(dat1$timestamp), dat1$SpCond, bg_SpCond, SpCond_mass)
## Summary
## Oak @ Lolomai: Q = 491.2901 L/sec (0.49 cms)


#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-15 10:00:00") #Lolomai 
peak_time <- dat1[which.max(dat1$SpCond),]$timestamp 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds # Oak @ Lolomai: 365.76 m upstream
v <- 365.76/time_diff_sec
v
## Example Summary
## Oak @ Lolomai: 0.346 m/s


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
## Enter average width measurement in m # NEED TO DOUBLE CHECK
w <- mean(20.5,10.5,16,9,9.5,18,16.5,10.5)
## Calculate effective depth
z <- (Q/100)/(w*v)
z


####################################
## Read in sample collection data ##
####################################
sample <- read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/Conduct_downloads/AssayHOBO/General_20210715_samplecollection.csv", skip=1, header=T)
sample <- sample[c(1:18),]
sample$date <- as.Date(sample$date, format="%m/%d/%y")

#sample$time <- as.POSIXct(as.character(sample$time), format="%H:%M:%S") ## modify the format to match your data
sample$timestamp <- as.POSIXct(paste(sample$date, sample$time), format="%Y-%m-%d %H:%M:%S")
summary(sample)

sample2 <- left_join(sample, dat1[c("timestamp", "Full.Range", "TempC", "SpCond")],
                     by = c("timestamp" = "timestamp"))
summary(sample2)

qplot(timestamp, SpCond, data = sample2, geom="point", ylab ="Flow", color = (Keep)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))



