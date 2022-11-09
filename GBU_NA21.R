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
GB_NA <- read.csv("./NA21_dat/samples/GBL20210623sample.csv")
GB_NA$datetime <- as.POSIXct(paste(GB_NA$date, GB_NA$time), format = "%m/%d/%y %H:%M:%S")
  
qplot(datetime, NO3, data = GB_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


# ## SPC data from HOBO
# GB_Hobo <-read.csv("./NA21_dat/GlenbrookNA20210722_20775520_4.csv", skip=1)
# summary(GB_Hobo)
# 
# 
# # modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
# GB_Hobo <- GB_Hobo[,c("Date.Time..GMT.07.00",
#                       "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
#                       "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]
# 
# colnames(GB_Hobo) <- c("DateTime","Cond","TempC")
# # Convert DateTime
# GB_Hobo$DateTime <- as.POSIXct(as.character(GB_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
# range(GB_Hobo$DateTime)
# 
# # 
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
## GBU 2021-06-23