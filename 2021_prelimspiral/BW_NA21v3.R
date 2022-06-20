library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)
library(lme4)
library(lmerTest)

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR/SummerResearch2021')){
  inputDir<- '/Users/kellyloria/Documents/UNR/SummerResearch2021/'
  outputDir<- '/Users/kellyloria/Documents/UNR/SummerResearch2021/DO_downloads/' 
}
## ---------------------------

# Read in the nitrogen uptake assay data:
BW_NA <- read.csv(paste0(inputDir,("/NA21/NA21_BW20210728.csv")))
BW_NA$datetime <- as.POSIXct(as.character(BW_NA$datetime), format="%Y-%m-%d %H:%M:%S") ## modify the format to match your data
# 

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))


## SPC data from HOBO
BW_Hobo <-read.csv("/Users/kellyloria/Documents/UNR/SummerResearch2021/NA21/GeneralandBW_20210727_20775520_4.csv", skip=1)
summary(BW_Hobo)

# modify the names to whatever names your sensor spits out # figure out the names after import by using names(dat) 
BW_Hobo <- BW_Hobo[,c("Date.Time..GMT.07.00",
                      "Full.Range..μS.cm..LGR.S.N..20775520..SEN.S.N..20775520.",
                      "Temp...C..LGR.S.N..20775520..SEN.S.N..20775520.")]

colnames(BW_Hobo) <- c("DateTime","Cond","TempC")
# Convert DateTime
BW_Hobo$DateTime <- as.POSIXct(as.character(BW_Hobo$DateTime), format="%y/%m/%d %H:%M:%S") 
range(BW_Hobo$DateTime)

# 
BW_Hobo$SpCond <- BW_Hobo$Cond/(1-(25-BW_Hobo$TempC)*0.021/100)

qplot(DateTime, Cond, data = BW_Hobo, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))

# Adjust the time range:
BW_Hobo <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 15:30:00") & DateTime <= as.POSIXct("2021-07-28 17:01:50"))

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
sub_bg <- subset(BW_Hobo, DateTime >= as.POSIXct("2021-07-28 15:31:00") & DateTime <= as.POSIXct("2021-07-28 15:39:00")) #Lolomai
bg_SpCond <- mean(sub_bg$SpCond)
## (2) Estimate conductivity slug based on mass of Cl added
SpCond_mass <- 2100*2000
## Calculate Q
## Units = L/sec
Q <- Qint(as.numeric(BW_Hobo$DateTime), BW_Hobo$SpCond, bg_SpCond, SpCond_mass)

#######################
## Estimate Velocity ##
#######################
inj_time <- as.POSIXct("2021-07-28 15:40:00") #Lolomai 
peak_time <- BW_Hobo[which.max(BW_Hobo$SpCond),]$DateTime 
time_diff_sec <- as.numeric(peak_time - inj_time)*60
## Velocity = distance in meters/time in seconds
v <- 250/time_diff_sec
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
w <- mean(5.49, 5.12, 3.98, 4.79, 6.1, 
          4.32, 4.45, 3.78, 6.02, 5.33)
## Calculate effective depth
z <- (Q/1000)/(w*v)
z


##
## Left join hobo conductivity estimates with samples:
##

BW_NA <- left_join(BW_NA, BW_Hobo[c("DateTime", "SpCond")],
                   by= c("datetime"="DateTime"))

summary(BW_NA)

qplot(datetime, Results, data = BW_NA, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"), 
                   breaks = date_breaks("15 min"))

##
## Calculations Notes:
##

# So lets start:
### 1. covert timestamp to seconds:
df <- BW_NA %>%
  arrange(Sample_ID, datetime) %>%
  group_by(Sample_ID) 
BW_NA$diff <- rep(0,nrow(BW_NA))
for(i in 1:(nrow(BW_NA)-1)){
  BW_NA$diff[i+1] <- BW_NA$datetime[i+1]- BW_NA$datetime[i]
}

df$diff_sec <- c(df$diff*60)

# want to add seconds...?
df$diff_sec_sum <- (Reduce("+", df$diff_sec, accumulate = TRUE))


### 2. Correct for background concentrations:
# correct for cl:
df$SPC <- c(df$SpCond-(bg_SpCond))
df$SPC[df$SPC <0]<-0
df$Cl_mg <- ((0.05/0.105)*df$SPC)

# correct for NO3:
df$NO3_mg <- c(df$Results-(0.005))

df$N_Cl <- df$NO3_mg/df$Cl_mg
df$N_Cl[df$N_Cl==Inf]<-0
df$N_Cl[df$N_Cl=NaN]<-0



# To integrate the breakthrough curve (BTC)
df$BCT <- rep(0,nrow(df))
for(i in 1:(nrow(df)-1)){
  df$BCT[i+1] <- (df$N_Cl[i+1]*df$diff_sec_sum[i+1])- (df$N_Cl[i]*df$diff_sec_sum[i])
}

# reach lenght = 250 m
L <- 250

#TMR?
df$TMR <- (df$BCT*Q) 

# inject concentration:
slug_mgCl <- ((2000*1000)) # weighed in grams
slug_mgN <- ((250*1000)) 
slug_mgN/slug_mgCl * Q

# inputs for kw plots?
inj<- log((slug_mgN/slug_mgCl*Q))
grab <-log(df$TMR) # TMR still looks higher than inj... should be lower... 

grab[grab==-Inf]<-0
grab<- replace(grab, is.na(grab), 0)


## Remove blank and head of reach .. so first 2 rows
grab<- grab[-c(1)]

BW_NA$L <- c(L)

#new df
distance <- c(0, 
              BW_NA$L[-c(1)])

TMR_con <- c(inj,
             grab)

sample_ID <- c("AAA",
               BW_NA$Sample_ID[-c(1)])

newdf <- data_frame(distance, TMR_con, sample_ID)

kw_glm <- glm(TMR_con~scale(distance)*sample_ID, data=newdf)
summary(kw_glm)

(summary(kw_glm)$coefficients)[2:13,]

# make new df for uptake summary valus
BW_uptake <- df[-c(1), c("Sample_ID","datetime","diff_sec","diff_sec_sum","Cl_mg","NO3_mg","N_Cl","TMR")]


BW_uptake$kw = c(-1.665239,
                 -7.869793,
                 -7.869793,
                 -1.063300,
                 -1.879801,
                 -7.869793,
                 -2.786330,
                 -7.869793,
                 -7.869793,
                 -1.146898,
                 -1.644562,
                 -1.350970)

# sw: Added nutrient uptake length calculated with the slug BTC-integrated approach (m) 
BW_uptake$sw <- -1/BW_uptake$kw 

# Uadd-int: Added nutrient areal uptake rate calculated with the slug BTC-integrated approach (mg m^-2 sec^-1)
BW_uptake$Uadd = ((Q * (BW_uptake$NO3_mg))/ (BW_uptake$sw * (w)))

plot(BW_uptake$NO3_mg, BW_uptake$Uadd)

# Vf-add-int: Added nutrient uptake velocity calculated with the dynamic TASCC approach (m sec^-1)
BW_uptake$Vf = BW_uptake$Uadd/ (BW_uptake$NO3_mg)
plot(BW_uptake$NO3_mg, BW_uptake$Vf)

# library(ggplot2)
# library(ggformula) # optional
# library(MASS)
# library(drc)
# library(gridExtra)
# 
# 
# BW_NA$site <- "Blackwood"
# BW_plot <- BW_NA[,c("datetime",
#                     "SpCond_C",
#                     "Cl_ugL",
#                     "NO3_C",
#                     "TMR_rat",
#                     "kw",
#                     "sw",
#                     "Uadd",
#                     "Vfaddint",
#                     "site")]
# 
# 
# #### Nutrient uptake rate for GB
# model.drm1 <- drm (Uadd ~ NO3_C, data = BW_plot, fct = MM.2())
# summary(model.drm1)
# 
# mm2 <- data.frame(NO3_C = seq(0, max(BW_plot$NO3_C), length.out = 100))
# mm2$Uadd <- predict(model.drm1, newdata = mm2)
# 
# summary(model.drm1)
# (model.drm1)
# 
# 
# Uadd_plotBW <- ggplot(BW_plot, aes(x=NO3_C, y=Uadd)) + 
#   #ylim(0,20) + xlim(0,1500) + 
#   geom_line(data = mm2, aes(x = NO3_C, y = Uadd), colour = "black") +
#   geom_point(size = 3, shape= 17, col = "#a67d17") +
#   labs(x = expression(paste("Nitrate ugL")),
#        y= expression(paste("Uadd (", mu, "g m^-2 min^-1)"))) +
#   theme_classic() #+ annotate("text", x = c(45,40,40), y = c(20, 18, 16), label = c("Umax = 954.43", "Km = 24.21", "p = 0.002 "))
# 
# max(BW_plot$sw)
# max(na.omit(BW_plot$Vfaddint))
# max(BW_plot$Uadd)
