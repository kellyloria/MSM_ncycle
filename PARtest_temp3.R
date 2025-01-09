# load in packages:
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)
library(lme4)
library(lmerTest)


## Test of light on the near shore of 
# Blackwood creek on 2022-06-07 and in the FA cold room on 2022-06-08
# HOBO pendant: 21116494

# Depth profile order (deepest to shallowest at each site):
# BWNS2 -- 3m at 10:27 PST, 
# BW 10m -- 1m 11:30 PST; 5m at 11:25 PST; 9m at 11:20 PST
# BW 15m -- 1m 12:05 PST; 7m at 12:00 PST; 14m at 11:55 PST
# BW 20m -- 1m 12:25 PST; 10m at 12:30 PST; 19m at 12:35 PST

# Read in data 
Fdat <- read.csv("./PAR_dat/21116494_20220608_BWVanDornReadings.csv", header = T)

Fdat.ts <- Fdat %>% select(Date.Time..PDT., Ch..1...Temperature.....C., Ch..2...Light....lux.) %>% 
  rename(datetime = Date.Time..PDT., 
         temp.obs = Ch..1...Temperature.....C., 
         light.obs = Ch..2...Light....lux.) %>%
  mutate(datetime = as.POSIXct((datetime), format ="%m/%d/%Y %H:%M:%S"),
         time= format(as.POSIXct(datetime), format = "%H:%M:%S"))

# limit time for deployment time cutoff first and last days
Fdat.ts <- subset(Fdat.ts, 
                datetime > "2022-06-07 9:00:00" & 
                  datetime < "2022-06-07 16:00:00")

ggplot(Fdat.ts, aes(x=datetime, y=temp.obs)) + 
  geom_line() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(Fdat.ts, aes(x=datetime, y=light.obs)) + 
  geom_line() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

## surface lux:

surface <- subset(Fdat.ts, 
                  datetime > "2022-06-07 9:35:00" & 
                    datetime < "2022-06-07 9:47:00")

ggplot(surface, aes(x=datetime, y=light.obs)) + 
  geom_line() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 

## summary values:
sur_lux <- mean(surface$light.obs)
sur_temp <- mean(surface$temp.obs)

surface$location <- "surface"
surface$depth <- 0

#
# Look at BWNS2
BWNS2 <- subset(Fdat.ts, 
                  datetime > "2022-06-07 10:26:20" & 
                    datetime < "2022-06-07 10:27:50")

ggplot(BWNS2, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BWNS2, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

## 
BWNS2_lux <- mean(BWNS2$light.obs)
BWNS2_temp <- mean(BWNS2$temp.obs)

BWNS2$location <- "BWNS2"
BWNS2$depth <- 3


# Look at 10m
BW10m_9 <- subset(Fdat.ts, 
                datetime > "2022-06-07 11:26:50" & 
                  datetime < "2022-06-07 11:28:20")

ggplot(BW10m_9, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW10m_9, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW10m_9$location <- "10m"
BW10m_9$depth <- 9

## 
BW10m_9_lux <- mean(BW10m_9$light.obs)
BW10m_9_temp <- mean(BW10m_9$temp.obs)

###
# Look at 5m
BW10m_5 <- subset(Fdat.ts, 
                  datetime > "2022-06-07 11:31:20" & 
                    datetime < "2022-06-07 11:32:30")

ggplot(BW10m_5, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW10m_5, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 


BW10m_5$location <- "10m"
BW10m_5$depth <- 5

## 
BW10m_5_lux <- mean(BW10m_5$light.obs)
BW10m_5_temp <- mean(BW10m_5$temp.obs)



###
# Look at 1m
BW10m_1 <- subset(Fdat.ts, 
                  datetime > "2022-06-07 11:34:50" & 
                    datetime < "2022-06-07 11:36:20")

ggplot(BW10m_1, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW10m_1, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW10m_1$location <- "10m"
BW10m_1$depth <- 1

## 
BW10m_1_lux <- mean(BW10m_1$light.obs)
BW10m_1_temp <- mean(BW10m_1$temp.obs)


##
###
# Look at 15m
BW15m_14 <- subset(Fdat.ts, 
                  datetime > "2022-06-07 11:58:10" & 
                    datetime < "2022-06-07 11:59:30")

ggplot(BW15m_14, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW15m_14, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW15m_14$location <- "15m"
BW15m_14$depth <- 14

## 
BW15m_14_lux <- mean(BW15m_14$light.obs)
BW15m_14_temp <- mean(BW15m_14$temp.obs)




###
##
# Look at 15m
BW15m_7 <- subset(Fdat.ts, 
                   datetime > "2022-06-07 12:02:00" & 
                     datetime < "2022-06-07 12:03:30")

ggplot(BW15m_7, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW15m_7, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW15m_7$location <- "15m"
BW15m_7$depth <- 7

## 
BW15m_7_lux <- mean(BW15m_7$light.obs)
BW15m_7_temp <- mean(BW15m_7$temp.obs)


##
##
BW15m_1 <- subset(Fdat.ts, 
                  datetime > "2022-06-07 12:05:10" & 
                    datetime < "2022-06-07 12:06:40")

ggplot(BW15m_1, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW15m_1, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 


BW15m_1$location <- "15m"
BW15m_1$depth <- 1



## 
BW15m_1_lux <- mean(BW15m_1$light.obs)
BW15m_1_temp <- mean(BW15m_1$temp.obs)



###
###
##
##
BW20m_1 <- subset(Fdat.ts, 
                  datetime > "2022-06-07 12:25:00" & 
                    datetime < "2022-06-07 12:26:20")

ggplot(BW20m_1, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW20m_1, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW20m_1$location <- "20m"
BW20m_1$depth <- 1


## 
BW20m_1_lux <- mean(BW20m_1$light.obs)
BW20m_1_temp <- mean(BW20m_1$temp.obs)


##
##
BW20m_10 <- subset(Fdat.ts, 
                   datetime > "2022-06-07 12:28:00" & 
                     datetime < "2022-06-07 12:28:50")

ggplot(BW20m_10, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW20m_10, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW20m_10$location <- "20m"
BW20m_10$depth <- 10

## 
BW20m_10_lux <- mean(na.omit(BW20m_10$light.obs))
BW20m_10_temp <- mean(na.omit(BW20m_10$temp.obs))


##
##
BW20m_19 <- subset(Fdat.ts, 
                   datetime > "2022-06-07 12:32:40" & 
                     datetime < "2022-06-07 12:34:50")

ggplot(BW20m_19, aes(x=datetime, y=temp.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Temp C") 


ggplot(BW20m_19, aes(x=datetime, y=light.obs)) + 
  geom_point() +  
  theme_classic() + 
  labs(x ="Date", y = "Light (lux)") 

BW20m_19$location <- "20m"
BW20m_19$depth <- 19
## 
BW20m_19_lux <- mean(na.omit(BW20m_19$light.obs))
BW20m_19_temp <- mean(na.omit(BW20m_19$temp.obs))

# data frame:
BW_hobo <- rbind(surface,BWNS2,BW10m_9, BW10m_5, BW10m_1, BW15m_14, BW15m_7, BW15m_1, BW20m_19, BW20m_10, BW20m_1)


plot_BW <- ggplot(BW_hobo, # call the dataframe with data 
                  aes(x=depth, y=light.obs)) + # name your x and y axis
  geom_point(aes(x=depth, y=light.obs, color =location), alpha = 0.8, size=2) +
  theme_bw() +
  facet_grid(location~.)


plot_BW <- ggplot(BW_hobo, # call the dataframe with data 
                  aes(x=depth, y=light.obs)) + # name your x and y axis
  geom_point(aes(x=depth, y=light.obs, color =location), alpha = 0.8, size=2) +
  theme_bw() 

plot_BW <- ggplot(BW_hobo, # call the dataframe with data 
                  aes(x=depth, y=temp.obs)) + # name your x and y axis
  geom_point(aes(x=depth, y=temp.obs, color =location), alpha = 0.8, size=2) +
  theme_bw() +
  facet_grid(location~.)


#write.csv(x = BW_hobo, file = "./PAR_dat/BWNS202206profile.csv", row.names = TRUE)


### Lux level in cold room:

Ldat <- read.csv("./PAR_dat/21116494_20220620_lablighttest.csv", header = T)


Ldat.ts <- Ldat %>% select(Date.Time..PDT., Ch..1...Temperature.....C., Ch..2...Light....lux.) %>% 
  rename(datetime = Date.Time..PDT., 
         temp.obs = Ch..1...Temperature.....C., 
         light.obs = Ch..2...Light....lux.) %>%
  mutate(datetime = as.POSIXct((datetime), format ="%m/%d/%Y %H:%M:%S"),
         time= format(as.POSIXct(datetime), format = "%H:%M:%S"))

Ldat.ts <- subset(Ldat.ts, 
                  datetime > "2022-06-08 6:00:00" & 
                    datetime < "2022-06-08 18:00:00")


ggplot(Ldat.ts, aes(x=datetime, y=temp.obs)) + 
  geom_line() +  
  theme_classic()


ggplot(Ldat.ts, aes(x=datetime, y=light.obs)) + 
  geom_line() +  
  theme_classic() 


Ldat.ts$datetime <- round_date(Ldat.ts$datetime, "1 minutes")
Ldat.ts1<- Ldat.ts %>% mutate(
  time2= format(as.POSIXct(datime2), format = "%H:%M:%S"))


## fix notes:
NfixNote <- read.csv("./N_Fixation/NfixBW_streamNS220608.csv", header = T)

NfixNote <- NfixNote %>%
  mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S"))
str(NfixNote)
str(Ldat.ts)

NfixNote$datetime <- round_date(NfixNote$datetime, "1 minutes")

NfixNote.ts1 <- left_join(NfixNote, Ldat.ts [c("datetime", "temp.obs", "light.obs")],
                         by = c("datetime"="datetime"))

names(Ldat.ts)
names(NfixNote)


#write.csv(x = NfixNote.ts1, file = "./N_Fixation/NfixBW_streamNS220608_HB.csv", row.names = TRUE)

# mean lux 818.4 at BWNS2
Fluxmean <- mean(7887.1, 11376.6, 15008.0, 9465.8, 12028.0, 13768.6, 2408.8, 7544.4, 5705.4, 818.4)
Lluxmean <- mean(2090.24, 2088.96, 2088.96, 2092.16, 2093.44, 2093.44, 2090.24, 2087.04, 2090.24, 2092.16, 2093.44, 2092.16)
Fluxmean - Lluxmean

# PAR verse lux?
OddPAR <- read.csv("./PAR_dat/PARHOBOtest_West_F37459D7C4C7_1655765991820.csv", header = T)
OddPAR.ts <- OddPAR %>% select(dateTime, data2, data1) %>% 
  rename(datetime = dateTime, 
         temp.obs = data2, 
         PAR = data1) %>%
  mutate(datetime = as.POSIXct((datetime), format ="%d-%m-%Y %H:%M:%S"))



HBlux <- read.csv("./PAR_dat/PARHOBOtest21116494.csv", header = T)

HBlux.ts <- HBlux %>% select(Date.Time..PDT., Ch..1...Temperature.....C., Ch..2...Light....lux.) %>% 
  rename(datetime = Date.Time..PDT., 
         temp.obsH = Ch..1...Temperature.....C., 
         light.obsH = Ch..2...Light....lux.) %>%
  mutate(datetime = as.POSIXct((datetime), format ="%m/%d/%Y %H:%M:%S"))

HBlux.ts$datetime <- round_date(HBlux.ts$datetime, "1 minutes")



Lightest <- left_join(OddPAR.ts, HBlux.ts [c("datetime", "temp.obsH", "light.obsH")],
                          by = c("datetime"="datetime"))


Lightest.ts1 <- subset(Lightest, 
                  datetime > "2022-06-20 15:00:00" & 
                    datetime < "2022-06-20 15:50:00")


ggplot(Lightest.ts1, aes(x=PAR, y=light.obsH)) + 
  geom_point() +  
  theme_classic() 

ggplot(Lightest.ts1, aes(x=temp.obs , y=temp.obsH)) + 
  geom_point() +  
  theme_classic() 

paroffset <- glm(PAR~light.obsH, data=Lightest.ts1)
summary(paroffset)

plot(residuals(paroffset))

mean(Lightest.ts1$light.obsH-Lightest.ts1$PAR)
# PAR values are often 47242.41 lower than lux... 
