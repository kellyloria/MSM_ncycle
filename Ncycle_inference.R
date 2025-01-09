lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2", "ggpubr"), require, character.only=T)
site_colors <- c(
  "BWL" = "#3283a8",
  "BWU" = "#3258a8",
  "GBL" = "#a67d17",
  "GBU" = "#a65d17"
)
# Create a sequence of dates
date_seq <- seq(from = as.Date("2021-03-20"), to = as.Date("2024-10-10"), by = "day")

# Convert the sequence to a dataframe
date_df <- data.frame(date = date_seq)
date_df$site <- "BWL"
date_df1 <- data.frame(date = date_seq)
date_df1$site <- "GBL"
date_df2 <- data.frame(date = date_seq)
date_df2$site <- "GBU"
date_df3 <- data.frame(date = date_seq)
date_df3$site <- "BWU"

date_datf <- rbind(date_df,date_df1,date_df2, date_df3)

water_year <- function(data) {
  # Ensure the `date` column exists
  if (!"date" %in% names(data)) {
    stop("The dataframe must contain a column named 'date'.")
  }
  
  # Convert the `date` column to Date format (if not already)
  data$date <- as.Date(data$date)
  
  # Add the water_year column
  data$water_year <- with(data, {
    year <- as.numeric(format(date, "%Y"))
    month <- as.numeric(format(date, "%m"))
    ifelse(month >= 10, year + 1, year)
  })
  
  return(data)
}


date_datf <- water_year(date_datf)

covariat_dat <- readRDS("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /data/CH1_covariate_dat.rds") %>%
  dplyr::select("date", "Site","bulk.density","AFDM_mgg","AFDM_mgcm2", "Chla_ugL_Q", "Pheo_ugL_Q")
summary(covariat_dat)

covariat_dat <- date_datf %>%
  left_join(covariat_dat, by=c("date", "site"="Site"))

bg_nuts <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/NS_chem_dat_nh4_24.rds") %>%
  #mutate(date = as.Date(date, format="%m/%d/%y")) %>%
  filter(location=="stream")

WQ_dat <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/Stream_Lake_YSI_WaterQuality.csv")%>%
  mutate(date = as.Date(date, format="%m/%d/%y")) %>%
  dplyr::select("site", "date","pH")
  

# bg_nuts <- read.csv("/Users/kellyloria/Downloads/NS_chem_dat_24.csv")%>%
#   mutate(date = as.Date(date, format="%m/%d/%y")) %>%
#   filter(location=="stream")

covariat_datq <- covariat_dat%>%
  left_join(bg_nuts[ , c("site", "substrate", "date", "NO3_mgL_dl", "NH4_mgL_dl", "PO4_ugL_dl", "DOC_mgL_dl")], by=c("date", "site"))

covariat_datq <- covariat_datq%>%
  left_join(WQ_dat,  by=c("date", "site"))


## bring in SPC data 
### missing 2024 BW SPC
SPC_BWL <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWL_SPCv2.rds") %>%
  filter(datetime>as.POSIXct("2021-04-29 24:00:00"))
SPC_BWL$site <- "BWL"
SPC_BWU <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWU_SPCv2.rds")
SPC_BWU$site <- "BWU"
SPC_GBL <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_GBL_SPCv2.rds")
SPC_GBL$site <- "GBL"
SPC_GBU <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_GBU_SPCv2.rds")
SPC_GBU$site <- "GBU"

SPC_dat <- rbind(SPC_BWL, SPC_BWU, SPC_GBL, SPC_GBU)

SPC_datD <- SPC_dat%>%
  mutate(date =as.Date(datetime))%>%
  group_by(site, date) %>%
  summarise(
    wt= mean(wtr, na.rm=T),
    wt_sd= sd(wtr, na.rm=T),
    SPC_m=mean(SPC, na.rm=T),
    SPC_sd=sd(SPC, na.rm=T))

hist(SPC_datD$wt_sd)
hist(SPC_datD$SPC_sd)

SPC_datD1 <- SPC_datD %>%
  filter(wt_sd<4.1 & SPC_sd <35)
  
hist(SPC_datD1$wt_sd)
hist(SPC_datD1$SPC_sd)
  
covariat_datq <- covariat_datq%>%
  left_join(SPC_datD1, by=c("date", "site"))


### 
## Bring in WQ data : 
BWL_DO <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/DO_dat/24_BWL_DO_flag_sat_light.rds")

BWL_DO_D <- BWL_DO%>%
  mutate(date =as.Date(datetime))%>%
  group_by(date) %>%
  summarise(
    wtemp= mean(wtr, na.rm=T),
    wtemp_sd= sd(wtr, na.rm=T),
    do_sat_m=mean(DO.sat, na.rm=T),
    do_sat_sd=sd(DO.sat, na.rm=T),
    Q_m=mean(dischargeCMS, na.rm=T),
    Q_sd=sd(dischargeCMS, na.rm=T))

BWL_DO_D$site <- "BWL"

BWU_DO <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/DO_dat/24_BWU_DO_flag_sat_light.rds")

BWU_DO_D <- BWU_DO%>%
  mutate(date =as.Date(datetime))%>%
  group_by(date) %>%
  summarise(
    wtemp= mean(wtr, na.rm=T),
    wtemp_sd= sd(wtr, na.rm=T),
    do_sat_m=mean(DO.sat, na.rm=T),
    do_sat_sd=sd(DO.sat, na.rm=T),
    Q_m=mean(dischargeCMS, na.rm=T),
    Q_sd=sd(dischargeCMS, na.rm=T))
BWU_DO_D$site <- "BWU"

GBL_DO <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/DO_dat/24_GBL_DO_flag_sat_light.rds")
GBL_DO_D <- GBL_DO%>%
  mutate(date =as.Date(datetime))%>%
  group_by(date) %>%
  summarise(
    wtemp= mean(wtr, na.rm=T),
    wtemp_sd= sd(wtr, na.rm=T),
    do_sat_m=mean(DO.sat, na.rm=T),
    do_sat_sd=sd(DO.sat, na.rm=T),
    Q_m=mean(dischargeCMS, na.rm=T),
    Q_sd=sd(dischargeCMS, na.rm=T))
GBL_DO_D$site <- "GBL"


GBU_DO <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/DO_dat/24_GBU_DO_flag_sat.rds")

GBU_DO_D <- GBU_DO%>%
  mutate(date =as.Date(datetime))%>%
  group_by(date) %>%
  summarise(
    wtemp= mean(wtr, na.rm=T),
    wtemp_sd= sd(wtr, na.rm=T),
    do_sat_m=mean(DO.sat, na.rm=T),
    do_sat_sd=sd(DO.sat, na.rm=T),
    Q_m=mean(dischargeCMS, na.rm=T),
    Q_sd=sd(dischargeCMS, na.rm=T))
GBU_DO_D$site <- "GBU"

DO_dat <- rbind(GBU_DO_D, GBL_DO_D, BWU_DO_D, BWL_DO_D)

### N-uptake metrics
n_up <- read.csv("/Users/kellyloria/Documents/UNR/Ncycle/BTC_N_Table_2024.csv")%>%
  mutate(date = as.Date(date, format="%m/%d/%y")) 

covariat_datq <- covariat_datq%>%
  left_join(n_up, by=c("date", "site"))

covariat_datq <- covariat_datq%>%
  left_join(DO_dat, by=c("date", "site"))


hist(covariat_datq$AFDM_mgcm2)
hist(covariat_datq$Chla_ugL_Q)
hist(covariat_datq$AFDM_mgg)

covariat_datq$year <- year(covariat_datq$date)
covariat_datq$yday <- yday(covariat_datq$date)
covariat_datq$month <- month(covariat_datq$date)


##################
##################
##################
hist(covariat_datq$AFDM_mgg)
hist(log(covariat_datq$AFDM_mgg)+1)
ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtemp, y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) + 
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

## biomass
hist(covariat_datq$AFDM_mgcm2)
hist(log(covariat_datq$AFDM_mgcm2+1))

ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtemp, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)+
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=NH4_mgL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=PO4_ugL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=DOC_mgL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

####
# chla
## biomass
hist(covariat_datq$Chla_ugL_Q)
hist(log(covariat_datq$Chla_ugL_Q+1))

ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtemp, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)+
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=NH4_mgL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=PO4_ugL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=DOC_mgL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

###

####
# chla
## biomass

covariat_datq$biom_chla <- (covariat_datq$Chla_ugL_Q)/(covariat_datq$AFDM_mgcm2)
hist(covariat_datq$biom_chla)
hist(log(covariat_datq$biom_chla+1))

ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtemp, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)+
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=NH4_mgL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=PO4_ugL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=DOC_mgL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

####################
## explore monthly values
WQ_month_dat <- covariat_datq%>%
  group_by(site, month, water_year) %>%
  summarise(
    bulk_density_m = mean(bulk.density, na.rm=T),
    bulk_density_sd= sd(bulk.density, na.rm=T),
    chla_m = mean(Chla_ugL_Q, na.rm=T),
    chla_sd= sd(Chla_ugL_Q, na.rm=T),
    OM_mgg_m=mean(AFDM_mgg, na.rm=T),
    OM_mgg_sd=sd(AFDM_mgg, na.rm=T),
    biomass_mgcm2_m=mean(AFDM_mgcm2, na.rm=T),
    biomass_mgcm2_sd=sd(AFDM_mgcm2, na.rm=T),
    wt_m=mean(wt, na.rm=T),
    wt_sd=sd(wt, na.rm=T),
    wtemp_m=mean(wtemp, na.rm=T),
    wtemp_sd=sd(wtemp, na.rm=T),
    Q_m_m=mean(Q_m, na.rm=T),
    Q_sd=sd(Q_m, na.rm=T),
    do_sat_m_m=mean(do_sat_m, na.rm=T),
    do_sat_sd=sd(do_sat_m, na.rm=T),
    SPC_m_m=mean(SPC_m, na.rm=T),
    SPC_sd=sd(SPC_m, na.rm=T))

nut_month_dat <- covariat_datq%>%
  group_by(site, month, water_year, substrate) %>%
  summarise(
    NO3_mgL = mean(NO3_mgL_dl, na.rm=T),
    NO3_mgL_sd= sd(NO3_mgL_dl, na.rm=T),
    NH4_mgL = mean(NH4_mgL_dl, na.rm=T),
    NH4_mgL_sd= sd(NH4_mgL_dl, na.rm=T),
    PO4_ugL=mean(PO4_ugL_dl, na.rm=T),
    PO4_ugL_sd=sd(PO4_ugL_dl, na.rm=T),
    DOC_mgL=mean(DOC_mgL_dl, na.rm=T),
    DOC_mgL_sd=sd(DOC_mgL_dl, na.rm=T))

##############
##############
## BOX plots
# Install ggpattern if not installed
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  install.packages("ggpattern")
}

library(ggpattern)
library(ggplot2)
library(dplyr)

# Define fill patterns for water years
water_year_patterns <- c("wave","weave","stripe","circle", "stripe", "circle")
# Updated plot code
temp_plot <- ggplot(
  covariat_datq %>% filter(month > 5 & month < 10 & water_year>2020),
    aes(x = site, y = wtemp, color = site, fill = site)
  ) +
    geom_boxplot_pattern(
      aes(pattern = as.factor(water_year)),
      pattern_fill = "black",            # Pattern color is black
      pattern_alpha = 0.5,               # Semi-transparent patterns
      pattern_density = 0.35,             # Adjust pattern density
      pattern_spacing = 0.05             # Adjust pattern spacing
    ) +
    scale_color_manual(values = site_colors) +
    scale_fill_manual(values = site_colors) + # Fill by site
    scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom"
    )

# Display plot
print(temp_plot)


# Updated plot code
Q_plot <- ggplot(
  covariat_datq %>% filter(month > 7 & month < 10 & water_year>2020),
  aes(x = site, y = Q_m, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
Q_plot


# Updated plot code
SPC_plot <- ggplot(
  covariat_datq %>% filter(month > 4 & month < 10 & water_year>2020),
  aes(x = site, y = SPC_m, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )

# Display plot
print(SPC_plot)

# Updated plot code
### need to update : 
pH_plot <- ggplot(
  covariat_datq %>% filter(month > 4 & month < 10 & water_year>2020),
  aes(x = site, y = pH, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(pH_plot)



# Updated plot code
NO3_plot <- ggplot(
  covariat_datq %>% filter(month > 4 & month < 10 & water_year>2020),
  aes(x = site, y = NO3_mgL_dl, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(NO3_plot)


# Updated plot code
NH4_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = NH4_mgL_dl, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(NH4_plot)


# Updated plot code
PO4_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = PO4_ugL_dl, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(PO4_plot)




# Updated plot code
DOC_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = DOC_mgL_dl, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(DOC_plot)



###
# Updated plot code
OM_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = AFDM_mgg, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(OM_plot)

biom_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = AFDM_mgcm2, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(biom_plot)



chla_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = Chla_ugL_Q, color = site, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(biom_plot)

### big grid 
up_grid3 <- ggarrange(
  Q_plot,
  temp_plot,
  SPC_plot,
  pH_plot, # not done
  NO3_plot,
  NH4_plot,
  PO4_plot,
  DOC_plot,
  OM_plot,
  biom_plot,
  chla_plot,
  chla_plot,
  ncol = 4, nrow = 3,
  common.legend = TRUE, 
  legend = "bottom")

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/Draft_figure3_CH1_grid.png", plot = up_grid3, width = 12, height = 12, units = "in")

























##############
##############



##################
### uptake plots:

Up_plot<- ggplot(covariat_datq, aes(x=AFDM_mgcm2, color =site, shape=method)) + 
  geom_point(aes(y=Uadd.in.ug.L.s), size=2, alpha=0.75) +
  theme_bw() +
  labs(y=expression(Areal~uptake~rate~Uadd~(ug~L^-1~s^-1)), x= expression(Epilithic~biomass~(mg~cm^-2))) + 
  geom_errorbar(aes(ymin = (Uadd.in.ug.L.s - Uadd_sd), ymax = (Uadd.in.ug.L.s + Uadd_sd)), width=0) + 
  scale_color_manual(values = site_colors) +
    theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") +facet_grid(.~method)


sw_plot2<- ggplot(covariat_datq, aes(y=sw.in.m, x=AFDM_mgcm2, color =site, shape=method)) + 
  geom_point(aes(y=sw.in.m), size=2, alpha=0.75) +
  theme_bw() + #geom_smooth(method="lm", se=F)+
  labs(y=expression(Nutrient~uptake~length~(Sw)~(m)), x= expression(Epilithic~biomass~(mg~cm^-2))) + 
  geom_errorbar(aes(ymin = (sw.in.m - sw_sd), ymax = (sw.in.m + sw_sd)), width=0) +  
  scale_color_manual(values = site_colors) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") 


Upv_plot2<- ggplot(covariat_datq, aes(y=vf, x=AFDM_mgcm2, color =site, shape=method)) + 
  geom_point(aes(y=vf), size=2, alpha=0.75) +
  theme_bw() + #geom_smooth(method="lm", se=F)+
  labs(y=expression(Nutrient~uptake~velocity~(vf)~(m~s^-1)), x= expression(Epilithic~biomass~(mg~cm^-2))) + 
  geom_errorbar(aes(ymin = (vf - vf.sd), ymax = (vf + vf.sd)), width=0) +  
  scale_color_manual(values=c("#3283a8", "#3258a8", "#a67d17","#a65d17")) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") 



up_grid <- ggarrange(Up_plot,
                     Up_plot2,
                     ncol = 2, nrow = 1,
                     common.legend = TRUE, 
                     legend = "bottom")


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/prelim_up_Biomass_coresites_CH1_grid.png", plot = up_grid, width = 6, height = 3, units = "in")



#### glms
library(lme4)
library(lmerTest)
### for uptake rates 
hist(covariat_dat_up$Uadd_m)
hist(log(covariat_dat_up$Uadd_m+1))
names(covariat_dat_up)
## co variance matrix 
library(PerformanceAnalytics)
chart.Correlation(covariat_dat_up[,c("AFDM_mgcm2", "Chla_ugL_Q", "AFDM_mgg", "wt", "SPC")])

uadd_model <- lmer(log(Uadd_m+1) ~ scale(AFDM_mgcm2) + scale(wt) + (1|Site), data=covariat_dat_up)
summary(uadd_model)


uadd_model <- lmer(sw ~ scale(AFDM_mgcm2) + scale(wt) + (1|Site), data=covariat_dat_up)
summary(uadd_model)


Up_plot_om<- ggplot(covariat_dat_up%>%filter(Uadd<0.0045), aes(x=AFDM_mgg, color =Site, shape=up_method)) + 
  geom_point(aes(y=Uadd*1000), size=2, alpha=0.75) +
  theme_bw() +
  labs(y=expression(Areal~uptake~rate~Uadd~(ug~L^-1~s^-1)), x= expression(Organic~matter~(g~mg^-1))) + 
  scale_color_manual(values=c("#3283a8", "#3258a8", "#a67d17","#a65d17")) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") 


Up_plot_om2<- ggplot(covariat_dat_up%>%filter(sw<1800), aes(y=sw, x=AFDM_mgg, color =Site, shape=up_method)) + 
  geom_point(aes(y=sw), size=2, alpha=0.75) +
  theme_bw() + #geom_smooth(method="lm", se=F)+
  labs(y=expression(Nutrient~uptake~length~(Sw)~(m)), x= expression(Organic~matter~(g~mg^-1))) + 
  scale_color_manual(values=c("#3283a8", "#3258a8", "#a67d17","#a65d17")) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") 




Up_plot_chla<- ggplot(covariat_dat_up%>%filter(Uadd<0.0045), aes(x=Chla_ugL_Q, color =Site, shape=up_method)) + 
  geom_point(aes(y=Uadd*1000), size=2, alpha=0.75) +
  theme_bw() +
  labs(y=expression(Areal~uptake~rate~Uadd~(ug~L^-1~s^-1)), x= expression(Epilithic~Chla~(ug~L^-1))) + 
  scale_color_manual(values=c("#3283a8", "#3258a8", "#a67d17","#a65d17")) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") 


Up_plot_chla2<- ggplot(covariat_dat_up%>%filter(sw<1800), aes(y=sw, x=Chla_ugL_Q, color =Site, shape=up_method)) + 
  geom_point(aes(y=sw), size=2, alpha=0.75) +
  theme_bw() + #geom_smooth(method="lm", se=F)+
  labs(y=expression(Nutrient~uptake~length~(Sw)~(m)), x= expression(Epilithic~Chla~(ug~L^-1))) + 
  scale_color_manual(values=c("#3283a8", "#3258a8", "#a67d17","#a65d17")) +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom") 



up_grid3 <- ggarrange(Up_plot,
                      Up_plot_om,
                      Up_plot_chla,
                      Up_plot2,
                      Up_plot_om2,
                      Up_plot_chla2,
                      ncol = 3, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/prelim_up_OM_bio_coresites_CH1_grid.png", plot = up_grid3, width = 9, height = 6, units = "in")


### some WQ data 
BWL_WQ <- readRDS("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /data/24_BWL_DO_flag_sat_light.rds")

BWL_WQ<- BWL_WQ%>%filter(datetime > as.POSIXct("2021-01-01 07:15:00 PDT"))
BWL_WQ1 <-BWL_WQ[!is.na(BWL_WQ$do.obs), ]

range(BWL_WQ1$datetime)


BWU_WQ <- readRDS("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /data/24_BWU_DO_flag_sat_light.rds")

BWU_WQ1 <-BWU_WQ[!is.na(BWU_WQ$do.obs), ]

range(BWU_WQ1$datetime)

















# #################################
# ##################################
# 
# GBL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output_v2.rds")
# names(GBL_nh3_datq)
# GBL_nh3_datq$up_method <- "NH4"
# str(GBL_nh3_datq)
# 
# GBL_nh3_day <- GBL_nh3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NH3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#          /(ln_injectate_ratio*-1))*100))
# 
# GBU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_output_v2.rds")
# names(GBU_nh3_datq)
# str(GBU_nh3_datq)
# GBU_nh3_datq$up_method <- "NH4"
# 
# GBU_nh3_day <- GBU_nh3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NH3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# BWL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output_v2.rds")
# names(BWL_nh3_datq)
# str(BWL_nh3_datq)
# BWL_nh3_datq$up_method <- "NH4"
# 
# 
# BWL_nh3_day <- BWL_nh3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NH3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# 
# BWU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output_v2.rds")
# names(BWU_nh3_datq)
# str(BWU_nh3_datq)
# BWU_nh3_datq$up_method <- "NH4"
# 
# BWU_nh3_day <- BWU_nh3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NH3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# 
# 
# ## example analysis GBL
# covariat_dat_GBL <- covariat_datq%>%
#   filter(Site=="GBL")%>%
#   left_join(GBL_nh3_day, by=c("date"))
# 
# GBU_nh3_day1 <- GBU_nh3_day %>%
#   filter(!date==as.Date("2022-06-23"))
# 
# GBU_nh3_day1 <- GBU_nh3_day1 %>%
#   filter(!date==as.Date("2022-07-22"))
# 
# covariat_dat_GBU <- covariat_datq%>%
#   filter(Site=="GBU")%>%
#   left_join(GBU_nh3_day1, by=c("date"))
# 
# BWL_nh3_day1 <- BWL_nh3_day %>%
#   filter(!date==as.Date("2022-05-26"))
# 
# covariat_dat_BWL <- covariat_datq%>%
#   filter(Site=="BWL")%>%
#   left_join(BWL_nh3_day1, by=c("date"))
# 
# BWU_nh3_day1 <- BWU_nh3_day %>%
#   filter(!date==as.Date("2023-07-18"))
# 
# BWU_nh3_day1 <- BWU_nh3_day1 %>%
#   filter(!date==as.Date("2023-08-10"))
# 
# covariat_dat_BWU <- covariat_datq%>%
#   filter(Site=="BWU")%>%
#   left_join(BWU_nh3_day1, by=c("date"))
# 
# #####
# ####
# 
# 
# GBL_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_output_v2.rds")
# names(GBL_NO3_datq)
# GBL_NO3_datq$up_method <- "NO3"
# str(GBL_NO3_datq)
# 
# GBL_NO3_day <- GBL_NO3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NO3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# GBU_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_output_v2.rds")
# names(GBU_NO3_datq)
# str(GBU_NO3_datq)
# GBU_NO3_datq$up_method <- "NO3"
# 
# GBU_NO3_day <- GBU_NO3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NO3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# BWL_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_output_v2.rds")
# names(BWL_NO3_datq)
# str(BWL_NO3_datq)
# BWL_NO3_datq$up_method <- "NO3"
# 
# BWL_NO3_day <- BWL_NO3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NO3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# 
# BWU_NO3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_output_v2.rds")
# names(BWU_NO3_datq)
# str(BWU_NO3_datq)
# BWU_NO3_datq$up_method <- "NO3"
# 
# BWU_NO3_day <- BWU_NO3_datq%>%
#   filter(sw>0 & Uadd>0)%>%
#   group_by(date, up_method, ln_injectate_ratio) %>%
#   summarise(
#     sw_m= mean(sw, na.rm=T),
#     sw_sd= sd(sw, na.rm=T),
#     Uadd_m=mean(Uadd, na.rm=T),
#     Uadd_sd=sd(Uadd, na.rm=T),
#     Uadd_int_m=mean(Uadd_int, na.rm=T),
#     Uadd_int_sd=sd(Uadd_int, na.rm=T),
#     Vf_add_int_m=mean(Vf_add_int, na.rm=T),
#     Vf_add_int_sd=sd(Vf_add_int, na.rm=T),
#     TMR_Cl = max(TMR_Cl, na.rm = T),
#     TMR_N = max(TMR_NO3, na.rm = T)) %>%
#   mutate(sample_recovery= c(log(TMR_N/TMR_Cl)),
#          percent_recovery = c(
#            (((sample_recovery*-1)- (ln_injectate_ratio*-1))
#             /(ln_injectate_ratio*-1))*100))
# 
# 
# 
# 
# ## example analysis GBL
# GBL_NO3_day1 <- GBL_NO3_day %>%
#   filter(!date==as.Date("2021-06-23"))
# 
# GBL_NO3_day1 <- GBL_NO3_day1 %>%
#   filter(!date==as.Date("2023-06-01"))
# 
# covariat_dat_GBL_n3 <- covariat_datq%>%
#   filter(Site=="GBL")%>%
#   left_join(GBL_NO3_day1, by=c("date"))
# 
# GBU_NO3_day1 <- GBU_NO3_day %>%
#   filter(!date==as.Date("2021-06-23"))
# 
# GBU_NO3_day1 <- GBU_NO3_day1 %>%
#   filter(!date==as.Date("2022-04-07"))
# 
# GBU_NO3_day1 <- GBU_NO3_day1 %>%
#   filter(!date==as.Date("2022-10-03"))
# 
# covariat_dat_GBU_n3 <- covariat_datq%>%
#   filter(Site=="GBU")%>%
#   left_join(GBU_NO3_day1, by=c("date"))
# 
# 
# BWL_NO3_day1 <- BWL_NO3_day %>%
#   filter(!date==as.Date("2022-08-24"))
# 
# covariat_dat_BWL_n3 <- covariat_datq%>%
#   filter(Site=="BWL")%>%
#   left_join(BWL_NO3_day1, by=c("date"))
# 
# BWU_NO3_day1 <- BWU_NO3_day %>%
#   filter(!date==as.Date("2023-07-18"))
# 
# covariat_dat_BWU_n3 <- covariat_datq%>%
#   filter(Site=="BWU")%>%
#   left_join(BWU_NO3_day, by=c("date"))
# 
# 
# ####
# ####
# 
# covariat_dat_up <- rbind(covariat_dat_GBL, covariat_dat_GBU, covariat_dat_BWL, covariat_dat_BWU,
#                          covariat_dat_GBL_n3, covariat_dat_GBU_n3, covariat_dat_BWL_n3, covariat_dat_BWU_n3)