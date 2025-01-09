lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2", "ggpubr", "ggpattern"), require, character.only=T)

site_colors <- c(
  "BWL" = "#054fb9",
  "BWU" = "#C0D6DF",
  "GBL" = "#DD6E42",
  "GBU" = "#E8DAB2"
)

# Define fill patterns for water years
water_year_patterns <- c("wave","weave","stripe","circle", "stripe", "circle")

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

### fxns:
## Fxn for water year:
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
## Function to perform ANOVA and post-hoc test within each site
comparisons <- list(
  c("2021", "2022"),
  c("2021", "2023"),
  c("2022", "2023"),
  c("2023", "2024")
)

####
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

### 
# fill in NA wtemp with temps from SPC 
covariat_datq <- covariat_datq  %>%
  mutate(wtemp = ifelse(is.na(wtemp), wt, wtemp))
  
str(covariat_datq)

##############
## BOX plots
temp_plot <- ggplot(
  covariat_datq %>% filter(month > 5 & month < 10 & water_year>2020),
    aes(x = site, y = wtemp, fill = site)
  ) +
    geom_boxplot_pattern(
      aes(pattern = as.factor(water_year)),
      pattern_fill = "black",            # Pattern color is black
      pattern_alpha = 0.5,               # Semi-transparent patterns
      pattern_density = 0.35,             # Adjust pattern density
      pattern_spacing = 0.05             # Adjust pattern spacing
    ) +
    #scale_color_manual(values = site_colors) +
    scale_fill_manual(values = site_colors) + # Fill by site
    scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
    theme_bw() +
  labs(y=expression(Water~temperature~(degree~C)), x=NULL) + 
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


#### ANOVA 
clean_data_p2 <- covariat_datq %>%
  filter(!is.na(wtemp), !is.nan(wtemp), month > 5, month < 10, water_year > 2020)
# Convert site and water_year to factors
clean_data_p2 <- clean_data_p2 %>%
  mutate(
    site = as.factor(site),
    water_year = as.factor(water_year)
  )
# Two-way ANOVA

get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov(wtemp ~ water_year, data = data)
  
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

anova_model_p2 <- aov(wtemp ~ site * water_year, data = clean_data_p2)
summary(anova_model_p2)
# Ensure that `site` and `water_year` are factors
covariat_dat_temp <- covariat_datq %>%
  mutate(
    site = as.factor(site),
    water_year = as.factor(water_year)
  )

# Apply the function to each site and combine the results
signif_pairs_all_sites_p2 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p2)  # Corrected object name here


# write.csv(x = signif_pairs_all_sites_p2, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_wtemp.csv", row.names = TRUE)

###
###
Q_plot <- ggplot(
  covariat_datq %>% filter(month > 5 & month < 10 & water_year>2020),
  aes(x = site, y = log(Q_m+1), fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(log(Q+1)~(m^3~s^-1)), x=NULL) + 
  
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
Q_plot

## ANOVA
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov(Q_m ~ water_year, data = data)
  
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p1 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p1)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p1, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_Q.csv", row.names = TRUE)


# Updated plot code
SPC_plot <- ggplot(
  covariat_datq %>% filter(month > 4 & month < 10 & water_year>2020),
  aes(x = site, y = SPC_m, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(SPC~(mu~S~cm^-1)), x=NULL) + 
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

## ANOVA
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov(SPC_m ~ water_year, data = data)
  
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p3 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p3)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p3, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_SPC.csv", row.names = TRUE)


# Updated plot code
### need to update : 
pH_plot <- ggplot(
  covariat_datq %>% filter(month > 4 & month < 10 & water_year>2020),
  aes(x = site, y = pH, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
 # scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(pH), x=NULL) + 
  
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(pH_plot)

# @#$!#%$#@%^
## ANOVA - NOT up dated 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov(pH ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p4 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p4)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p4, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_pH.csv", row.names = TRUE)
## *^$#%$#@$%#@

covariat_datq_no3 <- covariat_datq%>%
  filter(NO3_mgL_dl<0.778)

NO3_plot <- ggplot(
  covariat_datq_no3 %>% filter(water_year>2020),
  aes(x = site, y = NO3_mgL_dl*1000, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(N-NO[3]~(mu~g~L^-1)), x=NULL) + 
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(NO3_plot)

## ANOVA - NOT up dated 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((NO3_mgL_dl*1000) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p5 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p5)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p5, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_NO3.csv", row.names = TRUE)



# Updated plot code
NH4_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = NH4_mgL_dl*1000, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(N-NH[4]~(mu~g~L^-1)), x=NULL) + 
  
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(NH4_plot)


## ANOVA - NOT up dated 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((NH4_mgL_dl*1000) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p6 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p6)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p6, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_NH4.csv", row.names = TRUE)


# Updated plot code
PO4_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = PO4_ugL_dl, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(SRP~(mu~g~L^-1)), x=NULL) + 
  
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(PO4_plot)

## ANOVA - NOT up dated 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((PO4_ugL_dl) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p7 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p7)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p7, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_SRP.csv", row.names = TRUE)



# Updated plot code
DOC_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = DOC_mgL_dl, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(DOC~(mg~L^-1)), x=NULL) + 
  
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(DOC_plot)

## ANOVA - NOT up dated 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((DOC_mgL_dl) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p8 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p8)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p8, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_DOC.csv", row.names = TRUE)




### NEED to double check GBL sites 
# Updated plot code
OM_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = AFDM_mgg, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(OM~(mg~g^-1)), x=NULL) + 
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(OM_plot)

## ANOVA - NOT up dated 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((AFDM_mgg) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p9 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p9)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p9, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_OM.csv", row.names = TRUE)


biom_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = AFDM_mgcm2, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  labs(y=expression(Biomass~(mg~cm^-2)), x=NULL) + 
  
  #scale_color_manual(values = site_colors) +
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



## ANOVA 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((AFDM_mgcm2) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p10 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p10)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p10, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_biomass.csv", row.names = TRUE)


chla_plot <- ggplot(
  covariat_datq %>% filter(water_year>2020),
  aes(x = site, y = Chla_ugL_Q, fill = site)
) +
  geom_boxplot_pattern(
    aes(pattern = as.factor(water_year)),
    pattern_fill = "black",            # Pattern color is black
    pattern_alpha = 0.5,               # Semi-transparent patterns
    pattern_density = 0.35,             # Adjust pattern density
    pattern_spacing = 0.05             # Adjust pattern spacing
  ) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) + # Fill by site
  scale_pattern_manual(values = water_year_patterns) + # Patterns for water_year
  theme_bw() +
  labs(y=expression(Chla~(mu~g~L^-1)), x=NULL) + 
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )
print(chla_plot)



## ANOVA 
get_significant_pairs <- function(data) {
  # Perform the ANOVA for each site separately
  anova_model <- aov((Chla_ugL_Q) ~ water_year, data = data)
  # Only perform post-hoc if ANOVA is significant
  if (summary(anova_model)[[1]][[5]][1] < 0.05) {
    posthoc <- TukeyHSD(anova_model)$water_year
    signif_pairs <- as.data.frame(posthoc) %>%
      rownames_to_column("comparison") %>%
      filter(`p adj` < 0.05) %>%
      mutate(
        stars = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      separate(comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(site = unique(data$site))  # Add site info for each result
    return(signif_pairs)
  } else {
    return(NULL)  # No significant results
  }
}

# Apply the function to each site and combine the results
signif_pairs_all_sites_p11 <- covariat_dat_temp %>% 
  group_by(site) %>% 
  group_split() %>% 
  map_df(get_significant_pairs)

# View the significant pairs
print(signif_pairs_all_sites_p11)  # Corrected object name here
# write.csv(x = signif_pairs_all_sites_p11, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/WQ_Anova\'s/24_NCYC_ANOVA_chla.csv", row.names = TRUE)


### NEP ###



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
  labels=c("a", "b", "c", "d","e", "f", "g", "h", "i", "j", "k", "l"),
  common.legend = TRUE, 
  legend = "bottom")

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/Draft_figure3_CH1_grid.png", plot = up_grid3, width = 14, height = 10, units = "in")

