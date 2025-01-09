lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2", "ggpubr", "ggpattern",
         "lme4", "lmerTest", "MuMIn", "PerformanceAnalytics", "car"), require, character.only=T)
site_colors <- c(
  "BWL" = "#054fb9",
  "BWU" = "#C0D6DF",
  "GBL" = "#DD6E42",
  "GBU" = "#E8DAB2"
)


catch_colors <- c(
  "BW" = "#054fb9",
  "GB" = "#DD6E42"
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

## 
names(covariat_datq)
str(covariat_datq)
cor_data <- covariat_datq[, c(5, 10:14, 17, 33, 37)]
str(cor_data)

chart.Correlation(cor_data, histogram=TRUE, pch=19)

# ## composite varible of streamflow
# library(dplyr)
# 
# # Perform PCA for each site and extract PC1
# composite_pca <- covariat_datq %>%
#   group_by(site) %>%
#   filter(!is.na(Q_m) & !is.na(wtemp) & !is.na(SPC_m)) %>% # Remove rows with missing values
#   do({
#     # Perform PCA on the selected columns
#     pca_result <- prcomp(select(., Q_m, wtemp, SPC_m), scale. = TRUE)
#     
#     # Add PC1 as a column
#     mutate(., composite_var = pca_result$x[, 1])
#   }) %>%
#   ungroup()
# 
# # Merge the composite variable back into the original dataset
# covariat_datq <- left_join(covariat_datq, select(composite_pca, date, site, composite_var), 
#                            by = c("date", "site"))
# 
# 
# # Check variance explained by PC1
# pca_var_explained <- covariat_datq %>%
#   filter(!is.na(Q_m) & !is.na(wtemp) & !is.na(SPC_m)) %>% 
#   group_by(site) %>%
#   summarize(variance_explained = {
#     pca_result <- prcomp(select(., Q_m, wtemp, SPC_m), scale. = TRUE)
#     summary(pca_result)$importance[2, 1]
#   })
# 
# 
# ###
# # Compute the composite variable
# covariat_datq <- covariat_datq %>%
#   group_by(site) %>%
#   mutate(
#     Q_m_norm = (Q_m - mean(Q_m, na.rm = TRUE)) / sd(Q_m, na.rm = TRUE),
#     wtemp_norm = (wtemp - mean(wtemp, na.rm = TRUE)) / sd(wtemp, na.rm = TRUE),
#     SPC_m_norm = (SPC_m - mean(SPC_m, na.rm = TRUE)) / sd(SPC_m, na.rm = TRUE),
#     composite_var = Q_m_norm * wtemp_norm * SPC_m_norm
#   ) %>%
#   ungroup()

covariat_datq$TIN_mgL <- c(covariat_datq$NO3_mgL_dl + covariat_datq$NH4_mgL_dl)

glm_dat <-covariat_datq%>%
  filter(site=="BWL" | site=="BWU") %>%
  filter(!is.na(Q_m) & !is.na(wtemp) & !is.na(SPC_m) & !is.na(TIN_mgL) & !is.na(PO4_ugL_dl) &  !is.na(DOC_mgL_dl)) # Remove rows with missing values
  

glm_dat2 <-covariat_datq%>%
  filter(site=="BWL" | site=="BWU") %>%
  filter(!is.na(TIN_mgL) & !is.na(PO4_ugL_dl) &  !is.na(DOC_mgL_dl)) # Remove rows with missing values


####
#### GLMS
####

biomass_mod_gl <- lmer(log(AFDM_mgg+1)~ 
                         scale(NO3_mgL_dl) + 
                         scale(NH4_mgL_dl) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(SPC_m) + 
                         scale(wtemp) + 
                         scale(Q_m) + 
                         (1|water_year), data=glm_dat)

summary(biomass_mod_gl)
vif(biomass_mod_gl)
hist(residuals(biomass_mod_gl))


biomass_mod_glb <- lmer(log(AFDM_mgg+1)~ 
                         scale(NO3_mgL_dl) + 
                         scale(NH4_mgL_dl) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(SPC_sd) + 
                         scale(wtemp) + 
                         scale(Q_m) + 
                         (1|water_year), data=glm_dat)

summary(biomass_mod_glb)
vif(biomass_mod_glb)
hist(residuals(biomass_mod_glb))


biomass_mod_1 <- lmer(log(AFDM_mgg+1)~ 
                         scale(NO3_mgL_dl) + 
                         scale(NH4_mgL_dl) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         (1|water_year), data=glm_dat)

summary(biomass_mod_1)
vif(biomass_mod_1)
hist(residuals(biomass_mod_1))
r.squaredGLMM(biomass_mod_1)


biomass_mod_1a <- lmer(log(AFDM_mgg+1)~ 
                        scale(NO3_mgL_dl) + 
                        scale(NH4_mgL_dl) + 
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        (1|water_year), data=glm_dat2)

summary(biomass_mod_1a)
vif(biomass_mod_1a)
hist(residuals(biomass_mod_1a))
r.squaredGLMM(biomass_mod_1a)


biomass_mod_2 <- lmer(log(AFDM_mgg+1)~ 
                          scale(NO3_mgL_dl) + 
                          scale(NH4_mgL_dl) + 
                          scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                          scale(SPC_sd) + 
                          (1|water_year), data=glm_dat)

summary(biomass_mod_2)
vif(biomass_mod_2)
hist(residuals(biomass_mod_2))



biomass_mod_3 <- lmer(log(AFDM_mgg+1)~ 
                        scale(NO3_mgL_dl) + 
                        scale(NH4_mgL_dl) + 
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        scale(wtemp) + 
                        (1|water_year), data=glm_dat)

summary(biomass_mod_3)
vif(biomass_mod_3)
hist(residuals(biomass_mod_3))


biomass_mod_4 <- lmer(log(AFDM_mgg+1)~ 
                        scale(NO3_mgL_dl) + 
                        scale(NH4_mgL_dl) + 
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        scale(Q_m) + 
                        (1|water_year), data=glm_dat)

summary(biomass_mod_4)
vif(biomass_mod_4)
hist(residuals(biomass_mod_4))
AIC(biomass_mod_gl, biomass_mod_gl1, biomass_mod_2, biomass_mod_3, biomass_mod_4)


biomass_mod_5 <- lmer(log(AFDM_mgg+1)~ 
                        scale(TIN_mgL)+
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        scale(SPC_sd) + 
                        (1|water_year), data=glm_dat)

summary(biomass_mod_5)
vif(biomass_mod_5)
hist(residuals(biomass_mod_5))



biomass_mod_6 <- lmer(log(AFDM_mgg+1)~ 
                        scale(TIN_mgL)+
                        scale(PO4_ugL_dl)+ 
                        scale(SPC_m) + 
                        scale(DOC_mgL_dl)+
                        (1|water_year), data=glm_dat)

summary(biomass_mod_6)
vif(biomass_mod_6)
hist(residuals(biomass_mod_6))



biomass_mod_7 <- lmer(log(AFDM_mgg+1)~ 
                        scale(TIN_mgL)+
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        scale(Q_m) + 
                        (1|water_year), data=glm_dat)

summary(biomass_mod_7)
vif(biomass_mod_7)
hist(residuals(biomass_mod_7))




biomass_mod_8 <- lmer(log(AFDM_mgg+1)~ 
                        scale(TIN_mgL)+
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        scale(wtemp) + 
                        (1|water_year), data=glm_dat)

summary(biomass_mod_8)
vif(biomass_mod_8)
hist(residuals(biomass_mod_8))


biomass_mod_9 <- lmer(log(AFDM_mgg+1)~ 
                        scale(TIN_mgL)+
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        scale(wtemp) + 
                        scale(Q_m) + 
                        (1|water_year), data=glm_dat)

summary(biomass_mod_9)
vif(biomass_mod_9)
hist(residuals(biomass_mod_9))

biomass_mod_1b <- lmer(log(AFDM_mgg+1)~ 
                        scale(NO3_mgL_dl*1000) + 
                        scale(NH4_mgL_dl*1000) + 
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+ (1|substrate), data=glm_dat)
summary(biomass_mod_1b)
vif(biomass_mod_1b)
hist(residuals(biomass_mod_1b))
r.squaredGLMM(biomass_mod_1b)


AIC(biomass_mod_9, biomass_mod_8, biomass_mod_7, biomass_mod_6, biomass_mod_5, biomass_mod_4, biomass_mod_3, biomass_mod_2, biomass_mod_1, biomass_mod_glb, biomass_mod_gl, biomass_mod_1b)


#######
#######
# coeff plot 
biomass_mod_1 <- lmer(log(AFDM_mgg+1)~ 
                        scale(NO3_mgL_dl*1000) + 
                        scale(NH4_mgL_dl*1000) + 
                        scale(PO4_ugL_dl)+ 
                        scale(DOC_mgL_dl)+
                        (1|water_year), data=glm_dat)

summary(biomass_mod_1)
vif(biomass_mod_1)
hist(residuals(biomass_mod_1))
r.squaredGLMM(biomass_mod_1)


library(broom.mixed)

# Extract fixed effects and confidence intervals
fixed_effects <- broom.mixed::tidy(biomass_mod_1, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects1 <- fixed_effects[fixed_effects$term != "(Intercept)", ]
fixed_effects1 <- fixed_effects1[order(fixed_effects1$estimate), ]  # Order by estimate
fixed_effects1$Catchment <- "BW"
### gB
glm_dat_gb <-covariat_datq%>%
  filter(site=="GBL" | site=="GBU") %>%
  filter(!is.na(Q_m) & !is.na(wtemp) & !is.na(SPC_m) & !is.na(TIN_mgL) & !is.na(PO4_ugL_dl) &  !is.na(DOC_mgL_dl)) # Remove rows with missing values



biomass_mod_gl_g <- lmer(log(AFDM_mgg+1)~ 
                         scale(NO3_mgL_dl) + 
                         scale(NH4_mgL_dl) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(SPC_m) + 
                         scale(wtemp) + 
                         scale(Q_m) + 
                         (1|water_year), data=glm_dat_gb)

summary(biomass_mod_gl_g)
vif(biomass_mod_gl_g)
hist(residuals(biomass_mod_gl_g))





biomass_mod_1_g <- lmer(log(AFDM_mgg+1)~ 
                           scale(NO3_mgL_dl) + 
                           scale(NH4_mgL_dl) + 
                           scale(PO4_ugL_dl)+ 
                           scale(DOC_mgL_dl)+
                           (1|water_year), data=glm_dat_gb)

summary(biomass_mod_1_g)
vif(biomass_mod_1_g)
hist(residuals(biomass_mod_1_g))




biomass_mod_2_g <- lmer(log(AFDM_mgg+1)~ 
                          scale(NO3_mgL_dl) + 
                          scale(NH4_mgL_dl) + 
                          scale(PO4_ugL_dl)+ 
                          scale(DOC_mgL_dl)+
                          scale(SPC_m) + 
                          (1|water_year), data=glm_dat_gb)

summary(biomass_mod_2_g)
vif(biomass_mod_2_g)
hist(residuals(biomass_mod_2_g))





biomass_mod_3_g <- lmer(log(AFDM_mgg+1)~ 
                          scale(NO3_mgL_dl) + 
                          scale(NH4_mgL_dl) + 
                          scale(PO4_ugL_dl)+ 
                          scale(DOC_mgL_dl)+
                          scale(wtemp) + 
                          (1|water_year), data=glm_dat_gb)

summary(biomass_mod_3_g)
vif(biomass_mod_3_g)
hist(residuals(biomass_mod_3_g))





biomass_mod_4_g <- lmer(log(AFDM_mgg+1)~ 
                          scale(NO3_mgL_dl) + 
                          scale(NH4_mgL_dl) + 
                          scale(PO4_ugL_dl)+ 
                          scale(DOC_mgL_dl)+
                          scale(Q_m) + 
                          (1|water_year), data=glm_dat_gb)

summary(biomass_mod_4_g)
vif(biomass_mod_4_g)
hist(residuals(biomass_mod_4_g))


AIC(biomass_mod_gl_g, biomass_mod_1_g, biomass_mod_2_g, biomass_mod_3_g,biomass_mod_4_g)

biomass_mod_1_g <- lmer(log(AFDM_mgg+1)~ 
                          scale(NO3_mgL_dl*1000) + 
                          scale(NH4_mgL_dl*1000) + 
                          scale(PO4_ugL_dl)+ 
                          scale(DOC_mgL_dl)+
                          (1|water_year), data=glm_dat_gb)

summary(biomass_mod_1_g)
vif(biomass_mod_1_g)
hist(residuals(biomass_mod_1_g))
r.squaredGLMM(biomass_mod_1_g)


fixed_effects_gb <- broom.mixed::tidy(biomass_mod_1_g, effects = "fixed", conf.int = TRUE)

# Filter to exclude intercept and arrange predictors
fixed_effects_gb1 <- fixed_effects_gb[fixed_effects_gb$term != "(Intercept)", ]
fixed_effects_gb1 <- fixed_effects_gb1[order(fixed_effects_gb1$estimate), ]  # Order by estimate
fixed_effects_gb1$Catchment <- "GB"

fxd_eff <- rbind(fixed_effects_gb1, fixed_effects1)

fxd_eff1 <- fxd_eff %>%
  mutate(Significance = case_when(
    p.value <= 0.055 ~ "sig",
    p.value > 0.055 & p.value <= 0.09  ~ "marg.",
    p.value > 0.09  ~ "not sig."))


term_labels <- c(
  "scale(NO3_mgL_dl * 1000)" = expression(NO[3]~(µg~L^-1)),
  "scale(PO4_ugL_dl)" = expression(PO[4]~(µg~L^-1)),
  "scale(DOC_mgL_dl)" = expression(DOC~(mg~L^-1)),
  "scale(NH4_mgL_dl * 1000)" = expression(NH[4]~(µ~L^-1))
)

# Add an offset column based on Catchment
fxd_eff1 <- fxd_eff1 %>%
  mutate(y_position = as.numeric(reorder(term, estimate)) + ifelse(Catchment == "BW", 0.2, 0))


Biomass_coeff_plot <-ggplot(fxd_eff1, aes(x = estimate, y = reorder(term, estimate), color = Catchment, shape = Significance)) +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 0.85, position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_shape_manual(values = c(1, 4, 19)) +
  scale_color_manual(values = catch_colors) +
  geom_vline(xintercept = 0) +
  labs(
    x = "Scaled effect size",
    y = "Predictors",
    title = "Epilithic biomass"
  ) +
  scale_y_discrete(labels = term_labels) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )





##########
#########
########

chla_mod_gl <- lmer(log(Chla_ugL_Q+1)~ 
                         scale(NO3_mgL_dl) + 
                         scale(NH4_mgL_dl) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(SPC_m) + 
                         scale(wtemp) + 
                         scale(Q_m) + 
                         (1|water_year), data=glm_dat)

summary(chla_mod_gl)
vif(chla_mod_gl)
hist(residuals(chla_mod_gl))
r.squaredGLMM(chla_mod_gl)


chla_mod_g2 <- lmer(log(Chla_ugL_Q+1)~ 
                      scale(TIN_mgL) + 
                      scale(NH4_mgL_dl) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(SPC_m) + 
                      scale(wtemp) + 
                      scale(Q_m) + 
                      (1|water_year), data=glm_dat)

summary(chla_mod_g2)
vif(chla_mod_g2)
hist(residuals(chla_mod_g2))



chla_mod_1 <- lmer(log(Chla_ugL_Q+1)~ 
                      scale(NO3_mgL_dl) + 
                      scale(NH4_mgL_dl) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(SPC_m) + 
                      #scale(wtemp) + 
                      scale(Q_m) + 
                      (1|water_year), data=glm_dat)

summary(chla_mod_1)
vif(chla_mod_1)
hist(residuals(chla_mod_1))




chla_mod_2 <- lmer(log(Chla_ugL_Q+1)~ 
                     scale(NO3_mgL_dl) + 
                     scale(NH4_mgL_dl) + 
                     scale(PO4_ugL_dl)+ 
                     scale(DOC_mgL_dl)+
                     #scale(SPC_m) + 
                     scale(wtemp) + 
                     scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_2)
vif(chla_mod_2)
hist(residuals(chla_mod_2))





chla_mod_3 <- lmer(log(Chla_ugL_Q+1)~ 
                     scale(NO3_mgL_dl) + 
                     scale(NH4_mgL_dl) + 
                     scale(PO4_ugL_dl)+ 
                     scale(DOC_mgL_dl)+
                     scale(SPC_m) + 
                     #scale(wtemp) + 
                     scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_3)
vif(chla_mod_3)
hist(residuals(chla_mod_3))






chla_mod_4 <- lmer(log(Chla_ugL_Q+1)~ 
                     scale(NO3_mgL_dl) + 
                     scale(NH4_mgL_dl) + 
                     scale(PO4_ugL_dl)+ 
                     scale(DOC_mgL_dl)+
                     #scale(SPC_m) + 
                     scale(wtemp) + 
                     #scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_4)
vif(chla_mod_4)
hist(residuals(chla_mod_4))




chla_mod_5 <- lmer(log(Chla_ugL_Q+1)~ 
                     #scale(NO3_mgL_dl) + 
                     #scale(NH4_mgL_dl) + 
                     #scale(PO4_ugL_dl)+ 
                     #scale(DOC_mgL_dl)+
                     scale(SPC_m) + 
                     scale(wtemp) + 
                     scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_5)
vif(chla_mod_5)
hist(residuals(chla_mod_5))
r.squaredGLMM(chla_mod_5)




chla_mod_6 <- lmer(log(Chla_ugL_Q+1)~ 
                     #scale(NO3_mgL_dl) + 
                     #scale(NH4_mgL_dl) + 
                     #scale(PO4_ugL_dl)+ 
                     #scale(DOC_mgL_dl)+
                     scale(SPC_m) + 
                     scale(wtemp) + 
                    # scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_6)
vif(chla_mod_6)
hist(residuals(chla_mod_6))




chla_mod_7 <- lmer(log(Chla_ugL_Q+1)~ 
                     scale(NO3_mgL_dl) + 
                     scale(NH4_mgL_dl) + 
                     scale(PO4_ugL_dl)+ 
                     scale(DOC_mgL_dl)+
                     scale(SPC_m) + 
                     scale(wtemp) + 
                     #scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_7)
vif(chla_mod_7)
hist(residuals(chla_mod_7))


AIC(chla_mod_g2, chla_mod_gl, chla_mod_1, chla_mod_2, chla_mod_3, chla_mod_4, chla_mod_5, chla_mod_6, chla_mod_7)


chla_mod_glg <- lmer(log(Chla_ugL_Q+1)~ 
                      scale(NO3_mgL_dl) + 
                      scale(NH4_mgL_dl) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(SPC_m) + 
                      scale(wtemp) + 
                      scale(Q_m) + 
                      (1|water_year), data=glm_dat_gb)

summary(chla_mod_glg)
vif(chla_mod_glg)
hist(residuals(chla_mod_glg))


chla_mod_g2g <- lmer(log(Chla_ugL_Q+1)~ 
                      scale(TIN_mgL) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(SPC_m) + 
                      scale(wtemp) + 
                      scale(Q_m) + 
                      (1|water_year), data=glm_dat_gb)

summary(chla_mod_g2g)
vif(chla_mod_g2g)
hist(residuals(chla_mod_g2g))





chla_mod_4g <- lmer(log(Chla_ugL_Q+1)~ 
                      scale(NO3_mgL_dl) + 
                      #scale(NH4_mgL_dl) + 
                      #scale(PO4_ugL_dl)+ 
                      #scale(DOC_mgL_dl)+
                      scale(SPC_m) + 
                      scale(wtemp) + 
                      scale(Q_m) + 
                      (1|water_year), data=glm_dat_gb)

summary(chla_mod_4g)
vif(chla_mod_4g)
hist(residuals(chla_mod_4g))
r.squaredGLMM(chla_mod_4g)



chla_mod_5 <- lmer(log(Chla_ugL_Q+1)~ 
                     scale(wtemp) + 
                     scale(Q_m) + 
                     (1|water_year), data=glm_dat)

summary(chla_mod_5)
vif(chla_mod_5)
hist(residuals(chla_mod_5))
r.squaredGLMM(chla_mod_5)


# Extract fixed effects and confidence intervals
fixed_effects_chl <- broom.mixed::tidy(chla_mod_5, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects_chl1 <- fixed_effects_chl[fixed_effects_chl$term != "(Intercept)", ]
fixed_effects_chl1 <- fixed_effects_chl1[order(fixed_effects_chl1$estimate), ]  # Order by estimate
fixed_effects_chl1$Catchment <- "BW"

## GB
chla_mod_5g <- lmer(log(Chla_ugL_Q+1)~ 
                     scale(SPC_m) + 
                     scale(wtemp) + 
                     scale(Q_m) + 
                     (1|water_year), data=glm_dat_gb)

summary(chla_mod_5g)
vif(chla_mod_5g)
hist(residuals(chla_mod_5g))
r.squaredGLMM(chla_mod_5g)


fixed_effects_chlg <- broom.mixed::tidy(chla_mod_5g, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects_chlg1 <- fixed_effects_chlg[fixed_effects_chlg$term != "(Intercept)", ]
fixed_effects_chlg1 <- fixed_effects_chlg1[order(fixed_effects_chlg1$estimate), ]  # Order by estimate
fixed_effects_chlg1$Catchment <- "GB"


fxd_eff_chla <- rbind(fixed_effects_chlg1, fixed_effects_chl1)

fxd_eff_chla1 <- fxd_eff_chla %>%
  mutate(Significance = case_when(
    p.value <= 0.055 ~ "sig.",
    p.value > 0.055 & p.value <= 0.09  ~ "marg.",
    p.value > 0.09  ~ "not sig."))


term_labels_chla <- c(
 "scale(SPC_m)" = expression(SPC~(µS~cm^-1)),
  "scale(wtemp)" = expression(Temp~(degree~C)),
  "scale(Q_m)" = expression(Q~(m^3~s^-1))
)

# Add an offset column based on Catchment
fxd_eff_chla1 <- fxd_eff_chla1 %>%
  mutate(y_position = as.numeric(reorder(term, estimate)) + ifelse(Catchment == "BW", 0.2, 0))


chla_coeff_plot <-ggplot(fxd_eff_chla1, aes(x = estimate, y = reorder(term, estimate), color = Catchment, shape = Significance)) +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 0.85, position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_shape_manual(values = c(1, 4, 19)) +
  scale_color_manual(values = catch_colors) +
  geom_vline(xintercept = 0) +
  labs(
    x = "Scaled effect size",
    y = "Predictors",
    title = "Epilithic chl-a"
  ) +
  scale_y_discrete(labels = term_labels_chla) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )



chla_coeff_plot_temp <-ggplot(fxd_eff_chla1, aes(x = estimate, y = reorder(term, estimate), color = Catchment, shape = Significance)) +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 0.85, position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_shape_manual(values = c(1, 4, 19)) +
  scale_color_manual(values = catch_colors) +
  geom_vline(xintercept = 0) +
  labs(
    x = "Scaled effect size",
    y = "Predictors",
    title = "*Will be NEP"
  ) +
  scale_y_discrete(labels = term_labels_chla) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )


### big grid 
coeff_grid <- ggarrange(
  Biomass_coeff_plot,
  chla_coeff_plot,
  chla_coeff_plot_temp,
  ncol = 1, nrow = 3,
  common.legend = TRUE, 
  legend = "bottom")


coeff_grid <- ggarrange(
  Biomass_coeff_plot + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)),
  chla_coeff_plot + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)),
  chla_coeff_plot_temp + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)),
  labels=c("a", "b", "c"),
  
  ncol = 1, nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/Draft_Fig4_coeff_ecol_CH1_grid.png", plot = coeff_grid, width = 3.5, height = 8, units = "in")













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

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/Draft_figure3_CH1_grid.png", plot = up_grid3, width = 14, height = 10, units = "in")

