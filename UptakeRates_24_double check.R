
library(dplyr)
library(pracma)
## =====================
## GBL
## =====================
### GBL NO3 ### 24_BWU_NH3_BTC_TMRout.rds
GBL_no3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_TMRout.rds")
## check range in TMR
hist(GBL_no3_dat_TMR$TMR_Cl)
hist(GBL_no3_dat_TMR$TMR_NO3)

names(GBL_no3_dat_TMR)

str(GBL_no3_dat_TMR)
GBL_no3_dat_TMR<- as.data.frame(GBL_no3_dat_TMR)

# Step 2: Spiraling metrics calculation
GBL_no3_dat_spiraling <- GBL_no3_dat_TMR %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NO3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(GBL_no3_dat_spiraling$ln_injectate_ratio)
range(GBL_no3_dat_spiraling$ln_injectate_ratio) # negative 
hist(GBL_no3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(GBL_no3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate NO3-N_add-int grouped by date
GBL_no3_dat_spiraling <- GBL_no3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative NO3-N
    NO3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    NO3_N_add_obs = ifelse(NO3_corrected > 0, NO3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_NO3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_add_obs),
    int_NO3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute NO3_N_add_int (flow-weighted integration)
    NO3_N_add_int = sqrt(
      (int_NO3_N_add_obs / int_Q) * (int_NO3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  NO3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NO3 = numeric(), slope_sample = numeric(), kw = numeric(),
  NO3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(GBL_no3_dat_spiraling$date)) {
  temp_data <- GBL_no3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "GBL_NO3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        NO3 = temp_dat$NO3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NO3 = temp_dat$TMR_NO3[2],
        NO3_N_add_int = temp_dat$NO3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(GBL_no3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NO3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * NO3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * NO3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & NO3_N_add_int > 0,
                        Uadd_int / NO3_N_add_int, 
                        NA)
  )

### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=Uadd_int*1000)) +
#  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("GBL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=sw)) +
#  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  labs(x = expression(paste("TMR nitrate (ugL)")),
       y= expression(paste("Uadd (ug m^-2 min^-1)"))) +
  ggtitle("GBL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")

str(data_q1)

# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_output_v2.rds")


################
### GBL NH3 ###
GBL_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_TMRout.rds")

# Step 2: Spiraling metrics calculation
GBL_nh3_dat_spiraling <- GBL_nh3_dat_TMR %>%
  group_by(date) %>%
  filter(!sample_no=="NA_H1")%>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NH3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(GBL_nh3_dat_spiraling$ln_injectate_ratio)
range(GBL_nh3_dat_spiraling$ln_injectate_ratio) # negative 
hist(GBL_nh3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(GBL_nh3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate nh3-N_add-int grouped by date
GBL_nh3_dat_spiraling <- GBL_nh3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative nh3-N
    nh3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    nh3_N_add_obs = ifelse(NH3_corrected > 0, NH3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_nh3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_add_obs),
    int_nh3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute nh3_N_add_int (flow-weighted integration)
    nh3_N_add_int = sqrt(
      (int_nh3_N_add_obs / int_Q) * (int_nh3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  nh3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NH3 = numeric(), slope_sample = numeric(), kw = numeric(),
  nh3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(GBL_nh3_dat_spiraling$date)) {
  temp_data <- GBL_nh3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "GBL_nh3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        nh3 = temp_dat$NH3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NH3 = temp_dat$TMR_NH3[2],
        nh3_N_add_int = temp_dat$nh3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(GBL_nh3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NH3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = c(1/ kw), # Spiraling length
    Uadd = offset_Q * 1000 * nh3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * nh3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & nh3_N_add_int > 0,
                        Uadd_int / nh3_N_add_int, 
                        NA)
  )


### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  #ggtitle("GBL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  theme_bw() + facet_wrap(date~., scales = "free")



# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output_v2.rds")


## =====================
## GBU
## =====================
### GBU NO3 ###
GBU_no3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_TMRout.rds")

###
####

## check range in TMR
hist(GBU_no3_dat_TMR$TMR_Cl)
hist(GBU_no3_dat_TMR$TMR_NO3)

names(GBU_no3_dat_TMR)

str(GBU_no3_dat_TMR)
GBU_no3_dat_TMR<- as.data.frame(GBU_no3_dat_TMR)

# Step 2: Spiraling metrics calculation
GBU_no3_dat_spiraling <- GBU_no3_dat_TMR %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NO3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(GBU_no3_dat_spiraling$ln_injectate_ratio)
range(GBU_no3_dat_spiraling$ln_injectate_ratio) # negative 
hist(GBU_no3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(GBU_no3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate NO3-N_add-int grouped by date
GBU_no3_dat_spiraling <- GBU_no3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative NO3-N
    NO3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    NO3_N_add_obs = ifelse(NO3_corrected > 0, NO3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_NO3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_add_obs),
    int_NO3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute NO3_N_add_int (flow-weighted integration)
    NO3_N_add_int = sqrt(
      (int_NO3_N_add_obs / int_Q) * (int_NO3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  NO3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NO3 = numeric(), slope_sample = numeric(), kw = numeric(),
  NO3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(GBU_no3_dat_spiraling$date)) {
  temp_data <- GBU_no3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "GBU_NO3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        NO3 = temp_dat$NO3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NO3 = temp_dat$TMR_NO3[2],
        NO3_N_add_int = temp_dat$NO3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(GBU_no3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NO3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * NO3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * NO3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & NO3_N_add_int > 0,
                        Uadd_int / NO3_N_add_int, 
                        NA)
  )

### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  labs(x = expression(paste("TMR nitrate (ugL)")),
       y= expression(paste("Uadd (ug m^-2 min^-1)"))) +
  ggtitle("GBU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  labs(x = expression(paste("TMR nitrate (ugL)")),
       y= expression(paste("Uadd (ug m^-2 min^-1)"))) +
  ggtitle("GBU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")

str(data_q)

# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_output_v2.rds")


################
### GBU NH3 ###
GBU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_TMRout.rds")
## check range in TMR
hist(GBU_nh3_datq$TMR_Cl)
hist(GBU_nh3_datq$TMR_NH3)

# Step 2: Spiraling metrics calculation
GBU_nh3_dat_spiraling <- GBU_nh3_datq %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NH3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(GBU_nh3_dat_spiraling$ln_injectate_ratio)
range(GBU_nh3_dat_spiraling$ln_injectate_ratio) # negative 
hist(GBU_nh3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(GBU_nh3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate nh3-N_add-int grouped by date
GBU_nh3_dat_spiraling <- GBU_nh3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative nh3-N
    nh3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    nh3_N_add_obs = ifelse(NH3_corrected > 0, NH3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_nh3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_add_obs),
    int_nh3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute nh3_N_add_int (flow-weighted integration)
    nh3_N_add_int = sqrt(
      (int_nh3_N_add_obs / int_Q) * (int_nh3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  nh3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NH3 = numeric(), slope_sample = numeric(), kw = numeric(),
  nh3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(GBU_nh3_dat_spiraling$date)) {
  temp_data <- GBU_nh3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "GBU_nh3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        nh3 = temp_dat$NH3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NH3 = temp_dat$TMR_NH3[2],
        nh3_N_add_int = temp_dat$nh3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(GBU_nh3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NH3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * nh3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * nh3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & nh3_N_add_int > 0,
                        Uadd_int / nh3_N_add_int, 
                        NA)
  )


### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  #ggtitle("GBU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("GBU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")

str(data_q)

# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_output_v2.rds")


## =====================
## BWL
## =====================

### BWL NO3 ###
BWL_no3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_TMRout.rds")
## check range in TMR
hist(BWL_no3_dat_TMR$TMR_Cl)
hist(BWL_no3_dat_TMR$TMR_NO3)

names(BWL_no3_dat_TMR)

str(BWL_no3_dat_TMR)
BWL_no3_dat_TMR<- as.data.frame(BWL_no3_dat_TMR)

# Step 2: Spiraling metrics calculation
BWL_no3_dat_spiraling <- BWL_no3_dat_TMR %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NO3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(BWL_no3_dat_spiraling$ln_injectate_ratio)
range(BWL_no3_dat_spiraling$ln_injectate_ratio) # negative 
hist(BWL_no3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(BWL_no3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate NO3-N_add-int grouped by date
BWL_no3_dat_spiraling <- BWL_no3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative NO3-N
    NO3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    NO3_N_add_obs = ifelse(NO3_corrected > 0, NO3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_NO3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_add_obs),
    int_NO3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute NO3_N_add_int (flow-weighted integration)
    NO3_N_add_int = sqrt(
      (int_NO3_N_add_obs / int_Q) * (int_NO3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  NO3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NO3 = numeric(), slope_sample = numeric(), kw = numeric(),
  NO3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(BWL_no3_dat_spiraling$date)) {
  temp_data <- BWL_no3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "BWL_NO3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        NO3 = temp_dat$NO3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NO3 = temp_dat$TMR_NO3[2],
        NO3_N_add_int = temp_dat$NO3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(BWL_no3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NO3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * NO3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * NO3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & NO3_N_add_int > 0,
                        Uadd_int / NO3_N_add_int, 
                        NA)
  )

### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("BWL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("BWL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")

str(data_q1)

# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_output_v2.rds")


### BWL NH3 ###
BWL_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_TMRout.rds")

## check range in TMR
hist(BWL_nh3_dat_TMR$TMR_Cl)
hist(BWL_nh3_dat_TMR$TMR_NH3)

# Step 2: Spiraling metrics calculation
BWL_nh3_dat_spiraling <- BWL_nh3_dat_TMR %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NH3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(BWL_nh3_dat_spiraling$ln_injectate_ratio)
range(BWL_nh3_dat_spiraling$ln_injectate_ratio) # negative 
hist(BWL_nh3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(BWL_nh3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate nh3-N_add-int grouped by date
BWL_nh3_dat_spiraling <- BWL_nh3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative nh3-N
    nh3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    nh3_N_add_obs = ifelse(NH3_corrected > 0, NH3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_nh3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_add_obs),
    int_nh3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute nh3_N_add_int (flow-weighted integration)
    nh3_N_add_int = sqrt(
      (int_nh3_N_add_obs / int_Q) * (int_nh3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  nh3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NH3 = numeric(), slope_sample = numeric(), kw = numeric(),
  nh3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(BWL_nh3_dat_spiraling$date)) {
  temp_data <- BWL_nh3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "BWL_nh3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        nh3 = temp_dat$NH3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NH3 = temp_dat$TMR_NH3[2],
        nh3_N_add_int = temp_dat$nh3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(BWL_nh3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NH3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * nh3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * nh3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & nh3_N_add_int > 0,
                        Uadd_int / nh3_N_add_int, 
                        NA)
  )


### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  #ggtitle("BWL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("BWL nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")

# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output_v2.rds")

## =====================
## BWU
## =====================

### BWU NO3 ###
BWU_no3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_TMRout.rds")

## check range in TMR
hist(BWU_no3_dat_TMR$TMR_Cl)
hist(BWU_no3_dat_TMR$TMR_NO3)
names(BWU_no3_dat_TMR)

str(BWU_no3_dat_TMR)
BWU_no3_dat_TMR<- as.data.frame(BWU_no3_dat_TMR)

# Step 2: Spiraling metrics calculation
BWU_no3_dat_spiraling <- BWU_no3_dat_TMR %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g/ NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NO3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(BWU_no3_dat_spiraling$ln_injectate_ratio)
range(BWU_no3_dat_spiraling$ln_injectate_ratio) # negative 
hist(BWU_no3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(BWU_no3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate NO3-N_add-int grouped by date
BWU_no3_dat_spiraling <- BWU_no3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative NO3-N
    NO3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    NO3_N_add_obs = ifelse(NO3_corrected > 0, NO3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_NO3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_add_obs),
    int_NO3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * NO3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute NO3_N_add_int (flow-weighted integration)
    NO3_N_add_int = sqrt(
      (int_NO3_N_add_obs / int_Q) * (int_NO3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  NO3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NO3 = numeric(), slope_sample = numeric(), kw = numeric(),
  NO3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(BWU_no3_dat_spiraling$date)) {
  temp_data <- BWU_no3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "BWU_NO3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        NO3 = temp_dat$NO3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NO3 = temp_dat$TMR_NO3[2],
        NO3_N_add_int = temp_dat$NO3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(BWU_no3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NO3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * NO3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * NO3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & NO3_N_add_int > 0,
                        Uadd_int / NO3_N_add_int, 
                        NA)
  )

### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("BWU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NO3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  labs(x = expression(paste("TMR nitrate (ugL)")),
       y= expression(paste("Uadd (ug m^-2 min^-1)"))) +
  ggtitle("BWU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")
# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_output_v2.rds")

#########################
#########################
### BWU NH3 ###
BWU_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_TMRout.rds")

## check range in TMR
hist(BWU_nh3_dat_TMR$TMR_Cl)
hist(BWU_nh3_dat_TMR$TMR_NH3)

# Step 2: Spiraling metrics calculation
BWU_nh3_dat_spiraling <- BWU_nh3_dat_TMR %>%
  filter(!sample_no=="NA_H1")%>%
  group_by(date) %>%
  mutate(
    # Log ratio of injectate
    injectate_ratio = (carboy_actual_N_g / NaCl_carboy),
    ln_injectate_ratio = log(carboy_actual_N_g / NaCl_carboy),
    # Log ratio of TMR
    ln_TMR_ratio = log(TMR_NH3 / TMR_Cl),
    log_distance = log(reach_length)
  )

## check range in TMR
hist(BWU_nh3_dat_spiraling$ln_injectate_ratio)
range(BWU_nh3_dat_spiraling$ln_injectate_ratio) # negative 
hist(BWU_nh3_dat_spiraling$ln_TMR_ratio) 
range(na.omit(BWU_nh3_dat_spiraling$ln_TMR_ratio)) 


# Example: Calculate nh3-N_add-int grouped by date
BWU_nh3_dat_spiraling <- BWU_nh3_dat_spiraling %>%
  group_by(date) %>%
  arrange(datetime) %>%  # Ensure data is sorted by time
  mutate(
    # Conservative nh3-N
    nh3_N_cons = Cl_corrected * injectate_ratio,
    
    # Replace negatives with a small positive value (to avoid log-scale issues later)
    nh3_N_add_obs = ifelse(NH3_corrected > 0, NH3_corrected, 0.001),
    
    # Calculate time difference in seconds (replace NA with 0 for the first row)
    dt = as.numeric(difftime(datetime, lag(datetime, default = first(datetime)), units = "secs")),
    
    # Integrals using cumulative trapezoidal approximation
    int_nh3_N_add_obs = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_add_obs),
    int_nh3_N_cons = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt * nh3_N_cons),
    int_Q = cumsum(((lag(offset_Q, default = first(offset_Q)) + offset_Q) / 2) * dt)
  ) %>%
  ungroup() %>%
  mutate(
    # Compute nh3_N_add_int (flow-weighted integration)
    nh3_N_add_int = sqrt(
      (int_nh3_N_add_obs / int_Q) * (int_nh3_N_cons / int_Q)
    )
  )


# Step 3: Slope calculations for uptake kinetics (adapted TASCC method)
out <- data.frame(
  Site = character(), date = as.Date(character()), datetime = as.POSIXct(character()),
  nh3 = numeric(), Cl = numeric(), ln_injectate_ratio = numeric(), offset_Q = numeric(),
  TMR_Cl = numeric(), TMR_NH3 = numeric(), slope_sample = numeric(), kw = numeric(),
  nh3_N_add_int = numeric(),
  stringsAsFactors = FALSE
)

for (d in unique(BWU_nh3_dat_spiraling$date)) {
  temp_data <- BWU_nh3_dat_spiraling %>% filter(date == d)
  
  if (nrow(temp_data) > 1) {
    for (i in 2:nrow(temp_data)) {
      temp_dat <- temp_data[c(i - 1, i), ]
      slope_sample_time <- (temp_dat$ln_TMR_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
        as.numeric(difftime(temp_dat$datetime[2], temp_dat$datetime[1], units = "secs"))
      
      kw <- ((temp_dat$ln_injectate_ratio[2] - temp_dat$ln_TMR_ratio[1]) /
               (-temp_dat$reach_length[2]))
      
      temp_out <- data.frame(
        Site = "BWU_nh3", date = temp_dat$date[2], datetime = temp_dat$datetime[2],
        nh3 = temp_dat$NH3_corrected[2], Cl = temp_dat$Cl_corrected[2],
        ln_injectate_ratio = temp_dat$ln_injectate_ratio[2], offset_Q = temp_dat$offset_Q[2],
        TMR_Cl = temp_dat$TMR_Cl[2], TMR_NH3 = temp_dat$TMR_NH3[2],
        nh3_N_add_int = temp_dat$nh3_N_add_int[2],
        slope_sample = slope_sample_time, kw = kw
      )
      out <- rbind(out, temp_out)
    }
  }
}


# Step 4: Add additional spiraling metrics
data_q <- out %>%
  left_join(BWU_nh3_dat_spiraling %>%
              group_by(date) %>%
              summarize(Cadd = first(na.omit(NH3_corrected)),
                        w = first(na.omit(w))), by = "date") %>%
  mutate(
    sw = 1 / kw, # Spiraling length
    Uadd = offset_Q * 1000 * nh3_N_add_int / (sw * w) # Uptake rate
  )

data_q1 <- data_q %>%
  mutate(
    Uadd_int = ifelse(!is.na(sw) & sw > 0 & w > 0,
                      offset_Q * nh3_N_add_int / (sw * w), 
                      NA),
    Vf_add_int = ifelse(!is.na(Uadd_int) & nh3_N_add_int > 0,
                        Uadd_int / nh3_N_add_int, 
                        NA)
  )


### quick plot:
Uadd_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=Uadd_int*1000)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  #ggtitle("BWU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")


sw_plotGB <- ggplot(data_q1, aes(x=TMR_NH3*1000, y=sw)) +
  #  ylim(-5,30) + xlim(0,40000)+
  geom_point(col = "#a67d17") +
  ggtitle("BWU nitrate") +
  theme_bw() + facet_wrap(date~., scales = "free")

# saveRDS(data_q1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output_v2.rds")

