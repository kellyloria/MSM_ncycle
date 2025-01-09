# Lets start calculating uptake rates
## 24_GBL_NO3_BTC_output.rds

GBL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output_v2.rds")
names(GBL_nh3_datq)
str(GBL_nh3_datq)

GBU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_output.rds")
names(GBU_nh3_datq)
str(GBU_nh3_datq)

BWL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_output.rds")
names(BWL_nh3_datq)
str(BWL_nh3_datq)

BWU_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_output.rds")
names(BWU_nh3_datq)
str(BWU_nh3_datq)


#####################
## Bayes approach 
##=========================================== 
library(dplyr)
library(rstan)

GBL_nh3_dat <- as.data.frame(GBL_nh3_dat)
summary(GBL_nh3_dat)

GBL_nh3_dat_sm <- GBL_nh3_dat%>% filter(Uadd>0 & NH3 > 0.0035)
summary(GBL_nh3_dat_sm)

GBL_nh3_dat_sm <- GBL_nh3_dat_sm%>% filter(date != as.Date("2023-06-01"))
summary(GBL_nh3_dat_sm)

GBL_nh3_dat_sm <- GBL_nh3_datq %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NH3 = 0,  # Set NH3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(GBL_nh3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)           # Reorder the dataset by date and Uadd

# Check the result
summary(GBL_nh3_dat_sm)

#=========================================== 
GBU_nh3_datq <- as.data.frame(GBU_nh3_datq)
summary(GBU_nh3_datq)

GBU_nh3_dat_sm <- GBU_nh3_datq%>% filter(Uadd>0 & NH3 > 0.0035)
summary(GBU_nh3_dat_sm)

GBU_nh3_dat_sm <- GBU_nh3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NH3 = 0,  # Set NH3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(GBU_nh3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)     


#=========================================== 
BWL_nh3_datq <- as.data.frame(BWL_nh3_datq)
summary(BWL_nh3_datq)

BWL_nh3_dat_sm <- BWL_nh3_datq%>% filter(Uadd>0 & NH3 > 0.0035)
summary(BWL_nh3_dat_sm)

BWL_nh3_dat_sm <- BWL_nh3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NH3 = 0,  # Set NH3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(BWL_nh3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)     



#=========================================== 
BWU_nh3_datq <- as.data.frame(BWU_nh3_datq)
summary(BWU_nh3_datq)

BWU_nh3_dat_sm <- BWU_nh3_datq%>% filter(Uadd>0 & NH3 > 0.0035)
summary(BWU_nh3_dat_sm)

BWU_nh3_dat_sm <- BWU_nh3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NH3 = 0,  # Set NH3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(BWU_nh3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)     

GBL_nh3_datq$Uadd_int

GBL_nh3_datq$Uadd_int1 <- ifelse(is.na(GBL_nh3_datq$Uadd_int), GBL_nh3_datq$Uadd_int, 0.000001)
GBL_nh3_datq$Uadd_int1 <- ifelse(is.na(GBL_nh3_datq$Uadd_int) | is.nan(GBL_nh3_datq$Uadd_int), 0.000001, GBL_nh3_datq$Uadd_int)


# Prepare data for Stan
data_stan <- GBL_nh3_datq %>%
  mutate(date_idx = as.integer(as.factor(date))) %>%
  list(
    N = nrow(.),
    D = n_distinct(.$date_idx),
    date = .$date_idx,
    TMR_NH3 = .$TMR_NH3,
    Uadd = .$Uadd_int1
  )

library(rstan)
library(StanHeaders)
# Compile and fit the Stan model
MM_model <- "/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/MichaelisMentenCurve_nh3.stan"  # Save the Stan code in this file
fit <- stan(
  file = MM_model,
  data = data_stan,
  iter = 5000, chains = 3,
  warmup = 2500, thin = 1,
  control = list(adapt_delta = 0.95))

# Extract and visualize results
print(fit, pars = c("mu_Vmax", "mu_Km", "sigma", "Vmax", "Km"))

plot(fit)

#####################
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)  # For working with posterior samples
library(bayesplot)  # For posterior predictive checks

# Extract posterior samples
posterior <- as_draws_df(fit)

# Extract relevant parameter samples
param_samples <- posterior %>%
  select(contains("Vmax"), contains("Km"), sigma)

# Prepare a dataframe for predictions by date
data_stan_df <- as.data.frame(data_stan)

prediction_data <- data_stan_df %>%
  mutate(
    date_idx = as.factor(date),
    TMR_NH3_pred = NA_real_,
    TMR_NH3_lower = NA_real_,
    TMR_NH3_upper = NA_real_
  )

# Function to calculate TMR_NH3 predictions
predict_TMR_NH3 <- function(Uadd, Vmax, Km, sigma) {
  predicted <- (Vmax * Uadd) / (Km + Uadd)
  predicted + rnorm(length(predicted), mean = 0, sd = sigma)
}

# Loop over each date to calculate posterior predictions
posterior_predictions <- lapply(1:length(unique(data_stan$date)), function(d) {
  date_filter <- data_stan$date == d
  Uadd_values <- data_stan$Uadd[date_filter]
  
  # Extract posterior samples for the date
  Vmax_samples <- posterior %>% select(contains(paste0("Vmax[", d, "]"))) %>% unlist()
  Km_samples <- posterior %>% select(contains(paste0("Km[", d, "]"))) %>% unlist()
  sigma_samples <- posterior %>% pull(sigma)
  
  # Calculate predicted TMR_NH3 for each posterior draw
  predicted_samples <- sapply(1:length(Vmax_samples), function(i) {
    predict_TMR_NH3(Uadd_values, Vmax_samples[i], Km_samples[i], sigma_samples[i])
  })
  
  # Summarize predictions (mean and 95% credible intervals)
  tibble(
    Uadd = Uadd_values,
    TMR_NH3_mean = rowMeans(predicted_samples),
    TMR_NH3_lower = apply(predicted_samples, 1, quantile, probs = 0.025),
    TMR_NH3_upper = apply(predicted_samples, 1, quantile, probs = 0.975),
    date_idx = d
  )
}) %>% bind_rows()


# Calculate RMSE for each date
rmse_results <- posterior_predictions %>%
  left_join(data_stan_df, by = c("date_idx","Uadd"="Uadd_int1")) %>%
  group_by(date_idx) %>%
  summarize(
    RMSE = sqrt(mean((TMR_NH3 - TMR_NH3_mean)^2, na.rm = TRUE))
  )

# Combine data for plotting
all_results <- posterior_predictions %>%
  left_join(data_stan_df, by = c("date_idx", "Uadd"="Uadd_int1"))

# Plot
# all_results%>%filter(date_idx!="8")
plot <- ggplot(all_results, aes(y = Uadd*1000, x = TMR_NH3*1000)) +
  geom_point(aes(color = "Observed")) +
  geom_line(aes(color = "Observed")) +
  geom_point(aes(y = TMR_NH3_mean *1000, color = "Predicted")) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = TMR_NH3_mean*1000, color = "Predicted")) +
  geom_errorbar(aes(ymin = TMR_NH3_lower *1000, ymax = TMR_NH3_upper*1000), alpha = 0.3, color = "grey50") +
  #geom_ribbon(aes(ymin = TMR_NH3_lower*1000, ymax = TMR_NH3_upper*1000), alpha = 0.3, fill = "#grey50") +
  facet_wrap(~date_idx, scales = "free") +
  labs(
    title = "GBL KNH3",
    y = "Uadd (ug m^2 min^-1)",
    x = "TMR NH3 (ug L^-1)"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c("Observed" = "black", "Predicted" = "grey50")
  ) +
  geom_text(
    data = rmse_results,
    aes(x = Inf, y = Inf, label = paste("RMSE:", round(RMSE, 2))),
    hjust = 1.1, vjust = 1.1,
    size = 3,
    color = "black",
    inherit.aes = FALSE
  )

# Display plot
print(plot)

# ggsave(plot = plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBL_NH3_uadd_plotv2.png",sep=""),width=9.5,height=6.5,dpi=300)
# 
saveRDS(all_results, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/fits/24_GBL_NH3_BTC_TMRfit.rds")
write_csv(all_results, "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/fits/24_GBL_NH3_BTC_TMRfit.cvs")










