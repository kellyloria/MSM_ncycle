# Lets start calculating uptake rates
## 24_GBL_NO3_BTC_output.rds

GBL_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_output.rds")
names(GBL_no3_datq)
str(GBL_no3_datq)

GBU_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_output.rds")
names(GBU_no3_datq)
str(GBU_no3_datq)

BWL_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_output.rds")
names(BWL_no3_datq)
str(BWL_no3_datq)

BWU_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_output.rds")
names(BWU_no3_datq)
str(BWU_no3_datq)


GBL_nh3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_output.rds")
names(GBL_nh3_datq)
str(GBL_nh3_datq)

unique(GBL_no3_datq$date)

GBL_no3_dat <- na.omit(GBL_no3_datq)
GBL_test <- GBL_no3_dat%>%filter(date==as.Date("2023-07-10"))
summary(GBL_test)

# not working:
model.drmt <- drm(scale(Uadd) ~ scale(NO3), data = GBL_test, fct = MM.2(), 
                  start = c(Vmax = 0.0008, Km = 0.001))

# works 
model.nls <- nls(Uadd ~ (Vmax * NO3) / (Km + NO3), data = GBL_test, 
                 start = list(Vmax = 0.0008, Km = 0.001))

# quick visualization:
Uadd_plotGB <- ggplot(GBL_no3_datq, aes(x=TMR_NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 1, shape=15, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") + facet_wrap(.~date) 
  theme_bw()

#####################
#####################
## Bayes approach 
##=========================================== 
library(dplyr)
library(rstan)

GBL_no3_dat <- as.data.frame(GBL_no3_dat)
summary(GBL_no3_dat)

GBL_no3_dat_sm <- GBL_no3_dat%>% filter(Uadd>0 & NO3 > 0.0035)
summary(GBL_no3_dat_sm)

GBL_no3_dat_sm <- GBL_no3_dat_sm%>% filter(date != as.Date("2023-06-01"))
summary(GBL_no3_dat_sm)

GBL_no3_dat_sm <- GBL_no3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NO3 = 0,  # Set NO3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(GBL_no3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)           # Reorder the dataset by date and Uadd

# Check the result
summary(GBL_no3_dat_sm)

#=========================================== 
GBU_no3_datq <- as.data.frame(GBU_no3_datq)
summary(GBU_no3_datq)

GBU_no3_dat_sm <- GBU_no3_datq%>% filter(Uadd>0 & NO3 > 0.0035)
summary(GBU_no3_dat_sm)

GBU_no3_dat_sm <- GBU_no3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NO3 = 0,  # Set NO3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(GBU_no3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)     


#=========================================== 
BWL_no3_datq <- as.data.frame(BWL_no3_datq)
summary(BWL_no3_datq)

BWL_no3_dat_sm <- BWL_no3_datq%>% filter(Uadd>0 & NO3 > 0.0035)
summary(BWL_no3_dat_sm)

BWL_no3_dat_sm <- BWL_no3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NO3 = 0,  # Set NO3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(BWL_no3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)     



#=========================================== 
BWU_no3_datq <- as.data.frame(BWU_no3_datq)
summary(BWU_no3_datq)

BWU_no3_dat_sm <- BWU_no3_datq%>% filter(Uadd>0 & NO3 > 0.0035)
summary(BWU_no3_dat_sm)

BWU_no3_dat_sm <- BWU_no3_dat_sm %>%
  group_by(date) %>%
  arrange(date, Uadd) %>% # Arrange within each group by Uadd
  summarise(
    TMR_NO3 = 0,  # Set NO3_TMR to 0
    Uadd = 0,     # Set Uadd to 0
    .groups = "drop"
  ) %>%
  bind_rows(BWU_no3_dat_sm) %>% # Combine the original dataset with the new rows
  arrange(date, Uadd)     


# Prepare data for Stan
data_stan <- BWU_no3_dat_sm %>%
  mutate(date_idx = as.integer(as.factor(date))) %>%
  list(
    N = nrow(.),
    D = n_distinct(.$date_idx),
    date = .$date_idx,
    TMR_NO3 = .$TMR_NO3,
    Uadd = .$Uadd
  )

library(rstan)
library(StanHeaders)
# Compile and fit the Stan model
MM_model <- "/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/MichaelisMentenCurve.stan"  # Save the Stan code in this file
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
    TMR_NO3_pred = NA_real_,
    TMR_NO3_lower = NA_real_,
    TMR_NO3_upper = NA_real_
  )

# Function to calculate TMR_NO3 predictions
predict_TMR_NO3 <- function(Uadd, Vmax, Km, sigma) {
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
  
  # Calculate predicted TMR_NO3 for each posterior draw
  predicted_samples <- sapply(1:length(Vmax_samples), function(i) {
    predict_TMR_NO3(Uadd_values, Vmax_samples[i], Km_samples[i], sigma_samples[i])
  })
  
  # Summarize predictions (mean and 95% credible intervals)
  tibble(
    Uadd = Uadd_values,
    TMR_NO3_mean = rowMeans(predicted_samples),
    TMR_NO3_lower = apply(predicted_samples, 1, quantile, probs = 0.025),
    TMR_NO3_upper = apply(predicted_samples, 1, quantile, probs = 0.975),
    date_idx = d
  )
}) %>% bind_rows()


# Calculate RMSE for each date
rmse_results <- posterior_predictions %>%
  left_join(data_stan_df, by = c("date_idx", "Uadd")) %>%
  group_by(date_idx) %>%
  summarize(
    RMSE = sqrt(mean((TMR_NO3 - TMR_NO3_mean)^2, na.rm = TRUE))
  )

# Combine data for plotting
all_results <- posterior_predictions %>%
  left_join(data_stan_df, by = c("date_idx", "Uadd"))

# Plot
# all_results%>%filter(date_idx!="8")
plot <- ggplot(all_results, aes(y = Uadd*1000, x = TMR_NO3*1000)) +
  geom_point(aes(color = "Observed")) +
  geom_line(aes(color = "Observed")) +
  geom_point(aes(y = TMR_NO3_mean *1000, color = "Predicted")) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = TMR_NO3_mean*1000, color = "Predicted")) +
  geom_errorbar(aes(ymin = TMR_NO3_lower *1000, ymax = TMR_NO3_upper*1000), alpha = 0.3, color = "grey50") +
  #geom_ribbon(aes(ymin = TMR_NO3_lower*1000, ymax = TMR_NO3_upper*1000), alpha = 0.3, fill = "#grey50") +
  facet_wrap(~date_idx, scales = "free") +
  labs(
    title = "BWU KNO3",
    y = "Uadd (ug m^2 min^-1)",
    x = "TMR NO3 (ug L^-1)"
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

# ggsave(plot = plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBL_NO3_uadd_plotv2.png",sep=""),width=9.5,height=6.5,dpi=300)
# 
saveRDS(all_results, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/fits/24_BWU_NO3_BTC_TMRfit.rds")
write_csv(all_results, "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/fits/24_BWU_NO3_BTC_TMRfit.cvs")

















# GB
GBL_06 <- read.csv("./BTC_out/GBL_BTC_20210603.csv")

model.drm1 <- drm (Uadd ~ NO3, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()



# GB
GBL_06 <- read.csv("./BTC_out/GBL_NO3_BTC_220623.csv")

model.drm1 <- drm (Uadd ~ NO3, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()





# GB
GBL_06 <- read.csv("./BTC_out/GBL_BTC_NH4_20220623.csv")

model.drm1 <- drm (Uadd ~ NH4, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_06, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_07 <- read.csv("./BTC_out/GB_BTC_20210722.csv")
summary(GBL_07)

model.drm1 <- drm (Uadd ~ NO3, data = GBL_07a, fct = MM.2())
summary(model.drm1)



summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_07, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
 # geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()




GBL_12 <- read.csv("./BTC_out/GBL_BTC_221212.csv")
summary(GBL_12)
GBL_12a<-GBL_12[c(-2),]

model.drm1 <- drm (Uadd ~ NH4, data = GBL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_12a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_10 <- read.csv("./BTC_out/GBL_NO3_BTC_221003.csv")
summary(GBL_10)
GBL_10a<-GBL_10[c(-1),]

model.drm1 <- drm (Uadd ~ NO3, data = GBL_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_10a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_10a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_10 <- read.csv("./BTC_out/GBL_BTC_NH4_20221003.csv")
GBL_10a<-GBL_10[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NH4, data = GBL_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_10a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()

#############
GBL_11 <- read.csv("./BTC_out/GBL_BTC_NH4_20221104.csv")
GBL_11a<-GBL_11[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NH4, data = GBL_11a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_11a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()

################
GBL_12 <- read.csv("./BTC_out/GBL_NO3_BTC_221212.csv")
GBL_10a<-GBL_10[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NO3, data = GBL_12, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_12$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_12, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_04 <- read.csv("./BTC_out/GBL_NO3_BTC_220407.csv")
summary(GBL_04)

model.drm1 <- drm (Uadd ~ NO3, data = GBL_04, fct = MM.2())
summary(model.drm1)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_04, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_04 <- read.csv("./BTC_out/GBL_BTC_NH4_20220407.csv")
summary(GBL_04)

model.drm1 <- drm (Uadd ~ NH4, data = GBL_04, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_04, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()



GBL_12 <- read.csv("./BTC_out/GBL_NO3_BTC_221212.csv")
summary(GBL_12)

GBL_12a<-GBL_12[c(-2),]

model.drm1 <- drm (Uadd ~ NO3, data = GBL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_12a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_12a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_10 <- read.csv("./BTC_out/GBL_NO3_BTC_BWL221003.csv")
summary(GBL_10)
GBL_10a<-na.omit(GBL_10[c(-1,-2),])

model.drm1 <- drm (Uadd ~ NO3, data = GBL_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_10a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_10a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()




GBU_10 <- read.csv("./BTC_out/GBU_BTC_20221003_MW.csv")
summary(GBU_10)

GBU_10a<-na.omit(GBU_10[c(-11,-13),])

model.drm1 <- drm (Uadd ~ NH4, data = GBU_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_10a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()



GBU_11 <- read.csv("./BTC_out/GBU_BTC_20221104_MW.csv")
summary(GBU_11)

GBU_10a<-na.omit(GBU_10[c(-11,-13),])

model.drm1 <- drm (Uadd ~ NH4, data = GBU_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_11, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()





GBU_06 <- read.csv("./BTC_out/GBU_BTC_NH4_20220623.csv")
summary(GBU_06)
GBU_06a<-na.omit(GBU_06[c(-2),])

model.drm1 <- drm (Uadd ~ NH4, data = GBU_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_06$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


Uadd_plotGB <- ggplot(GBU_06, aes(x=log(NH4+1), y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()





GBU_04 <- read.csv("./BTC_out/GBU_NO3_BTC_220407.csv")
summary(GBU_04)

model.drm1 <- drm (Uadd ~ NH4, data = GBU_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_04, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()



GBU_06 <- read.csv("./BTC_out/GBU_NO3_BTC_220623.csv")
summary(GBU_06)

model.drm1 <- drm (Uadd ~ NO3, data = GBU_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()


# GB
GBL_03 <- read.csv("./BTC_out/GBL_BTC_NH4_20230327.csv")

model.drm1 <- drm (Uadd ~ NH4, data = GBL_03, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_03, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_03 <- read.csv("./BTC_out/GBL_NO3_BTC_230327.csv")
summary(GBL_03)
GBL_10a<-na.omit(GBL_10[c(-1,-2),])

model.drm1 <- drm (Uadd ~ NO3, data = GBL_03, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_10a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_03, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  ggtitle("Glenbrook L Creek") +
  theme_bw()












###
###
###

# BWL
BWL_05 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL220526.csv")

model.drm1 <- drm (Uadd ~ NO3, data = BWL_05, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_05$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_05, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

BWL_05 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL220526.csv")
BWL_05a<-BWL_05[c(-8),]

model.drm1 <- drm (Uadd ~ NH4, data = BWL_05a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_05a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_05a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_08 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL220824.csv")
#BWL_05a<-BWL_05[c(-8),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_08, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_08$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_08, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_10 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221012.csv")
#BWL_05a<-BWL_05[c(-8),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_10, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_10$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_10, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_10 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL221012.csv")
BWL_10a<-BWL_10[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_10, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_10$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_10, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

###

BWL_11 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221121v2.csv")
BWL_11a<-BWL_11[c(-1, -2, -3),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_11a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_11a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_11a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

BWL_11 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL221121.csv")
BWL_11a<-BWL_11[c(-1,-2,-15, -14),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_11a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_11a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_11a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

###
BWL_12 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221219.csv")
BWL_12a<-BWL_12[c(-1,-2,-3,-4,-18),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_12a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_12 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL221219.csv")
BWL_12a<-BWL_12[c(-1,-2,-14,-15),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_12a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_12a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

###

BWL_02 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL230215.csv")
BWL_12a<-BWL_12[c(-1,-2,-3,-4,-18),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_02, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_02, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_02 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL230215.csv")
BWL_12a<-BWL_12[c(-1,-2,-14,-15),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_02, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_02$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_02, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

####
BWL_04 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL230405.csv")
BWL_12a<-BWL_12[c(-1,-2,-3,-4,-18),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_04, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_04$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_04, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_04 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL230405.csv")

model.drm1 <- drm (Uadd ~ NO3, data = BWL_04, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_02$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_04, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

######
#####
BWL_07 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL230718.csv")
model.drm1 <- drm (Uadd ~ NH4, data = BWL_07, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_04$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_07, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

BWL_07 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL230718.csv")
BWL_07a<-BWL_07[c(-12),]


model.drm1 <- drm (Uadd ~ NO3, data = BWL_07a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_07a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_07a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

