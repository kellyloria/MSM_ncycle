###
###
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(unitted)
library(lubridate)

#############################
GBL_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_input.rds")
names(GBL_no3_datq)

GBL_no3_dat_TMR <- GBL_no3_datq %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
    integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
    integrated_NO3 = c(0, cumsum((NO3_corrected[-1] + NO3_corrected[-n()]) / 2 * diff(cumulative_time))),
    TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
    TMR_NO3 = offset_Q[1] * integrated_NO3) %>% ungroup()

names(GBL_no3_dat_TMR)

### estimate fractional mass recovery:
GBL_no3_dat_TMR1 <- GBL_no3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(TMR_NO3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NO3/TMR_Cl))
    
hist(GBL_no3_dat_TMR1$F_TMR_NO3)
hist(GBL_no3_dat_TMR1$F_TMR_Cl)
hist(GBL_no3_dat_TMR1$N_Cl)

###  Quick plot with a second y-axis
scale_factor = 1
GBL_no3_plot <- ggplot(GBL_no3_dat_TMR1, aes(x = cumulative_time)) +
  geom_point(aes(y = TMR_NO3), col = "red") +
  geom_line(aes(y = TMR_NO3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date) +
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
GBL_no3_plot

GBL_no3_dat_plot <- GBL_no3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(F_TMR_NO3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

hist(GBL_no3_dat_plot$inject_ratio)
# Display the plot
# Quick plot with a second y-axis
scale_factor = 1
GBL_no3_FTMR_plot <- ggplot(GBL_no3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NO3), col = "red") +
  geom_line(aes(y = F_TMR_NO3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "GBL KNO3",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
GBL_no3_FTMR_plot

# ggsave(plot = GBL_no3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBL_FTMR_NO3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 10       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
GBL_NO3_BTC_plot <- ggplot(GBL_no3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NO3_corrected), col = "red") + 
  geom_line(aes(y = NO3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NO3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "GBL KNO3",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
GBL_NO3_BTC_plot

# ggsave(plot = GBL_NO3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBL_BTC_NO3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)
## saveRDS(GBL_no3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NO3_BTC_TMRout.rds")



##############################
## read in data and plot TMR 
GBL_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_input.rds")%>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NH3 = c(0, cumsum((NH3_corrected[-1] + NH3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NH3 = offset_Q[1] * integrated_NH3) %>% ungroup()

names(GBL_nh3_dat_TMR)

### estimate fractional mass recovery:
GBL_nh3_dat_TMR1 <- GBL_nh3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(TMR_NH3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NH3/TMR_Cl))

hist(GBL_nh3_dat_TMR1$F_TMR_NH3)
hist(GBL_nh3_dat_TMR1$F_TMR_Cl)
hist(GBL_nh3_dat_TMR1$N_Cl)


GBL_nh3_dat_plot <- GBL_nh3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(F_TMR_NH3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

scale_factor = 1
GBL_nh3_FTMR_plot <- ggplot(GBL_nh3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NH3), col = "red") +
  geom_line(aes(y = F_TMR_NH3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "GBL Nh3Cl",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
GBL_nh3_FTMR_plot

# # ggsave(plot = GBL_nh3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBL_FTMR_NH3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 10       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
GBL_NH3_BTC_plot <- ggplot(GBL_nh3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NH3_corrected), col = "red") + 
  geom_line(aes(y = NH3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NH3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "GBL Nh3Cl",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
GBL_NH3_BTC_plot

# ggsave(plot = GBL_NH3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBL_BTC_NH3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)

## saveRDS(GBL_nh3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBL_NH3_BTC_TMRout.rds")

    
##############################

#############################
GBU_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_input.rds")
names(GBU_no3_datq)

GBU_no3_dat_TMR <- GBU_no3_datq %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NO3 = c(0, cumsum((NO3_corrected[-1] + NO3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NO3 = offset_Q[1] * integrated_NO3) %>% ungroup()

names(GBU_no3_dat_TMR)

### estimate fractional mass recovery:
GBU_no3_dat_TMR1 <- GBU_no3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(TMR_NO3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NO3/TMR_Cl))

hist(GBU_no3_dat_TMR1$F_TMR_NO3)
hist(GBU_no3_dat_TMR1$F_TMR_Cl)
hist(GBU_no3_dat_TMR1$N_Cl)

GBU_no3_dat_plot <- GBU_no3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(F_TMR_NO3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))
summary(GBU_no3_dat_plot)

hist(GBU_no3_dat_plot$inject_ratio)
# Display the plot
# Quick plot with a second y-axis
scale_factor = 1
GBU_no3_FTMR_plot <- ggplot(GBU_no3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NO3), col = "red") +
  geom_line(aes(y = F_TMR_NO3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "GBU KNO3",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
GBU_no3_FTMR_plot

# ggsave(plot = GBU_no3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBU_FTMR_NO3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 10       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
GBU_NO3_BTC_plot <- ggplot(GBU_no3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NO3_corrected), col = "red") + 
  geom_line(aes(y = NO3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NO3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "GBU KNO3",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
GBU_NO3_BTC_plot

# ggsave(plot = GBU_NO3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBU_BTC_NO3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)

## saveRDS(GBU_no3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NO3_BTC_TMRout.rds")



##############################
###############################
## read in data and plot TMR 
GBU_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_input.rds")%>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NH3 = c(0, cumsum((NH3_corrected[-1] + NH3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NH3 = offset_Q[1] * integrated_NH3) %>% ungroup()

names(GBU_nh3_dat_TMR)

### estimate fractional mass recovery:
GBU_nh3_dat_TMR1 <- GBU_nh3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(TMR_NH3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NH3/TMR_Cl))

hist(GBU_nh3_dat_TMR1$F_TMR_NH3)
hist(GBU_nh3_dat_TMR1$F_TMR_Cl)
hist(GBU_nh3_dat_TMR1$N_Cl)


GBU_nh3_dat_plot <- GBU_nh3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(F_TMR_NH3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

scale_factor = 1
GBU_nh3_FTMR_plot <- ggplot(GBU_nh3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NH3), col = "red") +
  geom_line(aes(y = F_TMR_NH3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "GBU Nh3Cl",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
GBU_nh3_FTMR_plot

## ggsave(plot = GBU_nh3_dat_TMR1, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBU_FTMR_NH3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 10       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
GBU_NH3_BTC_plot <- ggplot(GBU_nh3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NH3_corrected), col = "red") + 
  geom_line(aes(y = NH3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NH3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "GBU Nh3Cl",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
GBU_NH3_BTC_plot

# ggsave(plot = GBU_NH3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/GBU_BTC_NH3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)

## saveRDS(GBU_nh3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_GBU_NH3_BTC_TMRout.rds")


##############################
###############################


#### BWL 
#############################
BWL_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_input.rds")
names(BWL_no3_datq)

BWL_no3_dat_TMR <- BWL_no3_datq %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NO3 = c(0, cumsum((NO3_corrected[-1] + NO3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NO3 = offset_Q[1] * integrated_NO3) %>% ungroup()

names(BWL_no3_dat_TMR)

### estimate fractional mass recovery:
BWL_no3_dat_TMR1 <- BWL_no3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(TMR_NO3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NO3/TMR_Cl))

hist(BWL_no3_dat_TMR1$F_TMR_NO3)
hist(BWL_no3_dat_TMR1$F_TMR_Cl)
hist(BWL_no3_dat_TMR1$N_Cl)

BWL_no3_dat_plot <- BWL_no3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(F_TMR_NO3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

hist(BWL_no3_dat_plot$inject_ratio)
# Display the plot
# Quick plot with a second y-axis
scale_factor = 1
BWL_no3_FTMR_plot <- ggplot(BWL_no3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NO3), col = "red") +
  geom_line(aes(y = F_TMR_NO3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "BWL KNO3",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
BWL_no3_FTMR_plot

# ggsave(plot = BWL_no3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWL_FTMR_NO3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 1       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
BWL_NO3_BTC_plot <- ggplot(BWL_no3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NO3_corrected), col = "red") + 
  geom_line(aes(y = NO3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NO3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "BWL KNO3",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
BWL_NO3_BTC_plot
# ggsave(plot = BWL_NO3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWL_BTC_NO3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)

## saveRDS(BWL_no3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NO3_BTC_TMRout.rds")

#############################
###
##############################
## BWL 
BWL_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_input.rds")%>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NH3 = c(0, cumsum((NH3_corrected[-1] + NH3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NH3 = offset_Q[1] * integrated_NH3) %>% ungroup()

names(BWL_nh3_dat_TMR)

### estimate fractional mass recovery:
BWL_nh3_dat_TMR1 <- BWL_nh3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(TMR_NH3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NH3/TMR_Cl))

hist(BWL_nh3_dat_TMR1$F_TMR_NH3)
hist(BWL_nh3_dat_TMR1$F_TMR_Cl)
hist(BWL_nh3_dat_TMR1$N_Cl)


BWL_nh3_dat_plot <- BWL_nh3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(F_TMR_NH3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

scale_factor = 1
BWL_nh3_FTMR_plot <- ggplot(BWL_nh3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NH3), col = "red") +
  geom_line(aes(y = F_TMR_NH3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "BWL Nh3Cl",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
BWL_nh3_FTMR_plot

# # ggsave(plot = BWL_nh3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWL_FTMR_NH3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 1       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
BWL_NH3_BTC_plot <- ggplot(BWL_nh3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NH3_corrected), col = "red") + 
  geom_line(aes(y = NH3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NH3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "BWL Nh3Cl",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
BWL_NH3_BTC_plot

# ggsave(plot = BWL_NH3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWL_BTC_NH3_plot.png",sep=""),width=9.5,height=6.5,dpi=300)

## saveRDS(BWL_nh3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWL_NH3_BTC_TMRout.rds")

#############################
#### BWU 
#############################
BWU_no3_datq <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_input.rds")
names(BWU_no3_datq)

BWU_no3_dat_TMR <- BWU_no3_datq %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NO3 = c(0, cumsum((NO3_corrected[-1] + NO3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NO3 = offset_Q[1] * integrated_NO3) %>% ungroup()

names(BWU_no3_dat_TMR)

### estimate fractional mass recovery:
BWU_no3_dat_TMR1 <- BWU_no3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(TMR_NO3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NO3/TMR_Cl))

hist(BWU_no3_dat_TMR1$F_TMR_NO3)
hist(BWU_no3_dat_TMR1$F_TMR_Cl)
hist(BWU_no3_dat_TMR1$N_Cl)

BWU_no3_dat_plot <- BWU_no3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NO3 = c(F_TMR_NO3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

hist(BWU_no3_dat_plot$inject_ratio)
# Display the plot
# Quick plot with a second y-axis
scale_factor = 1
BWU_no3_FTMR_plot <- ggplot(BWU_no3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NO3), col = "red") +
  geom_line(aes(y = F_TMR_NO3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "BWU KNO3",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
BWU_no3_FTMR_plot

# ggsave(plot = BWU_no3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWU_FTMR_NO3_plot.png",sep=""),width=9.5,height=3.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 1       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
BWU_NO3_BTC_plot <- ggplot(BWU_no3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NO3_corrected), col = "red") + 
  geom_line(aes(y = NO3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NO3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "BWU KNO3",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
BWU_NO3_BTC_plot
# ggsave(plot = BWU_NO3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWU_BTC_NO3_plot.png",sep=""),width=9.5,height=3.5,dpi=300)

## saveRDS(BWU_no3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NO3_BTC_TMRout.rds")


###
##############################
## BWU 
BWU_nh3_dat_TMR <- readRDS("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_input.rds")%>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(cumulative_time = as.numeric(datetime - min(datetime, na.rm = TRUE), units = "secs"),
         integrated_Cl = c(0, cumsum((Cl_corrected[-1] + Cl_corrected[-n()]) / 2 * diff(cumulative_time))),
         integrated_NH3 = c(0, cumsum((NH3_corrected[-1] + NH3_corrected[-n()]) / 2 * diff(cumulative_time))),
         TMR_Cl = offset_Q[1] * integrated_Cl, # Assuming consistent discharge
         TMR_NH3 = offset_Q[1] * integrated_NH3) %>% ungroup()

names(BWU_nh3_dat_TMR)

### estimate fractional mass recovery:
BWU_nh3_dat_TMR1 <- BWU_nh3_dat_TMR %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(TMR_NH3/nitrogen_carboy),
         F_TMR_Cl = c(TMR_Cl/NaCl_carboy),
         N_Cl = c(TMR_NH3/TMR_Cl))

hist(BWU_nh3_dat_TMR1$F_TMR_NH3)
hist(BWU_nh3_dat_TMR1$F_TMR_Cl)
hist(BWU_nh3_dat_TMR1$N_Cl)


BWU_nh3_dat_plot <- BWU_nh3_dat_TMR1 %>%
  group_by(date) %>% arrange(datetime) %>% 
  mutate(F_TMR_NH3 = c(F_TMR_NH3*100),
         F_TMR_Cl = c(F_TMR_Cl*100),
         inject_ratio= c((nitrogen_carboy/NaCl_carboy)*100))

scale_factor = 1
BWU_nh3_FTMR_plot <- ggplot(BWU_nh3_dat_plot, aes(x = cumulative_time)) +
  geom_point(aes(y = F_TMR_NH3), col = "red") +
  geom_line(aes(y = F_TMR_NH3), col = "red") +
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  geom_point(aes(y = F_TMR_Cl / scale_factor), col = "grey25") + 
  scale_y_continuous(name = "Fractional N recovered", 
                     limits = c(0, 1),  # Set limits for both primary and secondary y-axes
                     sec.axis = sec_axis(~ . * scale_factor, name = "Fractional Cl recovered")) +
  theme_bw() + facet_wrap(~ date, scales = "free") +
  labs(
    title = "BWU Nh3Cl",
    x = "cumulative time (s)"
  ) + 
  theme(axis.title.y = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "grey25"))
# Display the plot
BWU_nh3_FTMR_plot

# # ggsave(plot = BWU_nh3_FTMR_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWU_FTMR_NH3_plot.png",sep=""),width=9.5,height=5.5,dpi=300)


###### FIG 4 in covino et al 2018
scale_factor_TMR_Cl <- 1       # Scaling factor for TMR_Cl
scale_factor_N_Cl <- 0.01      # Scaling factor for N_Cl

# Plot
BWU_NH3_BTC_plot <- ggplot(BWU_nh3_dat_plot, aes(x = cumulative_time)) + 
  # Primary y-axis: TMR_NH3
  geom_point(aes(y = NH3_corrected), col = "red") + 
  geom_line(aes(y = NH3_corrected), col = "red") + 
  # Secondary y-axis: TMR_Cl
  geom_line(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  geom_point(aes(y = Cl_corrected / scale_factor_TMR_Cl), col = "grey25") + 
  # Additional "third axis" for N_Cl
  geom_point(aes(y = N_Cl / scale_factor_N_Cl), col = "blue", alpha=0.5) + 
  # Add horizontal line for inject_ratio
  geom_hline(aes(yintercept = inject_ratio), linetype = "dashed", color = "blue") + 
  # Axis customization
  scale_y_continuous(
    name = "NH3 (mg/L)",                              # Primary y-axis label
    sec.axis = sec_axis(
      ~ . * scale_factor_TMR_Cl,                  # Transformation for TMR_Cl
      name = "estimated Cl (mg/L)")) + 
  # Manually annotate third axis in plot legend
  guides(color = guide_legend(title = NULL)) + 
  theme_bw() + facet_wrap(~ date, scales = "free") + 
  # Adjust axis and legend appearance
  theme(
    axis.title.y = element_text(color = "red"),       # Primary y-axis color
    axis.title.y.right = element_text(color = "grey25"), # Secondary y-axis color
    legend.position = "bottom") + 
  labs(
    title = "BWU Nh3Cl",
    x = "cumulative time (s)"
  )+
  # Add custom legend for N_Cl
  labs(color = NULL) + 
  scale_color_manual(
    values = c("TMR_NH3" = "red", "TMR_Cl" = "grey25", "N_Cl" = "blue"),
    labels = c("TMR_NH3", "TMR_Cl", "N_Cl (scaled)"))
# Display the plot
BWU_NH3_BTC_plot

# ggsave(plot = BWU_NH3_BTC_plot, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/figures/BWU_BTC_NH3_plot.png",sep=""),width=9.5,height=5.5,dpi=300)


## saveRDS(BWU_nh3_dat_TMR1, file = "/Users/kellyloria/Documents/UNR/Ncycle/2024_workflow/24_BWU_NH3_BTC_TMRout.rds")































rawdat = list.files(paste("./BTC_out/",sep=""), full.names = T)


# Lets start calculating uptake rates

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

