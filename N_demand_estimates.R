

# Load necessary library
library(dplyr)

# Example input dataframe with GPP and ER values (replace with actual data)
data <- data.frame(
  date = seq(as.Date("2024-01-01"), as.Date("2024-01-10"), by = "day"),
  GPP = runif(10, 1, 10),  # Daily gross primary production (e.g., mmol C m^-2 day^-1)
  ER = runif(10, -10, -1)  # Daily ecosystem respiration (e.g., mmol C m^-2 day^-1)
)

# Parameters
ra_coeff <- 0.5          # Autotrophic respiration coefficient (Hall and Tank 2003)
auto_CN_ratio <- 16      # Autotrophic molar C:N ratio for low N (Stelzer and Lamberti 2001)
het_growth_eff <- 0.05    # Heterotrophic growth efficiency for energy limited systems (Hall and Tank 2003)
het_CN_ratio <- 20       # Heterotrophic molar C:N ratio (Hall and Tank 2003)

# Calculate nitrate demand
data <- data %>%
  mutate(
    # Autotrophic respiration (raGPP)
    raGPP = ra_coeff * GPP,
    
    # Gross autotrophic nitrate assimilation (mmol N m^-2 day^-1)
    auto_NO3_assim = GPP / auto_CN_ratio,
    
    # Heterotrophic respiration (Rh)
    Rh = ER - raGPP,
    
    # Gross heterotrophic nitrate assimilation (mmol N m^-2 day^-1)
    het_NO3_assim = (Rh * het_growth_eff) / het_CN_ratio,
    
    # Total nitrate demand
    total_NO3_demand = auto_NO3_assim + het_NO3_assim
  )

# Display results
print(data)

# Optionally, write results to a CSV file
# write.csv(data, "nitrate_demand_results.csv", row.names = FALSE)
