# Load required packages
library(lavaan)
library(semPlot)

# Load required libraries
library(lavaan)
library(semPlot)

# Define the SEM model
model <- "
  # Regression equations
  Nup ~ nitrogen + NEPt + flow + SPC + wtemp
  NEPt ~ NEPt_1
"

# Simulate data
set.seed(123)

# Number of observations
n <- 200

# Simulated variables
sampletype <- factor(sample(c("Type1", "Type2"), n, replace = TRUE)) # Categorical predictor
nitrogen <- rnorm(n, mean = ifelse(sampletype == "Type1", 5, 7), sd = 1)  # Dependent on sampletype
NEPt_1 <- rnorm(n, mean = 3, sd = 1)  # Lagged NEP variable
NEPt <- 0.7 * NEPt_1 + rnorm(n, sd = 0.5)  # Current NEP based on lagged NEP
flow <- rnorm(n, mean = 10, sd = 2)  # Flow
SPC <- rnorm(n, mean = 100, sd = 15)  # Specific conductance
wtemp <- rnorm(n, mean = 20, sd = 2)  # Water temperature
Nup <- 2 + 0.5 * nitrogen + 0.3 * NEPt + 0.2 * flow + 0.1 * SPC - 0.2 * wtemp + rnorm(n)  # Nutrient uptake

# Combine into a data frame
data <- data.frame(sampletype, nitrogen, NEPt_1, NEPt, flow, SPC, wtemp, Nup)

# Fit the SEM model
fit <- sem(model, data = data)

# Summary of the SEM model
summary(fit, fit.measures = TRUE, standardized = TRUE)

# Visualize the SEM model
semPaths(
  fit,
  whatLabels = "std",      # Show standardized coefficients
  layout = "tree",         # Tree layout for clarity
  edge.label.cex = 1,      # Text size for edge labels
  sizeMan = 8,             # Size of manifest variables
  sizeLat = 10,            # Size of latent variables
  residuals = TRUE,        # Show residuals
  intercepts = FALSE,      # Hide intercepts for cleaner visualization
  title = TRUE             # Show title
)

lavaanPlot::lavaanPlot(model = fit, coefs = TRUE, stand = TRUE, covs = TRUE)

###########
##########

# Load required libraries
library(piecewiseSEM)

# Simulate data for pSEM
set.seed(123)

# Number of observations
n <- 200

# Simulated variables
nitrogen <- rnorm(n, mean = 5, sd = 1)  # Nitrogen concentration
flow <- rnorm(n, mean = 10, sd = 2)     # Flow
SPC <- rnorm(n, mean = 100, sd = 15)    # Specific conductance
wtemp <- rnorm(n, mean = 20, sd = 2)    # Water temperature
NEPt_1 <- rnorm(n, mean = 3, sd = 1)    # Lagged NEP

# Current NEP depends on lagged NEP, nitrogen, flow, and temperature
NEPt <- 0.5 * NEPt_1 + 0.3 * nitrogen + 0.2 * flow - 0.1 * wtemp + rnorm(n, sd = 0.5)

# Nutrient uptake (Nup) depends on nitrogen, NEPt, flow, SPC, and wtemp
Nup <- 2 + 0.4 * nitrogen + 0.5 * NEPt + 0.2 * flow + 0.1 * SPC - 0.2 * wtemp + rnorm(n)

# Combine into a data frame
data <- data.frame(nitrogen, flow, SPC, wtemp, NEPt_1, NEPt, Nup)

# Define piecewise SEM model
mod_list <- list(
  lm(NEPt ~ NEPt_1 + nitrogen + flow + wtemp, data = data),
  lm(Nup ~ nitrogen + NEPt + flow + SPC + wtemp, data = data)
)

# Evaluate the SEM
sem_fit <- psem(lm(NEPt ~ NEPt_1 + nitrogen + flow + wtemp, data = data),
                lm(Nup ~ nitrogen + NEPt + flow + wtemp, data = data))

# Summary of the SEM
summary(sem_fit)

# Plot the SEM
plot(sem_fit)

