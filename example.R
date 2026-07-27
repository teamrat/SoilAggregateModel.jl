# Load data in R:
# Clear environment
rm(list = ls(all = FALSE))
library(tidyverse)

pools    <- read.csv("output/timeseries_pools.csv")
fluxes   <- read.csv("output/timeseries_fluxes.csv")
microbes <- read.csv("output/timeseries_microbes.csv")
spatial  <- read.csv("output/spatial_profile_day365.csv")

# Plot carbon pools over time:
pools_long <- pools %>%
  pivot_longer(-time_days, names_to = "pool", values_to = "carbon")

ggplot(pools_long, aes(x = time_days, y = carbon, color = pool)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "Time (days)", y = "Carbon (μg-C, log scale)",
       title = "Carbon pool dynamics over 1 year")

# Plot CO2 flux:
ggplot(fluxes, aes(x = time_days, y = CO2_flux)) +
  geom_line() +
  labs(x = "Time (days)", y = "CO2 flux (μg-C/day)",
       title = "Respiration rate over 1 year")

# Plot final spatial profile:
ggplot(spatial, aes(x = radius_mm, y = B_bacteria)) +
  geom_line() +
  labs(x = "Radius (mm)", y = "Bacterial biomass (μg-C/mm³)",
       title = "Bacterial distribution at day 365")
