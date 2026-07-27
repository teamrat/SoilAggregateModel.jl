# diagnostics_plot.R
# Diagnostic plots for SoilAggregateModel output
# ================================================
# Assumes summary_data is already loaded as a data.frame
# with columns from result_to_dataframe() + diameter sweep metadata
# Load data in R:
# Clear environment
rm(list = ls(all = FALSE))
library(tidyverse)


library(ggplot2)
library(dplyr)
library(tidyr)

# ── Color palette ──
pool_colors <- c(
  "POM"  = "#8B4513",   # brown
  "DOC"  = "#1E90FF",   # blue
  "Bact" = "#FF4500",   # red-orange
  "Fungi"= "#228B22",   # green
  "EPS"  = "#DAA520",   # goldenrod
  "MAOC" = "#708090",   # slate gray
  "CO2"  = "#2F2F2F"    # near-black
)


summary_data = read.csv(file="output/summary.csv")

# ============================================================
# 1. Mass conservation check per diameter class
# ============================================================
# C_system_total should be constant for each diameter.
# Plot relative error over time.

p_conservation <- summary_data %>%
  ggplot(aes(x = time_days, y = C_balance_error, color = factor(diam_mm))) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(labels = scales::scientific) +
  labs(
    title = "Carbon Conservation Error",
    subtitle = "Relative error: (C_total(t) - C_total(0)) / C_total(0)",
    x = "Time (days)", y = "Relative error",
    color = "POM Ø (mm)"
  ) +
  theme_minimal(base_size = 12)

print(p_conservation)

# ============================================================
# 2. Carbon pool fate: per-mille stacked area (population-weighted)
# ============================================================
# Average across diameter classes (mass-weighted by P_0)
# or just show the middle diameter (1.25 mm) as representative

df_mid <- summary_data %>%
  filter(abs(diam_mm - 1.25) < 0.01)

df_fate <- df_mid %>%
  select(time_days, POM_permille, C_permille, B_permille,
         F_permille, E_permille, M_permille, CO2_permille) %>%
  pivot_longer(-time_days, names_to = "pool", values_to = "permille") %>%
  mutate(pool = recode(pool,
                       "POM_permille" = "POM",
                       "C_permille"   = "DOC",
                       "B_permille"   = "Bact",
                       "F_permille"   = "Fungi",
                       "E_permille"   = "EPS",
                       "M_permille"   = "MAOC",
                       "CO2_permille" = "CO2"
  )) %>%
  mutate(pool = factor(pool, levels = c("CO2", "MAOC", "EPS", "Fungi", "Bact", "DOC", "POM")))

p_fate_stacked <- ggplot(df_fate, aes(x = time_days, y = permille, fill = pool)) +
  geom_area(alpha = 0.85) +
  scale_fill_manual(values = pool_colors) +
  labs(
    title = "Carbon Fate — Stacked (D = 1.25 mm)",
    subtitle = "Per mille of initial total carbon",
    x = "Time (days)", y = "‰ of initial C",
    fill = "Pool"
  ) +
  theme_minimal(base_size = 12)

print(p_fate_stacked)

# ============================================================
# 3. Carbon pools as individual lines (log scale optional)
# ============================================================

p_fate_lines <- df_fate %>%
  ggplot(aes(x = time_days, y = permille, color = pool)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pool_colors) +
  labs(
    title = "Carbon Pool Dynamics — Lines (D = 1.25 mm)",
    x = "Time (days)", y = "‰ of initial C",
    color = "Pool"
  ) +
  theme_minimal(base_size = 12)

print(p_fate_lines)

# ============================================================
# 4. POM depletion across diameter classes
# ============================================================

p_pom <- summary_data %>%
  ggplot(aes(x = time_days, y = POM_permille, color = factor(diam_mm))) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "POM Depletion by Diameter Class",
    x = "Time (days)", y = "POM (‰ of initial C)",
    color = "POM Ø (mm)"
  ) +
  theme_minimal(base_size = 12)

print(p_pom)

p_agg <- summary_data %>%
  ggplot(aes(x = time_days, y = D_agg_mm, color = factor(diam_mm))) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Aggregation by Diameter Class",
    x = "Time (days)", y = "D (mm)",
    color = "POM Ø (mm)"
  ) +
  theme_minimal(base_size = 12)

print(p_agg)

# ============================================================
# 5. CO₂ flux (instantaneous rate) across diameters
# ============================================================

p_flux <- summary_data %>%
  ggplot(aes(x = time_days, y = CO2_flux, color = factor(diam_mm))) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "CO₂ Flux by Diameter Class",
    subtitle = "Instantaneous respiration rate per aggregate",
    x = "Time (days)", y = "CO₂ flux (µg-C/day)",
    color = "POM Ø (mm)"
  ) +
  theme_minimal(base_size = 12)

print(p_flux)

# ============================================================
# 6. Cumulative CO₂ across diameters
# ============================================================

p_co2_cum <- summary_data %>%
  ggplot(aes(x = time_days, y = CO2_cumulative, color = factor(diam_mm))) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Cumulative CO₂ by Diameter Class",
    x = "Time (days)", y = "Cumulative CO₂ (µg-C)",
    color = "POM Ø (mm)"
  ) +
  theme_minimal(base_size = 12)

print(p_co2_cum)

# ============================================================
# 7. Biomass dynamics (B + F) across diameters
# ============================================================

p_biomass <- summary_data %>%
  ggplot(aes(x = time_days, color = factor(diam_mm))) +
  geom_line(aes(y = B_permille, linetype = "Bacteria"), linewidth = 0.7) +
  geom_line(aes(y = F_permille, linetype = "Fungi"), linewidth = 0.7) +
  labs(
    title = "Microbial Biomass by Diameter Class",
    x = "Time (days)", y = "Biomass (‰ of initial C)",
    color = "POM Ø (mm)", linetype = "Group"
  ) +
  theme_minimal(base_size = 12)

print(p_biomass)

# ============================================================
# 8. Mass balance table at key times
# ============================================================

cat("\n══════════════════════════════════════════════════════════\n")
cat("Mass balance at key time points (D = 1.25 mm, ‰ of initial C)\n")
cat("══════════════════════════════════════════════════════════\n\n")

df_table <- df_mid %>%
  filter(time_days %in% c(0, 1, 5, 10, 15, 21, 30, 45)) %>%
  mutate(
    Total = POM_permille + C_permille + B_permille + F_permille +
      E_permille + M_permille + CO2_permille
  ) %>%
  select(time_days, POM_permille, C_permille, B_permille, F_permille,
         E_permille, M_permille, CO2_permille, Total)

print(df_table, digits = 3, row.names = FALSE)

# ============================================================
# 9. Save all plots to PDF
# ============================================================

pdf("diagnostics_plots.pdf", width = 10, height = 7)
print(p_conservation)
print(p_fate_stacked)
print(p_fate_lines)
print(p_pom)
print(p_flux)
print(p_co2_cum)
print(p_biomass)
dev.off()

cat("\n✓ Plots saved to diagnostics_plots.pdf\n")

summary_data %>% 
  filter(abs(diam_mm - 1.25) < 0.01, time_days == 0) %>% 
  select(POM, C_total, B_total, F_total, E_total, M_total, CO2_cumulative)



# NEW========================

# ============================================================
# Load population data and plot in µg-C/g-soil
# ============================================================

pop_data <- read.csv("output/population.csv")

# Constants (must match Julia script)
rho_bulk       <- 1300.0
soil_V         <- 100^3       # mm³ per liter
soil_mass_per_L <- soil_V * rho_bulk * 1e-6  # grams soil per liter

# Convert CO2_total (µg-C/L-soil) to µg-C/g-soil
pop_data$CO2_per_g <- pop_data$CO2_total / soil_mass_per_L
pop_data$CO2_flux_per_g <- pop_data$CO2_flux_total / soil_mass_per_L

# De Gryze Soil 3 data
data_times <- c(0, 5, 10, 15, 21)
data_CO2   <- c(0, 481.954, 1263.633, 1766.159, 2139.010)
df_data <- data.frame(time = data_times, CO2 = data_CO2)

# POM input for reference line
POM_input_per_g <- 4430  # µg-C/g-soil

p_co2 <- ggplot() +
  geom_line(data = pop_data, aes(x = time_days, y = CO2_per_g),
            linewidth = 0.8, color = "black") +
  geom_point(data = df_data, aes(x = time, y = CO2),
             color = "firebrick", size = 3) +
  geom_hline(yintercept = POM_input_per_g, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  annotate("text", x = 40, y = POM_input_per_g + 100,
           label = "Total POM input", hjust = 1, color = "gray50") +
  labs(
    title = "Cumulative CO₂ — Model vs Soil 3 Data",
    x = "Time (days)", y = "Cumulative CO₂ (µg-C/g-soil)"
  ) +
  theme_minimal(base_size = 12)

print(p_co2)

p_flux <- ggplot(pop_data, aes(x = time_days, y = CO2_flux_per_g)) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "CO₂ Flux",
    x = "Time (days)", y = "CO₂ flux (µg-C/g-soil/day)"
  ) +
  theme_minimal(base_size = 12)

print(p_flux)

# ============================================================
# Population-level carbon pools in µg-C/g-soil
# ============================================================

# Constants
rho_bulk        <- 1300.0
soil_V          <- 100^3
soil_mass_per_L <- soil_V * rho_bulk * 1e-6
rho_POM         <- 200.0
I_input         <- 4.43e-3
C_input_vol     <- I_input * rho_bulk
phi_POM         <- C_input_vol / rho_POM
f_pack          <- (1.0 / phi_POM)^(1/3)
f_domain_min    <- 10.0
f_domain        <- max(f_pack, f_domain_min)
omega           <- (f_domain / f_pack)^3

# POM size distribution
POM_mean  <- 1.25
POM_sigma <- 0.23
bin_edges <- seq(0.5, 2.0, by = 0.3)
diam_all  <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2
f_POM     <- diff(pnorm(bin_edges, POM_mean, POM_sigma))
V_pack    <- (4/3) * pi * (diam_all * f_pack / 2)^3
N_POM     <- f_POM * soil_V / V_pack

# Build lookup: diameter → N_POM
n_lookup <- data.frame(diam_mm = diam_all, N_POM = N_POM)

# Merge N_POM into summary_data
df <- merge(summary_data, n_lookup, by = "diam_mm")

# For each time, compute population-weighted pools:
#   pool_per_g = Σ_i (N_POM_i × pool_i) / ω / soil_mass_per_L

df_pop_pools <- df %>%
  group_by(time_days) %>%
  summarise(
    POM_per_g  = sum(N_POM * POM)             / omega / soil_mass_per_L,
    DOC_per_g  = sum(N_POM * C_total)         / omega / soil_mass_per_L,
    B_per_g    = sum(N_POM * B_total)         / omega / soil_mass_per_L,
    F_per_g    = sum(N_POM * F_total)         / omega / soil_mass_per_L,
    E_per_g    = sum(N_POM * E_total)         / omega / soil_mass_per_L,
    M_per_g    = sum(N_POM * M_total)         / omega / soil_mass_per_L,
    CO2_per_g  = sum(N_POM * CO2_cumulative)  / omega / soil_mass_per_L,
    .groups = "drop"
  ) %>%
  mutate(
    Total_per_g = POM_per_g + DOC_per_g + B_per_g + F_per_g +
      E_per_g + M_per_g + CO2_per_g
  )

# Print budget at key times
cat("\n══════════════════════════════════════════════════════════\n")
cat("Population carbon budget (µg-C/g-soil)\n")
cat("══════════════════════════════════════════════════════════\n\n")
df_pop_pools %>%
  filter(time_days %in% c(0, 5, 10, 15, 21, 45)) %>%
  print(width = Inf)

# Reference values
POM_input_per_g <- I_input * 1e6  # 4430 µg-C/g-soil
SOC_per_g       <- 0.0221 * 1e6   # 22100 µg-C/g-soil

cat("\nReference:\n")
cat("  POM input:  ", POM_input_per_g, " µg-C/g-soil\n")
cat("  Initial SOC:", SOC_per_g, " µg-C/g-soil\n")
cat("  Total t=0:  ", round(df_pop_pools$Total_per_g[1], 1), " µg-C/g-soil\n")
cat("  Expected:   ", round(POM_input_per_g + SOC_per_g, 1), " µg-C/g-soil\n")

# ── Plot: all pools + data ──
df_long <- df_pop_pools %>%
  select(time_days, POM_per_g, DOC_per_g, B_per_g, F_per_g,
         E_per_g, M_per_g, CO2_per_g) %>%
  pivot_longer(-time_days, names_to = "pool", values_to = "ug_per_g") %>%
  mutate(pool = recode(pool,
                       "POM_per_g" = "POM",
                       "DOC_per_g" = "DOC",
                       "B_per_g"   = "Bact",
                       "F_per_g"   = "Fungi",
                       "E_per_g"   = "EPS",
                       "M_per_g"   = "MAOC",
                       "CO2_per_g" = "CO2"
  ))

pool_colors <- c(
  "POM"  = "#8B4513", "DOC"  = "#1E90FF", "Bact" = "#FF4500",
  "Fungi"= "#228B22", "EPS"  = "#DAA520", "MAOC" = "#708090",
  "CO2"  = "#2F2F2F"
)

# De Gryze Soil 3 data
df_data <- data.frame(time = c(0, 5, 10, 15, 21),
                      CO2 = c(0, 481.954, 1263.633, 1766.159, 2139.010))

p_budget <- ggplot(df_long, aes(x = time_days, y = ug_per_g, color = pool)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = df_data, aes(x = time, y = CO2),
             color = "firebrick", size = 3, inherit.aes = FALSE) +
  geom_line(data = df_pop_pools, aes(x = time_days, y = Total_per_g),
            linetype = "dashed", color = "black", linewidth = 0.5,
            inherit.aes = FALSE) +
  scale_color_manual(values = pool_colors) +
  labs(
    title = "Population Carbon Budget (µg-C/g-soil)",
    subtitle = "Dashed = total C (should be constant). Red points = Soil 3 data.",
    x = "Time (days)", y = "µg-C/g-soil",
    color = "Pool"
  ) +
  theme_minimal(base_size = 12)

print(p_budget)
