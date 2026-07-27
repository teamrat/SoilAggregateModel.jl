"""
Example: 1-Year Aggregate Simulation with CSV Export

This script demonstrates how to:
1. Run a 1-year simulation with daily outputs
2. Extract time-series data for all carbon pools using package functions
3. Calculate CO2 flux using package functions
4. Export to CSV files for analysis in R or Python

Output files:
- timeseries_pools.csv: All carbon pools over time
- timeseries_fluxes.csv: CO2 flux
- spatial_profile_day365.csv: Final spatial distribution
"""

## Setup
using Pkg
Pkg.activate(@__DIR__() * "/..")

using SoilAggregateModel
using DataFrames
using CSV

println("="^70)
println("Example: 1-Year Aggregate Simulation")
println("="^70)
println()

## Configure simulation parameters
bio = BiologicalProperties()
soil = SoilProperties()

# Environmental conditions (constant for simplicity)
T_func = t -> 293.15  # 20°C in K
ψ_func = t -> -33.0   # Field capacity in kPa
O2_func = t -> 0.27   # 21% O₂ in μg/mm³ (NOT molar fraction)

# Time span and output frequency
t_start = 0.0
t_end = 365.0  # 1 year
output_interval = 1.0  # Daily outputs

# Generate output times (daily)
output_times = collect(t_start:output_interval:t_end)

println("Simulation configuration:")
println("  Duration: $(t_end) days")
println("  Output frequency: every $(output_interval) days")
println("  Number of outputs: $(length(output_times))")
println("  Grid points: 50 (default)")
println()

## Run simulation
println("Running simulation...")
println("(This may take ~1 minute)")
println()

result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (t_start, t_end);
                       output_times=output_times,
                       n_grid=50,
                       dt_max=0.1)

println("✓ Simulation complete!")
println("  Total steps: $(result.diagnostics["n_steps"])")
println("  Final time: $(result.diagnostics["final_time"]) days")
println()

## Extract time-series data using package functions

println("Extracting time-series data using package functions...")

# Use integrated_pools function from postprocessing
pools = integrated_pools(result)

# Use co2_flux function from postprocessing
flux = co2_flux(result)

println("✓ Time-series extracted")
println()

## Create DataFrames and export to CSV

# Create output directory if it doesn't exist
output_dir = joinpath(@__DIR__(), "..", "output")
if !isdir(output_dir)
    mkdir(output_dir)
    println("Created output directory: $output_dir")
end

println("Exporting to CSV...")
println()

# 1. Carbon pools time-series
df_pools = DataFrame(
    time_days = pools.t,
    C_dissolved = pools.C_total,
    B_bacteria = pools.B_total,
    F_insulated = pools.F_i_total,
    F_noninsulated = pools.F_n_total,
    F_mobile = pools.F_m_total,
    F_total = pools.F_i_total .+ pools.F_n_total .+ pools.F_m_total,
    E_EPS = pools.E_total,
    M_MAOC = pools.M_total,
    POM = pools.P,
    CO2_cumulative = pools.CO2,
    Total_carbon = pools.C_total .+ pools.B_total .+
                   pools.F_i_total .+ pools.F_n_total .+ pools.F_m_total .+
                   pools.E_total .+ pools.M_total .+ pools.P .+ pools.CO2
)

csv_pools = joinpath(output_dir, "timeseries_pools.csv")
CSV.write(csv_pools, df_pools)
println("✓ Saved: $(csv_pools)")

# 2. CO2 flux time-series
df_fluxes = DataFrame(
    time_days = pools.t,
    CO2_flux_ugC_per_day = flux
)

csv_fluxes = joinpath(output_dir, "timeseries_fluxes.csv")
CSV.write(csv_fluxes, df_fluxes)
println("✓ Saved: $(csv_fluxes)")

# 3. Spatial profile at final time (day 365)
final_output = result.outputs[end]
final_state = final_output.state

df_spatial = DataFrame(
    radius_mm = result.grid.r_grid,
    C_dissolved = final_state.C,
    B_bacteria = final_state.B,
    F_insulated = final_state.F_i,
    F_noninsulated = final_state.F_n,
    F_mobile = final_state.F_m,
    E_EPS = final_state.E,
    M_MAOC = final_state.M,
    O_oxygen = final_state.O
)

csv_spatial = joinpath(output_dir, "spatial_profile_day365.csv")
CSV.write(csv_spatial, df_spatial)
println("✓ Saved: $(csv_spatial)")
println()

## Print summary statistics
println("="^70)
println("Summary Statistics")
println("="^70)
println()

println("Initial state (day 0):")
println("  Total C: $(round(df_pools.Total_carbon[1], digits=2)) μg-C")
println("  POM: $(round(pools.P[1], digits=2)) μg-C")
println("  Dissolved C: $(round(pools.C_total[1], digits=2)) μg-C")
println("  Bacteria: $(round(pools.B_total[1], digits=4)) μg-C")
println("  Fungi (total): $(round(pools.F_i_total[1] + pools.F_n_total[1] + pools.F_m_total[1], digits=4)) μg-C")
println()

println("Final state (day 365):")
println("  Total C: $(round(df_pools.Total_carbon[end], digits=2)) μg-C")
println("  POM remaining: $(round(pools.P[end], digits=2)) μg-C ($(round(100*pools.P[end]/pools.P[1], digits=1))%)")
println("  Dissolved C: $(round(pools.C_total[end], digits=2)) μg-C")
println("  Bacteria: $(round(pools.B_total[end], digits=2)) μg-C")
println("  Fungi (total): $(round(pools.F_i_total[end] + pools.F_n_total[end] + pools.F_m_total[end], digits=2)) μg-C")
println("    - F_insulated: $(round(pools.F_i_total[end], digits=2)) μg-C")
println("    - F_noninsulated: $(round(pools.F_n_total[end], digits=2)) μg-C")
println("    - F_mobile: $(round(pools.F_m_total[end], digits=2)) μg-C")
println("  EPS: $(round(pools.E_total[end], digits=2)) μg-C")
println("  MAOC: $(round(pools.M_total[end], digits=2)) μg-C")
println("  Cumulative CO2: $(round(pools.CO2[end], digits=2)) μg-C")
println()

println("CO2 flux (final day):")
println("  $(round(flux[end], digits=3)) μg-C/day")
println()

println("Peak values:")
F_total = pools.F_i_total .+ pools.F_n_total .+ pools.F_m_total
println("  Max bacteria: $(round(maximum(pools.B_total), digits=2)) μg-C (day $(pools.t[argmax(pools.B_total)]))")
println("  Max fungi: $(round(maximum(F_total), digits=2)) μg-C (day $(pools.t[argmax(F_total)]))")
println("  Max CO2 flux: $(round(maximum(flux), digits=2)) μg-C/day (day $(pools.t[argmax(flux)]))")
println()

println("Conservation:")
balance = carbon_balance_table(result)
println("  Final mass balance error: $(balance.relative_error[end]) (relative)")
println()

println("="^70)
println("Next steps for analysis in R:")
println("="^70)
println("""
# Load data in R:
library(tidyverse)

pools <- read_csv("output/timeseries_pools.csv")
fluxes <- read_csv("output/timeseries_fluxes.csv")
spatial <- read_csv("output/spatial_profile_day365.csv")

# Plot carbon pools over time:
pools_long <- pools %>%
  pivot_longer(-time_days, names_to = "pool", values_to = "carbon")

ggplot(pools_long, aes(x = time_days, y = carbon, color = pool)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "Time (days)", y = "Carbon (μg-C, log scale)",
       title = "Carbon pool dynamics over 1 year")

# Plot CO2 flux:
ggplot(fluxes, aes(x = time_days, y = CO2_flux_ugC_per_day)) +
  geom_line() +
  labs(x = "Time (days)", y = "CO2 flux (μg-C/day)",
       title = "Respiration rate over 1 year")

# Plot final spatial profile:
ggplot(spatial, aes(x = radius_mm, y = B_bacteria)) +
  geom_line() +
  labs(x = "Radius (mm)", y = "Bacterial biomass (μg-C/mm³)",
       title = "Bacterial distribution at day 365")
""")

println("="^70)
println("✓ Script complete!")
println("="^70)
