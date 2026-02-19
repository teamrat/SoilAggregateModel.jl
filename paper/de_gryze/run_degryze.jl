"""
De Gryze et al. (2006) Forward Simulation
==========================================
European Journal of Soil Science, April 2006, 57, 235–246

Replicates the MATLAB `de_gryze_test.m` simulation:
  - Sweep over 5 POM diameter classes (0.5 → 2.0 mm, 0.3 mm bins)
  - Constant environment: ψ = -29 kPa (field capacity), T = 20°C, 21% O₂
  - 60-day simulation (data extends to day 21; extra time for post-peak dynamics)
  - Single particle diameter d_p = 30 µm

Each POM diameter produces one aggregate lifecycle. Together they represent a
population whose mean weight diameter (MWD) and cumulative CO₂ can be compared
with incubation data from 5 soils.

Output structure:
  output/
    summary.csv           — combined table: diameter × time → aggregate_D, cum_CO2, pools
    population.csv        — population-weighted MWD, CO₂, WAS

MATLAB correspondence:
  de_gryze_test.m         → this script (driver)
  param.m                 → BiologicalProperties() + SoilProperties() defaults
  single_aggregate_beta.m → run_aggregate()

Domain tessellation: see domain_tessellation_supplemental.md
"""

## ============================================================
## Setup
## ============================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))  # activate SoilAggregateModel root

using Revise
using SoilAggregateModel
using DataFrames
using CSV
using Distributions

# Postprocessing utilities (kept outside the package to avoid DataFrames dependency)
include(joinpath(@__DIR__, "postprocess_dataframe.jl"))

println("="^72)
println("De Gryze et al. (2006) — Forward Simulation")
println("="^72)
println()

## ============================================================
## POM Size Distribution
## ============================================================
#
# De Gryze: wheat stems cut to 0.5–2.0 mm pieces.
# Assume Normal distribution: µ = 1.25 mm, σ = 0.23 mm
# 5 bins of 0.3 mm width spanning 0.5–2.0 mm.
# Bin edges and midpoints all in [mm].
#
POM_mean  = 1.25     # mm
POM_sigma = 0.23     # mm

bin_edges = collect(0.5:0.3:2.0)   # [0.5, 0.8, 1.1, 1.4, 1.7, 2.0] mm
diam_all  = [(bin_edges[i] + bin_edges[i+1]) / 2.0 for i in 1:length(bin_edges)-1]
# diam_all = [0.65, 0.95, 1.25, 1.55, 1.85] mm

F_POM = cdf.(Normal(POM_mean, POM_sigma), bin_edges)
f_POM = diff(F_POM)  # relative frequency per bin (sums to ~1)

# Particle diameter for aggregation stability [µm]
particle_D = 30.0  # µm

## ============================================================
## Environmental Conditions (constant)
## ============================================================
#   T = 20°C = 293.15 K
#   ψ = -29 kPa (field capacity)
#   O₂: 21% molar fraction → concentration via ideal gas law

T_const    = 293.15   # K
ψ_const    = -29.0    # kPa
O2_frac    = 0.21     # molar fraction
M_O2       = 0.032    # kg/mol
P_atm      = 101000.0 # Pa
R_gas      = 8.314    # J/mol/K
O2_const   = O2_frac * P_atm * M_O2 / (R_gas * T_const)  # µg/mm³

T_func  = t -> T_const
ψ_func  = t -> ψ_const
O2_func = t -> O2_const

## ============================================================
## Time and Grid
## ============================================================

t_max      = 45.0    # days
dt_output  = 0.125   # output every 3 hours
n_grid     = 200     # radial grid points per aggregate

## ============================================================
## Domain Tessellation
## ============================================================
#
# All densities in µg/mm³ (= kg/m³).
#
# De Gryze experiment: 1.5 g wheat stems per 150 g soil = 10 g/kg
# Wheat stems: 44.3% C → 4.43 g-C/kg-soil = 4.43e-3 µg-C/µg-soil
#
# f_domain = (ρ_POM / (I × ρ_B))^(1/3)
#
ρ_POM        = 200.0         # POM carbon density [µg-C/mm³]
I_input      = 4.43e-3       # amendment rate [µg-C/µg-soil] (dimensionless)
ρ_bulk       = 1300.0        # soil bulk density [µg/mm³]
f_domain_min = 10.0          # minimum domain factor for numerical resolution

C_input_vol = I_input * ρ_bulk            # [µg-C/mm³]
φ_POM       = C_input_vol / ρ_POM        # POM volume fraction [-]
f_pack      = (1.0 / φ_POM)^(1/3)       # packing factor
f_domain    = max(f_pack, f_domain_min)  # model domain factor
ω           = (f_domain / f_pack)^3      # overlap correction
domain_factor = f_domain

## ============================================================
## Population: number of POM particles per liter soil
## ============================================================

soil_V = 100.0^3   # mm³ (1 liter)

# Packing cell volume for each size class [mm³]
V_pack = [(4.0/3.0) * π * (d * f_pack / 2)^3 for d in diam_all]

# Number of particles per liter soil
N_POM = f_POM .* soil_V ./ V_pack

# POM carbon per particle [µg-C]
P_0_per_particle = [(4.0/3.0) * π * (d/2)^3 * ρ_POM for d in diam_all]
total_POM_C = sum(N_POM .* P_0_per_particle)  # µg-C / L-soil

## ============================================================
## Print configuration
## ============================================================

println("Configuration:")
println("  POM diameters: $(length(diam_all)) bins, $(diam_all[1])–$(diam_all[end]) mm")
println("  POM distribution: N($(POM_mean), $(POM_sigma)) mm")
println("  Bin frequencies: $(round.(f_POM, sigdigits=4))")
println("  Environment: T=$(T_const)K, ψ=$(ψ_const) kPa, O₂=$(round(O2_const, digits=4)) µg/mm³")
println("  Duration: $(t_max) days, output every $(dt_output) days")
println("  Grid: $(n_grid) nodes per aggregate")
println()
println("Domain tessellation:")
println("  I_input: $(I_input * 1000) g-C/kg-soil ($(round(C_input_vol, sigdigits=4)) µg-C/mm³)")
println("  ρ_POM: $(ρ_POM) µg-C/mm³, ρ_bulk: $(ρ_bulk) µg/mm³")
println("  φ_POM: $(round(φ_POM * 100, sigdigits=3))%")
println("  Packing factor: $(round(f_pack, sigdigits=4)). Selected domain factor $(f_domain)")
println("  Overlap ω = $(round(ω, sigdigits=4))")
println("  Total POM particles: $(round(sum(N_POM), sigdigits=4)) per liter soil")
println("  Total POM carbon: $(round(total_POM_C, sigdigits=4)) µg-C per liter soil")
println("  Total POM carbon: $(round(total_POM_C / (soil_V * ρ_bulk * 1e-6) , sigdigits=4)) µg-C per g-soil")
println()

## ============================================================
## Parameters
## ============================================================

bio  = BiologicalProperties(
    #MAOC
    κ_s_ref = 0.01,
    κ_d_ref = 0.001,

    #FUNGI MINIMUMS
    F_i_min = 1e-6,
    F_n_min = 2e-4,
    F_m_min = 1e-6,

    #TRANSPORT
    D_Fn0   = 0.00001,
    D_Fm0   = 0.001,

    #MAXIMUM UPTAKE RATE 
    r_B_max = 8.0,
    r_F_max = 0.2,
    R_P_max = 5,
    Y_B_max = 0.4,
    B_S = 0.05,

    C_B = 5.0e-5, 

    #DEATH RATE
    μ_B = 0.0036,
    μ_F = 0.02
    
)

# De Gryze silt loam: 18% clay, 56% silt, 26% sand
soil = SoilProperties(
    k_L = 1000,
    D_B_rel = 0.00001,
    ρ_b = ρ_bulk,                    # µg/mm³ (consistent with domain tessellation)
    f_clay_silt = 0.74               # 18% clay + 56% silt
)

# Initial conditions: average of 5 soils (SOC 1.78–3.10%)
# Native POM removed (H₂O₂ + combustion), but fine fraction retains
# its microbial community and MAOC.
ic = InitialConditions(
    SOC   = 0.0220,       # 2.20% — midrange of 5 soils
    s_M   = 0.6,          # 40% MAOC saturation
    f_bact  = 0.01,       # 1% of SOC in bacteria
    f_fungi = 0.01,       # 1% of SOC in fungi
    f_eps   = 0.005,      # 0.5% of SOC in EPS
    T_0   = T_const,
    ψ_0   = ψ_const,
    O2_gas = O2_const
)

SOC_vol = ic.SOC * ρ_bulk
println("Initial conditions (SOC-partitioned, ω-diluted):")
println("  SOC: $(ic.SOC * 100)% → $(round(SOC_vol, sigdigits=4)) µg-C/mm³ (physical)")
println("  After ω dilution: $(round(SOC_vol/ω, sigdigits=4)) µg-C/mm³ (model)")
println("  B_0 (physical): $(round(ic.f_bact * SOC_vol, sigdigits=4)) µg-C/mm³")
println("  B_0 (model):    $(round(ic.f_bact * SOC_vol / ω, sigdigits=4)) µg-C/mm³")
println()

## ============================================================
## Run simulations
## ============================================================

output_times = collect(0.0:dt_output:t_max)

output_dir = joinpath(@__DIR__, "output")
isdir(output_dir) || mkpath(output_dir)

# Spatial snapshot times for diagnostics
snap_times = [0.0, 1.0, 5.0, 6.0, 28.0, 29.0, 30.0]

df_summary, df_snaps = run_diameter_sweep(diam_all, bio, soil, T_func, ψ_func, O2_func;
                                 t_max=t_max, output_times=output_times,
                                 n_grid=n_grid, domain_factor=domain_factor,
                                 ρ_POM=ρ_POM, ic=ic, ω=ω,
                                 snap_times=snap_times)

## ============================================================
## Save combined summary
## ============================================================

CSV.write(joinpath(output_dir, "summary.csv"), df_summary)
println("✓ Summary: $(nrow(df_summary)) rows ($(length(diam_all)) diameters × $(length(output_times)) times)")

CSV.write(joinpath(output_dir, "spatial_profiles.csv"), df_snaps)
println("✓ Spatial profiles: $(nrow(df_snaps)) rows at t = $(snap_times) days")
println()

## ============================================================
## Population-level outputs
## ============================================================

df_pop = population_outputs(df_summary, N_POM;
                            ω=ω, sieve_sizes=[0.25, 0.5, 1.0, 2.0])

CSV.write(joinpath(output_dir, "population.csv"), df_pop)
println("✓ Population outputs saved")
println()

## ============================================================
## Load experimental data
## ============================================================

data_dir = @__DIR__

df_mwd_data = nothing
if isfile(joinpath(data_dir, "degryze2006.csv"))
    df_mwd_data = CSV.read(joinpath(data_dir, "degryze2006.csv"), DataFrame;
                           header=2, missingstring="")
    println("MWD data loaded: $(nrow(df_mwd_data)) time points, $(ncol(df_mwd_data)-1) soils")
end

df_co2_data = nothing
if isfile(joinpath(data_dir, "degryze_CO2_2006.csv"))
    df_co2_data = CSV.read(joinpath(data_dir, "degryze_CO2_2006.csv"), DataFrame)
    if eltype(df_co2_data[:, 1]) <: AbstractString
        df_co2_data[!, 1] = parse.(Float64, df_co2_data[:, 1])
    end
    println("CO₂ data loaded: $(nrow(df_co2_data)) time points, $(ncol(df_co2_data)-1) soils")
end
println()

## ============================================================
## Summary printout
## ============================================================

println("="^72)
println("Population summary at key time points:")
println("="^72)

# Convert CO₂ to µg-C/g-soil for display
soil_mass_per_L = soil_V * ρ_bulk * 1e-6  # grams of soil per liter

for t_check in [0.0, 7.0, 14.0, 21.0, 60.0]
    row = df_pop[argmin(abs.(df_pop.time_days .- t_check)), :]
    co2_per_g = row.CO2_total / soil_mass_per_L
    println("  Day $(Int(t_check)): MWD = $(round(row.MWD_mm, digits=3)) mm, " *
            "CO₂ = $(round(co2_per_g, digits=1)) µg-C/g-soil, " *
            "WAS>0.25mm = $(round(row.WAS_0_25mm * 100, digits=1))%")
end
println()

## ============================================================
## Diagnostic plots: Model vs Data
## ============================================================

using Plots
gr()

soil_colors = [:royalblue, :darkorange, :seagreen, :firebrick, :purple]
soil_names  = ["Soil 1", "Soil 2", "Soil 3", "Soil 4", "Soil 5"]

# --- Panel 1: MWD ---
p1 = plot(df_pop.time_days, df_pop.MWD_mm,
          lw=2, color=:black, label="Model",
          xlabel="Time (days)", ylabel="MWD (mm)",
          title="Mean Weight Diameter",
          xlim=(0, 45), ylim=(0,5), legend=:topleft)

if df_mwd_data !== nothing
    t_mwd = df_mwd_data[:, 1]
    for (i, col) in enumerate(names(df_mwd_data)[2:end])
        vals = df_mwd_data[:, col]
        mask = .!ismissing.(vals)
        if any(mask)
            scatter!(p1, t_mwd[mask], Float64.(vals[mask]),
                     color=soil_colors[i], label=soil_names[i],
                     ms=5, markerstrokewidth=0.5)
        end
    end
end

# --- Panel 2: Cumulative CO₂ ---
# CO₂_total is µg-C/L-soil (overlap-corrected). Convert to µg-C/g-soil.
co2_model_ugC_per_gsoil = df_pop.CO2_total ./ soil_mass_per_L

p2 = plot(df_pop.time_days, co2_model_ugC_per_gsoil,
          lw=2, color=:black, label="Model",
          xlabel="Time (days)", ylabel="Cum. CO₂ (µg-C/g-soil)",
          title="Cumulative Respiration",
          xlim=(0, 45), legend=:topleft)

if df_co2_data !== nothing
    t_co2 = df_co2_data[:, 1]
    for (i, col) in enumerate(names(df_co2_data)[2:end])
        vals = df_co2_data[:, col]
        mask = .!ismissing.(vals)
        if any(mask)
            scatter!(p2, t_co2[mask], Float64.(vals[mask]),
                     color=soil_colors[i], label=soil_names[i],
                     ms=5, markerstrokewidth=0.5)
        end
    end
end

# --- Panel 3: CO₂ flux ---
p3 = plot(df_pop.time_days, df_pop.CO2_flux_total ./ soil_mass_per_L,
          lw=2, color=:black, label="Model",
          xlabel="Time (days)", ylabel="CO₂ flux (µg-C/g-soil/day)",
          title="Respiration Rate",
          xlim=(0, 45), legend=:topright)

# --- Panel 4: WAS ---
p4 = plot(df_pop.time_days, df_pop.WAS_0_25mm .* 100, lw=2, label=">0.25 mm",
          xlabel="Time (days)", ylabel="Aggregate volume (%)",
          title="Wet Aggregate Stability",
          xlim=(0, 45), ylim=(0, 105), legend=:bottomright)
plot!(p4, df_pop.time_days, df_pop.WAS_0_5mm  .* 100, lw=2, label=">0.5 mm")
plot!(p4, df_pop.time_days, df_pop.WAS_1_0mm  .* 100, lw=2, label=">1.0 mm")
plot!(p4, df_pop.time_days, df_pop.WAS_2_0mm  .* 100, lw=2, label=">2.0 mm")

# --- Combine ---
p_all = plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 700),
             plot_title="De Gryze (2006) — Model vs Data")

savefig(p_all, joinpath(output_dir, "degryze_model_vs_data.png"))
println("✓ Plot saved: $(joinpath(output_dir, "degryze_model_vs_data.png"))")
display(p_all)

println()
println("="^72)
println("✓ De Gryze forward simulation complete")
println("="^72)

## At day 21, for each diameter class:
println("="^72)

diams = sort(unique(df_summary.diam_mm))
for (i, d) in enumerate(diams)
    row = df_summary[(df_summary.diam_mm .== d) .& (abs.(df_summary.time_days .- 21.0) .< 0.01), :]
    co2 = row.CO2_cumulative[1]
    println("  d=$(d) mm  N=$(round(N_POM[i], digits=1))  CO2=$(round(co2, digits=2))  N×CO2=$(round(N_POM[i]*co2, digits=1))")
end

raw_sum = sum(N_POM[i] * df_summary[(df_summary.diam_mm .== diams[i]) .& (abs.(df_summary.time_days .- 21.0) .< 0.01), :CO2_cumulative][1] for i in 1:5)
println("\nΣ(N×CO2) = ", round(raw_sum, digits=1))
println("÷ ω      = ", round(raw_sum / ω, digits=1))
println("÷ soil_g  = ", round(raw_sum / ω / soil_mass_per_L, digits=1), " µg-C/g-soil")
println("\ndf_pop CO2 at day 21: ", round(df_pop[argmin(abs.(df_pop.time_days .- 21.0)), :CO2_total] / soil_mass_per_L, digits=1))


total_POM_per_g = sum(N_POM[i] * P_0_per_particle[i] for i in 1:5) / ω / soil_mass_per_L
println("Total POM input: ", round(total_POM_per_g, digits=1), " µg-C/g-soil")
println("Expected: 4430 µg-C/g-soil")

println("Σ(N_i × P_0_i) = ", round(sum(N_POM[i] * P_0_per_particle[i] for i in 1:5), digits=1))
println("÷ soil_mass_per_L = ", round(sum(N_POM[i] * P_0_per_particle[i] for i in 1:5) / soil_mass_per_L, digits=1))
println("÷ ω ÷ soil_mass_per_L = ", round(sum(N_POM[i] * P_0_per_particle[i] for i in 1:5) / ω / soil_mass_per_L, digits=1))