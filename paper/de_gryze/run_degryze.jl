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
include(joinpath(@__DIR__, "degryze_soils.jl"))

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
f_POM = diff(F_POM)
# RENORMALISE. The Normal is truncated at the 0.5 and 2.0 mm bin edges, so the
# raw bin fractions sum to 0.99889, not 1. That 0.11 % shortfall propagated
# straight through pom_population and was the entire residual in
# "Total POM input 4425.1 vs expected 4430" (4430 x 0.99889 = 4425.1).
# The amendment identity in src/physics/tessellation.jl is exact only for a
# normalised distribution.
f_POM ./= sum(f_POM)

# Particle diameter for aggregation stability [µm]
particle_D = 30.0  # µm

## ============================================================
## Soil and initial conditions  (must precede the tessellation)
## ============================================================
# Which of the five soils this run represents. Soil 3 (loam) is the mid-texture
# case; the paper's own conclusion is that soil 5 crosses a ~15 % clay threshold
# and behaves differently (p. 242).
SOIL_ID = 3

# Soil and initial conditions come from degryze_soils.jl, which overrides ONLY
# the fields De Gryze measures (ρ_b, f_clay_silt, θ_s; SOC, T_0, ψ_0) and leaves
# every other package default alone. Overrides below are run choices, not
# measurements — see docs/REFERENCE.md §5a for their standing.
soil = degryze_soil(SOIL_ID;
                    k_L     = 1000,
                    D_B_rel = 0.00001)   # §5a Group B: cited value is 0.01 (Wu 2006)

ic = degryze_ic(SOIL_ID, soil; s_M = 0.6)


## ============================================================
## Environmental Conditions (constant)
## ============================================================
# T and ψ are NOT free choices here: T is the incubation temperature and ψ is
# derived from 60 % WFPS for this soil's porosity (degryze_soils.jl).

T_const    = ic.T_0        # 298.15 K = 25 °C, De Gryze 2006 p.237
ψ_const    = ic.ψ_0        # from 60 % WFPS via the model's retention curve
O2_frac    = DEGRYZE_INCUBATION.O2_frac
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
## Domain tessellation and population — measured inputs only
## ============================================================
# De Gryze 2006 p.236-237: 1.5 g wheat stems per 150 g soil (= 10 g/kg),
# stems 44.3 % C  ->  4.43 g-C/kg-soil.  All densities µg/mm³ (= kg/m³).
# The geometry itself lives in src/physics/tessellation.jl.
ρ_POM   = 200.0      # POM carbon density [µg-C/mm³]
I_input = 4.43e-3    # amendment rate [µg-C/µg-soil]
ρ_bulk  = soil.ρ_b   # Table 1 bulk density for this soil [µg/mm³]
soil_V  = 100.0^3    # reference soil volume [mm³] = 1 litre

# Domain geometry, overlap correction and particle counts — src/physics/tessellation.jl
tess = domain_tessellation(ρ_POM=ρ_POM, I_input=I_input, ρ_b=ρ_bulk)
pop  = pom_population(diam_all, f_POM, tess; ρ_POM=ρ_POM, soil_volume_mm3=soil_V)

C_input_vol   = I_input * ρ_bulk
φ_POM         = tess.φ_POM
f_pack        = tess.f_pack
f_domain      = tess.f_domain
ω             = tess.ω
domain_factor = tess.f_domain

N_POM            = pop.N_POM
P_0_per_particle = pop.P_0_per_particle
total_POM_C      = pop.total_POM_C

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
# Stability threshold actually in force this run. G_c = τ_w·d_32/κ_b sets the
# binding concentration an aggregate must reach, so print it: a run whose
# threshold is not the one intended is otherwise indistinguishable from a
# correct one.
println("Aggregate stability:")
println("  d_32 (Sauter mean): $(round(soil.d_32*1000, digits=3)) µm")
println("  κ_b: $(soil.κ_b) Pa·mm/(µg/mm³),  w_E: $(soil.w_E)")
let L_s = 13.0, f_s = 34.0/60.0, μ = 1.002e-3, ν_w = 1.004
    v_s = π * f_s * L_s
    δ_s = sqrt(2.0 * ν_w / (2π * f_s))
    τ_w = sqrt(2.0) * μ * (v_s * 1e-3) / (δ_s * 1e-3)
    println("  τ_w: $(round(τ_w, sigdigits=5)) Pa,  δ_s: $(round(δ_s, digits=4)) mm")
    println("  G_c = τ_w·d_32/κ_b = $(round(τ_w * soil.d_32 / soil.κ_b, sigdigits=7)) µg/mm³")
    println("       (pre-2026-07-28 value for this soil: 0.0194084)")
end
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
    # Sensitivity setting, not a calibrated value: probing whether the CO2
    # overshoot and the too-fast MWD rise share a cause in POM supply rate.
    # optimized_params_soil3.txt was fitted against the broken POM
    # normalization and must be re-fitted, not reused.
    R_P_max = 2.5,
    Y_B_max = 0.4,
    B_S = 0.05,

    C_B = 5.0e-5, 

    #DEATH RATE
    μ_B = 0.0036,
    μ_F = 0.02
    
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

# Mineral mass distributed across the sieve classes. Required whenever
# aggregates are reported as a fraction of the whole sample: the model predicts
# aggregates, not the soil's own particle-size distribution.
#
# <53 µm  = silt + clay, measured (f_clay_silt).
# The sand is split across 53–250 and 250–2000 µm by assumption A1 of
# degryze_EJSS_2006_spec.md — equal mass per log interval between the bounding
# sieves, which contains no reference to any measured MWD:
#     f(53–250) = ln(250/53)/ln(2000/53) = 0.4272
# >2000 µm = 0: sand is ≤2000 µm by definition and the soil was crushed through
# 250 µm, so no primary particle can occupy that class.
# Mineral distribution across the sieve classes, and eq. (1) nominals —
# both defined in degryze_soils.jl. Assumptions A1 / A1b / A1c apply; see
# degryze_EJSS_2006_spec.md §0a before reporting anything derived from them.
# Assumption A1c (well-mixed matrix, absorption in proportion to abundance) —
# see degryze_EJSS_2006_spec.md. The `shell` column printed below is the check on
# it: where the shell is thinner than the coarse sieve class, the split between
# the three finer classes is indicative only. `pct_gt2000um` is unaffected,
# because its mineral fraction is zero.
mineral_classes = degryze_mineral_classes(SOIL_ID)
class_nominals  = DEGRYZE_CLASS_NOMINALS

# No ω here: the sums inside are built from physical particle counts and
# physical diameters, so they are already per-sample totals (REFERENCE.md §5b).
df_pop = population_outputs(df_summary, N_POM;
                            sieve_sizes=DEGRYZE_SIEVES,
                            mineral_class_fractions=mineral_classes,
                            class_nominal_mm=class_nominals,
                            class_labels=DEGRYZE_CLASS_LABELS,
                            cell_volume_mm3=pop.V_pack,
                            soil_volume_mm3=soil_V,
                            ρ_b=ρ_bulk,
                            f_C_POM=0.443)

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
    println("  Day $(Int(t_check)): D_agg(mean) = $(round(row.MWD_agg_only, digits=3)) mm, " *
            "CO₂ = $(round(co2_per_g, digits=1)) µg-C/g-soil, " *
            ">2mm = $(round(row.pct_gt2000um, digits=1))%, " *
            "0.25-2mm = $(round(row.pct_um250_2000, digits=1))%, " *
            "f_agg = $(round(row.f_agg, digits=3)), " *
            "shell = $(round(row.shell_mm*1000, digits=0)) µm, " *
            "cell_occ(max) = $(round(row.max_cell_occ, digits=3))")
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
p1 = plot(df_pop.time_days, df_pop.pct_gt2000um, lw=2, label="Model",
          xlabel="Time (days)", ylabel="% of sample mass > 2 mm",
          title="Large macroaggregates", legend=:bottomright)
# De Gryze Table 3: measured formation rates, % large macroaggregates per day.
for (lab, m) in (("Soil 1", 0.57), ("Soil 2", 0.59), ("Soil 3", 0.83),
                 ("Soil 4", 0.91), ("Soil 5", 2.02))
    plot!(p1, 0:21, (0:21) .* m, ls=:dot, lw=1.5, label=lab)
end

# Model CO₂ per gram of soil. Both sides of this comparison are TOTAL soil
# respiration: the model integrates respiration from all carbon it carries
# (residue plus the initial microbial, EPS and MAOC pools), and Figure 3 is bulk
# respiration with no unamended control. No partition factor is applied — see
# degryze_EJSS_2006_spec.md §0a A3.
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

# --- Panel 4: aggregate size classes, as NON-OVERLAPPING ranges ---
# These are De Gryze eq. (1)'s own four fractions [A%] [B%] [C%] [D%], so they
# sum to 100 and are directly comparable to the paper. The cumulative
# Cumulative "fraction above sieve X" is deliberately not reported: nested
# curves cannot show mass MOVING between classes, which is the whole signal.
# It is the running sum of these from the top if ever needed.
p4 = plot(xlabel="Time (days)", ylabel="% of sample mass",
          title="Aggregate size distribution", legend=:right)
plot!(p4, df_pop.time_days, df_pop.pct_gt2000um,   lw=2, label="> 2.00 mm")
plot!(p4, df_pop.time_days, df_pop.pct_um250_2000, lw=2, label="0.25 - 2.00 mm")
plot!(p4, df_pop.time_days, df_pop.pct_um53_250,   lw=2, label="0.05 - 0.25 mm")
plot!(p4, df_pop.time_days, df_pop.pct_lt53um,     lw=2, label="< 0.05 mm", ls=:dash)

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
println("÷ soil_g  = ", round(raw_sum / soil_mass_per_L, digits=1), " µg-C/g-soil")
println("\ndf_pop CO2 at day 21: ", round(df_pop[argmin(abs.(df_pop.time_days .- 21.0)), :CO2_total] / soil_mass_per_L, digits=1))


# POM is NOT ω-diluted: it is a lumped scalar at each domain centre, not
# background soil carbon spread over overlapping domains. Dividing by ω here
# understated the residue input by 28.8x and made CO2 look like 25x overshoot.
total_POM_per_g = sum(N_POM[i] * P_0_per_particle[i] for i in 1:5) / soil_mass_per_L
println("Total POM input: ", round(total_POM_per_g, digits=1), " µg-C/g-soil")
println("Expected: 4430 µg-C/g-soil")

println("Σ(N_i × P_0_i) = ", round(sum(N_POM[i] * P_0_per_particle[i] for i in 1:5), digits=1))
println("÷ soil_mass_per_L = ", round(sum(N_POM[i] * P_0_per_particle[i] for i in 1:5) / soil_mass_per_L, digits=1))
println("÷ ω ÷ soil_mass_per_L = ", round(sum(N_POM[i] * P_0_per_particle[i] for i in 1:5) / ω / soil_mass_per_L, digits=1))