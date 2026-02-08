# run_simulations.jl — Theme 1: Single Aggregate Physics
#
# Simulates one representative aggregate (0.5 mm POM, 60 months, standard conditions)
# Exports CSVs for R plotting: timeseries, radial profiles, CO2 flux, aqueous concentrations
#
# Output CSVs:
# - data/timeseries.csv: integrated pools over time (both domains)
# - data/radial_profiles.csv: spatial profiles at snapshot times
# - data/co2_flux.csv: CO2 flux time series
# - data/aqueous_*mo.csv: aqueous concentrations at each snapshot month
#
# Manuscript figures: Fig 3 (carbon fate, MAOC decoupling), Fig 4 (stability lag),
# Fig 5 (radial profiles)

include("../common.jl")

println("="^70)
println("Theme 1: Single Aggregate Physics")
println("="^70)
println("POM: 0.5 mm diameter")
println("Conditions: T=20°C, ψ=-33 kPa, O₂=21%")
println("Duration: 60 months ($(60*30.44) days)")
println()

# --- Run simulation ---
println("Running simulation...")
result = run_aggregate(
    STANDARD_BIO, STANDARD_SOIL,
    t -> STANDARD_T, t -> STANDARD_ψ, t -> STANDARD_O2,
    (0.0, 60 * 30.44);  # 60 months in days
    n_grid = STANDARD_n_grid,
    r_0 = STANDARD_r_0,
    r_max = STANDARD_r_max,
    output_times = TIMES_DENSE
)

println("  Completed: $(result.diagnostics["n_steps"]) timesteps")

# Check conservation
bal = carbon_balance_table(result)
max_error = maximum(abs.(bal.relative_error))
println("  Carbon conservation: max relative error = $(max_error)")
if max_error > 1e-12
    @warn "Conservation error exceeds tolerance! Expected < 1e-12, got $(max_error)"
end
println()

# --- CSV 1: Integrated pools timeseries (both domains) ---
println("Exporting timeseries...")
pools = integrated_pools(result)

timeseries_data = (
    t_days = pools.t,
    t_months = pools.t ./ 30.44,
    P = pools.P,
    CO2 = pools.CO2,
    r_agg = pools.r_agg,
    # Total domain
    C_total = pools.C_total,
    B_total = pools.B_total,
    F_i_total = pools.F_i_total,
    F_n_total = pools.F_n_total,
    F_m_total = pools.F_m_total,
    E_total = pools.E_total,
    M_total = pools.M_total,
    # Aggregate domain
    C_agg = pools.C_agg,
    B_agg = pools.B_agg,
    F_i_agg = pools.F_i_agg,
    F_n_agg = pools.F_n_agg,
    F_m_agg = pools.F_m_agg,
    E_agg = pools.E_agg,
    M_agg = pools.M_agg
)
write_csv("data/timeseries.csv", timeseries_data)

# --- CSV 2: Radial profiles at snapshot times ---
println("Exporting radial profiles...")
profs = radial_profiles(result; times = SNAPSHOT_TIMES)

n_snap = length(profs)
n_r = result.grid.n
total_rows = n_snap * n_r

t_days_vec = Float64[]
t_months_vec = Float64[]
r_vec = Float64[]
C_vec = Float64[]
B_vec = Float64[]
F_i_vec = Float64[]
F_n_vec = Float64[]
F_m_vec = Float64[]
E_vec = Float64[]
M_vec = Float64[]
O_vec = Float64[]

for p in profs
    for i in 1:n_r
        push!(t_days_vec, p.t)
        push!(t_months_vec, p.t / 30.44)
        push!(r_vec, p.r[i])
        push!(C_vec, p.C[i])
        push!(B_vec, p.B[i])
        push!(F_i_vec, p.F_i[i])
        push!(F_n_vec, p.F_n[i])
        push!(F_m_vec, p.F_m[i])
        push!(E_vec, p.E[i])
        push!(M_vec, p.M[i])
        push!(O_vec, p.O[i])
    end
end

profile_data = (
    t_days = t_days_vec,
    t_months = t_months_vec,
    r = r_vec,
    C = C_vec,
    B = B_vec,
    F_i = F_i_vec,
    F_n = F_n_vec,
    F_m = F_m_vec,
    E = E_vec,
    M = M_vec,
    O = O_vec,
)
write_csv("data/radial_profiles.csv", profile_data)

# --- CSV 3: CO₂ flux ---
println("Exporting CO2 flux...")
flux = co2_flux(result)
flux_data = (
    t_days = [rec.t for rec in result.outputs],
    t_months = [rec.t / 30.44 for rec in result.outputs],
    CO2_flux = flux,
    CO2_cumulative = [rec.state.CO2_cumulative for rec in result.outputs]
)
write_csv("data/co2_flux.csv", flux_data)

# --- CSV 4: Aqueous concentrations at snapshots ---
println("Exporting aqueous concentrations...")
for (month, t_target) in zip(SNAPSHOT_MONTHS, SNAPSHOT_TIMES)
    # Find nearest output
    idx_nearest = argmin(abs.([rec.t for rec in result.outputs] .- t_target))
    rec = result.outputs[idx_nearest]
    aq = aqueous_concentrations(rec, result.grid, result.params, result.env)
    aq_data = (
        r = result.grid.r_grid,
        C_aq = aq.C_aq,
        O_aq = aq.O_aq
    )
    write_csv("data/aqueous_$(month)mo.csv", aq_data)
end

println()
println("="^70)
println("Theme 1 Complete")
println("="^70)
println()

# --- Qualitative validation checks ---
println("Qualitative physics checks:")
println("-" ^70)

# Check 1: POM depletion
P_initial = pools.P[1]
P_final = pools.P[end]
P_ratio = P_final / P_initial
println("  POM: $(round(P_initial, digits=1)) → $(round(P_final, digits=1)) μg-C ($(round(100*P_ratio, digits=1))% remaining)")
if P_ratio > 0.5
    @warn "  ⚠  POM not depleted enough after 60 months (expected <50%, got $(round(100*P_ratio, digits=1))%)"
else
    println("  ✓  POM depletion looks reasonable")
end

# Check 2: r_agg growth and decline
r_agg_max = maximum(pools.r_agg)
r_agg_final = pools.r_agg[end]
println("  r_agg: peak = $(round(r_agg_max, digits=3)) mm, final = $(round(r_agg_final, digits=3)) mm")
if r_agg_max < 0.3
    @warn "  ⚠  r_agg peak very small (expected 1-3 mm for 0.5mm POM, got $(round(r_agg_max, digits=3)) mm)"
    @warn "     → Check G_c calculation, F_i/E concentrations, stability criterion"
elseif r_agg_max > pools.r_agg[1] * 1.5
    println("  ✓  r_agg growth observed")
else
    @warn "  ⚠  r_agg did not grow significantly above r_0"
end

# Check 3: MAOC accumulation
M_initial = pools.M_total[1]
M_final = pools.M_total[end]
println("  MAOC: $(round(M_initial, digits=1)) → $(round(M_final, digits=1)) μg-C ($(round(M_final-M_initial, digits=1)) net increase)")
if M_final < M_initial * 1.5
    @warn "  ⚠  MAOC not accumulating much (expected continuous increase)"
else
    println("  ✓  MAOC accumulation observed")
end

# Check 4: MAOC decoupling (M_total > M_agg after r_agg declines)
idx_peak_ragg = argmax(pools.r_agg)
if idx_peak_ragg < length(pools.r_agg) - 10
    M_total_post = pools.M_total[end]
    M_agg_post = pools.M_agg[end]
    if M_total_post > M_agg_post * 1.1
        println("  ✓  MAOC decoupling: M_total > M_agg after r_agg decline")
    else
        @warn "  ⚠  MAOC decoupling not clear (M_total ≈ M_agg at end)"
    end
end

# Check 5: CO2 flux peak timing
idx_peak_co2 = argmax(flux[2:end]) + 1  # skip first point
t_peak_co2_months = result.outputs[idx_peak_co2].t / 30.44
println("  CO₂ flux: peaks at $(round(t_peak_co2_months, digits=1)) months")
if t_peak_co2_months > 12
    @warn "  ⚠  CO₂ flux peaks late (expected 3-6 months, got $(round(t_peak_co2_months, digits=1)) months)"
else
    println("  ✓  CO₂ flux peaks reasonably early")
end

# Check 6: Radial gradients in profiles
prof_6mo = profs[findfirst(p -> abs(p.t - 6*30.44) < 15, profs)]
if !isnothing(prof_6mo)
    C_gradient = maximum(prof_6mo.C) - minimum(prof_6mo.C)
    if C_gradient > 1.0
        println("  ✓  Radial C gradient visible at 6 months ($(round(C_gradient, digits=1)) μg/mm³ range)")
    else
        @warn "  ⚠  Radial profiles nearly flat (expected sharp C gradient near POM)"
    end
end

println("-" ^70)
println()
println("All CSVs exported to data/")
println("Ready for R analysis")
