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

# Compute integrated pools (needed for diagnostics)
pools = integrated_pools(result)

# --- Fungal transition diagnostics at POM surface (node 0) ---
println("="^70)
println("Fungal Transition Diagnostics at POM Surface (r = r_0)")
println("="^70)
import SoilAggregateModel: Pi_protected, fungal_transitions

diagnostic_times = [0.0, 1.0, 7.0, 30.0, 6*30.44]  # 0, 1 day, 7 days, 30 days, 6 months
diagnostic_indices = Int[]
for t_target in diagnostic_times
    # Find closest output time
    idx = argmin(abs.(pools.t .- t_target))
    push!(diagnostic_indices, idx)
end

for (i, idx) in enumerate(diagnostic_indices)
    rec = result.outputs[idx]
    t = rec.t
    t_days = round(t, digits=2)
    t_months = round(t / 30.44, digits=2)

    # Extract state at node 0 (POM surface)
    F_i_0 = rec.state.F_i[1]
    F_n_0 = rec.state.F_n[1]
    F_m_0 = rec.state.F_m[1]

    # Compute protection ratio
    ε_F = STANDARD_BIO.ε_F
    Π = Pi_protected(F_m_0, F_i_0, F_n_0, ε_F)

    # Compute transition terms (temperature-corrected)
    env_vals = result.env.T(t), result.env.ψ(t), result.env.O2(t)
    T = env_vals[1]
    f_T_bac = exp(-STANDARD_BIO.Ea_B / 8.314 * (1/T - 1/STANDARD_BIO.T_ref))
    f_T_fun = exp(-STANDARD_BIO.Ea_F / 8.314 * (1/T - 1/STANDARD_BIO.T_ref))

    α_i_T = STANDARD_BIO.α_i * f_T_fun
    α_n_T = STANDARD_BIO.α_n * f_T_fun
    β_i_T = STANDARD_BIO.β_i * f_T_fun
    β_n_T = STANDARD_BIO.β_n * f_T_fun
    ζ_T = STANDARD_BIO.ζ * f_T_fun

    trans = fungal_transitions(F_i_0, F_n_0, F_m_0, Π, α_i_T, α_n_T, β_i_T, β_n_T,
                              ζ_T, STANDARD_BIO.delta, STANDARD_BIO.η_conv, ε_F)

    # Compute net rates manually for clarity
    Π_delta = Π^STANDARD_BIO.delta
    net_i_raw = β_i_T * Π - α_i_T * Π_delta
    net_n_raw = β_n_T * Π - α_n_T * Π_delta

    println()
    println("Time: $(t_days) days ($(t_months) months)")
    println("  Fungal pools [μg/mm³]:")
    println("    F_i = $(round(F_i_0, digits=6))")
    println("    F_n = $(round(F_n_0, digits=6))")
    println("    F_m = $(round(F_m_0, digits=6))")
    println("  Protection ratio: Π = $(round(Π, digits=6))")
    println("  Net transition rates [μg/mm³/day]:")
    println("    net_i (raw) = β_i·Π - α_i·Π^δ = $(round(net_i_raw * F_i_0, digits=8))")
    println("    net_n (raw) = β_n·Π - α_n·Π^δ = $(round(net_n_raw * F_n_0, digits=8))")
    println("  Actual fluxes:")
    println("    trans_i (η·net_i·F_i) = $(round(trans.trans_i, digits=8))")
    println("    trans_n (η·net_n·F_n) = $(round(trans.trans_n, digits=8))")
    println("    insulation (ζ·F_n → F_i) = $(round(trans.insulation, digits=8))")
    println("  Net change to F_n [μg/mm³/day]:")
    println("    dF_n/dt ≈ trans_n - insulation = $(round(trans.trans_n - trans.insulation, digits=8))")

    # === Bacterial budget diagnostics ===
    println()
    println("  Bacterial dynamics at node 0:")

    # Get state at node 0
    B_0 = rec.state.B[1]
    C_0 = rec.state.C[1]
    O_0 = rec.state.O[1]

    # Compute environment at node 0
    soil = STANDARD_SOIL
    ψ_val = STANDARD_ψ

    # Water content
    s_e = (1.0 + abs(soil.α_vg * ψ_val)^soil.n_vg)^(-(soil.n_vg - 1.0) / soil.n_vg)
    θ = soil.θ_r + (soil.θ_s - soil.θ_r) * s_e
    θ_a = soil.θ_s - θ

    # Aqueous concentrations
    C_aq = C_0 / (θ + soil.ρ_b * soil.k_d_eq)
    K_H_O = 0.032
    O_aq = O_0 * θ / (θ + K_H_O * θ_a)

    # Bacterial rates
    r_B_max_T = STANDARD_BIO.r_B_max * f_T_bac

    # R_B: gross uptake
    import SoilAggregateModel: R_B, R_Bb, Y_B_func, Gamma_B, R_rec_B, h_B
    R_B_val = R_B(C_aq, O_aq, B_0, ψ_val, r_B_max_T, STANDARD_BIO.K_B,
                  STANDARD_BIO.L_B, STANDARD_BIO.ν_B)

    # R_Bb: basal metabolic uptake
    R_Bb_val = R_Bb(STANDARD_BIO.C_B, O_aq, B_0, r_B_max_T, STANDARD_BIO.K_B,
                    STANDARD_BIO.L_B, STANDARD_BIO.B_min)

    # R_diff: uptake available for growth
    R_diff = R_B_val - R_Bb_val

    # Y_B: flexible yield (space-limited)
    Y_B_val = Y_B_func(R_diff, STANDARD_BIO.Y_B_max, STANDARD_BIO.K_Y, B_0, STANDARD_BIO.B_S, STANDARD_BIO.ε_Y)

    # Γ_B: bacterial growth
    Γ_B_val = Gamma_B(R_B_val, R_Bb_val, Y_B_val, STANDARD_BIO.γ, STANDARD_BIO.ε_Y)

    # R_rec_B: bacterial death
    μ_B_T = STANDARD_BIO.μ_B * f_T_bac
    R_rec_B_val = R_rec_B(μ_B_T, B_0, STANDARD_BIO.B_min)

    # Net dB/dt
    dB_dt = Γ_B_val - R_rec_B_val

    println("    B = $(round(B_0, digits=6)) μg/mm³")
    println("    R_B (gross uptake) = $(round(R_B_val, digits=8)) μg-C/mm³/day")
    println("    R_Bb (basal metabolism) = $(round(R_Bb_val, digits=8)) μg-C/mm³/day")
    println("    R_diff (for growth) = $(round(R_diff, digits=8)) μg-C/mm³/day")
    println("    Y_B (flexible yield) = $(round(Y_B_val, digits=6))")
    println("    Γ_B (growth) = $(round(Γ_B_val, digits=8)) μg-C/mm³/day")
    println("    R_rec_B (death) = $(round(R_rec_B_val, digits=8)) μg-C/mm³/day")
    println("    Net dB/dt = Γ_B - R_rec_B = $(round(dB_dt, digits=8)) μg-C/mm³/day")

    # === F_m budget diagnostics ===
    if i == 1  # Only for t=0
        println()
        println("  F_m budget analysis at t=0:")

        # Compute fungal uptake R_F at node 0
        soil = STANDARD_SOIL
        ψ_val = STANDARD_ψ  # -33 kPa (field capacity)

        # Van Genuchten water content
        s_e = (1.0 + abs(soil.α_vg * ψ_val)^soil.n_vg)^(-(soil.n_vg - 1.0) / soil.n_vg)
        θ = soil.θ_r + (soil.θ_s - soil.θ_r) * s_e
        θ_a = soil.θ_s - θ

        C = rec.state.C[1]
        O = rec.state.O[1]

        # Aqueous concentrations
        C_aq = C / (θ + soil.ρ_b * soil.k_d_eq)
        K_H_O = 0.032  # Approximate Henry's constant at 20°C
        O_aq = O * θ / (θ + K_H_O * θ_a)

        # Fungal uptake R_F
        r_F_max_T = STANDARD_BIO.r_F_max * f_T_fun
        R_F_val = r_F_max_T * (C_aq / (STANDARD_BIO.K_F + C_aq)) *
                  (O_aq / (STANDARD_BIO.L_F + O_aq)) *
                  (F_i_0 + STANDARD_BIO.λ * F_n_0) *
                  exp(STANDARD_BIO.ν_F * STANDARD_ψ)

        # Γ_F (growth entering F_m)
        Y_F_val = STANDARD_BIO.Y_F
        Γ_F_val = Y_F_val * R_F_val

        # Mobilization flux into F_m
        # When net_i < 0 or net_n < 0, mass flows into F_m
        # Total mobilization = -(trans_i + trans_n) when negative
        mobilization_source = 0.0
        if trans.trans_i < 0
            mobilization_source -= trans.trans_i
        end
        if trans.trans_n < 0
            mobilization_source -= trans.trans_n
        end

        # Diffusive flux out of node 0 (spherical geometry)
        # J_diff ≈ -D × ∂F_m/∂r at r_0
        # Use finite difference: (F_m[2] - F_m[1]) / h
        h = result.grid.h
        F_m_1 = rec.state.F_m[2]
        dFm_dr = (F_m_1 - F_m_0) / h

        # Diffusion coefficient for F_m
        D_Fm = STANDARD_BIO.D_Fm0  # mm²/day

        # Flux density [μg/mm²/day]
        J_diff = -D_Fm * dFm_dr

        # Total flux through spherical surface at r_0 [μg/day]
        r_0 = result.grid.r_0
        A_0 = 4.0 * π * r_0^2
        Flux_diff_total = J_diff * A_0

        # Convert to volumetric rate [μg/mm³/day] in node 0
        # Volume of shell 0: V_0 ≈ 4πr_0²h
        V_0 = A_0 * h
        diffusive_loss_rate = Flux_diff_total / V_0

        println("    Sources [μg/mm³/day]:")
        println("      Γ_F (fungal growth) = $(round(Γ_F_val, digits=8))")
        println("      Mobilization (from F_i+F_n) = $(round(mobilization_source, digits=8))")
        println("    Sinks [μg/mm³/day]:")
        println("      Diffusive loss = $(round(diffusive_loss_rate, digits=8))")
        println("      Immobilization (β·Π terms) = $(round(max(0, trans.trans_i + trans.trans_n), digits=8))")
        println()
        println("    C_aq at node 0 = $(round(C_aq, digits=8)) μg/mm³")
        println("    O_aq at node 0 = $(round(O_aq, digits=8)) μg/mm³")
        println("    Fungal uptake biomass (F_i + λ·F_n) = $(round(F_i_0 + STANDARD_BIO.λ * F_n_0, digits=8)) μg/mm³")
        println("    R_F = $(round(R_F_val, digits=8)) μg-C/mm³/day")
    end
end
println("="^70)
println()

# --- CSV 1: Integrated pools timeseries (both domains) ---
println("Exporting timeseries...")

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
