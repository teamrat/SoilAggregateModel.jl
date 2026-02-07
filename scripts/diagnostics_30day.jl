# diagnostics_30day.jl
# Run 30-day simulation and generate diagnostic table

using SoilAggregateModel
import SoilAggregateModel: Workspace, update_temperature_cache!, water_content, D_eff_DOC

using Printf
using Statistics
using Dates

# Setup parameters
bio = BiologicalProperties()
soil = SoilProperties()

# Constant environment
T(t) = 293.15  # 20°C
ψ(t) = -10.0   # kPa
O2(t) = 0.3    # μg/mm³

# Run 30-day simulation (default outputs every day or so)
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                      n_grid=50, dt_initial=0.01)

# === Compute diagnostics ===

# Grid setup
n_grid = 50
r_0 = 0.1
r_max = 2.0
h = (r_max - r_0) / (n_grid - 1)
r_grid = [r_0 + i*h for i in 0:n_grid-1]

# Find output indices closest to days 1, 10, 30
times = [r.t for r in result.outputs]
idx_1 = argmin(abs.(times .- 1.0))
idx_10 = argmin(abs.(times .- 10.0))
idx_30 = argmin(abs.(times .- 30.0))

# Create workspace for diagnostics
workspace = Workspace(n_grid)

# Function to compute C_aq and diagnostics for a given state
function compute_diagnostics(state, t_label)
    # Update temperature cache
    T_val = T(0.0)
    ψ_val = ψ(0.0)
    update_temperature_cache!(workspace.f_T, T_val, bio, soil)

    # Compute water content at each node
    θ = zeros(n_grid)
    for i in 1:n_grid
        θ[i] = water_content(ψ_val, state.E[i], state.F_i[i], soil)
    end

    # Compute C_aq at each node
    C_aq = zeros(n_grid)
    for i in 1:n_grid
        # R = θ + ρ_b·k_d (retardation factor denominator)
        R_denom = θ[i] + soil.ρ_b * soil.k_d_eq
        C_aq[i] = state.C[i] / R_denom
    end

    # Compute effective diffusion coefficient for DOC
    D_C_eff = zeros(n_grid)
    for i in 1:n_grid
        D_C_eff[i] = D_eff_DOC(workspace.f_T.D_DOC_w, θ[i],
                              soil.θ_s, soil.ρ_b, soil.k_d_eq)
    end

    # Retardation factor (average)
    R_avg = mean(θ .+ soil.ρ_b * soil.k_d_eq)
    retardation = R_avg / mean(θ)

    # Biomass (total across grid, volumetrically weighted)
    B_total = 0.0
    for i in 1:n_grid
        W_i = 4.0 * π * r_grid[i]^2 * h
        B_total += state.B[i] * W_i
    end

    # MAOC (total across grid)
    M_total = 0.0
    for i in 1:n_grid
        W_i = 4.0 * π * r_grid[i]^2 * h
        M_total += state.M[i] * W_i
    end

    return (
        t = t_label,
        C_aq_min = minimum(C_aq),
        C_aq_max = maximum(C_aq),
        C_aq_mean = mean(C_aq),
        D_C_eff_min = minimum(D_C_eff),
        D_C_eff_max = maximum(D_C_eff),
        D_C_eff_mean = mean(D_C_eff),
        retardation = retardation,
        B_total = B_total,
        M_total = M_total,
        CO2_cum = state.CO2_cumulative
    )
end

# Compute diagnostics at each timepoint (use actual times from output)
diag_1 = compute_diagnostics(result.outputs[idx_1].state, result.outputs[idx_1].t)
diag_10 = compute_diagnostics(result.outputs[idx_10].state, result.outputs[idx_10].t)
diag_30 = compute_diagnostics(result.outputs[idx_30].state, result.outputs[idx_30].t)

# === Generate diagnostic report ===

output_file = "diagnostics_30day_output.txt"
open(output_file, "w") do io
    println(io, "="^80)
    println(io, "30-DAY SIMULATION DIAGNOSTIC REPORT")
    println(io, "="^80)
    println(io)
    println(io, "Date: ", Dates.now())
    println(io, "Grid: ", n_grid, " points, r_0 = ", r_0, " mm, r_max = ", r_max, " mm")
    println(io, "Environment: T = 20°C (293.15 K), ψ = -10 kPa, O₂ = 0.3 μg/mm³")
    println(io)

    println(io, "-"^80)
    println(io, "CRITICAL PARAMETERS (after fixes)")
    println(io, "-"^80)
    @printf(io, "ρ_b         = %.1f μg/mm³ (FIXED: was 1.3e6)\n", soil.ρ_b)
    @printf(io, "ρ_POM       = %.1f μg-C/mm³ (FIXED: was 1.0e6)\n", bio.ρ_POM)
    @printf(io, "k_d_eq      = %.3f mm³/μg\n", soil.k_d_eq)
    @printf(io, "λ           = %.2f (FIXED: was 0.5)\n", bio.λ)
    println(io, "K_Y units   = [μg-C/mm³/day] (FIXED: was documented as [-])")
    println(io)

    println(io, "-"^80)
    println(io, "RETARDATION FACTOR & DIFFUSION")
    println(io, "-"^80)
    @printf(io, "ρ_b·k_d     = %.1f\n", soil.ρ_b * soil.k_d_eq)
    @printf(io, "Retardation = %.1f (R = (θ + ρ_b·k_d) / θ)\n", diag_1.retardation)
    @printf(io, "C_aq/C      ≈ 1/%.1f (aqueous fraction)\n", diag_1.retardation)
    println(io)
    @printf(io, "D_C_eff range: %.4e - %.4e mm²/day\n", diag_1.D_C_eff_min, diag_1.D_C_eff_max)
    @printf(io, "D_C_eff mean:  %.4e mm²/day\n", diag_1.D_C_eff_mean)
    println(io)

    println(io, "-"^80)
    println(io, "C_aq CONCENTRATIONS (Aqueous DOC)")
    println(io, "-"^80)
    println(io, "Time       Min [μg/mm³]  Mean [μg/mm³]  Max [μg/mm³]")
    println(io, "-"^60)
    @printf(io, "Day %.1f   %.4e      %.4e       %.4e\n",
            diag_1.t, diag_1.C_aq_min, diag_1.C_aq_mean, diag_1.C_aq_max)
    @printf(io, "Day %.1f  %.4e      %.4e       %.4e\n",
            diag_10.t, diag_10.C_aq_min, diag_10.C_aq_mean, diag_10.C_aq_max)
    @printf(io, "Day %.1f  %.4e      %.4e       %.4e\n",
            diag_30.t, diag_30.C_aq_min, diag_30.C_aq_mean, diag_30.C_aq_max)
    println(io)

    println(io, "-"^80)
    println(io, "BIOMASS DYNAMICS (Total B, volumetrically integrated)")
    println(io, "-"^80)
    println(io, "Time       B_total [μg-C]")
    println(io, "-"^35)
    @printf(io, "Day %.1f   %.4e\n", diag_1.t, diag_1.B_total)
    @printf(io, "Day %.1f  %.4e\n", diag_10.t, diag_10.B_total)
    @printf(io, "Day %.1f  %.4e\n", diag_30.t, diag_30.B_total)
    println(io)

    println(io, "-"^80)
    println(io, "MAOC ACCUMULATION (Total M, volumetrically integrated)")
    println(io, "-"^80)
    println(io, "Time       M_total [μg-C]")
    println(io, "-"^35)
    @printf(io, "Day %.1f   %.4e\n", diag_1.t, diag_1.M_total)
    @printf(io, "Day %.1f  %.4e\n", diag_10.t, diag_10.M_total)
    @printf(io, "Day %.1f  %.4e\n", diag_30.t, diag_30.M_total)
    println(io)

    println(io, "-"^80)
    println(io, "CUMULATIVE CO₂ RESPIRATION")
    println(io, "-"^80)
    println(io, "Time       CO₂ [μg-C]")
    println(io, "-"^30)
    @printf(io, "Day %.1f   %.4e\n", diag_1.t, diag_1.CO2_cum)
    @printf(io, "Day %.1f  %.4e\n", diag_10.t, diag_10.CO2_cum)
    @printf(io, "Day %.1f  %.4e\n", diag_30.t, diag_30.CO2_cum)
    println(io)

    println(io, "-"^80)
    println(io, "CARBON CONSERVATION")
    println(io, "-"^80)
    println(io, "Time       Mass Balance Error (relative)")
    println(io, "-"^45)
    @printf(io, "Day %.1f   %.4e\n", diag_1.t, result.outputs[idx_1].mass_balance_error)
    @printf(io, "Day %.1f  %.4e\n", diag_10.t, result.outputs[idx_10].mass_balance_error)
    @printf(io, "Day %.1f  %.4e\n", diag_30.t, result.outputs[idx_30].mass_balance_error)
    println(io)

    # Verification checks from PARAMETER_BUG_REPORT.md
    println(io, "-"^80)
    println(io, "VERIFICATION CHECKS (from PARAMETER_BUG_REPORT.md)")
    println(io, "-"^80)

    # Check 1: C_aq in reasonable range
    c_aq_ok = 0.1 <= diag_30.C_aq_mean <= 10.0
    println(io, "✓ C_aq range check: ", c_aq_ok ? "PASS" : "FAIL")
    @printf(io, "  Expected: 0.1-10 μg/mm³, Observed: %.4e μg/mm³\n", diag_30.C_aq_mean)
    println(io)

    # Check 2: D_C_eff in reasonable range
    d_eff_ok = 0.01 <= diag_30.D_C_eff_mean <= 0.1
    println(io, "✓ D_C_eff check: ", d_eff_ok ? "PASS" : "FAIL")
    @printf(io, "  Expected: 0.01-0.1 mm²/day, Observed: %.4e mm²/day\n", diag_30.D_C_eff_mean)
    println(io)

    # Check 3: Carbon conservation
    carbon_ok = abs(result.outputs[idx_30].mass_balance_error) < 1e-12
    println(io, "✓ Carbon conservation: ", carbon_ok ? "PASS" : "FAIL")
    @printf(io, "  Expected: < 1e-12, Observed: %.4e\n", result.outputs[idx_30].mass_balance_error)
    println(io)

    # Check 4: Retardation factor
    ret_note = diag_1.retardation > 1000 ? " (HIGH - see Bug #4)" : " (reasonable)"
    println(io, "✓ Retardation factor: ", @sprintf("%.1f", diag_1.retardation), ret_note)
    println(io)

    # Check 5: Biomass dynamics
    biomass_growing = diag_10.B_total > diag_1.B_total
    println(io, "✓ Biomass growth: ", biomass_growing ? "PASS" : "FAIL")
    @printf(io, "  Day 1→10: %.4e → %.4e\n", diag_1.B_total, diag_10.B_total)
    println(io)

    # Check 6: CO2 accumulation
    co2_positive = diag_30.CO2_cum > 0.0
    println(io, "✓ CO₂ accumulation: ", co2_positive ? "PASS" : "FAIL")
    @printf(io, "  Day 30: %.4e μg-C\n", diag_30.CO2_cum)
    println(io)

    println(io, "-"^80)
    println(io, "SIMULATION DIAGNOSTICS")
    println(io, "-"^80)
    @printf(io, "Total timesteps: %d\n", result.diagnostics["n_steps"])
    @printf(io, "Final time:      %.6f days\n", result.diagnostics["final_time"])
    @printf(io, "Output records:  %d\n", length(result.outputs))
    println(io)

    println(io, "="^80)
    println(io, "END OF REPORT")
    println(io, "="^80)
end

println("Diagnostic report saved to: ", output_file)
println()
println("Summary:")
@printf("  C_aq range:     %.4e - %.4e μg/mm³\n", diag_30.C_aq_min, diag_30.C_aq_max)
@printf("  D_C_eff mean:   %.4e mm²/day\n", diag_30.D_C_eff_mean)
@printf("  Retardation:    %.1f\n", diag_1.retardation)
@printf("  Biomass day 30: %.4e μg-C\n", diag_30.B_total)
@printf("  MAOC day 30:    %.4e μg-C\n", diag_30.M_total)
@printf("  CO₂ day 30:     %.4e μg-C\n", diag_30.CO2_cum)
@printf("  C balance err:  %.4e\n", result.outputs[idx_30].mass_balance_error)
