"""
High-Biomass Per-Step Strang Conservation Test

Run to day 60 with adaptive timestepper to get realistic high biomass state,
then run 100 fixed-dt Strang steps with full instrumentation.

This tests whether conservation breaks down at high biomass levels where
reaction rates are large and clipping events may occur.
"""

## Setup
using Pkg
Pkg.activate(@__DIR__() * "/..")

using Revise
using SoilAggregateModel
import SoilAggregateModel: diffusion_step!, reaction_step!, GridInfo
import SoilAggregateModel: J_P, R_P
import SoilAggregateModel: update_temperature_cache!, update_water_content!, update_effective_diffusion!
import SoilAggregateModel: Workspace, AggregateState

# Helper function to compute total carbon
function total_carbon(state, grid)
    total = state.P + state.CO2_cumulative
    for i in 1:grid.n
        total += (state.C[i] + state.B[i] + state.F_i[i] + state.F_n[i] +
                  state.F_m[i] + state.E[i] + state.M[i]) * grid.W[i]
    end
    return total
end

println("="^70)
println("High-Biomass Per-Step Strang Conservation Test")
println("="^70)
println()

## Standard setup
bio = BiologicalProperties()
soil = SoilProperties()

# Environmental conditions
T_func = t -> 293.15  # K
ψ_func = t -> -33.0   # kPa
O2_func = t -> 0.21   # 21% O₂

println("Phase 1: Running adaptive timestepper to day 60...")
println()

## Run adaptive timestepper to day 60
result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 60.0);
                       output_times=[0.0, 60.0])

println("Adaptive run complete: $(result.diagnostics["n_steps"]) timesteps")
println()

# Extract state at day 60 from outputs
state_60 = result.outputs[end].state
T_60 = 293.15
ψ_60 = -33.0
O2_amb = 0.21

# Report biomass levels
println("Day-60 State:")
println("  P = $(round(state_60.P, digits=4)) μg-C")
println("  B[1] = $(round(state_60.B[1], digits=6)) μg-C/mm³")
println("  C[1] = $(round(state_60.C[1], digits=4)) μg-C/mm³")
println("  F_i[1] = $(round(state_60.F_i[1], digits=6)) μg-C/mm³")
println("  Total C = $(round(total_carbon(state_60, result.grid), digits=4)) μg-C")
println()

# Create fresh grid and workspace for fixed-dt run
n = result.grid.n
grid = GridInfo(n, 0.25, 5.0)

# Deep copy the state
state = AggregateState(
    copy(state_60.C), copy(state_60.B),
    copy(state_60.F_n), copy(state_60.F_m), copy(state_60.O),
    copy(state_60.F_i), copy(state_60.E), copy(state_60.M),
    state_60.P,
    state_60.CO2_cumulative
)

# Initialize workspace at day-60 conditions
workspace = Workspace(n)
update_temperature_cache!(workspace.f_T, T_60, bio, soil)
update_water_content!(workspace.θ, workspace.θ_a, ψ_60, state, soil)
update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

println("="^70)
println("Phase 2: Running 100 instrumented Strang steps at dt=0.001")
println("="^70)
println()

dt = 0.001
n_steps = 100

# Track maximum errors
max_rxn_residual = 0.0
max_full_cycle_err = 0.0

# Sample steps to report
sample_steps = [1, 2, 5, 10, 20, 50, 100]

for step in 1:n_steps
    global max_rxn_residual, max_full_cycle_err

    # Compute POM flux ONCE per step
    B_0 = state.B[1]
    F_n_0 = state.F_n[1]
    θ_0 = workspace.θ[1]
    θ_a_0 = workspace.θ_a[1]
    O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
    R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
    J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                  bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
    R_P_val = R_P(J_P_val, grid.r_grid[1])

    # Expected flux per diffusion half-step
    expected_flux_half = J_P_val * 4.0 * π * grid.r_grid[1]^2 * (dt/2)

    # Initial carbon for this full Strang cycle
    C0 = total_carbon(state, grid)

    # === Diffusion half-step 1 ===
    diffusion_step!(state, workspace, dt/2, grid.r_grid, grid.h, J_P_val, O2_amb)
    C1 = total_carbon(state, grid)
    diff1_err = (C1 - C0) - expected_flux_half

    # === Reaction step ===
    reaction_step!(state, workspace, dt, grid.r_grid, grid.h, bio, soil, ψ_60, R_P_val)
    C2 = total_carbon(state, grid)

    # The reaction step should conserve exactly EXCEPT for the POM decrease
    # C2 - C1 should equal -R_P_val * dt (the POM that disappeared)
    rxn_change = C2 - C1
    expected_pom_decrease = -R_P_val * dt
    rxn_residual = rxn_change - expected_pom_decrease

    # === Diffusion half-step 2 ===
    diffusion_step!(state, workspace, dt/2, grid.r_grid, grid.h, J_P_val, O2_amb)
    C3 = total_carbon(state, grid)
    diff2_err = (C3 - C2) - expected_flux_half

    # Full Strang cycle error
    # Expected: C3 - C0 = 2*flux - R_P*dt = 2*flux - POM_decrease
    expected_increase = 2.0 * expected_flux_half + expected_pom_decrease
    full_cycle_err = (C3 - C0) - expected_increase

    # Update workspace for next step
    update_water_content!(workspace.θ, workspace.θ_a, ψ_60, state, soil)
    update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

    # Track maximum errors
    max_rxn_residual = max(max_rxn_residual, abs(rxn_residual))
    max_full_cycle_err = max(max_full_cycle_err, abs(full_cycle_err))

    # Report sample steps
    if step in sample_steps
        println("Step $step:")
        println("  Diff1 error:      $(round(diff1_err, sigdigits=4))")
        println("  Rxn change:       $(round(rxn_change, sigdigits=6))")
        println("  Expected POM dec: $(round(expected_pom_decrease, sigdigits=6))")
        println("  Rxn residual:     $(round(rxn_residual, sigdigits=4))  ← should be ~0")
        println("  Diff2 error:      $(round(diff2_err, sigdigits=4))")
        println("  Full cycle error: $(round(full_cycle_err, sigdigits=4))  ← should be ~0")
        println()
    end
end

println("="^70)
println("Maximum Errors Across 100 Steps at High Biomass:")
println("="^70)
println()
println("Max reaction residual:   $(round(max_rxn_residual, sigdigits=6))")
println("Max full-cycle error:    $(round(max_full_cycle_err, sigdigits=6))")
println()

# Final conservation check
C_final = total_carbon(state, grid)
C_initial_60 = total_carbon(state_60, grid)
rel_err = (C_final - C_initial_60) / C_initial_60

println("="^70)
println("Overall 100-Step Run from Day 60:")
println("="^70)
println()
println("Initial (day 60): $(round(C_initial_60, digits=4)) μg-C")
println("Final (day 60.1): $(round(C_final, digits=4)) μg-C")
println("Net change:       $(round(C_final - C_initial_60, digits=6)) μg-C")
println("Relative error:   $(round(rel_err * 100, sigdigits=6))%")
println()

# Interpretation
println("="^70)
println("Interpretation:")
println("="^70)
if max_rxn_residual < 1e-10
    println("✓ Reaction step conserves at high biomass (< 1e-10 residual)")
else
    println("✗ LEAK DETECTED: Reaction step has $(round(max_rxn_residual, sigdigits=4)) residual")
end

if max_full_cycle_err < 1e-10
    println("✓ Full Strang cycle conserves at high biomass (< 1e-10 error)")
else
    println("✗ LEAK DETECTED: Strang cycle has $(round(max_full_cycle_err, sigdigits=4)) error per step")
end
println("="^70)
