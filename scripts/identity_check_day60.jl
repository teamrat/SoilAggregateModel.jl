"""
Direct Conservation Identity Check at Day 60

Evaluate compute_source_terms at high biomass and verify that:
S_C + S_B + S_Fn + S_Fm + S_Fi + S_E + S_M + Resp_total = 0

If this identity has non-zero error, there's an algebraic bug in reactions.jl.
"""

## Setup
using Pkg
Pkg.activate(@__DIR__() * "/..")

using SoilAggregateModel
import SoilAggregateModel: compute_source_terms

println("="^70)
println("Conservation Identity Check at Day 60 (High Biomass)")
println("="^70)
println()

# Standard setup
bio = BiologicalProperties()
soil = SoilProperties()

# Environmental conditions
T_func = t -> 293.15  # K
ψ_func = t -> -33.0   # kPa
O2_func = t -> 0.21   # 21% O₂

println("Phase 1: Running to day 60...")
println()

# Run adaptive timestepper to day 60
result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 60.0);
                       output_times=[0.0, 60.0])

println("Complete: $(result.diagnostics["n_steps"]) timesteps")
println()

# Extract state at day 60
state_60 = result.outputs[end].state
ψ = -33.0

# Initialize workspace at day-60 conditions
import SoilAggregateModel: Workspace, update_temperature_cache!, update_water_content!, update_effective_diffusion!

n = result.grid.n
workspace = Workspace(n)
update_temperature_cache!(workspace.f_T, 293.15, bio, soil)
update_water_content!(workspace.θ, workspace.θ_a, ψ, state_60, soil)
update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

println("="^70)
println("Phase 2: Checking conservation identity at nodes 1, 50, 100")
println("="^70)
println()

# Check nodes 1, 25, 50 (grid has n nodes)
test_nodes = [1, 25, n]

for i in test_nodes
    println("Node $i:")
    println("  C  = $(state_60.C[i])")
    println("  B  = $(state_60.B[i])")
    println("  F_n = $(state_60.F_n[i])")
    println("  F_m = $(state_60.F_m[i])")
    println("  F_i = $(state_60.F_i[i])")
    println("  E  = $(state_60.E[i])")
    println("  M  = $(state_60.M[i])")
    println("  O  = $(state_60.O[i])")
    println("  θ  = $(workspace.θ[i])")
    println("  θ_a = $(workspace.θ_a[i])")
    println()

    # Compute source terms
    sources = compute_source_terms(
        state_60.C[i], state_60.B[i], state_60.F_n[i], state_60.F_m[i],
        state_60.F_i[i], state_60.E[i], state_60.M[i], state_60.O[i],
        workspace.θ[i], workspace.θ_a[i], ψ, bio, soil, workspace.f_T
    )

    # Sum all source terms
    sum_S = sources.S_C + sources.S_B + sources.S_Fn + sources.S_Fm +
            sources.S_Fi + sources.S_E + sources.S_M

    # Conservation identity: sum_S + Resp_total should be zero
    identity_error = sum_S + sources.Resp_total

    println("  Source Terms:")
    println("    S_C  = $(sources.S_C)")
    println("    S_B  = $(sources.S_B)")
    println("    S_Fn = $(sources.S_Fn)")
    println("    S_Fm = $(sources.S_Fm)")
    println("    S_Fi = $(sources.S_Fi)")
    println("    S_E  = $(sources.S_E)")
    println("    S_M  = $(sources.S_M)")
    println("    Resp = $(sources.Resp_total)")
    println()
    println("  Sum of S_X terms = $(sum_S)")
    println("  Identity error   = $(identity_error)")
    println()

    # Compute what the per-step contribution would be
    r_i = result.grid.r_grid[i]
    h = result.grid.h
    volume_i = 4.0 * π * r_i^2 * h
    dt = 0.001

    # The per-step error in total carbon from this node would be:
    per_step_error_from_node = identity_error * volume_i * dt

    println("  Volume at node $i = $(volume_i) mm³")
    println("  Per-step carbon error from node $i = $(per_step_error_from_node) μg-C")
    println("  (at dt=0.001)")
    println()

    # Relative error
    total_rate = abs(sum_S)
    if total_rate > 0
        rel_error = abs(identity_error / total_rate)
        println("  Relative identity error = $(rel_error)")
    end
    println()
    println("-"^70)
    println()
end

# Sum up total per-step error from all nodes
println("="^70)
println("Total Per-Step Error Estimate")
println("="^70)
println()

dt = 0.001
total_per_step_error = 0.0

for i in 1:n
    global total_per_step_error
    sources = compute_source_terms(
        state_60.C[i], state_60.B[i], state_60.F_n[i], state_60.F_m[i],
        state_60.F_i[i], state_60.E[i], state_60.M[i], state_60.O[i],
        workspace.θ[i], workspace.θ_a[i], ψ, bio, soil, workspace.f_T
    )

    sum_S = sources.S_C + sources.S_B + sources.S_Fn + sources.S_Fm +
            sources.S_Fi + sources.S_E + sources.S_M
    identity_error = sum_S + sources.Resp_total

    r_i = result.grid.r_grid[i]
    h = result.grid.h
    volume_i = 4.0 * π * r_i^2 * h

    total_per_step_error += identity_error * volume_i * dt
end

println("Sum over all $(n) nodes:")
println("  Total per-step error = $(total_per_step_error) μg-C")
println()
println("Observed per-step error from Strang test: 2.49e-5 μg-C")
println("Ratio: $(total_per_step_error / 2.49e-5)")
println()

if abs(total_per_step_error - 2.49e-5) < 1e-6
    println("✓ MATCH! The identity error explains the observed leak.")
else
    println("✗ NO MATCH. The leak is elsewhere (clipping, diffusion, POM accounting).")
end
println("="^70)
