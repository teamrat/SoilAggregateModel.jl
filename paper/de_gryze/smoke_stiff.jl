"""
smoke_stiff.jl — first run of the stiff solver, deliberately small
==================================================================

Purpose is to find out whether it LOADS, COMPILES and RUNS, in that order, and
only then whether it agrees. Everything here is short on purpose: if the RHS has
a typo, better to learn it in 30 s than in 40.

Reports, in order:
  0. package load and first-solve latency (the numbers we owe you on cost)
  1. Jacobian sparsity — states, nonzero fraction
  2. a 0.5-day stiff solve: steps, RHS evaluations, Jacobians, wall time
  3. the same 0.5 days through the existing split solver
  4. pool-by-pool agreement between them

Step 4 is NOT expected to be exact. mol.jl lists three deliberate differences:
the split solver lags theta and D by one step, it clips negatives into CO2, and
it accumulates CO2 as a running sum. Disagreement at the 1e-3 level is the
lagging; disagreement at the 1e-1 level is a bug.

Usage:
    julia --project=. paper/de_gryze/smoke_stiff.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

print("loading SoilAggregateModel ... ")
t_load = @elapsed using SoilAggregateModel
println("$(round(t_load, digits=1)) s")

using Printf
const SAM = SoilAggregateModel

include(joinpath(@__DIR__, "degryze_soils.jl"))

# ── Setup: identical to run_degryze.jl / diagnose_speed.jl ───────────────────

const SOIL_ID = 3
const DIAM    = 1.25
const N_GRID  = 200

soil = degryze_soil(SOIL_ID; k_L = 1000, D_B_rel = 0.00001)
ic   = degryze_ic(SOIL_ID, soil; s_M = 0.6)

bio = BiologicalProperties(
    κ_s_ref = 0.01, κ_d_ref = 0.001,
    F_i_min = 1e-6, F_n_min = 2e-4, F_m_min = 1e-6,
    D_Fn0   = 0.00001, D_Fm0 = 0.001,
    r_B_max = 8.0, r_F_max = 0.2, R_P_max = 2.5,
    Y_B_max = 0.4, B_S = 0.05, C_B = 5.0e-5,
    μ_B = 0.0036, μ_F = 0.02,
)

T_const  = ic.T_0
ψ_const  = ic.ψ_0
O2_const = DEGRYZE_INCUBATION.O2_frac * 101000.0 * 0.032 / (8.314 * T_const)

# NOTE: const, deliberately. A closure over a non-const global makes T, ψ and O2
# type-unstable, which turns every call in update_temperature_cache! into a
# dynamic dispatch. run_degryze.jl has this issue; it is one candidate for the
# 2048 bytes/step the allocation probe found.
const TF  = t -> T_const
const PSI = t -> ψ_const
const O2F = t -> O2_const

tess  = SAM.domain_tessellation(ρ_POM = 200.0, I_input = 4.43e-3, ρ_b = soil.ρ_b)
r_0   = DIAM / 2.0
r_max = DIAM * tess.f_domain / 2.0
P_0   = (4.0/3.0) * π * r_0^3 * 200.0

T_END = 0.5    # days — deliberately short

println("="^76)
println("STIFF SOLVER SMOKE TEST — soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID)")
println("="^76)
@printf("  states = %d   t_end = %.2f d\n\n", 8*N_GRID + 2, T_END)

# ── 1-2. stiff solve ─────────────────────────────────────────────────────────

println("── stiff solver (FBDF, sparse Jacobian, KLU) ──")
t_first = @elapsed res_stiff = run_aggregate_stiff(
    bio, soil, TF, PSI, O2F, (0.0, T_END);
    n_grid = N_GRID, r_0 = r_0, r_max = r_max,
    ic = ic, P_0 = P_0, ω = tess.ω,
    output_times = [0.0, T_END], verbose = true)
@printf("  first call (includes compilation): %.1f s\n", t_first)

t_warm = @elapsed res_stiff = run_aggregate_stiff(
    bio, soil, TF, PSI, O2F, (0.0, T_END);
    n_grid = N_GRID, r_0 = r_0, r_max = r_max,
    ic = ic, P_0 = P_0, ω = tess.ω,
    output_times = [0.0, T_END])
d = res_stiff.diagnostics
@printf("  warm call: %.2f s   accepted %d   rejected %d   f-evals %d   jacobians %d\n",
        t_warm, d["n_accept"], d["n_reject"], d["n_f"], d["n_jac"])
@printf("  retcode %s   mean Δt %.4g d\n\n",
        d["retcode"], T_END / max(d["n_accept"], 1))

# ── 3. split solver, same span ───────────────────────────────────────────────

println("── split solver (Strang + Crank-Nicolson + Forward Euler) ──")
run_aggregate(bio, soil, TF, PSI, O2F, (0.0, 0.02);
              n_grid = N_GRID, r_0 = r_0, r_max = r_max,
              ic = ic, P_0 = P_0, ω = tess.ω, output_times = [0.0, 0.02])
t_split = @elapsed res_split = run_aggregate(
    bio, soil, TF, PSI, O2F, (0.0, T_END);
    n_grid = N_GRID, r_0 = r_0, r_max = r_max,
    ic = ic, P_0 = P_0, ω = tess.ω,
    output_times = [0.0, T_END])
@printf("  warm call: %.2f s   steps %d   mean Δt %.4g d\n\n",
        t_split, res_split.diagnostics["n_steps"],
        T_END / res_split.diagnostics["n_steps"])

@printf("  SPEEDUP: %.1fx  (%.2f s -> %.2f s)\n\n", t_split / t_warm, t_split, t_warm)

# ── 4. agreement ─────────────────────────────────────────────────────────────

println("── agreement at t = $(T_END) d ──")
a = res_stiff.outputs[end].state
b = res_split.outputs[end].state

function relerr(x, y)
    d = maximum(abs.(x .- y))
    s = max(maximum(abs.(y)), eps())
    return d / s
end

@printf("  %-6s %14s %14s %12s\n", "pool", "stiff max", "split max", "rel diff")
for (nm, fa, fb) in (("C", a.C, b.C), ("B", a.B, b.B), ("F_n", a.F_n, b.F_n),
                     ("F_m", a.F_m, b.F_m), ("O", a.O, b.O), ("F_i", a.F_i, b.F_i),
                     ("E", a.E, b.E), ("M", a.M, b.M))
    @printf("  %-6s %14.6g %14.6g %12.3g\n", nm, maximum(fa), maximum(fb), relerr(fa, fb))
end
@printf("  %-6s %14.6g %14.6g %12.3g\n", "P", a.P, b.P, abs(a.P - b.P) / max(abs(b.P), eps()))
@printf("  %-6s %14.6g %14.6g %12.3g\n", "CO2", a.CO2_cumulative, b.CO2_cumulative,
        abs(a.CO2_cumulative - b.CO2_cumulative) / max(abs(b.CO2_cumulative), eps()))

@printf("\n  carbon balance error:  stiff %.3g   split %.3g\n",
        res_stiff.outputs[end].mass_balance_error,
        res_split.outputs[end].mass_balance_error)

println("\n" * "="^76)
println("done")
println("="^76)
