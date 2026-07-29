"""
diagnose_speed.jl — why the De Gryze run takes 40 s per aggregate
=================================================================

Read-only diagnostic. It changes nothing in `src/`; it calls the same functions
the solver calls and reports what they say. Run it, paste the output back.

The run costs about 390 000 steps for 45 days, i.e. Δt pinned at the 1e-4 floor
for the entire simulation (`output/summary.csv`: n_steps = 391 773, wall = 39 s
for D = 1.25 mm). Two independent things could be responsible and they need
separate evidence:

  PART 1  step count      which species and which node set max_rel_change, and
                          what Δt each species would allow on its own
  PART 2  cost per step   where the ~100 µs per step actually goes
  PART 3  allocation      whether the "zero allocation" claim in the solver
                          docstrings still holds

Usage:
    julia --project=. paper/de_gryze/diagnose_speed.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using SoilAggregateModel
using Printf
using Profile

const SAM = SoilAggregateModel

include(joinpath(@__DIR__, "degryze_soils.jl"))

# ── Setup: identical to run_degryze.jl for soil 3, middle diameter ────────────

const SOIL_ID = 3
const DIAM    = 1.25          # mm — the modal POM class
const N_GRID  = 200

soil = degryze_soil(SOIL_ID; k_L = 1000, D_B_rel = 0.00001)
ic   = degryze_ic(SOIL_ID, soil; s_M = 0.6)
# These MUST match run_degryze.jl exactly — a diagnostic of a different
# parameter set answers a different question.
bio  = BiologicalProperties(
    κ_s_ref = 0.01,
    κ_d_ref = 0.001,
    F_i_min = 1e-6,
    F_n_min = 2e-4,
    F_m_min = 1e-6,
    D_Fn0   = 0.00001,
    D_Fm0   = 0.001,
    r_B_max = 8.0,
    r_F_max = 0.2,
    R_P_max = 2.5,
    Y_B_max = 0.4,
    B_S     = 0.05,
    C_B     = 5.0e-5,
    μ_B     = 0.0036,
    μ_F     = 0.02,
)

T_const  = ic.T_0
ψ_const  = ic.ψ_0
O2_const = DEGRYZE_INCUBATION.O2_frac * 101000.0 * 0.032 / (8.314 * T_const)

T_func  = t -> T_const
ψ_func  = t -> ψ_const
O2_func = t -> O2_const

ρ_POM   = 200.0
I_input = 4.43e-3
tess    = SAM.domain_tessellation(ρ_POM = ρ_POM, I_input = I_input, ρ_b = soil.ρ_b)
ω       = tess.ω
domain_factor = tess.f_domain

r_0   = DIAM / 2.0
r_max = DIAM * domain_factor / 2.0
P_0   = (4.0/3.0) * π * r_0^3 * ρ_POM

const SPECIES = ("C", "B", "F_n", "F_m", "F_i", "E", "M", "O")
const DT_PROBE = 1e-4          # the floor the run is pinned at
const THRESHOLD = 1e-6         # reaction_step.jl:49

println("="^76)
println("SPEED DIAGNOSTIC — soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID)")
println("="^76)
@printf("  r_0 = %.4f mm   r_max = %.4f mm   ω = %.4g\n", r_0, r_max, ω)
@printf("  T = %.2f K   ψ = %.2f kPa   O2_gas = %.5g µg/mm³\n\n",
        T_const, ψ_const, O2_const)

# ── PART 1: what limits the timestep ─────────────────────────────────────────
#
# reaction_step.jl:70-81 forms rel_X = |S_X · Δt| / max(u_X, 1e-6) for the 8
# species at every node and takes the max. adapt_timestep halves Δt whenever
# that max exceeds 0.10. So the species with the largest |S|/u sets Δt for the
# whole system, and the Δt it permits is 0.10·u/|S|.
#
# This block reconstructs that quantity from saved states, using the same
# workspace-update sequence the timestepper uses (timestepper.jl:156-162).

probe_times = [0.0, 0.25, 1.0, 2.0, 5.0, 10.0, 21.0, 45.0]

println("Running the reference simulation (45 days)...")
t_run = @elapsed result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 45.0);
                                        n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                                        ic = ic, P_0 = P_0, ω = ω,
                                        output_times = probe_times)
@printf("  wall %.1f s   n_steps %d   mean Δt %.4g d   n_rejected %d\n\n",
        t_run, result.diagnostics["n_steps"],
        45.0 / result.diagnostics["n_steps"], result.diagnostics["n_rejected"])

grid    = result.grid
r_grid  = grid.r_grid
ws      = SAM.Workspace(N_GRID)
out_t   = [o.t for o in result.outputs]

"""Per-species max of |S|/u over nodes, at one saved state."""
function limiter_scan(state, ws, r_grid, bio, soil, T, ψ)
    n = length(state.C)
    SAM.update_temperature_cache!(ws.f_T, T, bio, soil)
    SAM.update_water_content!(ws.θ, ws.θ_a, ψ, state, soil)
    SAM.update_effective_diffusion!(ws, state, soil, bio, ws.f_T)

    worst     = fill(-Inf, 8)      # max |S|/u per species  [1/day]
    worst_i   = zeros(Int, 8)      # node where it occurs
    worst_u   = fill(NaN, 8)       # the pool value there
    worst_S   = fill(NaN, 8)       # the source term there

    for i in 1:n
        s = SAM.compute_source_terms(state.C[i], state.B[i], state.F_n[i],
                                     state.F_m[i], state.F_i[i], state.E[i],
                                     state.M[i], state.O[i],
                                     ws.θ[i], ws.θ_a[i], ψ, bio, soil, ws.f_T)
        u = (state.C[i], state.B[i], state.F_n[i], state.F_m[i],
             state.F_i[i], state.E[i], state.M[i], state.O[i])
        S = (s.S_C, s.S_B, s.S_Fn, s.S_Fm, s.S_Fi, s.S_E, s.S_M, s.S_O)
        for k in 1:8
            rate = abs(S[k]) / max(u[k], THRESHOLD)
            if rate > worst[k]
                worst[k]   = rate
                worst_i[k] = i
                worst_u[k] = u[k]
                worst_S[k] = S[k]
            end
        end
    end
    return (worst = worst, node = worst_i, u = worst_u, S = worst_S)
end

println("="^76)
println("PART 1 — what sets max_rel_change")
println("="^76)
println("""
For each species: the largest |S|/u over the 200 nodes, the node and pool value
where it occurs, and the timestep that species alone would permit under the
0.10 criterion:  Δt_allowed = 0.10 · u / |S|.  dt_max is 0.1 d, dt_min 1e-4 d.
""")

for (k_t, t) in enumerate(out_t)
    st  = result.outputs[k_t].state
    sc  = limiter_scan(st, ws, r_grid, bio, soil, T_func(t), ψ_func(t))
    ord = sortperm(sc.worst, rev = true)

    @printf("\n── t = %6.2f d ───────────────────────────────────────────────────\n", t)
    @printf("  %-5s %12s %6s %12s %12s %12s\n",
            "sp", "|S|/u [1/d]", "node", "u", "S", "Δt_allowed")
    for k in ord
        flag = (k == ord[1]) ? "  <== LIMITER" : ""
        @printf("  %-5s %12.4g %6d %12.4g %12.4g %12.4g%s\n",
                SPECIES[k], sc.worst[k], sc.node[k], sc.u[k], sc.S[k],
                0.10 / sc.worst[k], flag)
    end
    @printf("  max_rel_change at Δt=%.0e would be %.4g (criterion: halve if > 0.10)\n",
            DT_PROBE, sc.worst[ord[1]] * DT_PROBE)
    # what Δt the run could take if the top species were handled differently
    if length(ord) >= 2
        @printf("  without %s the limiter would be %s, allowing Δt = %.4g d (%.0fx)\n",
                SPECIES[ord[1]], SPECIES[ord[2]],
                min(0.1, 0.10 / sc.worst[ord[2]]),
                min(0.1, 0.10 / sc.worst[ord[2]]) / DT_PROBE)
    end
end

# radial profile of the limiter at the worst probe time
println("\n── radial profile of |S|/u for every species at t = 5 d ──")
let k5 = findmin(abs.(out_t .- 5.0))[2]
    st = result.outputs[k5].state
    SAM.update_temperature_cache!(ws.f_T, T_func(out_t[k5]), bio, soil)
    SAM.update_water_content!(ws.θ, ws.θ_a, ψ_func(out_t[k5]), st, soil)
    SAM.update_effective_diffusion!(ws, st, soil, bio, ws.f_T)
    @printf("  %6s %8s", "node", "r [mm]")
    for s in SPECIES; @printf(" %10s", s); end
    println()
    for i in 1:10:N_GRID
        s = SAM.compute_source_terms(st.C[i], st.B[i], st.F_n[i], st.F_m[i],
                                     st.F_i[i], st.E[i], st.M[i], st.O[i],
                                     ws.θ[i], ws.θ_a[i], ψ_func(out_t[k5]),
                                     bio, soil, ws.f_T)
        u = (st.C[i], st.B[i], st.F_n[i], st.F_m[i], st.F_i[i], st.E[i], st.M[i], st.O[i])
        S = (s.S_C, s.S_B, s.S_Fn, s.S_Fm, s.S_Fi, s.S_E, s.S_M, s.S_O)
        @printf("  %6d %8.3f", i, r_grid[i])
        for k in 1:8
            @printf(" %10.3g", abs(S[k]) / max(u[k], THRESHOLD))
        end
        println()
    end
end

# ── PART 2: where the wall time goes ─────────────────────────────────────────

println("\n" * "="^76)
println("PART 2 — where the wall time goes (sampling profiler, 3-day run)")
println("="^76)

# warm up so compilation is not profiled
run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 0.05);
              n_grid = N_GRID, r_0 = r_0, r_max = r_max,
              ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 0.05])

Profile.clear()
Profile.init(n = 10^7, delay = 0.0005)
@profile run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 3.0);
                       n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                       ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 3.0])

prof_path = joinpath(@__DIR__, "output", "profile_flat.txt")
mkpath(dirname(prof_path))
open(prof_path, "w") do io
    Profile.print(IOContext(io, :displaysize => (10_000, 200));
                  format = :flat, sortedby = :count, mincount = 5)
end
println("  full flat profile written to $(prof_path)")
println("  top entries (count, then function):\n")
let lines = readlines(prof_path)
    # flat profile is ascending by count; the tail is the hot code
    for l in Iterators.reverse(lines[max(1, end-45):end])
        println("  ", l)
    end
end

# ── PART 3: allocation per step ──────────────────────────────────────────────

println("\n" * "="^76)
println("PART 3 — allocation")
println("="^76)

let
    run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 0.05);
                  n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                  ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 0.05])
    a1 = @allocated run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 0.5);
                                  n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                                  ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 0.5])
    a2 = @allocated run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 1.0);
                                  n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                                  ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 1.0])
    # the difference removes the fixed setup cost, leaving per-step allocation
    r1 = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 0.5);
                       n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                       ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 0.5])
    r2 = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 1.0);
                       n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                       ic = ic, P_0 = P_0, ω = ω, output_times = [0.0, 1.0])
    ds = r2.diagnostics["n_steps"] - r1.diagnostics["n_steps"]
    @printf("  0.5 d run: %12d bytes over %d steps\n", a1, r1.diagnostics["n_steps"])
    @printf("  1.0 d run: %12d bytes over %d steps\n", a2, r2.diagnostics["n_steps"])
    @printf("  marginal:  %12.1f bytes per step  (docstrings claim zero)\n",
            ds > 0 ? (a2 - a1) / ds : NaN)
end

println("\n" * "="^76)
println("done")
println("="^76)
