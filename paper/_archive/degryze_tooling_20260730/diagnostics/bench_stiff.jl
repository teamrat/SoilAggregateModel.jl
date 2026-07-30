"""
bench_stiff.jl — does the implicit step actually grow?
======================================================

smoke_stiff.jl compared the two solvers over 0.5 days and the stiff solver lost
by 18x. That window is the split solver's best case: it runs at Δt ≈ 1.3e-3
there, nowhere near its 1e-4 floor. The entire case for switching rests on days
5-45, where the split solver collapses to the floor and its cost per simulated
day grows.

This script answers three questions and nothing else.

  PART 1  Does FBDF's step grow with time, or is it flat?
          Run to horizons 1, 5, 15, 45 d and difference the step counts to get
          steps-per-day within each interval. Flat steps/day means the implicit
          method is accuracy-limited by something real and the switch does not
          pay. Falling steps/day is the whole argument.

  PART 2  Is the tolerance the binding constraint?
          abstol currently defaults to max(1e-10, 1e-8·|u0|). F_m starts at
          3.3e-8 (F_m_min/ω) and grows to ~3.5e-4, so it is being controlled to
          1e-10 absolute — four orders tighter than its own scale. Basing
          tolerance on the INITIAL value is wrong for any pool that starts near
          zero and grows. Sweep it.

  PART 3  Does the answer survive the loosening?
          Any tolerance that changes day-45 pools by more than the 1e-3 the two
          solvers already disagree by is too loose, regardless of speed.

Usage:
    julia --project=. paper/de_gryze/bench_stiff.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using SoilAggregateModel
using Printf
const SAM = SoilAggregateModel

include(joinpath(@__DIR__, "..", "degryze_soils.jl"))

const SOIL_ID = 3
const DIAM    = 1.25
const N_GRID  = 200

soil = degryze_soil(SOIL_ID; k_L = 1000, D_B_rel = 0.00001)
ic   = degryze_ic(SOIL_ID, soil)

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
const TF  = t -> T_const
const PSI = t -> ψ_const
const O2F = t -> O2_const

tess  = SAM.domain_tessellation(ρ_POM = 200.0, I_input = 4.43e-3, ρ_b = soil.ρ_b)
r_0   = DIAM / 2.0
r_max = DIAM * tess.f_domain / 2.0
P_0   = (4.0/3.0) * π * r_0^3 * 200.0

stiff(t_end; kw...) = run_aggregate_stiff(
    bio, soil, TF, PSI, O2F, (0.0, t_end);
    n_grid = N_GRID, r_0 = r_0, r_max = r_max,
    ic = ic, P_0 = P_0, ω = tess.ω,
    output_times = [0.0, t_end], kw...)

split(t_end; kw...) = run_aggregate(
    bio, soil, TF, PSI, O2F, (0.0, t_end);
    n_grid = N_GRID, r_0 = r_0, r_max = r_max,
    ic = ic, P_0 = P_0, ω = tess.ω,
    output_times = [0.0, t_end], kw...)

println("="^80)
println("STIFF SOLVER BENCHMARK — soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID)")
println("="^80)

print("warming up ... ")
stiff(0.05); split(0.05)
println("done\n")

# ── PART 1 ───────────────────────────────────────────────────────────────────

println("PART 1 — does the implicit step grow?")
println("-"^80)
horizons = [1.0, 5.0, 15.0, 45.0]
rows = NamedTuple[]
for T in horizons
    t = @elapsed r = stiff(T)
    d = r.diagnostics
    push!(rows, (T = T, wall = t, acc = d["n_accept"], rej = d["n_reject"],
                 nf = d["n_f"], njac = d["n_jac"], rc = d["retcode"],
                 res = r))
    @printf("  t_end %5.1f d   wall %7.2f s   accepted %7d   rejected %5d   f %8d   jac %4d   %s\n",
            T, t, d["n_accept"], d["n_reject"], d["n_f"], d["n_jac"], d["retcode"])
end

println("\n  steps per day, by interval — the number that decides this:")
@printf("  %-14s %12s %12s %14s\n", "interval [d]", "steps", "steps/day", "mean Δt [d]")
prev_T = 0.0; prev_acc = 0
for r in rows
    ds = r.acc - prev_acc
    dT = r.T - prev_T
    @printf("  %-14s %12d %12.1f %14.3g\n",
            "$(prev_T)–$(r.T)", ds, ds/dT, dT/max(ds,1))
    global prev_T = r.T; global prev_acc = r.acc
end
println("""
  For reference the split solver runs the same 45 days in 391,773 steps —
  8,706 steps/day, essentially flat at its floor, with the REQUIRED step still
  falling (6.7e-5 d at day 21, 4.9e-5 d at day 45).
""")

t_split45 = @elapsed res_split45 = split(45.0)
@printf("  split solver, 45 d: wall %6.2f s   steps %d\n", t_split45, res_split45.diagnostics["n_steps"])
@printf("  stiff  solver, 45 d: wall %6.2f s   steps %d\n", rows[end].wall, rows[end].acc)
@printf("  SPEEDUP: %.2fx\n\n", t_split45 / rows[end].wall)

# ── PART 2 ───────────────────────────────────────────────────────────────────

println("PART 2 — tolerance sensitivity at 45 d")
println("-"^80)
@printf("  %-28s %10s %10s %10s %10s\n", "setting", "wall [s]", "accepted", "rejected", "retcode")

configs = [
    ("reltol 1e-6, abstol(u0)",  (reltol = 1e-6,)),
    ("reltol 1e-5, abstol(u0)",  (reltol = 1e-5,)),
    ("reltol 1e-4, abstol(u0)",  (reltol = 1e-4,)),
    ("reltol 1e-6, abstol 1e-8", (reltol = 1e-6, abstol = 1e-8, abstol_scale = nothing)),
    ("reltol 1e-5, abstol 1e-8", (reltol = 1e-5, abstol = 1e-8, abstol_scale = nothing)),
    ("reltol 1e-4, abstol 1e-7", (reltol = 1e-4, abstol = 1e-7, abstol_scale = nothing)),
]

results = Any[]
for (label, kw) in configs
    local t, r
    try
        t = @elapsed r = stiff(45.0; kw...)
        @printf("  %-28s %10.2f %10d %10d %10s\n",
                label, t, r.diagnostics["n_accept"], r.diagnostics["n_reject"],
                r.diagnostics["retcode"])
        push!(results, (label = label, wall = t, res = r))
    catch e
        @printf("  %-28s   FAILED: %s\n", label, sprint(showerror, e)[1:min(end,60)])
        push!(results, (label = label, wall = NaN, res = nothing))
    end
end

# ── PART 3 ───────────────────────────────────────────────────────────────────

println("\nPART 3 — does loosening change the answer?")
println("-"^80)
println("""
  Reference is the tightest stiff run. The two solvers already disagree by
  ~1.3e-3 on CO2 from coefficient lagging alone, so a tolerance that moves
  day-45 pools by more than that is buying speed with accuracy.
""")

ref = results[1].res
if ref !== nothing
    a = ref.outputs[end].state
    @printf("  %-28s %11s %11s %11s %11s\n", "setting", "C", "B", "CO2", "P")
    for rr in results
        rr.res === nothing && continue
        b = rr.res.outputs[end].state
        rel(x, y) = maximum(abs.(x .- y)) / max(maximum(abs.(x)), eps())
        @printf("  %-28s %11.3g %11.3g %11.3g %11.3g\n", rr.label,
                rel(a.C, b.C), rel(a.B, b.B),
                abs(a.CO2_cumulative - b.CO2_cumulative) / max(abs(a.CO2_cumulative), eps()),
                abs(a.P - b.P) / max(abs(a.P), eps()))
    end

    println("\n  and against the split solver at 45 d (its production settings):")
    c = res_split45.outputs[end].state
    @printf("    C %.3g   B %.3g   F_i %.3g   E %.3g   M %.3g   CO2 %.3g   P %.3g\n",
            maximum(abs.(a.C .- c.C))/maximum(abs.(c.C)),
            maximum(abs.(a.B .- c.B))/maximum(abs.(c.B)),
            maximum(abs.(a.F_i .- c.F_i))/maximum(abs.(c.F_i)),
            maximum(abs.(a.E .- c.E))/maximum(abs.(c.E)),
            maximum(abs.(a.M .- c.M))/maximum(abs.(c.M)),
            abs(a.CO2_cumulative - c.CO2_cumulative)/abs(c.CO2_cumulative),
            abs(a.P - c.P)/max(abs(c.P), eps()))
    @printf("    carbon balance: stiff %.3g   split %.3g\n",
            ref.outputs[end].mass_balance_error,
            res_split45.outputs[end].mass_balance_error)
end

println("\n" * "="^80)
println("done")
println("="^80)
