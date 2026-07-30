"""
verify_solver_agreement.jl — do the two solvers describe the same model?
========================================================================

`run_aggregate_stiff` is the production solver. `run_aggregate` (Strang splitting
+ Crank–Nicolson) is the reference implementation, kept because it is an
independent discretisation of the same equations and because it tracks respired
carbon separately, so it can close the carbon balance as a measurement rather
than by construction. REFERENCE.md §17a and §20a rest on the two agreeing. This
script is the evidence for that claim, and there is nothing else.

It lives here rather than in `test/` because it needs a configuration to run and
`degryze_config.jl` is the one we have; the package must not depend on it. It is
not part of the forward run — it verifies the solver the forward run uses.

WHAT AGREEMENT MEANS HERE

`mol.jl` documents three deliberate differences, so exact agreement is not the
expectation and would in fact be suspicious:

  1. the split solver lags θ and the effective diffusivities by one step
  2. it clips negative pools to zero and redirects the clipped carbon into CO₂
  3. it accumulates CO₂ as a running sum, where the stiff path recovers it by
     difference from the conserved total

Disagreement at 1e-3 is (1). Disagreement at 1e-1 is a bug. The quantity to
watch is DOC at the POM-adjacent nodes, where the lag bites hardest, and `r_agg`,
because that is the output that feeds MWD.

PART 2 IS THE PART THAT SETTLES IT

If they disagree, agreeing on who is wrong needs one more fact: whether the split
answer is converged. Its controller halves Δt only until max|S·Δt|/u ≤ 0.10, so
lowering `dt_min` does not force small steps everywhere — it only lets Δt fall
where the criterion asks. If refining the floor moves the split answer toward the
stiff one, the split run was under-resolved. If it stays put, the stiff path is
the suspect.

Replaces `diagnostics/converge_doc.jl`, archived 2026-07-30: that script forked
the configuration and its stated premise ("9.6 % on C at day 45") was measured at
`k_L = 1000` against the current 25000 — the one parameter that moves DOC.

Usage:
    julia --project=. paper/de_gryze/verify_solver_agreement.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using SoilAggregateModel
using Printf

# THE configuration. Nothing is redefined here.
include(joinpath(@__DIR__, "degryze_config.jl"))

# Modal class, geometry formed exactly as run_diameter_sweep forms it.
const DIAM  = DIAM_ALL[argmin(abs.(DIAM_ALL .- POM_MEAN))]
const R_0   = DIAM / 2.0
const R_MAX = DIAM * DOMAIN_FACTOR / 2.0
const P0    = (4.0 / 3.0) * π * R_0^3 * ρ_POM
# Horizon. Defaults to the De Gryze window; pass one on the command line to test
# a longer one. 45 d is the interesting case: the split solver's required step is
# still falling there (4.9e-5 d against its 1e-4 floor), so it is where the two
# solvers were first seen to disagree.
#
#     julia --project=. paper/de_gryze/verify_solver_agreement.jl 45
#
# At 45 d the refined split run takes several minutes — its cost per simulated
# day grows while the stiff solver's falls.
const T_END = isempty(ARGS) ? T_MAX : parse(Float64, ARGS[1])
const TIMES = collect(0.0:DT_OUTPUT:T_END)

common = (n_grid = N_GRID, r_0 = R_0, r_max = R_MAX,
          ic = IC, P_0 = P0, ω = OMEGA, output_times = TIMES)

run_stiff() = run_aggregate_stiff(BIO, SOIL, T_FUNC, PSI_FUNC, O2_FUNC,
                                  (0.0, T_END); common...)
run_split(; dt_min = 1e-4) = run_aggregate(BIO, SOIL, T_FUNC, PSI_FUNC, O2_FUNC,
                                  (0.0, T_END); dt_min = dt_min, dt_max = 0.1, common...)

# Relative difference on a profile: scaled by the reference's own magnitude, so a
# pool that is uniformly tiny does not read as a large disagreement.
reldiff(a::AbstractVector, b::AbstractVector) =
    maximum(abs.(a .- b)) / max(maximum(abs.(b)), eps())
reldiff(a::Real, b::Real) = abs(a - b) / max(abs(b), eps())

const FIELDS = (:C, :B, :F_n, :F_m, :F_i, :E, :M, :O)

# Per-field gaps, so the verdict can test DIRECTION and not merely magnitude.
gaps(ra, rb) = Dict(f => reldiff(getfield(ra.outputs[end].state, f),
                                 getfield(rb.outputs[end].state, f)) for f in FIELDS)

function compare(label, ra, rb)
    a, b = ra.outputs[end].state, rb.outputs[end].state
    println("  $(label)")
    @printf("    %-6s", "field")
    for f in FIELDS; @printf(" %9s", f); end
    @printf(" %9s %9s %9s\n", "P", "CO2", "r_agg")
    ragg(r) = compute_r_agg(r.outputs[end], r.grid, r.params)
    @printf("    %-6s", "rel Δ")
    for f in FIELDS
        @printf(" %9.2e", reldiff(getfield(a, f), getfield(b, f)))
    end
    @printf(" %9.2e %9.2e %9.2e\n",
            reldiff(a.P, b.P),
            reldiff(a.CO2_cumulative, b.CO2_cumulative),
            reldiff(ragg(ra), ragg(rb)))
    return maximum(reldiff(getfield(a, f), getfield(b, f)) for f in FIELDS)
end

println("="^78)
println("SOLVER AGREEMENT — soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID), $(T_END) d")
println("  config: degryze_config.jl        reference: run_aggregate (split)")
println("="^78)

# Warm-up first. Without it the first timing is dominated by compilation: an
# earlier version of this script reported the stiff run at 11.7 s against the
# 1.16 s the profiler measures for the same 22 days, and printed a 1.6x speedup
# where the documented figure is 24x.
print("\nwarming up ... ")
run_aggregate_stiff(BIO, SOIL, T_FUNC, PSI_FUNC, O2_FUNC, (0.0, 0.05);
                    n_grid = N_GRID, r_0 = R_0, r_max = R_MAX,
                    ic = IC, P_0 = P0, ω = OMEGA, output_times = [0.0, 0.05])
run_aggregate(BIO, SOIL, T_FUNC, PSI_FUNC, O2_FUNC, (0.0, 0.05);
              n_grid = N_GRID, r_0 = R_0, r_max = R_MAX, dt_min = 1e-4, dt_max = 0.1,
              ic = IC, P_0 = P0, ω = OMEGA, output_times = [0.0, 0.05])
println("done")

print("\nstiff ... ");   t_stiff = @elapsed r_stiff = run_stiff();  @printf("%.1f s\n", t_stiff)
print("split ... ");     t_split = @elapsed r_split = run_split();  @printf("%.1f s\n", t_split)
@printf("  speedup %.1fx     split steps %d\n\n",
        t_split / t_stiff, r_split.diagnostics["n_steps"])

println("PART 1 — agreement at day $(T_END), split (production floor) vs stiff")
println("-"^78)
worst = compare("dt_min 1e-4 vs stiff", r_split, r_stiff)

# The split solver's own carbon closure — the reason it is kept. The stiff path
# reports NaN here by design, because it recovers CO₂ by difference and the
# balance is then identically zero.
@printf("\n  split-solver carbon balance error: %.3e   (stiff reports %s by design)\n",
        r_split.outputs[end].mass_balance_error,
        string(r_stiff.outputs[end].mass_balance_error))

println("\nPART 2 — is the split answer converged? refine its floor")
println("-"^78)
print("split, dt_min 1e-5 ... "); t_ref = @elapsed r_ref = run_split(dt_min = 1e-5)
@printf("%.1f s   steps %d\n", t_ref, r_ref.diagnostics["n_steps"])
compare("dt_min 1e-5 vs stiff", r_ref, r_stiff)
moved = compare("dt_min 1e-5 vs dt_min 1e-4", r_ref, r_split)

# ── Verdict ──────────────────────────────────────────────────────────────────
# Movement alone says nothing: the split answer could move AWAY from stiff on
# refinement, which is the opposite conclusion. What is tested is whether every
# field's gap SHRANK, and by how much. A gap falling in proportion to the floor
# is first order in Δt, which is the signature of the documented one-step lag in
# θ and D rather than of a defect.

g_prod, g_ref = gaps(r_split, r_stiff), gaps(r_ref, r_stiff)
shrank   = count(f -> g_ref[f] < g_prod[f], FIELDS)
ratios   = [g_prod[f] / max(g_ref[f], eps()) for f in FIELDS]
med_fac  = sort(ratios)[cld(length(ratios), 2)]

# The observables the paper reports, kept apart from the interior pools: they
# agree far better, and a single worst-field number hides that.
obs(ra, rb) = (P = reldiff(ra.outputs[end].state.P, rb.outputs[end].state.P),
               CO2 = reldiff(ra.outputs[end].state.CO2_cumulative,
                             rb.outputs[end].state.CO2_cumulative),
               r_agg = reldiff(compute_r_agg(ra.outputs[end], ra.grid, ra.params),
                               compute_r_agg(rb.outputs[end], rb.grid, rb.params)))
o_prod = obs(r_split, r_stiff)

println("\n" * "="^78)
@printf("VERDICT\n\n  worst interior field, split at the production floor: %.2e\n", worst)
@printf("  reported observables at that floor: P %.1e, CO2 %.1e, r_agg %.1e\n",
        o_prod.P, o_prod.CO2, o_prod.r_agg)
@printf("  refining dt_min 10x: %d of %d field gaps shrank, median factor %.1fx\n\n",
        shrank, length(FIELDS), med_fac)

if worst < 1e-2
    println("  The two agree to better than 1e-2 at the production floor, which is")
    println("  consistent with the documented one-step lag. Supported.")
elseif shrank == length(FIELDS) && med_fac > 3.0
    println("  Every gap shrank and they fall roughly in proportion to the floor,")
    println("  i.e. first order in dt. That is the documented lag in theta and D,")
    println("  not a defect: the split solver CONVERGES TO the stiff answer, and its")
    println("  production floor is simply not converged. The stiff answer is the")
    println("  trusted one, and any comparison made at dt_min = 1e-4 understates")
    println("  the agreement.")
elseif shrank < length(FIELDS)
    println("  $(length(FIELDS) - shrank) field(s) moved AWAY from the stiff answer on")
    println("  refinement. Convergence is not toward stiff. Do not quote either")
    println("  solver until this is explained.")
else
    println("  Gaps shrank but only by $(round(med_fac, sigdigits=2))x for a 10x floor")
    println("  refinement, so the residual is not purely a step-size effect.")
    println("  Something else contributes. Explain before quoting.")
end
# ── PART 3, only when something moved the wrong way ──────────────────────────
#
# A field that diverges on refinement is either a real structural disagreement or
# a tolerance artefact, and the two are told apart by WHERE the gap lives. If the
# worst node is one where the field has decayed far below its own initial value,
# suspect the stiff solver's absolute tolerance: `abstol_scale` sets it from |u0|
# (abstol = max(abstol, abstol_scale·|u0|)), so a pool that starts at 1e-2 and
# falls to 1e-9 is being controlled to ~10 % of itself by the end of a long run.
# `bench_stiff.jl` flagged the mirror image of this for F_m, which starts tiny and
# grows.

diverged = [f for f in FIELDS if g_ref[f] >= g_prod[f]]
if !isempty(diverged)
    println("\nPART 3 — where the diverging gap lives")
    println("-"^78)
    sr, ss = r_ref.outputs[end].state, r_stiff.outputs[end].state
    @printf("  %-5s %6s %10s %12s %12s %12s %10s\n",
            "field", "node", "r [mm]", "split", "stiff", "profile max", "value/u0")
    for f in diverged
        x, y = getfield(sr, f), getfield(ss, f)
        i = argmax(abs.(x .- y))
        u0 = getfield(r_stiff.outputs[1].state, f)[i]
        @printf("  %-5s %6d %10.4f %12.4e %12.4e %12.4e %10.2e\n",
                f, i, r_stiff.grid.r_grid[i], x[i], y[i],
                maximum(abs, y), u0 == 0 ? NaN : abs(y[i]) / abs(u0))
    end
    println("""
  A `value/u0` far below 1 with the worst node in a depleted region means the
  stiff solver's absolute tolerance, not a structural disagreement. Test it by
  re-running with a tolerance floor tied to the field's own scale rather than its
  initial value:

      run_aggregate_stiff(...; reltol = 1e-8, abstol = 1e-14, abstol_scale = nothing)

  If the gap closes, the tolerance is the cause and `abstol_scale` needs to stop
  being anchored to u0. If it does not, the disagreement is structural and
  neither solver should be quoted at this horizon.""")
end

println()
println("  This is a CONSISTENCY result between two discretisations, not a")
println("  verification against an analytic solution. Record it in REFERENCE §17a")
println("  with today's date and the parameters — the claim is only as current as")
println("  the last run of this script.")
println("="^78)
