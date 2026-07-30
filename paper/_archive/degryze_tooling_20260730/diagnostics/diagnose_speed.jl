"""
diagnose_speed.jl — where does the PRODUCTION run spend its time?
=================================================================

Read-only. It changes nothing in `src/`; it calls what the solver calls and
reports what comes back.

WHY THIS WAS REWRITTEN (2026-07-30)

The previous version profiled `run_aggregate` — the split solver — against a
hand-copied parameter set, under a header that said "These MUST match
run_degryze.jl exactly — a diagnostic of a different parameter set answers a
different question." It had stopped matching: `κ_d_ref = 0.001`, `R_P_max = 2.5`,
`k_L = 1000`, and `101000.0` for one atmosphere. Its output is the only evidence
behind the audit's cost arm (§B), which therefore rests on a profile of the
solver that is being archived, taken at parameters no longer in use.

This version includes `degryze_config.jl`, so there is one configuration, and
profiles `run_aggregate_stiff`, which is what actually runs.

WHAT IT ANSWERS

  PART 1  What does one aggregate cost, and how many steps does it take?

  PART 2  Where does the time go? A flat profile, plus a bucketed summary
          aimed at the three specific hoists the cost arm proposes:

            pow / ^        — θ_s^(2/3) and θ_s^2 are soil constants recomputed
                             per node per step in D_eff_oxygen and the
                             Millington-Quirk tortuosity. A non-integer power is
                             a log and an exp, so this bucket is the test of
                             whether hoisting them is worth anything.
            reactions.jl   — compute_source_terms forms ~10 rate×Arrhenius
                             products per NODE from a cache that is per STEP.
            linear algebra — Jacobian construction and the KLU solve. If this
                             dominates, arithmetic hoisting in the RHS cannot
                             pay and the cost arm should be closed on the
                             measurement rather than acted on.

  PART 3  Does the run allocate per step?

Usage:
    julia --project=. paper/de_gryze/diagnostics/diagnose_speed.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using SoilAggregateModel
using Printf
using Profile

const SAM = SoilAggregateModel

# THE configuration — the same file run_degryze.jl and fitting/setup.jl include.
# Nothing is redefined here. If a parameter needs changing it changes there.
include(joinpath(@__DIR__, "..", "degryze_config.jl"))

# Geometry for the modal class, formed exactly as run_diameter_sweep forms it.
# The class nearest the distribution mean. POM_MEAN itself is not a bin
# midpoint, so profiling it would profile an aggregate the sweep never runs.
const DIAM  = DIAM_ALL[argmin(abs.(DIAM_ALL .- POM_MEAN))]
const R_0   = DIAM / 2.0
const R_MAX = DIAM * DOMAIN_FACTOR / 2.0
const P0    = (4.0 / 3.0) * π * R_0^3 * ρ_POM

stiff(t_end; output_times = [0.0, t_end]) = run_aggregate_stiff(
    BIO, SOIL, T_FUNC, PSI_FUNC, O2_FUNC, (0.0, t_end);
    n_grid = N_GRID, r_0 = R_0, r_max = R_MAX,
    ic = IC, P_0 = P0, ω = OMEGA, output_times = output_times)

println("="^78)
println("PRODUCTION SPEED DIAGNOSTIC — soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID)")
println("  solver: run_aggregate_stiff        config: degryze_config.jl")
println("="^78)

print("\nwarming up ... ")
stiff(0.05)
println("done")

# ── PART 1 ───────────────────────────────────────────────────────────────────

println("\nPART 1 — cost of one aggregate")
println("-"^78)
t_full = @elapsed res = stiff(T_MAX; output_times = collect(0.0:DT_OUTPUT:T_MAX))
d = res.diagnostics
@printf("  %.1f days      wall %7.2f s\n", T_MAX, t_full)
@printf("  accepted %7d   rejected %6d   f evals %8d   jacobians %5d   %s\n",
        d["n_accept"], d["n_reject"], d["n_f"], d["n_jac"], d["retcode"])
@printf("  mean Δt %.4g d      wall per accepted step %6.1f µs\n",
        T_MAX / max(d["n_accept"], 1), 1e6 * t_full / max(d["n_accept"], 1))
@printf("  f evals per accepted step: %.1f     (each one is %d nodes)\n",
        d["n_f"] / max(d["n_accept"], 1), N_GRID)
@printf("  full sweep of %d classes would be about %.0f s\n",
        N_POM_BINS, N_POM_BINS * t_full)

# ── PART 2 ───────────────────────────────────────────────────────────────────

println("\nPART 2 — where the time goes")
println("-"^78)

# A shorter horizon: the profile only needs enough samples, and the mix of work
# per step does not depend on how many steps there are.
const T_PROF = min(5.0, T_MAX)

Profile.clear()
Profile.init(n = 10^7, delay = 0.0005)
@profile stiff(T_PROF)

flat = sprint() do io
    Profile.print(IOContext(io, :displaysize => (10_000, 240));
                  format = :flat, sortedby = :count, mincount = 1)
end

outdir = joinpath(@__DIR__, "..", "output")
isdir(outdir) || mkpath(outdir)
flatfile = joinpath(outdir, "profile_stiff_flat.txt")
open(flatfile, "w") do io
    println(io, "run_aggregate_stiff, soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID), $(T_PROF) d")
    println(io, "config: degryze_config.jl")
    print(io, flat)
end
println("  flat profile written to output/profile_stiff_flat.txt")

# ── Bucketed self-time ───────────────────────────────────────────────────────
# The `Overhead` column is self-time: samples where that line was executing
# rather than waiting on a callee. Summing it by bucket is what says whether a
# hoist can pay. Total-time columns cannot be summed — they nest.

const BUCKETS = [
    ("pow / ^ (θ_s^(2/3), θ_a^(10/3))", ["math.jl"],                    ["pow", "^"]),
    ("exp",                             ["special/exp.jl"],             String[]),
    ("log",                             ["special/log.jl"],             String[]),
    ("reactions.jl  (source terms)",    ["solver/reactions.jl"],        String[]),
    ("biology/*.jl  (rate functions)",  ["biology/"],                   String[]),
    ("mol.jl  (RHS assembly)",          ["solver/mol.jl"],              String[]),
    ("effective_diffusion.jl",          ["effective_diffusion.jl"],     String[]),
    ("water_retention.jl",              ["water_retention.jl"],         String[]),
    ("linear algebra / KLU",            ["KLU", "LinearSolve", "SparseArrays", "SparseMatrix"], String[]),
    ("Jacobian / integrator",           ["OrdinaryDiffEq", "FiniteDiff", "SparseDiffTools",
                                         "DiffEqBase", "SciMLBase"],    String[]),
]

# Parsed inside a function: a top-level `for` opens a soft scope, so
# accumulating into globals from one silently creates locals instead.
function parse_self_time(flat::AbstractString)
    total = 0
    tally = Dict{String,Int}(b[1] => 0 for b in BUCKETS)
    other = 0
    rows  = NamedTuple[]
    for line in split(flat, '\n')
        m = match(r"^\s*(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(.*)$", line)
        m === nothing && continue
        self = parse(Int, m.captures[2])
        self == 0 && continue
        file, fun = m.captures[3], m.captures[5]
        total += self
        push!(rows, (self = self, file = file, line = m.captures[4], fun = fun))
        hit = false
        for (name, files, funs) in BUCKETS
            if any(f -> occursin(f, file), files) &&
               (isempty(funs) || any(f -> occursin(f, fun), funs))
                tally[name] += self
                hit = true
                break
            end
        end
        hit || (other += self)
    end
    sort!(rows, by = r -> -r.self)
    return (total = total, tally = tally, other = other, rows = rows)
end

const PT = parse_self_time(flat)
const self_total = PT.total

println()
if self_total == 0
    println("  no self-time samples — raise the horizon or lower the delay")
else
    @printf("  self-time, %d samples total (%.2f s of wall at 0.5 ms per sample)\n",
            self_total, self_total * 0.0005)
    @printf("  %-36s %8s %8s\n", "bucket", "samples", "% self")
    for (name, _, _) in BUCKETS
        v = PT.tally[name]
        v == 0 && continue
        @printf("  %-36s %8d %7.1f%%\n", name, v, 100 * v / self_total)
    end
    @printf("  %-36s %8d %7.1f%%\n", "everything else", PT.other,
            100 * PT.other / self_total)

    println("\n  top self-time lines:")
    @printf("  %8s  %-44s %6s  %s\n", "% self", "file", "line", "function")
    for r in first(PT.rows, 18)
        @printf("  %7.1f%%  %-44s %6s  %s\n",
                100 * r.self / self_total, first(r.file, 44), r.line, first(r.fun, 46))
    end

    println("""
  Percentages are of samples carrying self-time, which is what a hoist can
  remove. Total-time columns in the flat file cannot be summed — they nest.
  A transcendental called inside a biology function shows up under exp/log/pow,
  not under that function, so those buckets read low by design.

  Reading it. The three hoists in the audit's cost arm are worth doing only in
  proportion to their bucket. A `pow / ^` bucket in the low single digits means
  hoisting θ_s^(2/3) buys that much at most, and less in practice because
  θ_a^(10/3) is genuinely per-node and stays. If linear algebra plus Jacobian is
  most of the time, no arithmetic hoisting in the RHS can matter and the cost arm
  closes on this measurement.""")
end

println("\nPART 3 — allocation")
println("-"^78)
stiff(0.05)                                  # compile
a1 = @allocated stiff(0.5)
a2 = @allocated stiff(1.0)
@printf("  0.5 d: %10.3f MB\n", a1 / 2^20)
@printf("  1.0 d: %10.3f MB\n", a2 / 2^20)
@printf("  ratio: %.2f   (≈2 means allocation scales with step count, i.e. per step)\n",
        a2 / max(a1, 1))
println("""
  The stiff solver allocates by design — it deep-copies a state per output
  record and OrdinaryDiffEq keeps its own history. What matters is whether the
  RHS allocates per call, which shows up as a ratio near 2 with a large
  absolute size.""")

println("\n" * "="^78)
println("done — paste the PART 2 table back")
println("="^78)
