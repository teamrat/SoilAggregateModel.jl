"""
converge_doc.jl — which solver is right about DOC at day 45?
=============================================================

bench_stiff.jl found the two solvers disagree by 9.6% on C at day 45, while
everything else agrees to under 0.4%. Six stiff runs spanning a 100x tolerance
range agree with each other to 1.5e-4, so the stiff answer is converged. This
asks whether the SPLIT solver is converged, by the only test that settles it:
refine its floor and see which way it moves.

The split controller halves dt only until max|S·dt|/u <= 0.10, so lowering
dt_min does not force tiny steps everywhere — it only lets dt fall where the
criterion demands it. At day 45 the criterion demands 4.9e-5, below the
production floor of 1e-4. So dt_min = 1e-5 and 1e-6 should cost about the same
and both should land near 4.9e-5. If they agree with each other AND with the
stiff answer, the production run is under-resolved on DOC. If they stay put,
the discrepancy is something else and the stiff path is the suspect.

Reports the radial DOC profile, not just its maximum, because the disagreement
is expected at the POM-adjacent nodes and a max-only comparison would hide
where it lives. Also reports r_agg, which is the output that actually feeds MWD.

Usage:
    julia --project=. paper/de_gryze/converge_doc.jl
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
const T_END   = 45.0

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

common = (n_grid = N_GRID, r_0 = r_0, r_max = r_max,
          ic = ic, P_0 = P_0, ω = tess.ω,
          output_times = [0.0, T_END])

println("="^80)
println("DOC CONVERGENCE — soil $(SOIL_ID), D = $(DIAM) mm, n = $(N_GRID), t = $(T_END) d")
println("="^80)

print("warm-up ... ")
run_aggregate_stiff(bio, soil, TF, PSI, O2F, (0.0, 0.05); common...)
run_aggregate(bio, soil, TF, PSI, O2F, (0.0, 0.05); common...)
println("done\n")

runs = Pair{String,Any}[]

print("stiff (FBDF, reltol 1e-6) ... ")
t = @elapsed r_stiff = run_aggregate_stiff(bio, soil, TF, PSI, O2F, (0.0, T_END);
                                           common..., reltol = 1e-6)
@printf("%.1f s, %d steps\n", t, r_stiff.diagnostics["n_accept"])
push!(runs, "stiff reltol 1e-6" => r_stiff)

# A second stiff run at a different tolerance AND a different error norm, so
# "the stiff answers agree" is not merely the same run twice.
print("stiff (FBDF, reltol 1e-7) ... ")
t = @elapsed r_stiff7 = run_aggregate_stiff(bio, soil, TF, PSI, O2F, (0.0, T_END);
                                            common..., reltol = 1e-7,
                                            abstol = 1e-11, abstol_scale = nothing)
@printf("%.1f s, %d steps\n", t, r_stiff7.diagnostics["n_accept"])
push!(runs, "stiff reltol 1e-7" => r_stiff7)

for dtm in (1e-4, 1e-5, 1e-6)
    @printf("split (dt_min = %.0e) ... ", dtm)
    local t, r
    t = @elapsed r = run_aggregate(bio, soil, TF, PSI, O2F, (0.0, T_END);
                                   common..., dt_min = dtm)
    @printf("%.1f s, %d steps\n", t, r.diagnostics["n_steps"])
    push!(runs, "split dt_min $(dtm)" => r)
end

# ── comparison ───────────────────────────────────────────────────────────────

ref = r_stiff                       # converged reference
a   = ref.outputs[end].state
grid = ref.grid

println("\n" * "-"^80)
println("Relative difference from the stiff reference at t = $(T_END) d")
println("  max_i |x_i - a_i| / max_i |a_i|")
println("-"^80)
@printf("  %-22s %10s %10s %10s %10s %10s %10s\n",
        "run", "C", "B", "F_i", "E", "M", "CO2")
for (label, r) in runs
    b = r.outputs[end].state
    rel(x, y) = maximum(abs.(x .- y)) / max(maximum(abs.(y)), eps())
    @printf("  %-22s %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g\n", label,
            rel(b.C, a.C), rel(b.B, a.B), rel(b.F_i, a.F_i),
            rel(b.E, a.E), rel(b.M, a.M),
            abs(b.CO2_cumulative - a.CO2_cumulative) / max(abs(a.CO2_cumulative), eps()))
end

println("""

  Reading: if the split rows move monotonically toward zero as dt_min falls,
  the production run (dt_min = 1e-4) is under-resolved and the stiff answer is
  right. If they sit still, the disagreement is not a timestep effect and the
  stiff path is the one to doubt.
""")

# ── where in the domain ──────────────────────────────────────────────────────

println("-"^80)
println("Radial DOC profile — the disagreement should live near the POM surface")
println("-"^80)
@printf("  %5s %8s", "node", "r [mm]")
for (label, _) in runs
    @printf(" %14s", first(label, 14))
end
println()
for i in [1, 2, 3, 5, 8, 12, 20, 40, 80, 140, 200]
    i > N_GRID && continue
    @printf("  %5d %8.3f", i, grid.r_grid[i])
    for (_, r) in runs
        @printf(" %14.6g", r.outputs[end].state.C[i])
    end
    println()
end

# ── the output that actually matters ─────────────────────────────────────────

println("\n" * "-"^80)
println("r_agg — the quantity that feeds MWD")
println("-"^80)
@printf("  %-22s %12s %12s %12s\n", "run", "r_agg [mm]", "P [µg-C]", "CO2 [µg-C]")
for (label, r) in runs
    ra = compute_r_agg(r.outputs[end], r.grid, r.params)
    st = r.outputs[end].state
    @printf("  %-22s %12.6f %12.4f %12.4f\n", label, ra, st.P, st.CO2_cumulative)
end

println("\n" * "="^80)
println("done")
println("="^80)
