"""
sensitivity.jl — which parameters move MWD and CO₂, and in what ratio?
======================================================================

Soil 3, five POM diameters, days 0-21 (the data range). Runs BEFORE any fitting,
for three reasons:

1. The model overshoots BOTH observables in the same direction — MWD 1389 µm at
   day 21 against 1014 measured, CO₂ ~3600 against 2139 µg-C/g. That is one
   cause or two. If a single parameter moves both in the observed ratio, the fit
   is nearly over before it starts. If MWD and CO₂ pull in opposite directions,
   that is a structural finding worth more than any fit.

2. ~57 data points cannot identify ~12 free parameters. This says which ones the
   data can actually see.

3. It is the empirical check on the by-inspection analysis in the session notes:
   with T and ψ held constant, the six activation energies, ν_B and ν_F enter
   only as fixed multipliers and are structurally unidentifiable HERE. They are
   swept anyway and marked :collinear. Their elasticities are NOT zero — a
   multiplier changes the answer even when it cannot be separated from what it
   multiplies. The 2026-07-29 run settled that: Ea_B came out at -0.356. What is
   checked is that each one's elasticity vector is parallel to the rate it
   multiplies; see the COLLINEAR section printed at the end.

WHAT IT REPORTS

For each parameter p, the log-log elasticity from a two-sided ±10 % perturbation

    E_y = ln(y₊/y₋) / ln(1.1/0.9)

for y = MWD at day 21 and y = cumulative CO₂ at day 21. Elasticity is
dimensionless, so parameters with wildly different units are comparable.

Then the decisive column: the log-change in p that would close the MWD gap on
its own, and what that same change does to the CO₂ gap. A parameter that closes
both is the one-cause candidate.

Usage:
    julia --project=. paper/de_gryze/fitting/sensitivity.jl
"""

include(joinpath(@__DIR__, "setup.jl"))

const PERTURB = 0.10

# create_initial_state emits an @info per run. 79 evaluations makes that noise,
# and it is identical every time. Warnings and errors still get through — a
# solver failure during the sweep must not be silenced.
Logging.disable_logging(Logging.Info)

# ── One model evaluation ─────────────────────────────────────────────────────
# Day-21 values only. This measures the endpoint, not the trajectory: a
# parameter that changes the SHAPE while leaving day 21 alone reads as zero
# here. Elasticities against the full loss are a separate sweep.

function evaluate(bio, soil)
    df = run_model(bio, soil)
    return (MWD = at_time(df, :MWD_fixed_weight_mm, 21.0) * 1000.0,   # µm, as measured
            CO2 = at_time(df, :CO2_total, 21.0) / SOIL_MASS_PER_L)     # µg-C/g-soil
end

# ── Measured targets, soil 3, day 21 ─────────────────────────────────────────

MWD_obs = obs_at(MWD_DAYS, MWD_OBS, 21.0)
CO2_obs = obs_at(CO2_DAYS, CO2_OBS, 21.0)

# ── Baseline ─────────────────────────────────────────────────────────────────

println("="^100)
println("SENSITIVITY — soil $(SOIL_ID), days 0-$(Int(T_END)), n=$(N_GRID), ±$(Int(PERTURB*100))%")
println("="^100)

t0 = @elapsed base = evaluate(BIO0, SOIL0)
@printf("  baseline: MWD %.1f µm   CO₂ %.1f µg-C/g     (%.1f s per evaluation)\n",
        base.MWD, base.CO2, t0)
@printf("  measured: MWD %.1f µm   CO₂ %.1f µg-C/g\n", MWD_obs, CO2_obs)

need_MWD = log(MWD_obs / base.MWD)
need_CO2 = log(CO2_obs / base.CO2)
@printf("  gap to close (log): MWD %+.3f  (×%.2f)    CO₂ %+.3f  (×%.2f)\n\n",
        need_MWD, exp(need_MWD), need_CO2, exp(need_CO2))

# ── The sweep ────────────────────────────────────────────────────────────────
#
# :collinear entries are multipliers on other rates when T and ψ are constant.
# They are redundant with those rates, not inert — see the note printed below.

PARAMS = [
    # (owner, field, role)
    (:bio,  :R_P_max,        :free),
    (:soil, :κ_b,            :free),
    (:soil, :w_E,            :free),
    (:bio,  :γ,              :free),
    (:bio,  :Y_B_max,        :free),
    (:bio,  :Y_F,            :free),
    (:bio,  :r_B_max,        :free),
    (:bio,  :r_F_max,        :free),
    (:bio,  :μ_B,            :free),
    (:bio,  :μ_F,            :free),
    (:bio,  :μ_E_max,        :free),
    (:bio,  :K_E,            :free),
    (:bio,  :K_B,            :free),
    (:bio,  :K_F,            :free),
    (:bio,  :D_Fn0,          :free),
    (:bio,  :D_Fm0,          :free),
    (:bio,  :α_i,            :free),
    (:bio,  :α_n,            :free),
    (:bio,  :ζ,              :free),
    (:bio,  :λ,              :free),
    (:bio,  :η_conv,         :free),
    (:bio,  :B_S,            :free),
    (:bio,  :F_S,            :free),
    (:bio,  :ε_F,            :free),
    (:bio,  :κ_s_ref,        :free),
    (:bio,  :κ_d_ref,        :free),
    (:soil, :k_ma,           :free),   # sets M_max = k_ma·f_clay_silt·ρ_b
    (:soil, :ω_E,            :free),
    (:soil, :ω_F,            :free),
    (:bio,  :Ea_B,           :collinear),
    (:bio,  :Ea_F,           :collinear),
    (:bio,  :Ea_EPS,         :collinear),
    (:bio,  :Ea_POM,         :collinear),
    (:bio,  :Ea_MAOC_sorb,   :collinear),
    (:bio,  :Ea_MAOC_desorb, :collinear),
    (:bio,  :ν_B,            :collinear),
    (:bio,  :ν_F,            :collinear),
    (:bio,  :L_B,            :collinear),
    (:bio,  :L_F,            :collinear),
]

const DLNP = log((1 + PERTURB) / (1 - PERTURB))

results = NamedTuple[]
t_sweep = @elapsed for (owner, field, role) in PARAMS
    base_val = owner === :bio ? getfield(BIO0, field) : getfield(SOIL0, field)
    if base_val == 0.0
        println("  skip $(field): base value is 0, multiplicative perturbation undefined")
        continue
    end
    hi, lo = try
        (owner === :bio ?
            evaluate(with_field(BIO0, field, base_val * (1 + PERTURB)), SOIL0) :
            evaluate(BIO0, with_field(SOIL0, field, base_val * (1 + PERTURB)))),
        (owner === :bio ?
            evaluate(with_field(BIO0, field, base_val * (1 - PERTURB)), SOIL0) :
            evaluate(BIO0, with_field(SOIL0, field, base_val * (1 - PERTURB))))
    catch e
        println("  FAILED $(field): $(sprint(showerror, e)[1:min(end,70)])")
        continue
    end
    E_MWD = log(hi.MWD / lo.MWD) / DLNP
    E_CO2 = log(hi.CO2 / lo.CO2) / DLNP
    push!(results, (field = field, role = role, base = base_val,
                    E_MWD = E_MWD, E_CO2 = E_CO2))
    @printf("  %-16s E_MWD %+8.4f   E_CO₂ %+8.4f\n", field, E_MWD, E_CO2)
end
@printf("\n  sweep: %.1f s for %d parameters\n\n", t_sweep, length(results))
# ── Collinear parameters: what the sweep actually showed ─────────────────────

println("="^100)
println("COLLINEAR PARAMETERS — redundant, NOT inert")
println("="^100)
println("""
  These were first labelled 'controls' expected to show ZERO elasticity. That
  was wrong, and the 2026-07-29 run proved it: Ea_B came out at -0.356.

  With T and psi constant these parameters enter only as fixed multipliers on
  other rates. That makes them COLLINEAR with those rates — not independently
  identifiable — but a multiplier still changes the answer, so the elasticity is
  nonzero. "Cannot be separated from r_B_max" and "has no effect" are different
  statements and only the first is true.

  The correct check is whether the elasticity VECTOR is parallel to that of the
  rate it multiplies, scaled by dln(f)/dln(Ea) = Ea*(1/T_ref - 1/T)/R. For Ea_B
  that predicts E_CO2 = 0.413 * (E[r_B_max] + E[mu_B]) = +0.0171 against +0.0163
  observed.

  Hold them fixed during fitting. If a fit moves one, it is trading against the
  rate it multiplies and the pair is unidentifiable.
""")
@printf("  %-16s %12s %12s\n", "parameter", "E_MWD", "E_CO2")
for r in results
    r.role === :collinear || continue
    @printf("  %-16s %+12.4f %+12.4f\n", r.field, r.E_MWD, r.E_CO2)
end
println()

# ── The decisive table ───────────────────────────────────────────────────────

println("="^100)
println("CAN ONE PARAMETER CLOSE BOTH GAPS?")
println("="^100)
println("""
  Δln p    = the log-change in p that closes the MWD gap by itself
  CO₂ after = the CO₂ gap remaining once that change is made (0 = closed)
  Ratio    = E_MWD / E_CO₂. The gap ratio to match is $(round(need_MWD/need_CO2, digits=3)).
""")
@printf("  %-16s %10s %10s %9s %10s %12s\n",
        "parameter", "E_MWD", "E_CO₂", "ratio", "Δln p", "CO₂ after")
free = filter(r -> r.role === :free && abs(r.E_MWD) > 1e-4, results)
sort!(free, by = r -> -abs(r.E_MWD))
for r in free
    Δlnp = need_MWD / r.E_MWD
    co2_after = need_CO2 - r.E_CO2 * Δlnp
    ratio = abs(r.E_CO2) > 1e-9 ? r.E_MWD / r.E_CO2 : NaN
    flag = abs(co2_after) < 0.1 * abs(need_CO2) ? "  <== closes BOTH" : ""
    @printf("  %-16s %+10.4f %+10.4f %9.3f %+10.3f %+12.3f%s\n",
            r.field, r.E_MWD, r.E_CO2, ratio, Δlnp, co2_after, flag)
end

println("""

  Reading it. A parameter with 'CO₂ after' near zero moves both observables in
  the ratio the misfit demands — one cause, one knob. If nothing is near zero,
  no single parameter explains the joint overshoot and at least two are needed;
  the spread of the ratio column then says which pair is most orthogonal, which
  is what makes them separately identifiable.

  Δln p is a LOCAL extrapolation from a ±$(Int(PERTURB*100))% perturbation. A
  |Δln p| much beyond ~0.7 (a factor of 2) is outside the range this linearisation
  can be trusted over — treat it as a direction, not a value.
""")

println("="^100)
println("done")
println("="^100)
