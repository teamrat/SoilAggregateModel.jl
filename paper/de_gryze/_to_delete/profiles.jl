"""
profiles.jl — every state variable against radius, at many times
================================================================

Diagnostic. The population outputs say the aggregate radius jumps on day 3 and
is then flat for 42 days; this shows the fields that produce it.

Nine panels: the eight pools, plus the binder `F_i + w_E·E` against the
threshold it has to cross.

On the binder panel, three lines:

  * `G_c` as the model currently computes it — τ_w·d_32/κ_b, CONSTANT in r.
  * `G_c(r)` CANDIDATE — the same threshold scaled by r/δ_s, i.e. a disruptive
    stress that grows with aggregate size. **Not implemented in the model.**
    Drawn only to show where it would cut the binder profile. It is the
    MATLAB precursor's commented-out `strength ./ x` written as a threshold, and
    it is the deferred sieve-confinement item in BACKLOG.
  * `r_agg` at each time, as a marker on the binder curve.

Why the candidate is worth looking at: the binder rises ~6-fold near the POM
between days 5 and 30 but is FLAT OR FALLING at the radius where the constant
threshold actually cuts. A threshold that rises with r moves the crossing inward,
into the region that is still accumulating, which is the only way this model can
produce a radius that grows for 21 days.

Writes profiles.csv as well — the CSV is the deliverable, the figure is
diagnostic.

Usage:
    julia --project=. paper/de_gryze/fitting/profiles.jl
"""

include(joinpath(@__DIR__, "loss.jl"))

using Plots
gr()

Logging.disable_logging(Logging.Info)

const DIAM      = 1.25                       # modal POM class
const R_P_SCAN  = 0.9                        # scan.jl minimiser; BIO0 uses 2.5
const SNAP      = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 14.0, 21.0]

bio  = with_field(BIO0, :R_P_max, R_P_SCAN)
soil = SOIL0

r_0   = DIAM / 2.0
r_max = DIAM * TESS.f_domain / 2.0
P_0   = (4.0/3.0) * π * r_0^3 * ρ_POM

res = run_aggregate_stiff(bio, soil, TF, PSI, O2F, (0.0, T_END);
                          n_grid = N_GRID, r_0 = r_0, r_max = r_max,
                          ic = IC, P_0 = P_0, ω = TESS.ω,
                          output_times = SNAP)

r  = res.grid.r_grid
hy = wet_sieving_stress()
G_c = critical_binding(soil)

# CANDIDATE ONLY — not in the model. Threshold scaled by r/δ_s so it equals the
# current constant value at r = δ_s and rises linearly outside it.
G_c_r = G_c .* (r ./ hy.δ_s)

@printf("τ_w = %.5f Pa   δ_s = %.4f mm   d_32 = %.5f mm   κ_b = %.5f\n",
        hy.τ_w, hy.δ_s, soil.d_32, soil.κ_b)
@printf("G_c = %.5f µg/mm³ (constant)   r_0 = %.3f   r_max = %.3f mm\n\n",
        G_c, r_0, r_max)

# ── How far can DOC actually get? ────────────────────────────────────────────
#
# D_eff_DOC = D_DOC_w · τ · θ/(θ + ρ_b·k_d). The last factor is sorption
# retardation and it is the dominant one here.
let st0 = res.outputs[1].state
    tc = SoilAggregateModel.TemperatureCache()
    SoilAggregateModel.update_temperature_cache!(tc, T_CONST, bio, soil)
    α_eff = SoilAggregateModel.alpha_effective(st0.E[end], st0.F_i[end], soil)
    θ0    = van_genuchten(PSI_CONST, α_eff, soil.n_vg, soil.θ_r, soil.θ_s)
    τ     = θ0^2 / soil.θ_s^(2/3)
    Rf    = θ0 / (θ0 + soil.ρ_b * soil.k_d_eq)
    D_C   = SoilAggregateModel.D_eff_DOC(tc.D_DOC_w, θ0, soil.θ_s, soil.ρ_b, soil.k_d_eq)
    D_nr  = tc.D_DOC_w * τ                       # same, with no sorption

    println("DOC TRANSPORT")
    @printf("  θ = %.4f   θ_s = %.4f   τ = %.4f\n", θ0, soil.θ_s, τ)
    @printf("  k_d = %.4g mm³/µg = %.0f L/kg   ρ_b·k_d = %.1f\n",
            soil.k_d_eq, soil.k_d_eq * 1000, soil.ρ_b * soil.k_d_eq)
    @printf("  retardation factor θ/(θ+ρ_b·k_d) = %.5f   → transport slowed %.0fx\n",
            Rf, 1/Rf)
    @printf("  D_DOC_w = %.2f   D_eff = %.4f mm²/day   (unsorbed would be %.2f)\n",
            tc.D_DOC_w, D_C, D_nr)
    println("  diffusion length √(D·t) from the POM surface at r_0 = $(round(r_0,digits=3)) mm:")
    @printf("  %8s %12s %12s %12s\n", "day", "√(D·t)", "reaches r =", "unsorbed r =")
    for t in (1.0, 3.0, 5.0, 21.0)
        @printf("  %8.1f %12.3f %12.3f %12.3f\n",
                t, sqrt(D_C*t), r_0 + sqrt(D_C*t), r_0 + sqrt(D_nr*t))
    end
    @printf("  r_max = %.3f mm\n\n", r_max)
end

# ── CSV ──────────────────────────────────────────────────────────────────────

rows = DataFrame(time_days = Float64[], node = Int[], radius_mm = Float64[],
                 C = Float64[], B = Float64[], F_n = Float64[], F_m = Float64[],
                 O = Float64[], F_i = Float64[], E = Float64[], M = Float64[],
                 binding = Float64[], G_c = Float64[], G_c_r_candidate = Float64[])
for o in res.outputs, i in 1:length(r)
    st = o.state
    push!(rows, (o.t, i, r[i], st.C[i], st.B[i], st.F_n[i], st.F_m[i], st.O[i],
                 st.F_i[i], st.E[i], st.M[i],
                 st.F_i[i] + soil.w_E * st.E[i], G_c, G_c_r[i]))
end
out_dir = joinpath(@__DIR__, "..", "output")
mkpath(out_dir)
CSV.write(joinpath(out_dir, "profiles.csv"), rows)
println("✓ profiles.csv: $(nrow(rows)) rows")

# ── Figure ───────────────────────────────────────────────────────────────────

# Log y throughout: the pools span six orders of magnitude and a linear axis
# shows only the largest.
function panel(field, title)
    p = plot(xlabel = "r [mm]", ylabel = title, title = title,
             yscale = :log10, legend = false, xlim = (r_0, r_max))
    for (k, o) in enumerate(res.outputs)
        y = max.(getfield(o.state, field), 1e-12)
        plot!(p, r, y, lw = 1.6,
              color = cgrad(:viridis)[(k - 1) / max(length(res.outputs) - 1, 1)])
    end
    return p
end

panels = [panel(:C, "C  DOC"), panel(:B, "B  bacteria"),
          panel(:F_n, "F_n"), panel(:F_m, "F_m"), panel(:F_i, "F_i insulated"),
          panel(:E, "E  EPS"), panel(:M, "M  MAOC"), panel(:O, "O  oxygen")]

# Binder panel, with both thresholds and the resulting r_agg
pb = plot(xlabel = "r [mm]", ylabel = "F_i + w_E·E  [µg/mm³]",
          title = "BINDER vs THRESHOLD", yscale = :log10,
          xlim = (r_0, r_max), legend = :topright, legendfontsize = 6)
for (k, o) in enumerate(res.outputs)
    b = max.(o.state.F_i .+ soil.w_E .* o.state.E, 1e-12)
    col = cgrad(:viridis)[(k - 1) / max(length(res.outputs) - 1, 1)]
    plot!(pb, r, b, lw = 1.8, color = col, label = "t=$(o.t)")
    ra = compute_r_agg(o, res.grid, res.params)
    scatter!(pb, [ra], [max(G_c, 1e-12)], color = col, ms = 4,
             markerstrokewidth = 0.5, label = "")
end
plot!(pb, r, fill(G_c, length(r)), lw = 3, color = :red, ls = :solid,
      label = "G_c (model, constant)")
plot!(pb, r, G_c_r, lw = 3, color = :red, ls = :dash,
      label = "G_c·r/δ_s (candidate, NOT in model)")
vline!(pb, [hy.δ_s], lw = 1, color = :gray, ls = :dot, label = "δ_s")

fig = plot(panels..., pb, layout = (3, 3), size = (1400, 1100),
           left_margin = 5Plots.mm, bottom_margin = 4Plots.mm,
           plot_title = "soil $(SOIL_ID), D = $(DIAM) mm, R_P_max = $(R_P_SCAN) — dark→light = day $(SNAP[1])→$(SNAP[end])")
savefig(fig, joinpath(out_dir, "profiles.png"))
println("✓ profiles.png")

# ── Loss, on the same configuration ──────────────────────────────────────────
# The panels above are ONE diameter; the loss is the five-diameter population.
# Both come from the same setup.jl, so the figure and the score always describe
# the same parameter set — they cannot drift apart.

println("\n" * "="^78)
println("LOSS on this configuration (5-diameter population)")
println("="^78)
let res_l = model_loss(bio, soil)
    println("  MWD [µm]")
    @printf("  %6s %12s %12s %10s\n", "day", "model", "measured", "ratio")
    for (i, d) in enumerate(MWD_DAYS)
        @printf("  %6.0f %12.1f %12.1f %10.3f\n",
                d, res_l.mwd_mod[i], MWD_OBS[i], res_l.mwd_mod[i]/MWD_OBS[i])
    end
    println("\n  CO₂ [µg-C/g-soil]")
    @printf("  %6s %12s %12s %10s\n", "day", "model", "measured", "ratio")
    for (i, d) in enumerate(CO2_DAYS)
        @printf("  %6.0f %12.1f %12.1f %10.3f\n",
                d, res_l.co2_mod[i], CO2_OBS[i], res_l.co2_mod[i]/CO2_OBS[i])
    end
    mM = sum(res_l.r_MWD)/length(res_l.r_MWD)
    sM = sqrt(max(sum(abs2, res_l.r_MWD .- mM)/length(res_l.r_MWD), 0.0))
    @printf("\n  L %.5f = L_MWD %.5f + L_CO₂ %.5f\n", res_l.L, res_l.L_MWD, res_l.L_CO2)
    @printf("  MWD mean ln residual %+.3f (×%.2f), spread %.3f\n", mM, exp(mM), sM)
end
println("="^78)

# ── Where each threshold would put the radius ────────────────────────────────

println("\n  r_agg under each threshold [mm]")
@printf("  %8s %12s %12s\n", "day", "constant", "candidate r/δ_s")
for o in res.outputs
    b = o.state.F_i .+ soil.w_E .* o.state.E
    ra_c = compute_r_agg(o, res.grid, res.params)
    ra_r = r[1]
    for i in length(r):-1:1
        if b[i] >= G_c_r[i]
            ra_r = r[i]; break
        end
    end
    @printf("  %8.1f %12.4f %12.4f\n", o.t, ra_c, ra_r)
end
println("""
  The candidate column is arithmetic on the same fields, not a model run — the
  binder profile it is applied to was produced with the constant threshold. If
  the candidate radius grows across the 21 days where the constant one does not,
  that is the case for implementing it; the profile would then have to be
  recomputed with the new threshold in force.
""")
