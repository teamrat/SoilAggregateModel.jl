"""
loss.jl — the objective
=======================

    L = (1/N_MWD) Σ ln(MWD_mod/MWD_obs)²  +  (1/N_CO2) Σ ln(CO2_mod/CO2_obs)²

Log residuals: MWD is ~1000 µm and CO₂ ~2000 µg-C/g, so absolute residuals
would weight them by an accident of units. Logs are scale-free and need no
conversion factor. The misfits are also multiplicative (×1.37 on MWD, ×1.68 on
CO₂), which is the case logs are for.

Weighted per OBSERVABLE, not per point. Soil 3 has 7 usable MWD points and 4
CO₂ points; an unweighted sum would let MWD outvote CO₂ by count alone. Replace
1/N with 1/σ² per point if Figure 4's error bars are ever digitised — the paper
plots them but does not tabulate them.

Points used: MWD days 1, 4, 7, 10, 13, 16, 21 (day 0 excluded, see setup.jl);
CO₂ days 5, 10, 15, 21 (day 0 is identically zero on both sides).

Run directly to print the baseline decomposition:
    julia --project=. paper/de_gryze/fitting/loss.jl
"""

include(joinpath(@__DIR__, "setup.jl"))

"""
    model_loss(bio, soil) -> NamedTuple

`(L, L_MWD, L_CO2, r_MWD, r_CO2, mwd_mod, co2_mod)`.

`L_MWD` and `L_CO2` are the mean squared log residual for each observable, so
each is the square of a typical relative error: 0.01 means the model is off by
about 10 % on that observable. `L` is their sum.
"""
function model_loss(bio::BiologicalProperties, soil::SoilProperties)
    df = run_model(bio, soil)

    mwd_mod = [at_time(df, :MWD_fixed_weight_mm, t) * 1000.0 for t in MWD_DAYS]
    co2_mod = [at_time(df, :CO2_total, t) / SOIL_MASS_PER_L for t in CO2_DAYS]

    all(>(0.0), mwd_mod) || error("non-positive model MWD; log residual undefined")
    all(>(0.0), co2_mod) || error("non-positive model CO₂; log residual undefined")

    r_MWD = log.(mwd_mod ./ MWD_OBS)
    r_CO2 = log.(co2_mod ./ CO2_OBS)

    L_MWD = sum(abs2, r_MWD) / length(r_MWD)
    L_CO2 = sum(abs2, r_CO2) / length(r_CO2)

    return (L = L_MWD + L_CO2, L_MWD = L_MWD, L_CO2 = L_CO2,
            r_MWD = r_MWD, r_CO2 = r_CO2, mwd_mod = mwd_mod, co2_mod = co2_mod)
end

# ── Baseline report ──────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    Logging.disable_logging(Logging.Info)

    println("="^78)
    println("LOSS — soil $(SOIL_ID), days 0-$(Int(T_END)), n=$(N_GRID)")
    println("="^78)

    t = @elapsed res = model_loss(BIO0, SOIL0)
    @printf("  evaluation: %.1f s (first call includes compilation)\n\n", t)

    println("  MWD [µm]")
    @printf("  %6s %12s %12s %10s %10s\n", "day", "model", "measured", "ratio", "ln ratio")
    for (i, d) in enumerate(MWD_DAYS)
        @printf("  %6.0f %12.1f %12.1f %10.3f %+10.3f\n",
                d, res.mwd_mod[i], MWD_OBS[i], res.mwd_mod[i]/MWD_OBS[i], res.r_MWD[i])
    end

    println("\n  CO₂ [µg-C/g-soil]")
    @printf("  %6s %12s %12s %10s %10s\n", "day", "model", "measured", "ratio", "ln ratio")
    for (i, d) in enumerate(CO2_DAYS)
        @printf("  %6.0f %12.1f %12.1f %10.3f %+10.3f\n",
                d, res.co2_mod[i], CO2_OBS[i], res.co2_mod[i]/CO2_OBS[i], res.r_CO2[i])
    end

    @printf("\n  L_MWD %.5f   (typical relative error %.1f %%)\n",
            res.L_MWD, 100 * (exp(sqrt(res.L_MWD)) - 1))
    @printf("  L_CO₂ %.5f   (typical relative error %.1f %%)\n",
            res.L_CO2, 100 * (exp(sqrt(res.L_CO2)) - 1))
    @printf("  L     %.5f\n", res.L)

    # A single systematic offset and a shape error are different problems and
    # want different parameters. Mean log residual is the offset; the spread
    # about it is the shape.
    @printf("\n  MWD: mean ln residual %+.3f (×%.2f), spread %.3f\n",
            sum(res.r_MWD)/length(res.r_MWD), exp(sum(res.r_MWD)/length(res.r_MWD)),
            sqrt(max(sum(abs2, res.r_MWD .- sum(res.r_MWD)/length(res.r_MWD)) /
                     length(res.r_MWD), 0.0)))
    @printf("  CO₂: mean ln residual %+.3f (×%.2f), spread %.3f\n",
            sum(res.r_CO2)/length(res.r_CO2), exp(sum(res.r_CO2)/length(res.r_CO2)),
            sqrt(max(sum(abs2, res.r_CO2 .- sum(res.r_CO2)/length(res.r_CO2)) /
                     length(res.r_CO2), 0.0)))
    println("""

  If the spread is small next to the mean, the misfit is one scale factor and a
  single parameter can close it. If the spread is comparable, the model has the
  wrong shape in time and no rescaling will fix it.
""")
    println("="^78)
end
