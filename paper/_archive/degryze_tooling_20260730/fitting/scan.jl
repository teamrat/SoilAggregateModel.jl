"""
scan.jl — one parameter, wide range, full loss
==============================================

The ±10 % elasticities in sensitivity.jl are local. The baseline misfit needs
changes of a factor of ~3, far outside where a linearisation holds, so the
response has to be measured over the actual range.

Scans `R_P_max` by default. Edit SCAN below for others.

Why R_P_max first: the day-resolved loss (loss.jl) shows the model is too FAST,
not merely too high — MWD peaks at day 4 and then declines while the data rises
through day 21, and the CO₂ ratio falls from ×2.58 at day 5 to ×1.68 at day 21.
R_P_max sets both the timescale of substrate release and the total released, so
it is the only single parameter that can move shape and scale together. POM is
92 % consumed by day 21 at the current value of 2.5.

Reports for each value: the loss decomposition, and the day-4 and day-21 MWD
ratios so the SHAPE can be watched separately from the total. A value that
lowers L by flattening everything is not the same as one that fixes the timing.

Usage:
    julia --project=. paper/de_gryze/fitting/scan.jl
"""

include(joinpath(@__DIR__, "loss.jl"))

Logging.disable_logging(Logging.Info)

# (owner, field, values)
const SCAN = (:bio, :R_P_max, [0.15, 0.25, 0.4, 0.6, 0.9, 1.3, 1.8, 2.5])

owner, field, values = SCAN

println("="^104)
println("SCAN — $(field) over $(values), soil $(SOIL_ID), days 0-$(Int(T_END))")
println("="^104)
@printf("  %8s %10s %10s %10s %11s %10s %11s %10s %10s\n",
        field, "L", "L_MWD", "L_CO₂", "mean lnMWD", "spread", "mean lnCO₂",
        "MWD d4", "MWD d21")
println("  " * "-"^100)

const I_D4  = findfirst(==(4.0),  MWD_DAYS)
const I_D21 = findfirst(==(21.0), MWD_DAYS)

rows = NamedTuple[]
t_all = @elapsed for v in values
    bio  = owner === :bio  ? with_field(BIO0, field, v)  : BIO0
    soil = owner === :soil ? with_field(SOIL0, field, v) : SOIL0
    local r
    try
        r = model_loss(bio, soil)
    catch e
        @printf("  %8.3g   FAILED: %s\n", v, sprint(showerror, e)[1:min(end, 60)])
        continue
    end
    mM = sum(r.r_MWD) / length(r.r_MWD)
    sM = sqrt(max(sum(abs2, r.r_MWD .- mM) / length(r.r_MWD), 0.0))
    mC = sum(r.r_CO2) / length(r.r_CO2)
    push!(rows, (v = v, L = r.L, L_MWD = r.L_MWD, L_CO2 = r.L_CO2,
                 mM = mM, sM = sM, mC = mC,
                 d4 = exp(r.r_MWD[I_D4]), d21 = exp(r.r_MWD[I_D21])))
    @printf("  %8.3g %10.5f %10.5f %10.5f %+11.3f %10.3f %+11.3f %10.2f %10.2f\n",
            v, r.L, r.L_MWD, r.L_CO2, mM, sM, mC, exp(r.r_MWD[I_D4]), exp(r.r_MWD[I_D21]))
end

@printf("\n  %.1f s for %d evaluations\n", t_all, length(rows))

if !isempty(rows)
    best = rows[argmin([r.L for r in rows])]
    @printf("\n  lowest L at %s = %.3g   (L %.5f, was %.5f at %.3g)\n",
            field, best.v, best.L, rows[end].L, rows[end].v)
    println("""
  Watch the two right-hand columns, not L. Both should approach 1.00 together.
  If d21 improves while d4 stays large, the parameter is rescaling the level and
  leaving the timing wrong — the misfit is then a shape error that needs a
  different mechanism, not a smaller rate. The spread column says the same thing:
  it must fall, not just the mean.
""")
end
println("="^104)
