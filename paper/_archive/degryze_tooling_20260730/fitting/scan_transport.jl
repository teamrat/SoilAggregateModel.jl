"""
scan_transport.jl — does binder transport unpin the aggregate front?
====================================================================

Established so far:

  * R_P_max is a level knob. Over a factor of 17 the day-4 / day-21 MWD ratio
    stays at 2.8-2.9. It cannot fix the shape. (scan.jl)
  * Every POM diameter is a STEP in time: r_agg = r_0 through day 2, jumps at
    day 3, then flat for 42 days. All five step together, so a finer POM size
    distribution cannot produce a 21-day rise out of five simultaneous events.
  * E and F_i are IMMOBILE — only C, B, F_n, F_m, O diffuse. The binder's
    spatial extent is set entirely by its precursors.

Which makes the precursor transport the suspect. All three are running at about
a thousandth of their reference values:

    D_B_rel   1e-5   package default 0.001, cited 0.01 (Wu 2006)
                     — Group B anchor, 1000x off; BACKLOG item 5
    D_Fn0     1e-5   package default 0.01   — no anchor (Group C)
    D_Fm0     1e-3   package default 1.0    — no anchor (Group C)

The joint scan multiplies all three by a common factor; f = 1000 restores all
three to reference simultaneously (D_B_rel to the cited 0.01, D_Fn0 to 0.01,
D_Fm0 to 1.0). That is the direct test — raising one alone
may do nothing if the other two still pin the front.

Run at R_P_max = 0.9, the minimiser from scan.jl, so the level is roughly right
and the shape is what is being watched.

READ THE RATIO ROW, NOT L. Model/measured MWD at each measured day. The failure
mode is all seven ratios moving together, which is a level change. Success is
the early ones falling toward 1 while the late ones hold — the front propagating
over 21 days instead of jumping on day 3.

Usage:
    julia --project=. paper/de_gryze/fitting/scan_transport.jl
"""

include(joinpath(@__DIR__, "loss.jl"))

Logging.disable_logging(Logging.Info)

"""Copy of `x` with several fields replaced."""
function with_fields(x::T, pairs...) where {T}
    d  = Dict(pairs)
    kw = NamedTuple(n => get(d, n, getfield(x, n)) for n in fieldnames(T))
    return T(; kw...)
end

# Level roughly corrected first, so the shape is what moves.
const BASE_BIO  = with_field(BIO0, :R_P_max, 0.9)
const BASE_SOIL = SOIL0

const D_B_REL_0 = BASE_SOIL.D_B_rel     # 1e-5
const D_FN0_0   = BASE_BIO.D_Fn0        # 1e-5
const D_FM0_0   = BASE_BIO.D_Fm0        # 1e-3

# (label, values, apply) — apply(v) -> (bio, soil)
const SCANS = [
    ("ALL THREE ×f", [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0],
        v -> (with_fields(BASE_BIO, :D_Fn0 => D_FN0_0 * v, :D_Fm0 => D_FM0_0 * v),
              with_field(BASE_SOIL, :D_B_rel, D_B_REL_0 * v))),
    ("D_B_rel", [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2],
        v -> (BASE_BIO, with_field(BASE_SOIL, :D_B_rel, v))),
    ("D_Fn0", [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2],
        v -> (with_field(BASE_BIO, :D_Fn0, v), BASE_SOIL)),
    ("D_Fm0", [1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0],
        v -> (with_field(BASE_BIO, :D_Fm0, v), BASE_SOIL)),
]

println("="^108)
println("TRANSPORT SCAN — soil $(SOIL_ID), R_P_max = 0.9, days 0-$(Int(T_END))")
println("="^108)
println("  measured MWD [µm]: " * join((@sprintf("d%d=%.0f", d, o)
                                        for (d, o) in zip(MWD_DAYS, MWD_OBS)), "  "))
println()

t_all = @elapsed for (label, values, apply) in SCANS
    println("── $(label) " * "─"^(96 - length(label)))
    @printf("  %10s %8s %8s", "value", "L", "spread")
    for d in MWD_DAYS
        @printf(" %7s", "d$(Int(d))")
    end
    @printf(" %8s\n", "CO₂ d21")

    for v in values
        bio, soil = apply(v)
        local r
        try
            r = model_loss(bio, soil)
        catch e
            @printf("  %10.3g   FAILED: %s\n", v, sprint(showerror, e)[1:min(end, 60)])
            continue
        end
        m = sum(r.r_MWD) / length(r.r_MWD)
        s = sqrt(max(sum(abs2, r.r_MWD .- m) / length(r.r_MWD), 0.0))
        @printf("  %10.3g %8.4f %8.3f", v, r.L, s)
        for x in r.r_MWD
            @printf(" %7.2f", exp(x))
        end
        @printf(" %8.2f\n", exp(r.r_CO2[end]))
    end
    println()
end

@printf("  %.1f s total\n\n", t_all)
println("""
  The measured MWD rises 199 -> 1014 µm over days 1-21, a factor of 5.1. The
  model at baseline rises 334 -> 1420, a factor of 4.2, but reaches it by day 4
  and then decays. A ratio row that reads high-to-low across the days (4, 3, 2,
  ... 1.4) is that jump. A flat ratio row near 1 is a model that grows on the
  right timescale.

  If no value of any of these flattens the row, transport is not the limiter and
  the remaining candidate is that G_c does not depend on aggregate radius — a
  larger aggregate faces the same threshold, so nothing grades its growth. That
  is a model change, not a parameter, and needs a decision.
""")
println("="^108)
