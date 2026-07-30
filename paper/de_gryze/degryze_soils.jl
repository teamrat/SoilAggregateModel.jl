"""
    degryze_soils.jl

The five soils of De Gryze et al. (2006) as overrides on the model's defaults.

Source: De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R. (2006)
*Water repellence and soil aggregate dynamics in a loamy grassland soil as
affected by texture.* European Journal of Soil Science 57, 235–246.
Table 1 (p. 236), Table 2 (p. 240), Table 3 (p. 242).

This file defines **no types**. `SoilProperties` and `InitialConditions` already
carry a default for every field; each soil overrides only what the paper
measures and leaves the rest alone. Model inputs and measured observations are
kept apart on purpose — `degryze_soil()` and `degryze_ic()` return things the
solver consumes, `DEGRYZE_OBSERVED` holds the numbers it is compared against.

Assumptions used here are recorded in `degryze_EJSS_2006_spec.md` §0a.

> Table 1's "Loam /%" column means SILT (Belgian convention, Fr. *limon*); the
> three texture columns sum to 100, so the triple is sand / silt / clay.
"""

# ---------------------------------------------------------------------------
# Measured, Table 1 (p. 236). Bulk density converted g/cm³ -> µg/mm³ (×1000).
# ---------------------------------------------------------------------------
# Mineralogy. All five soils are one Belgian loess-derived grassland profile
# (p. 236), so low-activity minerals apply to all of them and k_ma is not a
# per-soil field here. `maoc_capacity` then gives M_max = k_ma·f_clay_silt·ρ_b
# = 30.9, 44.1, 36.8, 38.2, 48.5 µg-C/mm³ for soils 1-5. Georgiou et al. (2022)
# via K_MA_LOW_ACTIVITY; use K_MA_HIGH_ACTIVITY for a high-activity mineralogy.
const DEGRYZE_TEXTURE = (
    (id = 1, name = "sandy loam", SOC = 0.0178, sand = 0.53, silt = 0.40, clay = 0.07, ρ_b = 1370.0),
    (id = 2, name = "silt loam",  SOC = 0.0172, sand = 0.33, silt = 0.57, clay = 0.10, ρ_b = 1370.0),
    (id = 3, name = "loam",       SOC = 0.0214, sand = 0.44, silt = 0.45, clay = 0.11, ρ_b = 1370.0),
    (id = 4, name = "loam",       SOC = 0.0233, sand = 0.44, silt = 0.43, clay = 0.13, ρ_b = 1420.0),
    (id = 5, name = "clay loam",  SOC = 0.0310, sand = 0.21, silt = 0.52, clay = 0.27, ρ_b = 1280.0),
)

# ---------------------------------------------------------------------------
# Measured observations — what the model is compared AGAINST. Never a model input.
# Table 2 (p. 240) MWD in µm; Table 3 (p. 242) large-macroaggregate formation
# rate in % per day with its R².
# ---------------------------------------------------------------------------
const DEGRYZE_OBSERVED = (
    (id = 1, MWD_day1 = 284.0, MWD_day21 = 1006.0, MWD_field = 2238.0, slope_gt2000 = 0.57, R2 = 0.85),
    (id = 2, MWD_day1 = 210.0, MWD_day21 =  809.0, MWD_field = 3002.0, slope_gt2000 = 0.59, R2 = 0.93),
    (id = 3, MWD_day1 = 198.0, MWD_day21 = 1008.0, MWD_field = 3259.0, slope_gt2000 = 0.83, R2 = 0.94),
    (id = 4, MWD_day1 = 391.0, MWD_day21 = 1269.0, MWD_field = 2707.0, slope_gt2000 = 0.91, R2 = 0.91),
    (id = 5, MWD_day1 = 310.0, MWD_day21 = 2410.0, MWD_field = 3905.0, slope_gt2000 = 2.02, R2 = 0.94),
)

"Incubation conditions, identical across all five soils (p. 237)."
const DEGRYZE_INCUBATION = (
    T_K      = 298.15,   # 25 °C
    duration = 21.0,     # days
    O2_frac  = 0.21,     # headspace flushed daily -> ambient
    WFPS     = 0.60,     # water-filled pore space
)

"""
Residue: wheat stems, identical across all five soils (p. 236-237).
1.5 g per 150 g soil; 44.3 % C; cut to 0.5-2 mm. Fibre fractions sum to 72 %,
not 100 — the balance is unaccounted for in the paper.
"""
const DEGRYZE_RESIDUE = (
    mass_fraction   = 0.01,      # 1.5 g / 150 g
    f_C             = 0.443,
    C_per_g_soil    = 4430.0,    # µg-C/g-soil
    size_range_mm   = (0.5, 2.0),
    f_soluble       = 0.11,
    f_hemicellulose = 0.03,
    f_cellulose     = 0.46,
    f_lignin        = 0.11,
)

"Mineral particle density. NOT stated in the paper; standard assumption (spec §0a A2)."
const ρ_PARTICLE = 2.65e3   # µg/mm³

"""
    degryze_soil(id; kwargs...)

`SoilProperties` for De Gryze soil `id`, overriding **only** the three fields
the paper constrains:

- `ρ_b`        — Table 1 bulk density
- `f_clay_silt` — silt + clay from Table 1
- `θ_s`        — porosity, `1 - ρ_b/ρ_PARTICLE`
- `d_32`       — Sauter mean particle diameter from the Table 1 texture triple,
                 via the package's `sauter_from_texture`. This is the length
                 that carries the texture dependence of the stability criterion
                 (`docs/REFERENCE.md` §4.4); the paper supplies the texture, the
                 package supplies the method.

Everything else keeps its package default. Extra keyword arguments are passed
through, so a run can override further fields without editing this file.
"""
function degryze_soil(id::Int; kwargs...)
    t = DEGRYZE_TEXTURE[findfirst(x -> x.id == id, DEGRYZE_TEXTURE)]
    SoilProperties(; ρ_b = t.ρ_b,
                     f_clay_silt = t.silt + t.clay,
                     θ_s = 1.0 - t.ρ_b / ρ_PARTICLE,
                     d_32 = sauter_from_texture(t.sand, t.silt, t.clay),
                     k_ma = K_MA_LOW_ACTIVITY,
                     kwargs...)
end

"""
    degryze_d32(id)

Sauter mean particle diameter [mm] for De Gryze soil `id`, from the Table 1
texture triple. Exposed separately so the values can be inspected and reported;
`degryze_soil` sets the same number on `SoilProperties`.

With the package's default class midpoints the five soils give

| soil | sand/silt/clay | d₃₂ [µm] |
|---|---|---|
| 1 | 0.53 / 0.40 / 0.07 | 8.95 |
| 2 | 0.33 / 0.57 / 0.10 | 6.33 |
| 3 | 0.44 / 0.45 / 0.11 | 6.39 |
| 4 | 0.44 / 0.43 / 0.13 | 5.73 |
| 5 | 0.21 / 0.52 / 0.27 | 3.10 |

The ORDERING is insensitive to the clay class midpoint; the SPREAD is not (see
`sauter_mean_diameter`). Halving the clay midpoint to 0.5 µm compresses the
soil-5/soil-1 ratio from 0.35 to 0.31, leaving the ranking unchanged.
"""
function degryze_d32(id::Int)
    t = DEGRYZE_TEXTURE[findfirst(x -> x.id == id, DEGRYZE_TEXTURE)]
    sauter_from_texture(t.sand, t.silt, t.clay)
end

"""
    degryze_ic(id, soil; kwargs...)

`InitialConditions` for De Gryze soil `id`, overriding only:

- `SOC` — Table 1 organic carbon
- `T_0` — incubation temperature (25 °C)
- `ψ_0` — the potential matching 60 % WFPS **for this soil**, obtained from the
  model's own retention curve rather than stored as a constant

The paper specifies water as 60 % WFPS, not as a potential. Because porosity
differs across soils, ψ must be derived per soil; it lands between −28.4 and
−29.3 kPa for all five, so water content is not what separates them.

`soil` must be the `SoilProperties` this run will use — ψ depends on its
`θ_s`, `α_vg`, `n_vg` and `θ_r`.
"""
function degryze_ic(id::Int, soil::SoilProperties; kwargs...)
    t = DEGRYZE_TEXTURE[findfirst(x -> x.id == id, DEGRYZE_TEXTURE)]
    θ_target = DEGRYZE_INCUBATION.WFPS * soil.θ_s
    ψ = van_genuchten_inverse(θ_target, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
    InitialConditions(; SOC = t.SOC,
                        T_0 = DEGRYZE_INCUBATION.T_K,
                        ψ_0 = ψ,
                        kwargs...)
end

# ---------------------------------------------------------------------------
# Day-0 MWD, µm. Digitised, and identical to Table 2's day-0 column.
# This is a MEASUREMENT OF THE CLASS DISTRIBUTION, not of aggregation: at day 0
# the sample is primary particles plus unaggregated residue, so eq. (1) can be
# inverted for the one unknown in that distribution — the sand split.
# ---------------------------------------------------------------------------
const DEGRYZE_MWD_DAY0 = (452.0, 240.0, 199.0, 406.0, 358.0)

"""
Sand coarse fraction for soils the day-0 inversion cannot reach.

**Soil 5 only.** Its measured day-0 MWD of 358 µm exceeds what its texture can
produce with no aggregates: 21 % sand, all of it in 250–2000 µm, gives 266 µm
including the residue. No sand split reaches 358, so mass must already sit in
the >2 mm class at day 0 — water-stable aggregates survived the pretreatment,
and assumption A1b is false for this soil, not approximately.

0.9 is a stated choice, not a solve. The residual −112 µm is carried and
reported rather than absorbed.
"""
const DEGRYZE_SAND_COARSE_OVERRIDE = Dict(5 => 0.9)

"""
    sand_coarse_fraction(id; f_POM_mass) -> x

Fraction of a soil's SAND that sits in the 250–2000 µm class, obtained by
inverting De Gryze eq. (1) against the measured day-0 MWD.

At day 0 the model's sample is the mineral matrix plus bare residue particles,
all of which are 0.5–2.0 mm and therefore in class B (250–2000 µm). Writing `p`
for the residue mass fraction and `w` for the eq. (1) nominals,

    MWD(0) = (1 - p)·[w_D·(silt+clay) + w_C·s·(1-x) + w_B·s·x]  +  p·w_B

which is linear in `x` and inverts directly.

# Why this replaces assumption A1

A1 assumed cumulative sand mass linear in log(diameter) — equal mass per log
interval, `x = 0.5728` for every soil — chosen with no reference to any MWD. It
is a guess about a distribution the paper does not report, and the class weights
differ by 973.5 µm, so at 44 % sand the unknown swings MWD by ±428 µm. That is
comparable to the entire day-0 signal.

The day-0 MWD measures exactly that distribution. Using it is not fitting the
model: no biological parameter reaches day 0, the aggregation machinery has not
run, and one measurement per soil determines one unknown per soil. Net degrees
of freedom in the fit: zero.

# What it assumes, and where that shows

It assumes A1b — no water-stable aggregates at day 0 — because a miss cannot
otherwise be attributed between the sand split and day-0 aggregation. A1b is
known to be false in general. Its failure is VISIBLE rather than absorbed: any
soil carrying day-0 aggregates demands `x > 1`, which this function refuses. Only
soil 5 does, and its override records that.

Soils 2 and 3 land WELL BELOW the A1 value (0.51 and 0.25 against 0.573). Since
aggregation can only add mass to coarse classes, their measured MWD being below
A1's prediction is one-directional evidence that A1 assigned them too much
coarse sand — independent of A1b, and the reason A1 could not simply be kept.
"""
function sand_coarse_fraction(id::Int; f_POM_mass::Real)
    haskey(DEGRYZE_SAND_COARSE_OVERRIDE, id) &&
        return DEGRYZE_SAND_COARSE_OVERRIDE[id]

    t = DEGRYZE_TEXTURE[findfirst(x -> x.id == id, DEGRYZE_TEXTURE)]
    w = DEGRYZE_CLASS_NOMINALS .* 1000.0          # mm -> µm, eq. (1) weights
    s = t.sand
    p = f_POM_mass

    # Strip the residue's own contribution, then solve the mineral part for x.
    mwd_mineral = (DEGRYZE_MWD_DAY0[id] - p * w[3]) / (1.0 - p)
    x = (mwd_mineral - w[1] * (t.silt + t.clay) - w[2] * s) / (s * (w[3] - w[2]))

    0.0 <= x <= 1.0 || error(
        "soil $(id): day-0 MWD $(DEGRYZE_MWD_DAY0[id]) µm needs sand coarse " *
        "fraction $(round(x, digits=3)), outside [0,1]. No sand split produces " *
        "it, so the soil carries day-0 aggregates (A1b false). Add an entry to " *
        "DEGRYZE_SAND_COARSE_OVERRIDE and report the residual.")
    return x
end

"""
    degryze_mineral_classes(id; f_POM_mass)

Mineral mass fraction in the four wet-sieve classes
`[<53, 53–250, 250–2000, >2000]` µm, for use with `population_outputs`.

- `<53 µm` = silt + clay, measured (Table 1).
- The sand split between 53–250 and 250–2000 µm is **not reported by the paper**
  and is obtained from the measured day-0 MWD by `sand_coarse_fraction` — spec
  §0a assumption **A1′**, which replaced the log-linear guess A1 on 2026-07-29.
- `>2000 µm` = 0: sand is ≤2000 µm by definition and the soil was crushed
  through a 250 µm sieve, so no primary particle can occupy that class.

`f_POM_mass` is the residue mass as a fraction of the sample at day 0 —
`I_input / f_C_POM` for the amendment being simulated. It is required, not
defaulted: the inversion is only exact for the amendment actually run.

Assumptions A1b (no day-0 aggregation) and A1c (aggregates draw mineral material
proportionally) still apply. See §0a before reporting any statistic that depends
on them.
"""
function degryze_mineral_classes(id::Int; f_POM_mass::Real)
    t = DEGRYZE_TEXTURE[findfirst(x -> x.id == id, DEGRYZE_TEXTURE)]
    x = sand_coarse_fraction(id; f_POM_mass = f_POM_mass)
    [t.silt + t.clay, t.sand * (1.0 - x), t.sand * x, 0.0]
end

"De Gryze wet-sieve series [mm] — three sieves, four classes (p. 237)."
const DEGRYZE_SIEVES = [0.053, 0.25, 2.0]

"""
Column labels for the four classes, ascending. The assay names its own classes;
`population_outputs` takes these via `class_labels` and emits `pct_<label>`.
Matches the paper's own [D%] [C%] [B%] [A%] of eq. (1).
"""
const DEGRYZE_CLASS_LABELS = ["lt53um", "um53_250", "um250_2000", "gt2000um"]

"""
Nominal diameters per class [mm] for eq. (1), p. 237:

    MWD = (5000[A%] + 1125[B%] + 151.5[C%] + 26.5[D%]) / 100   [µm]

Three are class midpoints; 5000 µm is the midpoint of 2000 µm and the 8 mm
sieve used during destructive sampling. MWD therefore saturates at 5000 µm and
is bounded below by 26.5 µm.
"""
const DEGRYZE_CLASS_NOMINALS = [0.0265, 0.1515, 1.125, 5.0]
