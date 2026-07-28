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

"""
    degryze_mineral_classes(id)

Mineral mass fraction in the four wet-sieve classes
`[<53, 53–250, 250–2000, >2000]` µm, for use with `population_outputs`.

- `<53 µm` = silt + clay, measured (Table 1).
- The sand split between 53–250 and 250–2000 µm is **not reported**. It is
  obtained from `log_interpolate_fraction(53, 250, 2000)` — equal mass per log
  interval — which is spec §0a assumption A1, made with no reference to any
  measured MWD.
- `>2000 µm` = 0: sand is ≤2000 µm by definition and the soil was crushed
  through a 250 µm sieve, so no primary particle can occupy that class.

Two further assumptions apply to anything derived from this — A1b (no day-0
aggregation, known false) and A1c (aggregates draw mineral material
proportionally). See §0a before reporting any statistic that depends on it.
"""
function degryze_mineral_classes(id::Int)
    t = DEGRYZE_TEXTURE[findfirst(x -> x.id == id, DEGRYZE_TEXTURE)]
    f_sand = t.sand
    f_fine = log_interpolate_fraction(53, 250, 2000)
    [t.silt + t.clay, f_sand * f_fine, f_sand * (1.0 - f_fine), 0.0]
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
