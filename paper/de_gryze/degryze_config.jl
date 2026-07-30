# degryze_config.jl — THE configuration for this problem
#
# One definition of the soil, the initial condition, the parameters, the POM
# distribution and the domain tessellation. `run_degryze.jl` includes this and
# is currently its only consumer.
#
# Two copies existed until 2026-07-29 and had diverged: the fitting copy carried
# five corrected values the forward run did not, so the two scored different
# models (CLAUDE.md §8). Five further consumers kept their own forks; those were
# archived on 2026-07-30 (`paper/_archive/degryze_tooling_20260730/`). Anything
# added later — a fitting routine, a diagnostic — includes this file rather than
# restating any of it.
#
# Everything here is either MEASURED (from degryze_soils.jl, which reads the
# paper's tables) or a RUN CHOICE stated with its reason. No model logic.

using Distributions

include(joinpath(@__DIR__, "degryze_soils.jl"))

# ── Which soil ───────────────────────────────────────────────────────────────
# Soil 3 (loam) is the mid-texture case. The paper's own conclusion is that
# soil 5 crosses a ~15 % clay threshold and behaves differently (p. 242).
const SOIL_ID = 3

# ── Grid and time ────────────────────────────────────────────────────────────
const N_GRID    = 200        # radial nodes per aggregate
const T_MAX     = 22.0       # days; data ends at 21
const DT_OUTPUT = 0.125      # output every 3 hours

# ── POM size distribution ────────────────────────────────────────────────────
# Wheat stems cut to 0.5-2.0 mm (p. 236). Normal, µ = 1.25, σ = 0.23, in
# N_POM_BINS equal-width bins. RENORMALISED: the Normal is truncated at the 0.5
# and 2.0 mm edges, so the raw bin fractions do not sum to 1 and the shortfall
# used to appear as "total POM input 4425.1 vs expected 4430".
const POM_MEAN  = 1.25       # mm
const POM_SIGMA = 0.23       # mm
const N_POM_BINS = 10        # was 5 (0.3 mm bins) until 2026-07-29
const BIN_EDGES = collect(range(0.5, 2.0, length = N_POM_BINS + 1))
const DIAM_ALL  = [(BIN_EDGES[i] + BIN_EDGES[i+1]) / 2.0 for i in 1:N_POM_BINS]
const F_POM     = let f = diff(cdf.(Normal(POM_MEAN, POM_SIGMA), BIN_EDGES)); f ./ sum(f) end

# The amendment identity is exact for ANY bin count, and this is why:
#
#   total_POM_C = Σ N_i·P_0,i
#               = Σ [f_i·V_soil / ((4/3)π(d_i·f_pack/2)³)] · (4/3)π(d_i/2)³·ρ_POM
#               = V_soil·ρ_POM/f_pack³ · Σ f_i
#
# Every d_i cancels. The total depends on Σf_POM and nothing else, so refining
# the distribution redistributes carbon between classes without changing the
# amount. The precondition is the normalisation, asserted here rather than
# trusted — it is the same 0.11 % truncation shortfall that used to show up as
# "total POM input 4425.1 vs expected 4430".
@assert abs(sum(F_POM) - 1.0) < 1e-12 "F_POM must be normalised — the amendment identity depends on it"
@assert length(DIAM_ALL) == N_POM_BINS

# ── Soil ─────────────────────────────────────────────────────────────────────
# degryze_soil supplies only what the paper measures (ρ_b, texture, θ_s, SOC).
# The overrides below are run choices.
#
# k_L 10 -> 25000 mm³/µg   (Langmuir affinity)
#   Set by ONE physical constraint, not by any fitted observable: pore-water DOC
#   must fall in the range soils are observed to have, 10-100 mg/L. De Gryze
#   removed native POM, so background SOC is mineral-associated plus microbial —
#   there is no particulate background pool for it to sit in, and any carbon the
#   isotherm cannot hold is declared dissolved by default. At k_L = 1000 that
#   left 573 mg/L. At 25000 it is 48 mg/L and the minerals hold 28.7 of the 29.0.
#   k_L is Group C, no anchor; this is the first constraint placed on it.
#
# k_d_eq 0.05 -> 0.005 mm³/µg  (50 -> 5 L/kg)
#   Whole-soil DOC partition coefficients are typically 0.5-5 L/kg. At 50 the
#   retardation factor θ/(θ+ρ_b·k_d) is 0.0042, so DOC moves 229x slower than
#   tortuosity alone gives and the day-21 diffusion length is 1.02 mm from a POM
#   surface at r_0 = 0.625 — while r_agg sits at 1.275. The aggregate edge WAS
#   the DOC front. At 5 L/kg the factor is 0.041 and the front reaches ~3.2 mm.
#
# D_B_rel 0.001 -> 1e-5
#   §5a Group B: cited value is 0.01 (Wu 2006). Retained pending BACKLOG item 5;
#   the 2026-07-29 transport scan showed the model is insensitive to it over
#   three orders of magnitude, so it is not currently load-bearing.
#
# κ_b — NOT overridden here. The package default is 0.16 (src/parameters.jl),
#   changed from 0.0143869 on 2026-07-29 because that value only reproduced a
#   legacy G_c and had nothing measured behind it. REFERENCE.md §4.4.
#
#   G_c = τ_w·d_32/κ_b, so κ_b and the threshold move OPPOSITE ways. Every step
#   lowers w_E and raises κ_b by more than enough to offset the binder it costs,
#   because the simulated MWD sits below the measured.
#     κ_b       0.008   0.024   0.072   0.100    0.150    0.180    0.160
#     w_E       0.25    0.125   0.0625  0.05     0.05     0.05     0.05
#     G_c(δ_s)  0.0349  0.0116  0.00388 0.00279  0.00186  0.00155  0.00175  (soil 3)
#   0.2 overshot. 0.18 with A1' put day-21 MWD at 1.171 mm against 1.014
#   measured; 0.16 backs that off. Note 0.15 was measured under the OLD A1, so
#   the two are not directly comparable — A1' lowered the mineral baseline at
#   every time, not only at day 0.
#
# w_E 0.5 -> 0.25 -> 0.125 -> 0.0625 -> 0.05
#   Binder = F_i + w_E·E. w_E is the weight of EPS carbon against insulated
#   hyphal carbon and is FITTED (REFERENCE.md §4.4) — it was a hard-coded 0.5 in
#   the MATLAB precursor. Each halving halves EPS's share of the binder, so F_i
#   carries proportionally more.
#
#   IT ALSO LOWERS THE BINDER. Where EPS dominates, halving w_E roughly halves
#   F_i + w_E·E, which on its own shrinks the aggregate. κ_b is raised alongside
#   to more than offset that — see below.
#
# p_Gc 0.0 -> 1.0  (threshold now rises linearly with aggregate radius)
#   G_c(r) = G_c(δ_s)·(r/δ_s)^p_Gc. p_Gc = 0 is the package default and is the
#   model as it stood before 2026-07-29: one threshold for every aggregate size.
#   Implemented in src/postprocessing/aggregate_radius.jl; the derivation, the
#   two mechanisms that motivate it and the statement that the exponent is
#   FITTED are in critical_binding's docstring and REFERENCE.md §4.4a.
#
#   Why it is needed here, measured 2026-07-29 at p_Gc = 0: r_agg is a step —
#   r_0 through day 2, a jump on day 3, then flat for 42 days — because the flat
#   threshold cuts the binder profile out in the far field where the profile is
#   flat, while the binder accumulates 6-fold near the POM. Nothing grades
#   growth. A threshold rising with r moves the crossing inward into the region
#   still accumulating, which is the only mechanism in this model that can give
#   a radius that grows for 21 days.
#
#   p_Gc = 1 IS A STARTING POINT, not a fit. It is the form the MATLAB precursor
#   had commented out (`strength ./ x`) and the top of the range Dexter and Watts
#   report for aggregate tensile strength falling with diameter (p 0.5-1).
#   To scan it, sweep [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]. The scan harness
#   that did this is archived; see `paper/_archive/degryze_tooling_20260730/`.
#
#   κ_b sets the LEVEL of the threshold (its value at r = δ_s) and p_Gc sets its
#   SHAPE. To isolate:
#   Isolate it by setting p_Gc = 0, which is the package default and the model
#   as it stood before this change.
const SOIL = degryze_soil(SOIL_ID;
                          k_L     = 25000,
                          D_B_rel = 0.00001,
                          k_d_eq  = 0.005,
                          w_E     = 0.05,
                          p_Gc    = 1.0)

# ── Initial condition ────────────────────────────────────────────────────────
# The microbial fractions are NOT package defaults, for one reason.
#
# De Gryze p. 236: soils were air-dried, crushed through a 250 µm sieve, the
# >250 µm fraction dispersed in sodium hexametaphosphate with its organic matter
# destroyed by H₂O₂ and combustion at 500 °C, then reconstituted. Day-0 MWD is
# 199 µm. There are no macroaggregates and no native POM at t = 0.
#
#   f_insulated 0.5 -> 0.0
#     F_i is network-embedded hyphae — the aggregate-binding form. Crushing to
#     <250 µm destroys hyphal networks. At 0.5 the model opened with
#     F_i = 0.00483 = 25 % of G_c, uniform across the domain and inert for 21
#     days: half the glue at the aggregate edge was inherited, not built.
#
#   f_bact 0.01 -> 1e-4
#     Measured 2026-07-29: B collapses from 0.00966 to B_min = 1e-4, a factor of
#     97, within two days everywhere DOC does not reach. The condition is exact
#     and needs no run — growth requires C_aq > C_B, and the background gives
#     C_aq = 4.30e-5 against C_B = 5.0e-5, so the population was in net decay
#     from t = 0 and stripped 72 % of the background DOC in three hours.
#
#   f_eps 0.005 -> 5e-4
#     Same pretreatment argument, same order of reduction.
#
# Binder at t = 0 falls from 0.00725 (37.3 % of G_c) to 2.4e-4 (1.2 %), so the
# aggregate must be built by the incubation.
#
# NOT FIXED by this: nothing sustains background biomass. The far field has no
# carbon source and MAOC desorption is the only route back to DOC, so any
# biomass initialised there starves. f_bact is pinned by that, not chosen.
const IC = degryze_ic(SOIL_ID, SOIL;
                      f_insulated = 0.0,
                      f_bact      = 1.0e-4,
                      f_eps       = 5.0e-4)

# ── Biological parameters ────────────────────────────────────────────────────
# Group C (no anchor) unless noted. See REFERENCE.md §5a.
#
# κ_s_ref, κ_d_ref — working values. The ratio is the hysteresis strength.
# R_P_max, r_B_max, r_F_max, μ_B, μ_F, Y_B_max, B_S, C_B, D_Fn0, D_Fm0,
# F_*_min — working values.
const BIO = BiologicalProperties(
    κ_s_ref = 0.01, κ_d_ref = 1.0e-4,
    F_i_min = 1e-6, F_n_min = 2e-4, F_m_min = 1e-6,
    D_Fn0   = 0.00001, D_Fm0 = 0.001,
    r_B_max = 8.0, r_F_max = 0.2, R_P_max = 1.5,
    Y_B_max = 0.4, B_S = 0.05, C_B = 5.0e-5,
    μ_B = 0.0036, μ_F = 0.02,
)

# ── Environment ──────────────────────────────────────────────────────────────
# T is the incubation temperature and ψ is derived from 60 % WFPS at this soil's
# porosity. Neither is a free choice. Headspace flushed daily, so O₂ ≈ ambient.
const T_CONST   = IC.T_0
const PSI_CONST = IC.ψ_0
const O2_CONST  = DEGRYZE_INCUBATION.O2_frac * P_ATM * M_O2 / (R_GAS * T_CONST)
const T_FUNC    = t -> T_CONST
const PSI_FUNC  = t -> PSI_CONST
const O2_FUNC   = t -> O2_CONST

# ── Domain tessellation ──────────────────────────────────────────────────────
# 1.5 g wheat stems per 150 g soil (= 10 g/kg), stems 44.3 % C -> 4.43 g-C/kg.
# Geometry in src/physics/tessellation.jl.
const ρ_POM   = 200.0        # POM carbon density [µg-C/mm³]
const I_INPUT = 4.43e-3      # amendment rate [µg-C/µg-soil]
const ρ_BULK  = SOIL.ρ_b     # Table 1 bulk density
const SOIL_V  = 100.0^3      # reference volume [mm³] = 1 litre

const F_C_POM = 0.443        # wheat stem carbon fraction (p. 236)

# Residue mass as a fraction of the sample at day 0. population_statistics uses
# soil mass as the denominator for every class, so this is I_input/f_C_POM = 1 %.
# It is what the day-0 MWD inversion has to strip off before solving the mineral
# part, and it is the same number `f_agg` reports at t = 0.
const F_POM_MASS = I_INPUT / F_C_POM

# Mineral mass across the four sieve classes. The sand split is INVERTED from
# the measured day-0 MWD (spec §0a A1'), replacing the log-linear guess A1 on
# 2026-07-29. One call.
const MINERAL_CLASSES = degryze_mineral_classes(SOIL_ID; f_POM_mass = F_POM_MASS)

const TESS = domain_tessellation(ρ_POM = ρ_POM, I_input = I_INPUT, ρ_b = ρ_BULK)
const POP  = pom_population(DIAM_ALL, F_POM, TESS; ρ_POM = ρ_POM, soil_volume_mm3 = SOIL_V)

const OMEGA         = TESS.ω
const DOMAIN_FACTOR = TESS.f_domain
const SOIL_MASS_PER_L = SOIL_V * ρ_BULK * 1e-6   # grams of soil per litre

# Initial POM carbon per particle and in total, per litre of soil. Computed by
# pom_population, not recomputed here (CLAUDE.md §8).
const P_0_PER_PARTICLE = POP.P_0_per_particle
const TOTAL_POM_C      = POP.total_POM_C

# ── Solver ───────────────────────────────────────────────────────────────────
# :stiff is the default workhorse, :split the reference implementation.
# REFERENCE.md §20a.
const SOLVER_USED = :stiff   # one solver since 2026-07-30; kept for the figure label
