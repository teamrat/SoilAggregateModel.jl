# De Gryze EJSS 2006 — incubation replication spec

**Role in this project: the VALIDATION set** (out-of-sample; different soils and a
different aggregate statistic from the fitting set).

**Extracted:** 2026-07-27
**Source:** De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R. (2006)
*Water repellence and soil aggregate dynamics in a loamy grassland soil as
affected by texture.* European Journal of Soil Science **57**, 235–246.
doi:10.1111/j.1365-2389.2005.00733.x

**Companion document:** `degryze_SBB_2005_spec.md` — covers the other incubation
(the fitting set) and untangles the three-way citation confusion between these
papers. Read its §1 first if there is any doubt about which paper is which.

Everything in §1–§4 is from the paper. My own arithmetic and judgement is in
blocks marked **[DERIVED]** or **[MY ANALYSIS]**. Where the paper is silent I
write **not stated**.

---

## 0. Read this first — the four things that will bite the comparison

*(Items 1–3 were code defects when this spec was written; all three are fixed.
They are kept because they describe the paper, and because the fix is only as
durable as the reason for it. Item 4 is a property of the data and is permanent.)*

1. **MWD is not a mean diameter.** It is a fixed-weight sum over four sieve
   classes with *assigned* nominal diameters (5000, 1125, 151.5, 26.5 µm). A
   model that computes a true volume-weighted mean diameter is **not computing
   the same quantity** and will disagree systematically. See §3.2.
   *Fixed: `population_statistics` reports both, separately — `MWD_agg_only`
   (a true mean over aggregates) and `MWD_fixed_weight` (eq. 1's form). They
   are not interchangeable and must not be plotted on one axis.*
2. **Incubation is 21 days at 25 °C.** *Fixed: `T_const = ic.T_0` = 298.15 K,
   from `DEGRYZE_INCUBATION`. `run_degryze.jl` previously ran at 20 °C — a
   ~1.4× rate error at Q₁₀ ≈ 2, in the direction of under-predicting.*
3. **Water is specified as 60 % WFPS, not a matric potential.** Because bulk
   density differs across soils (1.28–1.42 g cm⁻³), a single ψ cannot represent
   all five. See §3.1. *Fixed: `degryze_ic` derives ψ per soil through the
   model's own retention curve.*
4. **The published CO₂ is bulk soil respiration with no unamended control.**
   Soil 5 respired 119 % of the residue's carbohydrate-C by day 21, so a
   substantial and unquantifiable share of the signal is native SOC. Calibrating
   a residue-only model against these curves will over-attribute. See §5.3.
   **Not fixable in code** — it is a limitation of the data.

---

## 0a. Assumptions we make in using this data

Every assumption required to compare the model against this paper, in one place.
Each states what is measured, what is assumed, and what the assumption is *not*
allowed to touch. Nothing here may be set by inverting an observable that the
comparison then reports as agreement.

### A1'. Sand particle-size distribution — inverted from the measured day-0 MWD

*Replaced the log-linear assumption A1 on 2026-07-29. A1 is kept below as
superseded, because the argument for dropping it is part of the result.*

**Problem.** Eq. (1) needs mass in four classes. Silt+clay (Table 1) gives the
<53 µm class directly. But sand spans *two* classes — 53–250 and 250–2000 µm —
and the paper reports only the total. The weights differ by 1125 - 151.5 = 973.5
µm, so at 44 % sand the unknown split swings MWD by ±428 µm, comparable to the
entire day-0 signal.

**What is done now.** At day 0 the sample is primary particles plus unaggregated
residue, and all residue is 0.5-2.0 mm, hence in class B. Writing `p` for the
residue mass fraction (= `I_input/f_C_POM` = 0.01) and `w` for the eq. (1)
nominals,

```
MWD(0) = (1 - p)·[w_D·(silt+clay) + w_C·s·(1-x) + w_B·s·x]  +  p·w_B
```

linear in `x`, inverted per soil in `degryze_soils.jl :: sand_coarse_fraction`.

| soil | sand % | day-0 MWD | `x` = sand in 250-2000 | A1 said |
|---|---|---|---|---|
| 1 | 53 | 452 | 0.683 | 0.573 |
| 2 | 33 | 240 | 0.508 | 0.573 |
| 3 | 44 | 199 | **0.253** | 0.573 |
| 4 | 44 | 406 | 0.741 | 0.573 |
| 5 | 21 | 358 | **1.455 — impossible** | 0.573 |

**This is not fitting the model.** No biological parameter reaches day 0; the
aggregation machinery has not run; one measurement per soil determines one
unknown per soil. Net degrees of freedom in the fit: zero. What it does is
replace a guess about a distribution the paper never reports with the one
measurement that observes that distribution directly.

**It assumes A1b.** A miss cannot otherwise be attributed between the sand split
and day-0 aggregation, so the inversion charges everything to the split. A1b is
known false in general (see below), and the inversion does not hide that — a
soil carrying day-0 aggregates demands `x > 1`, which `sand_coarse_fraction`
refuses with an error rather than clamping.

**Soil 5 is that case.** 21 % sand, all of it coarse, gives 266 µm including the
residue against 358 measured. No sand split reaches it. Mass must already sit in
the >2 mm class at day 0, so water-stable aggregates survived the pretreatment.
`x = 0.9` is set by hand for soil 5 and the residual, **-112 µm**, is carried and
reported. It is a stated mismatch, not an absorbed one.

**What the numbers say about A1.** Soils 2 and 3 come out well below 0.573 —
0.508 and 0.253. That direction is established independently of A1b: aggregation
can only add mass to coarse classes, so measured MWD below A1's prediction means
A1 gave those soils too much coarse sand. Soil 3 needs less than half of it. And
the five values span 0.25 to >1, so no single sand-shape rule fits five soils
from one meadow — A1's premise, not merely its constant, was wrong.

---

### A1 (SUPERSEDED). Sand particle-size distribution — log-linear

**Assumption.** Cumulative sand mass linear in log(diameter) across the 53-2000
µm span, i.e. equal mass per log interval:

```
f(53-250)   = ln(250/53) / ln(2000/53) = 0.4272
f(250-2000) = 1 - 0.4272                = 0.5728
```

Applied to each soil's measured sand percentage. **This was a statement about
particle size distributions, made without reference to any measured MWD** — that
independence was its whole merit, and it is what A1' gives up in exchange for
using the measurement.

**A1 describes DISPERSED soil.** Texture in Table 1 is measured after chemical
dispersion; it is the primary-particle distribution. MWD is measured by wet
sieving of UNDISPERSED soil, where primary particles may be bound together.
These are two different material states and A1 alone cannot predict the second.
A1' inherits this problem: the inversion produces the undispersed class
distribution directly, which is the right target, but it can no longer be
checked against Table 1 independently.

### A1b. No aggregation at day 0 — known to be false

To turn a sand split into a day-0 prediction, a second assumption is required:
that the reconstituted soil contains no water-stable aggregates, so the wet-sieve
classes are populated by primary particles only.

**This is wrong.** The 5-minute submersion and sieving act on undispersed
material, and re-aggregation on wetting is not prevented by the pretreatment.
Soil 5 proves it outright: its day-0 MWD exceeds what its texture can produce.

It is adopted for the other four soils because nothing in the paper constrains
day-0 aggregate content — only MWD is reported, never the size distribution.
Under A1' the cost is explicit: whatever day-0 aggregation those soils carry is
charged to their sand split instead.

A related consequence: **the <53 µm class equals silt+clay only in dispersed
soil.** In undispersed soil, clay bound into aggregates is retained on coarser
sieves, so assigning Table 1's silt+clay to the <53 µm class is itself part of
A1b, not an independent measurement.

### What the day-0 comparison used to be

Under A1 it was a plausibility check and explicitly not a target: a miss could
not be attributed between A1 and A1b, so fitting to it would tune one assumption
to absorb error in the other. Three of five soils landed within ±20 %:

| Soil | sand % | 53–250 | 250–2000 | <53 | MWD pred | MWD obs | obs − pred |
|---|---|---|---|---|---|---|---|
| 1 | 53 | 22.64 | 30.36 | 47.0 | 388 | 452 | +64 |
| 2 | 33 | 14.10 | 18.90 | 67.0 | 252 | 240 | **−12** |
| 3 | 44 | 18.80 | 25.20 | 56.0 | 327 | 199 | **−128** |
| 4 | 44 | 18.80 | 25.20 | 56.0 | 327 | 406 | +79 |
| 5 | 21 | 8.97 | 12.03 | 79.0 | 170 | 358 | +188 |

**Under A1' this check is spent.** Days 1 onward remain independent of it, and
soil 5's −112 µm residual is the one part of day 0 that still tests anything.

### A1c. Well-mixed matrix: material is absorbed in proportion to abundance (t > 0)

**Restated 2026-07-28.** Previously worded as "aggregates draw mineral material
proportionally", which reads as a claim about uptake physics. It is not one.

**What is measured.** Nothing about aggregate composition. De Gryze reports
sieve-class mass fractions of the whole sample; the paper never disperses a
fraction to ask what it is made of.

**What is assumed.** Given no information about which particles a growing
aggregate takes up, the matrix is taken to be well mixed and material
incorporated **in proportion to its relative abundance**. The unaggregated
remainder therefore keeps the same class distribution and only its total falls.
This is the least-assuming closure available — a statement about what the model
does not track.

**The classes are this assay's sieve series** — 53 / 250 / 2000 µm — not
sand / silt / clay. The algorithm in `population_statistics` is universal; only
the cutoffs are problem-specific, and mapping Table 1's texture onto these three
sieves is `degryze_mineral_classes`'s job.

**What this is NOT: a claim that the composition is knowable.** The model is a
continuum — a homogeneous shell of bulk-density matrix around a POM core — while
the class accounting is discrete. The two are compatible only while the shell is
much thicker than the particles it is said to contain. **A 0.6 mm aggregate
cannot engulf a 0.5 mm grain.**

At this experiment's operating point that limit is breached, not approached.
Day 21, soil 3, from the model's own output:

| Quantity | Value |
|---|---|
| mass-weighted shell thickness | 485–711 µm (`shell_mm` column) |
| coarse mineral class | 250–2000 µm |
| share of that class coarser than the shell | ~50–68 % |
| **share of the whole sample coarser than its own shell** | **~14 %** |
| `f_agg` from day 5 onward | 0.25 |

So the closure assigns roughly a seventh of the sample to aggregate interiors it
could not physically occupy, while a quarter of the sample is claimed as
aggregate. **Total mass balance is unaffected** — the four columns still sum to
100, and `pct_gt2000um` (whose mineral fraction is zero) is untouched by this
entirely. What is not defensible in detail is the *composition* of the
unaggregated remainder, and therefore the split between the three lower classes.

**Deliberately not fixed by a size rule.** Restricting absorption to particles
finer than the aggregate would splice a discrete particle model into a continuum
one — the same incoherence, not a repair. Resolving it means choosing one
representation, which is a modelling decision, not a bookkeeping one.

**Conventions fixed 2026-07-28.** Two choices in this bookkeeping were settled
explicitly rather than left implicit:

- **Shell mass uses BULK density, not particle density.** ρ_b is dry solid mass
  per total volume, so `ρ_b·V_shell` already *is* the solid mass in the shell —
  pores contribute volume, not mass. Using ρ_s = 2650 would assert the shell is
  pore-free, contradicting θ_s = 1 − ρ_b/ρ_s = 0.483, which the same model uses
  for water inside that same shell. (Factor 1.93 on every aggregate mass.)
- **The remainder is unaggregated MATRIX, not unaggregated mineral.** It is
  `sample − (shell + POM core)` — what the sieve weighs as not-retained,
  including residue-derived carbon dispersed into the matrix. This keeps the
  four columns summing to 100 % of the real sample. The mineral-only reading
  differs by the core's share: 1.00 % (day 0), 0.31 % (day 21), 0.08 % (day 45).

**Practical guidance.** Three diagnostic columns are emitted in
`population.csv` for every run:

| Column | Read it against | Failure signal |
|---|---|---|
| `shell_mm` | the coarse end of `DEGRYZE_SIEVES` | shell thinner than the coarse class ⇒ the three finer class shares are indicative only; `pct_gt2000um` is unaffected |
| `f_agg` | 1 | > 1 means more aggregate than sample |
| `max_cell_occ` | 1 | ≥ 1 means some POM size class consumed more soil than it owns |

**`f_POM` must be normalised.** `run_degryze.jl` truncates `Normal(1.25, 0.23)`
at the 0.5 and 2.0 mm bin edges, so the raw bin fractions summed to 0.99889.
That 0.11 % was the entire residual in "Total POM input 4425.1 vs expected
4430" — 4430 × 0.99889 = 4425.1. Renormalised as of 2026-07-28; the run should
now recover 4430.0.

### A2b. Bulk density is taken to be that of the AMENDED soil

**Not stated.** Table 1 gives ρ_b = 1.28–1.42 g cm⁻³ without saying whether it
was measured before or after the residue was mixed in.

**Assumed:** that it is the amended value, so `ρ_b · V` is the total sample mass
— soil plus residue. This removes any need to track the amendment separately
when forming `f_agg` and the sieve-class shares.

**Magnitude:** the amendment is 1.5 g per 150 g soil = 1 % of soil mass, so the
amended and bare-soil readings differ by ~1 %. The same ρ_b sets the aggregate
shell density, the POM volume fraction in the tessellation, and porosity through
θ_s = 1 − ρ_b/2650, so the convention propagates into all three at that level.

**Update this if a measured amended bulk density becomes available.**

### A2. Particle density ρ_s = 2.65 g cm⁻³

Not stated in the paper. Needed for porosity → WFPS → θ, and for the aggregate
shell mass. Standard mineral value.

### A3. CO₂ is compared total-to-total

The model integrates respiration from **all** carbon it carries — the residue
*and* the initial microbial, EPS and MAOC pools — into `CO2_cumulative`. It is
not a residue-only quantity. The published Figure 3 is likewise bulk soil
respiration with no unamended control.

**Both sides therefore include background mineralization, and the comparison is
made directly, total to total. No partition factor is applied.** The model's
carbon closure (`C_balance_error` ~1e-13) means counting the respiration flux
and subtracting current pools from initial total give the same number, so there
is no ambiguity about what is being compared.

### A4. What the model does not produce

The model carries no particle-size distribution for unaggregated mineral matrix.
It predicts the aggregate size classes only. Whole-sample eq. (1) MWD therefore
requires A1'; the class percentages and `pct_gt2000` do not, and `pct_gt2000` is
free of A1' entirely — sand is ≤2000 µm by definition and the soil was crushed
through 250 µm, so no primary particle can occupy that class. That makes
`pct_gt2000` the one comparison untouched by the day-0 inversion.

---

## 0b. How this specification is wired into the code  *(2026-07-27)*

`paper/de_gryze/degryze_soils.jl` holds everything in §1–§3 that the model
consumes, and nothing else. It defines **no types**: `degryze_soil(id)` and
`degryze_ic(id, soil)` return the package's own `SoilProperties` and
`InitialConditions`, overriding only the fields this paper measures —

| Override | From | Note |
|---|---|---|
| `ρ_b` | Table 1 | per soil, 1280–1420 µg/mm³ |
| `f_clay_silt` | Table 1 | silt + clay |
| `θ_s` | derived | `1 − ρ_b/2650`, assumption A2 |
| `SOC` | Table 1 | per soil |
| `T_0` | p. 237 | 298.15 K |
| `ψ_0` | derived | 60 % WFPS through the model's own retention curve, per soil |

— and leaving every other package default untouched. Measured **observations**
(Table 2 MWD, Table 3 slopes) are held separately in `DEGRYZE_OBSERVED`; they
are never model inputs.

`ψ_0` is computed rather than stored: the paper specifies water as 60 % WFPS,
and porosity differs across soils, so ψ must be derived per soil. It lands
between −28.4 and −29.3 kPa for all five — water content is not what
distinguishes them.

Domain geometry (`ω`, `f_pack`, `N_POM`) moved to `src/physics/tessellation.jl`
and is documented in `docs/REFERENCE.md` §5b. It was previously duplicated by
copy-paste in `run_degryze.jl` and `optimize_soil3.jl`, which had diverged —
`optimize_soil3.jl` was still running at 20 °C against the paper's 25 °C.
(`optimize_soil3.jl` was archived on 2026-07-30 to
`paper/_archive/degryze_tooling_20260730/`. It is named below only as history;
nothing in it is live.)

**Sieve-class outputs.** `population_outputs` takes the sieve series and the
class labels from the caller, because the sieve series is a property of the
assay, not of the model: `DEGRYZE_SIEVES = [0.053, 0.25, 2.0]` and
`DEGRYZE_CLASS_LABELS = ["lt53um", "um53_250", "um250_2000", "gt2000um"]`,
giving columns `pct_lt53um` … `pct_gt2000um`. These are eq. (1)'s own
[D%] [C%] [B%] [A%]. Because `run_degryze.jl` also supplies
`mineral_class_fractions` (from `degryze_mineral_classes`), the unaggregated
remainder is distributed across the same classes and the four columns sum to
100. Without that input they would carry aggregate mass alone and sum to
`f_agg·100`.

Cumulative "fraction above sieve X" is no longer reported. Nested curves cannot
show mass moving between classes, which is the signal; it is the running sum of
the class columns from the top if ever wanted.

`pct_gt2000um` is the comparison target for Table 3, and the only class column
free of assumption A1' — sand is ≤2000 µm by definition and the soil was crushed
through 250 µm, so no primary particle can occupy that class.

Assumption A1's sand split is applied by `degryze_mineral_classes(id;
f_POM_mass)`. Since 2026-07-29 it is **A1'**: the split comes from
`sand_coarse_fraction(id; f_POM_mass)`, which inverts eq. (1) against the
measured day-0 MWD in `DEGRYZE_MWD_DAY0`, and errors rather than clamping when
no split can reach the measurement. Soil 5 is the one such case and is set by
`DEGRYZE_SAND_COARSE_OVERRIDE`. The superseded log-linear A1 used
`log_interpolate_fraction` from `src/`; that helper is untouched and still
general. A1b and A1c still govern anything derived from this.

**Population statistics moved to `src/`** *(2026-07-28)*. Aggregate mass, sieve
binning, the fixed-weight MWD and the population sums are now
`population_statistics` / `aggregate_mass` / `sieve_class` in
`src/postprocessing/population.jl`, documented in `docs/REFERENCE.md` §5c. None
of it was De Gryze-specific. What remains in
`paper/de_gryze/postprocess_dataframe.jl` is packaging: column extraction,
DataFrame construction, class-column naming.

`population_outputs` **no longer takes an `ω` argument.** It never used it: the
sums are over physical particle counts and physical diameters and are already
per-sample totals, so ω would have been applied twice (§5b). The argument was
accepted and silently ignored, which is a trap, so it is gone. `run_degryze.jl`
and both call sites in `optimize_soil3.jl` were updated; the two optimizer
calls now pass `ρ_b` and `f_C_POM` explicitly instead of relying on the
wrapper's defaults — the values are unchanged (1300.0 and 0.443), so the
optimizer's numbers are unaffected.

`ρ_b` and `f_C_POM` have **no defaults** in the `src/` function. `f_C_POM` is
the residue carbon fraction (0.443 here, p. 236) and belongs to the experiment,
not to the model.


## 0c. Texture enters the model through d₃₂  *(2026-07-28)*

The aggregate stability criterion carries a texture dependence:

    F_i(r) + w_E·E(r)  ≥  G_c = τ_w · d₃₂ / κ_b

so a **coarser soil needs a higher binder concentration to hold**. The form comes
from bond counting plus Albalasmeh & Ghezzehei (2014); it is derived in
`docs/REFERENCE.md` §4.4 and is **not** fitted to this paper. De Gryze is the
test.

`degryze_soil(id)` sets `d_32` from the Table 1 texture triple using the
package's `sauter_from_texture`. The paper supplies the texture; the package
supplies the method.

`G_c` below is at `κ_b = 0.16` — the value at `r = δ_s` once the threshold
became size-dependent (REFERENCE.md §4.4a).

| soil | sand / silt / clay | d₃₂ [µm] | G_c(δ_s) [µg/mm³] | G_c ÷ soil 3 |
|---|---|---|---|---|
| 1 | 0.53 / 0.40 / 0.07 | 8.954 | 0.00244 | 1.400 |
| 2 | 0.33 / 0.57 / 0.10 | 6.327 | 0.00173 | 0.990 |
| 3 | 0.44 / 0.45 / 0.11 | 6.394 | 0.00175 | 1.000 |
| 4 | 0.44 / 0.43 / 0.13 | 5.734 | 0.00157 | 0.897 |
| 5 | 0.21 / 0.52 / 0.27 | 3.099 | 0.00085 | 0.485 |

**The last column is the content of this table; the fourth is not.** `κ_b` is a
single fitted constant multiplying every row, so it sets the level and cancels
from the ratios. The texture ordering and spacing come from `d₃₂` alone, out of
Table 1, with no fitting.

`κ_b` was 0.0143869 until 2026-07-29, chosen as `2.25 × d₃₂(soil 3)` so that
soil 3 reproduced a legacy `G_c = 0.0194` exactly. That equality preserved a
number with a dimensionally inconsistent derivation behind it (REFERENCE.md §26
erratum 11) and constrained nothing, so it was dropped.

### What this predicts, and how it fares

Against Table 3's large-macroaggregate formation rates (the cleaner target — see
§3.2, since day-21 MWD still contains the day-0 texture baseline):

| soil | G_c ÷ soil 3 | slope /% day⁻¹ | MWD day 21 /µm |
|---|---|---|---|
| 1 | 1.40 (highest threshold) | 0.57 (slowest) | 1006 |
| 2 | 0.99 | 0.59 | 809 |
| 3 | 1.00 | 0.83 | 1008 |
| 4 | 0.90 | 0.91 | 1269 |
| 5 | 0.48 (lowest threshold) | 2.02 (fastest) | 2410 |

Rank correlation ρ = 0.90 against the slopes, 0.70 against day-21 MWD. **With
n = 5 that is not a statistical result and must not be reported as one.** What
it is: the coarsest soil has the highest threshold and the slowest measured
formation, the finest has the lowest threshold and a rate 3.4× the rest, and
soils 1–4 cluster within a factor of 1.6 while soil 5 sits apart — which is the
paper's own "threshold at about 15 % clay" (p. 242) emerging from a smooth
relation rather than being imposed.

**Known shortfall, to be reported rather than tuned away.** Soil 5 gets a 2×
threshold advantage while its measured formation rate is 3.4× the pack, so the
model will under-predict how different soil 5 is. A geometric-mean diameter
would give it a 3.2× advantage but puts soil 2 in the wrong place; the Sauter
mean is the correct average for a bond count and is kept for that reason, not
for the fit.

**Assumption.** `d₃₂` is computed from a texture *triple*, not a measured
particle-size distribution, so the class midpoints in `TEXTURE_CLASS_DIAMETERS`
(clay 1 µm, silt 10 µm, sand 316 µm) are assumed. `d₃₂` is dominated by the
finest class, so the clay midpoint matters most: halving it to 0.5 µm changes
the soil-5/soil-1 ratio from 0.35 to 0.31 and leaves the ranking unchanged.


## 1. Soils

Sampled **autumn 2003**, fertilized sheep meadow, K.U. Leuven Centre for Animal
Husbandry, **Lovenjoel, Belgium (50°51′N, 4°47′E)**. Gleyic and Haplic Luvisols
(FAO 1998). Depth **0–10 cm**. Five soils selected along a natural texture
gradient.

**Table 1 (p. 236) — general characteristics of the five crushed (< 250 µm) soils.**
Values are means ± standard error.

| Soil | pH¹ | CEC² /cmol_c kg⁻¹ | Base sat. /% | Organic C /% | C/N | Sand /% | Loam /% | Clay /% | Texture | BD³ /g cm⁻³ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 4.4 | 5.6 ± 0.1 | 72 ± 0.7 | 1.78 ± 0.04 | 11.5 | 53 | 40 | 7 | Sandy loam | 1.37 |
| 2 | 4.6 | 6.1 ± 0.1 | 72 ± 0.4 | 1.72 ± 0.02 | 11.9 | 33 | 57 | 10 | Silt loam | 1.37 |
| 3 | 4.8 | 6.8 ± 0.2 | 90 ± 0.3 | 2.14 ± 0.01 | 9.5 | 44 | 45 | 11 | Loam | 1.37 |
| 4 | 4.2 | 6.6 ± 0.1 | 63 ± 1.5 | 2.33 ± 0.05 | 11.1 | 44 | 43 | 13 | Loam | 1.42 |
| 5 | 5.3 | 17.4 ± 0.3 | 84 ± 2.3 | 3.10 ± 0.01 | 9.7 | 21 | 52 | 27 | Clay loam | 1.28 |

¹ 10 g soil in 25 ml of 0.01 mol l⁻¹ CaCl₂. ² Cation exchange capacity.
³ Bulk density, calculated on a stone-free basis.

> **⚠ "Loam /%" in Table 1 means SILT.** This is the Belgian/European texture
> convention (Fr. *limon*). The three columns sum to 100 % for every soil, so the
> triple is **sand / silt / clay**. Do not read the "Loam" column as a loam
> *class* fraction.

**[DERIVED] Silt + clay fraction** (`f_clay_silt` in your `SoilProperties`):

| Soil | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| `f_clay_silt` | **0.47** | **0.67** | **0.56** | **0.56** | **0.79** |

Note the paper's own framing: pH is **4.2–5.3** (acidic), and soil 5 is distinct
on *every* axis — highest clay (27 %), highest CEC (17.4, ~2.6× the others),
highest organic C (3.10 %), highest pH, highest respiration, and 2–3× the MWD
response. Soils 1–4 form a tight cluster.

### 1.1 Soil pretreatment (p. 236, verbatim)

> "Briefly, soil was taken from the 0–10 cm depth layer, air-dried and crushed
> through a **250-µm aperture sieve**. Material larger than 2 mm (plant residues
> and stones) was discarded. Material larger than 250 µm – native particulate
> organic matter (POM) and sand – was dispersed using an aqueous solution of
> sodium hexametaphosphate 50 g l⁻¹. **To ensure that all five soils contained the
> same amount and quality of particulate organic matter during the incubation,
> the organic matter was removed from this fraction using an aqueous solution of
> 10 % hydrogen peroxide and a combustion (500 °C) treatment.** Sand and material
> smaller than 250 µm were then reconstituted according to the original
> proportion of these fractions and homogenized by end-over-end shaking."

**Four consequences for initialization:**

1. **Initial POM ≈ 0.** Native particulate OM was destroyed. All native C at t=0
   is mineral-associated / microbial. Your `run_degryze.jl:210–211` comment
   already states this correctly.
2. **Initial aggregates ≈ 0.** The soil was crushed through 250 µm and
   mechanically dispersed. There is no pre-existing macroaggregate structure.
3. **But day-0 MWD is 199–452 µm, not < 250 µm** — because the **sand fraction
   (> 250 µm) was added back**. At t = 0 the > 250 µm material is *primary sand
   grains*, not aggregates. Any model that treats the day-0 > 250 µm mass as
   aggregate is mis-attributing it.
4. Soil sampling/preparation is described more fully in the **companion paper**,
   De Gryze et al. (2005) — see §6.

---

## 2. The incubation

### 2.1 Setup (p. 237, verbatim)

> "From each of the five reconstituted soils, **14 samples** were prepared by
> mixing **150 g of soil with 1.5 g wheat stems** (De Gryze *et al.*, 2005). Water
> was added to obtain **60 % water-filled pore space**. To account for water uptake
> by the wheat residue, **an extra 1.5 ml of water was added**. This extra volume
> was determined in a preliminary experiment. After mixing, the soil was
> transferred to **stainless-steel cylinders (5 cm diameter)**, closed at the bottom
> by a PVC disk glued to the cylinder. Several times during transfer, **a plunger
> was used to compress the soil material evenly and bring all soils to the bulk
> densities of the untreated field samples (Table 1)**. Soils were incubated for
> **21 days at 25 °C** in glass canning jars. During the incubation, respiration
> (CO₂) was measured **every 3 hours** using a GC-8A gas chromatograph coupled to
> an automated gas measurement system (Swerts *et al.*, 1995). **Each day, the
> canning jars were flushed with compressed air**, and the moisture content of the
> samples was checked and, if necessary, adjusted. Samples were destructively
> analysed at **days 0 (control), 1, 4, 7, 13 and 21**."

**Destructive sampling procedure (same paragraph):**

> "This was done by removing the PVC disks from the cylinders, and gently pushing
> the soil material out of the cylinders using a plunger with the same diameter as
> the cylinders to minimize disturbance. The soil was then gently crumbled along
> the natural planes of weakness, until all material passed **an 8-mm aperture
> sieve**. After air-drying to constant weight, subsamples of 20 g were taken for
> water repellence measurements, and the remaining soil was used for aggregate
> separation."

### 2.2 Residue (p. 236, verbatim)

> "Leaves were separated from the stems of wheat residue (*Triticum aestivum* L.),
> which were then cut into pieces of **0.5–2 mm**. This residue contained **44.3 %
> carbon (C) and 0.17 % nitrogen (N)** and had **a C:N ratio of 261**. Residue
> composition according to the acid detergent fibre method (Van Soest & Wine,
> 1967) was **11 % soluble fraction, 3 % hemicellulose, 46 % cellulose, 11 %
> lignin and 1 % inorganic fraction**."

> ⚠ The material incubated is **wheat *stems*** — leaves were separated out and
> discarded. The fibre fractions sum to **72 %**, not 100 %; the balance is
> unaccounted for in the paper.

### 2.3 Consolidated setup table

| Quantity | Value | Source |
|---|---|---|
| Soil mass per replicate | **150 g** | p. 237 |
| Residue mass per replicate | **1.5 g wheat stems** (= 1 % w/w) | p. 237; fig. captions |
| Residue particle size | **0.5–2 mm** | p. 236 |
| Residue C | **44.3 %** | p. 236 |
| Residue N | **0.17 %** | p. 236 |
| Residue C:N | **261** | p. 236 |
| Residue fibre | 11 % soluble, 3 % hemicellulose, 46 % cellulose, 11 % lignin, 1 % inorganic | p. 236 |
| Water content | **60 % WFPS + 1.5 ml** | p. 237 |
| Bulk density | **per-soil, Table 1** (1.28–1.42) — compacted to field BD | p. 237 |
| Temperature | **25 °C**, constant | p. 237 |
| Duration | **21 days** | p. 237 |
| Vessel | stainless-steel cylinder, **5 cm diameter**, PVC-sealed base, in glass canning jar | p. 237 |
| Headspace | **flushed daily with compressed air** | p. 237 |
| Respiration cadence | **every 3 h** (GC-8A) | p. 237 |
| Destructive sampling | days **0, 1, 4, 7, 13, 21** *(but see §5.1)* | p. 237 |
| Replicates | 14 samples prepared per soil; **6 replicates at day 0** | p. 237 |
| Moisture control | checked daily, adjusted if necessary | p. 237 |

### 2.4 [DERIVED] Quantities you actually need in the model

Assuming particle density **ρ_s = 2.65 g cm⁻³** (**not stated in the paper** —
this is my assumption; porosity = 1 − BD/ρ_s).

| Soil | BD /g cm⁻³ | Porosity | θ at 60 % WFPS /cm³ cm⁻³ | Gravimetric w /g g⁻¹ | Soil vol /cm³ | Cylinder height /cm | Water added /ml |
|---|---|---|---|---|---|---|---|
| 1 | 1.37 | 0.4830 | **0.2898** | 0.2115 | 109.5 | 5.58 | 31.7 + 1.5 = **33.2** |
| 2 | 1.37 | 0.4830 | **0.2898** | 0.2115 | 109.5 | 5.58 | 31.7 + 1.5 = **33.2** |
| 3 | 1.37 | 0.4830 | **0.2898** | 0.2115 | 109.5 | 5.58 | 31.7 + 1.5 = **33.2** |
| 4 | 1.42 | 0.4642 | **0.2785** | 0.1961 | 105.6 | 5.38 | 29.4 + 1.5 = **30.9** |
| 5 | 1.28 | 0.5170 | **0.3102** | 0.2423 | 117.2 | 5.97 | 36.4 + 1.5 = **37.9** |

*Cylinder height is inferred from 150 g ÷ BD ÷ π(2.5 cm)²; the paper gives only
the 5 cm diameter.*

**Residue inputs per gram of soil** — 1.5 g × 44.3 % C ÷ 150 g:

| Quantity | Value |
|---|---|
| Total residue C | **4430 µg C g⁻¹ soil** |
| Total residue N | **17.0 µg N g⁻¹ soil** |
| Carbohydrate-C (soluble + hemicellulose + cellulose = 60 %) | **2658 µg C g⁻¹ soil** |
| Lignin-C (11 %) | **487 µg C g⁻¹ soil** |

> ✅ **Cross-check passed.** The paper independently states (p. 241): "we added
> only **2.6 mg carbohydrate-C g⁻¹** oven-dry soil." My 2658 µg C g⁻¹ reproduces
> this to three figures, which confirms both the 4430 µg C g⁻¹ total and the
> reading that "carbohydrate" = soluble + hemicellulose + cellulose.

---

## 3. The two observables

### 3.1 Respiration

- Measured **every 3 h** by GC; reported as **cumulative respiration, µg C g⁻¹
  soil** (Figure 3, p. 239).
- **Headspace flushed daily with compressed air** → CO₂ cannot accumulate to
  inhibitory levels and O₂ stays at ambient. ✅ Your constant-21 % O₂ assumption
  (`run_degryze.jl:80`) is well justified by this.
- **Rank order (p. 239, verbatim):** "The total cumulative respiration increased
  in the following order: soil 1 = soil 2 = soil 4 < soil 3 < soil 5."
- **Peak rate (p. 239, verbatim):** "The largest respiration rate was observed
  **between days 3 and 6** of the incubation."

**Digitized values in `degryze_CO2_2006.csv` (µg C g⁻¹ soil), which I verified
against Figure 3 — they are consistent:**

| Day | Soil 1 | Soil 2 | Soil 3 | Soil 4 | Soil 5 |
|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 | 0 |
| 5 | 571.8 | 408.9 | 482.0 | 356.2 | 648.2 |
| 10 | 1051.9 | 800.5 | 1263.6 | 726.9 | 1557.3 |
| 15 | 1339.0 | 1170.8 | 1766.2 | 1074.9 | 2423.7 |
| 21 | 1592.4 | 1490.4 | 2139.0 | 1448.2 | 3156.6 |

### 3.2 MWD — **the definition you must match exactly**

**Aggregate separation (p. 237, verbatim):**

> "We separated aggregates by wet sieving (Elliott, 1986). Before sieving, **50 g
> soil was submerged for 5 minutes in deionized water**. A series of three sieves
> was used consecutively to obtain four fractions: **> 2000 µm** (= large
> macro-aggregates), **250–2000 µm** (= small macro-aggregates), **53–250 µm** (=
> micro-aggregates), and **< 53 µm** (= silt and clay fraction). The large and
> small macro-aggregates were obtained by **slowly moving the sieve 2 cm in and out
> of the water for 50 cycles in 2 minutes**, and then pouring the remaining water
> and soil onto the next sieve… The micro-aggregate fraction was obtained by
> pouring the remaining water and soil on the sieve and **washing this fraction
> vigorously until the washing water was completely clear**. This adjustment to the
> original method was necessary to remove clay material completely from these
> loamy soils. All fractions were **oven-dried (50 °C)**. **Six replicates** of each
> soil were sieved at the start of the experiment, to obtain a measurement of
> aggregate stability at day 0."

**Equation (1), p. 237, verbatim:**

```
        5000[A%] + 1125[B%] + 151.5[C%] + 26.5[D%]
MWD = ──────────────────────────────────────────────
                        100
```

> "where [A%] is the fraction > 2000 µm in wt %, [B%] is the fraction 250–2000 µm
> in wt %, [C%] is the fraction 53–250 µm in wt % and [D%] is the fraction < 53 µm
> in wt %."

**[MY ANALYSIS] Why the coefficients are what they are, and why this matters.**

Three of the four are simple class midpoints:

- `1125` = (250 + 2000)/2 ✓
- `151.5` = (53 + 250)/2 ✓
- `26.5` = 53/2 ✓
- `5000` = **(2000 + 8000)/2** ✓ — the upper bound is the **8-mm sieve** used
  during destructive sampling (§2.1). Everything > 2000 µm is assigned a nominal
  5000 µm regardless of its actual size.

**Implications for the model comparison:**

1. You must bin model aggregates into **exactly these four classes** and apply
   **exactly these four weights**. Your `population_outputs` call
   (`run_degryze.jl:269`) uses `sieve_sizes=[0.25, 0.5, 1.0, 2.0]` — **four**
   sieves and **five** classes, with different boundaries. That is a different
   statistic. Change to `[0.053, 0.25, 2.0]` and hard-code the weights.
2. **MWD saturates at 5000 µm.** No matter how large the model's aggregates grow,
   MWD cannot exceed 5000 µm. Conversely, a model aggregate of 2100 µm and one of
   7000 µm contribute *identically*.
3. **MWD is bounded below by 26.5 µm**, not 0.
4. The **< 53 µm class is dispersed clay**, deliberately washed out. If your model
   has no explicit silt+clay pool, [D%] is whatever mass is not in aggregates —
   and at day 0 that is most of the sample.
5. **Pre-wetting is by 5-min submersion, i.e. slaking is included.** This is a
   destructive, slaking-inclusive stability test, not capillary rewetting. If your
   WAS metric assumes a gentler disruption, model MWD will be biased high.

**Measured MWD, Table 2 (p. 240)** — day 1 and day 21, plus the field values.
Values followed by a different uppercase letter differ significantly at P < 0.05.

| Soil | Day 1 MWD /µm | Day 21 MWD /µm | Field MWD /µm |
|---|---|---|---|
| 1 | 284 B | 1006 B | 2238 |
| 2 | 210 C | 809 B | 3002 |
| 3 | 198 C | 1008 B | 3259 |
| 4 | 391 A | 1269 B | 2707 |
| 5 | 310 B | **2410 A** | 3905 |

✅ Your `degryze2006.csv` day-1 and day-21 rows reproduce these exactly
(0.288/0.216/0.199/0.398/0.309 and 1.007/0.810/1.014/1.273/2.406 mm). The
digitization is sound.

**Onset of significance (p. 239, verbatim):** "The day onward from which the MWD
was significantly larger compared with the start of the incubation was **day 16
for soil 1, day 10 for soil 2** (apart from the smaller value at day 13), **day 13
for soil 3, day 7 for soil 4, and day 4 for soil 5**."

**Large-macroaggregate formation rate, Table 3 (p. 242)** — linear regression of
% > 2000 µm against time, n = 16:

| Soil | Slope /% large macro day⁻¹ | R² | P |
|---|---|---|---|
| 1 | 0.57 ± 0.15 | 0.85 | < 0.0001 |
| 2 | 0.59 ± 0.09 | 0.93 | < 0.0001 |
| 3 | 0.83 ± 0.12 | 0.94 | < 0.0001 |
| 4 | 0.91 ± 0.20 | 0.91 | < 0.0001 |
| 5 | **2.02 ± 0.37** | 0.94 | < 0.0001 |

**[MY ANALYSIS]** This table is a better calibration target than MWD itself — it
is a single scalar per soil, it is linear (so easy to match), and it isolates the
> 2000 µm class where the 5000-µm weight dominates MWD. Fit to these five slopes
before fitting the MWD trajectories.

**The paper's own interpretation (p. 242, verbatim):** "the aggregate formation
rate was similar among the sandy loam, silt loam and loam soils (< 14 % clay),
whereas it was much larger in the clay loam soil (clay content 27 %) … **clay
content does not seem to affect aggregate formation in soils where the clay
content is below a threshold value of about 15 %.**"

---

## 4. Not stated — you must assume, and should say so

| # | Missing | Why it matters | Suggested default |
|---|---|---|---|
| 1 | **Particle density ρ_s** | needed for porosity → WFPS → θ | 2.65 g cm⁻³ (used in §2.4) |
| 2 | **Matric potential** | paper gives WFPS only | derive per soil from a published WRC, or run θ-driven |
| 3 | **Cylinder height** | domain geometry | ~5.4–6.0 cm (derived, §2.4) |
| 4 | **Headspace volume of canning jars** | CO₂ back-pressure (moot — flushed daily) | — |
| 5 | **Residue size distribution within 0.5–2 mm** | POM binning | uniform is the neutral choice; see §5.4 |
| 6 | **Initial microbial biomass** | model IC | not measured |
| 7 | **Full day-0 aggregate size distribution** | model IC | only MWD is reported |
| 8 | **Whether Table 1 organic C is pre- or post-POM-removal** | initial C partition | Table 1 caption says "crushed (< 250 µm) soils", so presumably post |
| 9 | **Unamended control respiration** | residue vs native CO₂ attribution | **none exists** — see §5.3 |
| 10 | **N amendment** | C:N 261 implies severe N limitation | none reported |
| 11 | **Balance of the residue fibre fractions** (sum = 72 %) | POM pool partition | unallocated 28 % |

---

## 5. Internal inconsistencies in the paper — flagging, not resolving

### 5.1 Sampling days don't match the figures

Methods (p. 237) say destructive analysis at **days 0, 1, 4, 7, 13 and 21** — six
dates. But the Results text (p. 239) refers to **day 10 and day 16**, Figure 4
plots roughly **eight** points per soil, and Table 3 regressions report **n = 16**.

**[MY ANALYSIS]** 8 dates × 2 replicates = 16 reconciles Table 3. And 14 samples
prepared = 7 post-day-0 dates × 2 replicates, with day 0 handled by the separate
6-replicate batch. So the actual schedule was most likely **0, 1, 2–3, 4, 7, 10,
13, 16, 21** with n = 2 per date. Your `degryze2006.csv` (days 0, 1, 2, 4, 7, 10,
13, 16, 21, with gaps) is consistent with this and is the better guide than the
Methods sentence.

### 5.2 Which paper is which — your README cites the wrong one

`paper/de_gryze/README.md` currently reads:

> "De Gryze, S., Six, J., Brits, C., & Merckx, R. (2006). *A quantification of
> short-term macroaggregate dynamics: influences of wheat residue input and
> texture.* European Journal of Soil Science, 57, 235–246."

This **conflates two different papers**:

- **EJSS 57, 235–246 (2006)** = *Water repellence and soil aggregate dynamics in
  a loamy grassland soil as affected by texture*, by **De Gryze, Jassogne,
  Bossuyt, Six & Merckx** — the paper this spec describes.
- *A quantification of short-term macro-aggregate dynamics: influences of wheat
  residue input and texture* = **De Gryze, Six, Brits & Merckx (2005), Soil
  Biology and Biochemistry** — a separate companion study, cited in the 2006
  reference list. *(The 2006 reference list prints the volume as 35, 55–66; SBB
  vol. 35 is 2003, so the volume number appears to be a typo in the source.)*

The README also says "5 soil textures (sandy loam → clay)" — there is **no clay
soil**; the range is sandy loam → clay **loam**.

### 5.3 The carbon budget does not close

Discussion, p. 241, verbatim: "we added only 2.6 mg carbohydrate-C g⁻¹ oven-dry
soil. **In both their study and ours, about 20 % of this C was respired.**"

**[MY ANALYSIS] That statement is not consistent with Figure 3.** 20 % of
2.6 mg C g⁻¹ is ~520 µg C g⁻¹. Measured day-21 cumulative respiration:

| Soil | Day-21 CO₂ /µg C g⁻¹ | as % of **total** residue C (4430) | as % of **carbohydrate**-C (2658) |
|---|---|---|---|
| 1 | 1592 | 36 % | 60 % |
| 2 | 1490 | 34 % | 56 % |
| 3 | 2139 | 48 % | 80 % |
| 4 | 1448 | 33 % | 55 % |
| 5 | 3157 | **71 %** | **119 %** |

Soil 5 respired **more C than the entire carbohydrate fraction contained**. Either
the "20 %" figure is wrong, or — far more likely — **Figure 3 is bulk soil
respiration including native SOC mineralization, with no unamended control
subtracted.** The Methods describe no unamended treatment ("day 0 (control)"
refers to the day-0 destructive sampling, not an unamended jar).

**This is the most consequential caveat for your calibration.** If you fit a
residue-driven respiration model to these curves, the model must either (a) also
carry native SOC mineralization, or (b) be fitted only to the *shape* and
*rank order*, not the absolute magnitude. Soil 5's 119 % makes the second option
safer. Note also that soil 5 has the highest native organic C (3.10 %) — exactly
the soil where a native-SOC contribution would be largest, which supports
interpretation (a).

---

## 6. Gap analysis against `paper/de_gryze/run_degryze.jl`

**Line numbers are checked against the file as of 2026-07-28.** They move on
every edit — treat the symbol, not the number, as the reference.

| # | Item | Code | Paper | Status |
|---|---|---|---|---|
| 1 | **Temperature** | `T_const = ic.T_0` = 298.15 K (L97), from `DEGRYZE_INCUBATION` | **25 °C** | ✅ resolved |
| 2 | **MWD sieve classes** | `sieve_sizes=DEGRYZE_SIEVES` = `[0.053, 0.25, 2.0]` (L265), nominals `DEGRYZE_CLASS_NOMINALS` | 3 sieves → 4 classes at **53 / 250 / 2000 µm**, fixed weights 5000/1125/151.5/26.5 | ✅ resolved |
| 3 | **Soil texture** | `degryze_soil(SOIL_ID)` (L84) sets `f_clay_silt` from Table 1 | per soil 0.47/0.67/0.56/0.56/0.79 | ✅ resolved — the phantom 0.74 soil is gone |
| 4 | **Bulk density** | `ρ_bulk = soil.ρ_b` (L125) | per-soil 1280–1420 | ✅ resolved |
| 5 | **Water content** | `ψ_const = ic.ψ_0` (L98), derived by `van_genuchten_inverse` from 60 % WFPS | **60 % WFPS**; θ = 0.279–0.310 by soil | ✅ resolved |
| 6 | **Duration** | `t_max = 45.0` (L113) | **21 days** | ✅ over-run is deliberate; the READMEs no longer claim 28 or 60 |
| 7 | **Residue C fraction** | `f_C_POM = 0.443` (L271), consistent with `I_input = 4.43e-3` | **0.443** | ✅ resolved |
| 8 | **POM size distribution** | `Normal(1.25, 0.23)` over 5 bins (L59–66) | only "**0.5–2 mm**" — no distribution given | 🟡 **open** — unconstrained assumption, recorded in the de_gryze README |
| 9 | **POM bin count** | 5 bins (L62) | — | ✅ READMEs corrected from 23 to 5 |
| 10 | **Residue chemistry** | unused | 11 % soluble / 3 % hemi / 46 % cellulose / 11 % lignin available | 🟡 **open** — free constraint left on the table |
| 11 | ~~**`s_M` comment**~~ | — | — | ✅ **moot 2026-07-29** — `s_M` deleted; MAOC saturation is now an output. REFERENCE §5d, §26 erratum 14 |
| 12 | **Citation** | `paper/de_gryze/README.md` | see §5.2 | ✅ resolved |
| 13 | **O₂ = 21 % constant** | `DEGRYZE_INCUBATION.O2_frac` (L99) | ✅ justified by daily air flushing | ✅ correct |
| 14 | **Output cadence 3 h** | `dt_output = 0.125` (L114) | ✅ matches GC cadence | ✅ correct |
| 15 | **Initial SOC** | `degryze_ic` takes Table 1 SOC per soil | 1.78–3.10 % | ✅ resolved — no longer a five-soil mean |
| 16 | **POM removed at t=0** | matches pretreatment | ✅ | ✅ correct |
| 17 | **`D_background`, `bg_class_fractions`** | **removed** | were set by inverting the measured day-0 MWD | ✅ scrubbed — an observable must not re-enter as an input (CLAUDE.md rule 3) |

### 6.0 Open after this pass

Rows 8 and 10 above, plus, from outside this table:

- **`pct_gt2000um` is a step function.** `r_agg` snaps by day 2 and then
  freezes. Until this is understood, the shape comparison against Table 3 is
  not meaningful.
- **CO₂ overshoots ~1.7×** at day 21 for soil 3 (3595.8 vs 2139 µg-C/g). Both
  sides are total respiration — see A3 and §5.3 — so the discrepancy is real,
  not a partitioning artefact.
- **Group B parameters** sit off their cited anchors by 3.3× to 1500×
  (`docs/REFERENCE.md` §5a). Fitting over Group C while Group B is off its
  anchors will absorb the Group B error invisibly. Fix Group B first.
- ~~**`optimize_soil3.jl` has not been migrated** to `degryze_soils.jl`~~ —
  settled 2026-07-30 by archiving it rather than migrating it
  (`paper/_archive/degryze_tooling_20260730/`). It hard-coded `ρ_bulk = 1300`,
  `SOC = 0.0221`, `ψ = −29`, `f_clay_silt = 0.74` for soil 3 where
  `degryze_soil(3)` gives 1370 and `degryze_ic` derives ψ, and its fitted
  parameters were separately invalid (fitted against the broken POM
  normalization). Calibration is a separate future project: hand-tune to close,
  then one reusable fitting routine for all examples, tested against De Gryze.

### 6.1 [MY ANALYSIS] Suggested priority order

1. **Fix the MWD estimator** (#2). Until model and data compute the same
   statistic, nothing else can be diagnosed. Bin at 53/250/2000 µm, cap the top
   class at a nominal 5000 µm, apply eq. (1) verbatim.
2. **Set T = 298.15 K** (#1). Your open note at `run_degryze.jl:181–186` records a
   "~1.9× CO₂ overshoot" — worth checking whether raising T to 25 °C makes that
   *worse* before you spend more on `R_P_max`. It will: correcting T pushes
   respiration up ~40 %.
3. **Decide the CO₂ attribution question** (§5.3) before any further respiration
   calibration. If the data include native SOC, a 1.9× "overshoot" against a
   residue-only model may be partly real signal, and tuning `R_P_max` down to
   match would be fitting the wrong thing.
4. **Go per-soil** (#3, #4, #5). One parameter set cannot span soils 1–4 *and*
   soil 5 — the paper's own conclusion (§3.2, p. 242) is that soil 5 crosses a
   ~15 % clay threshold and behaves differently. Calibrate to the soils 1–4
   cluster; treat soil 5 as a validation case.
5. **Use Table 3 slopes as the primary MWD target** (§3.2) rather than the full
   trajectories — five scalars, all with R² ≥ 0.85.
6. Then the minor items (#7–#12).

---

## 7. Companion paper — OBTAINED, see `degryze_SBB_2005_spec.md`

> **Update 2026-07-27:** this paper has since been added to `references/` as
> `De Gryze SBB 2005.pdf` and is fully extracted in **`degryze_SBB_2005_spec.md`**.
> Two results there supersede open questions in this document:
> - **§5.3 above (the carbon budget) is resolved.** SBB 2005 has an unamended
>   control; only **~36–46 % of cumulative CO₂ at this loading is residue-derived**.
> - The correct citation is **Soil Biology & Biochemistry 37, 55–66** (not 35 —
>   verified from the PDF's own running header).
>
> The original note is kept below for the record.

**De Gryze, S., Six, J., Brits, C. & Merckx, R. (2005).** *A quantification of
short-term macro-aggregate dynamics: influences of wheat residue input and
texture.* **Soil Biology and Biochemistry** *(the 2006 reference list prints
"35, 55–66"; the volume number looks like a typo)*.

The 2006 paper defers to it **twice** for things you need:

- p. 236: "Soil sampling and preparation are **described in detail in De Gryze
  *et al.* (2005)**."
- p. 237: the 150 g / 1.5 g mixing protocol is cited to it.

It is **not in `references/`**. Given that it is the same soils, same residue,
and explicitly a residue-input × texture aggregate-dynamics study, it is the most
likely source for the day-0 aggregate size distribution, a fuller parameter
table, and possibly an unamended control. *(That it contains these is my
inference from how it is cited — I have not read it.)*
