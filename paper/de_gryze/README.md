# De Gryze et al. (2006) — Aggregate Formation Incubation

**Reference:** De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R.
(2006). *Water repellence and soil aggregate dynamics in a loamy grassland soil
as affected by texture.* **European Journal of Soil Science** 57, 235–246.

> Earlier versions of this file cited *"A quantification of short-term
> macroaggregate dynamics: influences of wheat residue input and texture"* for
> this volume. That is a **different paper** — De Gryze, Six, Brits & Merckx,
> *Soil Biology and Biochemistry*, the companion study specified in
> `degryze_SBB_2005_spec.md`. Our data are from the EJSS 2006 paper above.

**The authoritative note for this directory is
[`degryze_EJSS_2006_spec.md`](degryze_EJSS_2006_spec.md).** Read §0 and §0a
before changing anything here or reporting any number derived from it. This
README is an index, not a specification.

## Experiment (p. 236–237)

- **21-day** incubation at **25 °C**, wheat-residue amendment
- **5 soils**, sandy loam → clay **loam** (not clay; clay ranges 7–27 %)
- 1.5 g wheat stems per 150 g soil, 44.3 % C → 4430 µg-C per g soil, cut to
  0.5–2 mm
- Water held at **60 % WFPS**; headspace flushed daily, so O₂ ≈ ambient
- Observables: wet-sieved size distribution (3 sieves at 53 / 250 / 2000 µm,
  reported as MWD by eq. 1) and cumulative CO₂

Sampling dates in the Methods do not match the figures or the regression `n`;
see spec §5.1. `degryze2006.csv` is the better guide.

## Model setup — as the code actually stands

| | Value | Source |
|---|---|---|
| POM size classes | **5**, midpoints 0.65–1.85 mm | `Normal(1.25, 0.23)` over 0.5:0.3:2.0 mm bins — an **assumption**, the paper gives only the 0.5–2 mm range (spec §6 item 8) |
| Soil | `SOIL_ID = 3` (loam) | `degryze_soil(3)`, `degryze_ic(3, soil)` |
| Temperature | 298.15 K | measured, p. 237 |
| Water potential | derived per soil from 60 % WFPS | `van_genuchten_inverse`, −28.4 to −29.3 kPa across the five |
| Duration | 45 days | over-run past the 21-day data window |
| Sieves | `[0.053, 0.25, 2.0]` mm | measured, p. 237 |
| Comparison target | `pct_gt2000um` vs Table 3 | the only class column free of assumption A1 |

Each POM size class is one `run_aggregate()` call; the population is assembled
by `domain_tessellation` / `pom_population` (`src/physics/tessellation.jl`) and
summarized by `population_statistics` (`src/postprocessing/population.jl`).

## Files

| File | Description |
|---|---|
| `degryze_EJSS_2006_spec.md` | **Authoritative.** What the paper says, what it does not say, what we assume, and the gap analysis |
| `degryze_SBB_2005_spec.md` | Companion 2005 paper (has the unamended control the 2006 paper lacks) |
| `degryze_soils.jl` | The five soils as **overrides** on package defaults. Defines no types. Observations kept separate from inputs |
| `run_degryze.jl` | Forward run for one soil |
| `optimize_soil3.jl` | Parameter fit for soil 3. **Outputs currently invalid** — fitted against the broken POM normalization, and not yet migrated to `degryze_soils.jl` |
| `postprocess_dataframe.jl` | DataFrame packaging only. No physics — see its header |
| `degryze2006.csv` | MWD data (µm) for 5 soils |
| `degryze_CO2_2006.csv` | Cumulative CO₂ data for 5 soils |
| `degryze_*_plots.R` | R plotting of the CSV outputs |
| `output/` | Generated outputs (git-ignored) |

## MATLAB correspondence

| MATLAB | Julia |
|---|---|
| `de_gryze_test.m` | `run_degryze.jl` |
| `param.m` | `BiologicalProperties()` + `SoilProperties()` defaults |
| `single_aggregate_beta.m` | `run_aggregate()` |

## Status

- [x] Forward simulation script
- [x] Soils driven from measured Table 1 values, not a phantom texture
- [x] Sieve classes and MWD estimator matched to eq. (1)
- [x] POM normalization corrected (input was understated 28.8× by an ω division)
- [x] Population statistics moved into the package and tested
- [ ] **MWD / `pct_gt2000um` is a step function** — `r_agg` snaps by day 2 and
      then freezes. Open.
- [ ] **CO₂ overshoots** — 3595.8 vs 2139 µg-C/g at day 21 for soil 3 (~1.7×).
      Note the measured curve is bulk respiration with no unamended control
      (spec §5.3), so this comparison is total-to-total by construction.
- [ ] Group B parameters returned to their cited anchors (REFERENCE.md §5a)
      **before** any fitting
- [ ] Calibration across all five soils
- [ ] Figures for manuscript

## What this directory may not do

Per `CLAUDE.md` rule 4: anything a second problem could use belongs in `src/`.
Files here override defaults, hold measured data, make the call, and plot. They
define no types and no reusable logic. If something here would be useful to
another experiment, it is in the wrong place.
