# De Gryze et al. (2006) — aggregate formation incubation

De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R. (2006). *Water
repellence and soil aggregate dynamics in a loamy grassland soil as affected by
texture.* **European Journal of Soil Science** 57, 235–246.
PDF: `references/De Gryze et al EJSS 2006.pdf`.

> An earlier version of this file cited *"A quantification of short-term
> macroaggregate dynamics…"* for this volume. That is a **different paper** —
> De Gryze, Six, Brits & Merckx, *Soil Biology and Biochemistry*, the companion
> study specified in `degryze_SBB_2005_spec.md`.

This README is an index. **The authoritative note is
[`degryze_EJSS_2006_spec.md`](degryze_EJSS_2006_spec.md)** — read §0 and §0a
before changing anything here or quoting any number derived from it.

---

## What this problem is for

Five soils from one meadow, differing in texture and almost nothing else, given
identical residue and identical incubation. That isolates texture, which makes
it a **test** of the model's texture mechanism rather than a fitting target that
happens to have five curves. Texture enters the model through one length, the
Sauter mean diameter `d₃₂` (REFERENCE.md §4.4); it is not fitted per soil.

The soils were sieved to < 250 µm and stripped of native POM before incubation,
so every macroaggregate measured was formed during the 21 days.

---

## Layout

```
degryze_config.jl           THE configuration — soil, IC, parameters, POM
                            distribution, environment, tessellation
run_degryze.jl              the driver — run this
degryze_soils.jl            measured inputs per soil + DEGRYZE_OBSERVED
postprocess_dataframe.jl    DataFrame layer; kept outside src/ so the package
                            does not depend on DataFrames/CSV

data/                       measured data + its own README
output/                     generated, gitignored
                            baseline_split_20260728/ = pre-stiff-solver reference

degryze_EJSS_2006_spec.md   authoritative: assumptions, and what each is worth
degryze_SBB_2005_spec.md    companion 2005 study
```

Model physics lives in `src/`. Files here supply **only** what De Gryze
measures, plus the run choices — never model logic (CLAUDE.md §4).

**One configuration.** `degryze_config.jl` is included by `run_degryze.jl`, and
is currently its only consumer. Three copies existed on 2026-07-29 and had
diverged — the forward run and the objective were scoring different models — and
five further consumers kept their own forks until 2026-07-30. Change a parameter
in `degryze_config.jl` or nowhere.

**What is here is what the run needs.** On 2026-07-30 everything outside `src/`
that did not contribute to this run was archived to
`paper/_archive/degryze_tooling_20260730/`: the parameter fit, `fitting/`,
`diagnostics/`, and the superseded R plotting. That archive has its own README
saying what each was for and why none of it should be quoted. Nothing was
patched into currency — a script whose premise had changed is cheaper to rewrite
than to correct.

---

## Running it

```julia
julia --project=. paper/de_gryze/run_degryze.jl
```

Sweeps five POM diameter classes (0.65–1.85 mm) for one soil, set by `SOIL_ID`
in `degryze_config.jl`. Writes:

| output | what |
|---|---|
| `output/summary.csv` | diameter × time → `r_agg_mm`, pools, CO₂ |
| `output/spatial_profiles.csv` | diameter × time × node → all 8 fields, binder, thresholds |
| `output/population.csv` | population MWD, CO₂, size classes |
| `output/degryze_model_vs_data.png` | 5 panels, model against measurements |
| `output/degryze_profiles_d<D>mm.png` | 9 panels per POM class, every field against radius |

The CSVs are the deliverable — publication figures are built in R from them.
The PNGs are diagnostic, produced from the same single run, so a figure and a
score can never describe different parameter sets.

There is one solver, `run_aggregate_stiff`. The Strang/Crank–Nicolson reference
was archived 2026-07-30 — `_archive/split_solver_20260730/README.md` — so
`SOLVER_USED` and the `solver=` argument are gone.

Cost, measured 07-30 on this configuration: 1.16 s per aggregate at `n_grid=200`
over 22 days (2341 accepted steps, 39 Jacobians), so ~12 s for all ten classes.
The sparse linear solve is 33 % of it. `docs/_archive/AUDIT_20260729.md` §B has the breakdown.

### Fitting — a separate future project

Not done here, and no fitting code lives in this folder. The plan: get the fit
visually close by hand first, then build **one** reusable fitting routine,
central enough for every example to use, and test it against De Gryze. The
harness that used to sit in `fitting/` is archived, and `setup.jl` there is the
piece worth reading when that project starts — `run_model`, `with_field` with an
identity self-check, `at_time`, and the data readers are all reusable.

The objective it was built around, kept here because it is a design decision and
not a description of a file:

The objective is

```
L = (1/N_MWD) Σ ln(MWD_mod/MWD_obs)²  +  (1/N_CO2) Σ ln(CO2_mod/CO2_obs)²
```

Log residuals, so the two observables need no unit conversion to be comparable;
weighted per observable, not per point, so 7 MWD points cannot outvote 4 CO₂
points by count. MWD day 0 is excluded (no fitted parameter reaches it); CO₂ day
0 is excluded (identically zero on both sides).

---

## What the model is compared against

| observable | source | where |
|---|---|---|
| MWD time series, 9 dates × 5 soils, mm | digitised | `data/degryze2006.csv` |
| Cumulative CO₂, 5 dates × 5 soils, µg-C/g | digitised from Fig. 3 | `data/degryze_CO2_2006.csv` |
| MWD day 1 / 21, formation rates, R² | Tables 2 and 3 | `degryze_soils.jl` — `DEGRYZE_OBSERVED` |

`data/README.md` describes each and its caveats. Two to know before plotting
anything: the Table 3 **formation rates are fitted slopes, not measurements** —
never draw them as data; and the CO₂ figures are **total** soil respiration with
no unamended control, so no partition factor is applied on either side.

---

## Status

Reproduces the qualitative behaviour. Quantitatively it overshoots, in one
direction, on both observables — soil 3 MWD 334 µm at day 1 against 198
measured and 1389 at day 21 against 1014; CO₂ 3600 against 2139 µg-C/g at day
21. That the two overshoots share a direction and rough magnitude suggests one
cause rather than two.

The overshoot is not the main problem. The **shape** is: the model peaks on day
4 and declines, where the data rises to day 21. Measured, not assumed —
`R_P_max` moves the level and leaves the day-4/day-21 ratio at 2.8–2.9 over a
factor of 17, and multiplying all three binder-precursor diffusivities by 1000
moves L from 0.2000 to 0.1998. The remaining candidate was that `G_c` did not
depend on aggregate radius, so nothing graded growth.

Three changes are in place and **none has been run yet**, so every number above
still describes the model before them.

1. Five configuration values corrected, documented with their justification in
   `degryze_config.jl`: `f_insulated` 0.5→0.0, `f_bact` 0.01→1e-4, `f_eps`
   0.005→5e-4, `k_d_eq` 0.05→0.005, `κ_d_ref` 0.001→1e-4.
2. `w_E` 0.5→0.25, halving EPS's weight in the binder `F_i + w_E·E` so that
   `F_i` carries proportionally more. `w_E` is fitted (REFERENCE.md §4.4).
   `κ_b` is not overridden here; the package default moved 0.0143869→0.012.
3. The mineral sand split is **inverted from the measured day-0 MWD** per soil
   (spec §0a A1′), replacing the log-linear guess. Soil 3 goes 0.573→0.253 of
   sand in 250–2000 µm, dropping model day-0 MWD 335→199 µm. Soil 5 admits no
   solution — its day-0 MWD exceeds what its texture can produce without
   aggregates — so it is set to 0.9 and the −112 µm residual is reported.
4. `p_Gc` 0→1, making the threshold rise with radius:
   `G_c(r) = G_c(δ_s)·(r/δ_s)^p_Gc`. This is a **model change**, in
   `src/postprocessing/aggregate_radius.jl`; `p_Gc = 0` is the package default
   so nothing outside this problem moved. REFERENCE.md §4.4a has the two
   mechanisms that motivate it, the statement that the exponent is fitted, and
   why `p_Gc < 0` would reinstate a form already rejected once.

`w_E` and `p_Gc` both reach `r_agg` and nothing else — `r_agg` does not feed back
into the state equations — so CO₂ and every pool profile are unchanged by them.
Only MWD and the size classes move. `p_Gc = 0` recovers the pre-2026-07-29
threshold.

(3) adds no degrees of freedom to the fit: one measurement per soil determines
one unknown per soil, and no fitted parameter reaches day 0. It does spend the
day-0 point, which was already excluded from the loss.

Not yet fitted. Blockers are in `docs/BACKLOG.md`: `M_max` has two values 29×
apart, `delta` and `beta` are undecided, and five Group B parameters sit 1000×,
1500× and 20× off their cited anchors. Fitting free parameters while anchored
ones are wrong absorbs the error invisibly, which is why the backlog puts
anchors before fit.

`M_max` is settled since (`maoc_capacity`, REFERENCE §26 erratum 12), as is the
initial condition (erratum 14, one free parameter removed). The Group B anchors
remain.

The archived `optimize_soil3.jl` hard-coded soil 3's properties *and* the
measured CO₂ series, and its stored fit was made against the broken POM
normalisation. Not a starting point — do not restart from it.
