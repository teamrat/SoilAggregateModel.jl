# Structural audit — 2026-07-27

Findings from a systematic pass over the repository after a long working session.
Grouped by cause. Every item carries a file:line reference and states whether it
was introduced in this session or inherited.

Sections A–C and E–F are mine. Section D predates this session. Section G is the
underlying reason A–C keep recurring.

---

## A. Things placed in the wrong file

### A1. `paper/de_gryze/soils.jl` — invented a parallel type *(mine)*

`soils.jl:16` defines `struct DeGryzeSoil` with fields

```
id, name, SOC, sand, silt, clay, BD, MWD_day0, MWD_day21,
slope_gt2000, f_clay_silt, θ_s, θ_60WFPS, ψ
```

Three problems, compounding:

1. **It duplicates existing model types.** `SoilProperties` (`src/parameters.jl:215`)
   already has `θ_s`, `ρ_b`, `f_clay_silt`. `InitialConditions`
   (`src/physics/initial_conditions.jl:46`) already has `SOC`. The model cannot
   consume `DeGryzeSoil`, so any caller must unpack it and rebuild the real
   types by hand — which is the same pattern that put model logic in scripts.
2. **It welds inputs to observations.** `SOC`, `BD`, `f_clay_silt`, `θ_s`, `ψ`
   are model inputs. `MWD_day0`, `MWD_day21`, `slope_gt2000` are data to compare
   *against*. One type carrying both makes it impossible to pass the inputs
   without also passing the answers.
3. **Generic name, specific content.** The only struct in the repository defined
   outside `src/`.

**Correct shape:** no new type. Functions returning the existing
`SoilProperties` / `InitialConditions` with only the fields that differ from
default overridden; observations as a data table beside `degryze2006.csv`.

### A2. Mineral size distribution added to the run script *(mine)*

`run_degryze.jl` now computes `mineral_classes` and `class_nominals` inline. The
log-linear interpolation is a general method and belongs in `src/`; only the
sieve series and De Gryze's eq. (1) nominals are paper-specific.

### A3. Repository-root clutter *(mine)*

| Path | Status |
|---|---|
| `_backup_20260727/` | 11 pre-patch source copies. Gitignored. **Restoring from it silently reverts every fix made this session.** |
| `degryze_run.txt`, `test_output.txt`, `diag_dfm.txt` | run logs at repo root, untracked |
| `scripts/diag_dfm.jl` | one-off diagnostic; its question is answered (D_Fm was fine) |

---

## B. Contradictions created inside a single document

### B1. `docs/REFERENCE.md` now defines α and β two opposite ways *(mine)*

After the α/β swap in `src/parameters.jl`, §3.3 and §9 were not updated:

| Location | States |
|---|---|
| §3.3 table | `α_i` = **Mobilization** rate; `β_i` = **Immobilization** rate |
| §9 formula | `net_i = η·(β_i·Π − α_i·Π^δ)·F_i` — exponent on the **loss** term |
| §5a (added this session) | `α` = **immobilization** (gain, carries the exponent) |

**One file, three sections, two incompatible definitions.** §5a was written to be
the provenance record and contradicts the parameter table twelve sections above it.

### B2. §9 describes code deleted in February *(inherited, worsened)*

```
Insulation: trans_ni = ζ · F_n (one-way, F_n → F_i)
S_Fi gets: +net_i + trans_ni
S_Fn gets: +net_n − trans_ni
```

This is the independent-drain formulation removed in February when ζ became a
splitting fraction. `src/solver/reactions.jl` has `S_Fn = trans.immobil_n`.

### B3. §3.3 values that do not match the code *(inherited)*

| Symbol | §3.3 | `parameters.jl` |
|---|---|---|
| `δ` | 2.0 | 1.0 |
| `ε_F` | 1e-10 | 1e-4 |
| `ζ` | "day⁻¹, Insulation rate" | dimensionless splitting fraction |
| `D_Fn0` | *(blank)* | 0.01 |

---

## C. Divergence created between sibling scripts *(mine)*

`paper/de_gryze/optimize_soil3.jl` was never updated alongside `run_degryze.jl`:

| | `run_degryze.jl` | `optimize_soil3.jl` |
|---|---|---|
| Temperature | 298.15 K (25 °C) | **293.15 K (20 °C)** (`:52`) |
| Tessellation | lines 110–135 | **copy-pasted at lines 76–88** |
| `f_C_POM` | 0.443 | not passed |

Same experiment, two temperatures. The ω derivation exists in two scripts, so any
correction must be made twice or they drift further. `optimize_soil3.jl:267` calls
`population_outputs(...; sieve_sizes=Float64[])`, which after this session's
changes emits a single meaningless class column.

The optimizer's outputs are already invalid for a separate reason (fitted against
the broken POM normalization), so this compounds rather than causes.

---

## D. Inherited structural problems

### D1. `ρ_POM` has two values 95× apart

| Location | Value |
|---|---|
| `src/parameters.jl:183` | 200.0 |
| `paper/de_gryze/run_degryze.jl:109`, `optimize_soil3.jl:73` | 200.0 |
| `paper/simulations/common.jl:32` | **19098.6** |

Same name, same units, different physical claim.

### D2. 646 lines of model logic outside the package — **resolved 2026-07-28**

`paper/de_gryze/postprocess_dataframe.jl` defined `result_to_dataframe`,
`spatial_snapshots`, `run_diameter_sweep`, `population_outputs`. None are
De Gryze-specific. Reached by `include`, not `import`, from two scripts.

**What moved.** The physics inside `population_outputs` — aggregate mass, sieve
binning, the fixed-weight MWD, the population sums — is now
`src/postprocessing/population.jl` (`aggregate_mass`, `sieve_class`,
`population_statistics`), operating on plain arrays. Documented in
`docs/REFERENCE.md` §5c; tested in `test/test_population.jl`.

**What stayed, and why.** `result_to_dataframe`, `spatial_snapshots`,
`run_diameter_sweep` and a thin `population_outputs` wrapper remain in
`paper/de_gryze/`. They define no quantity: they tabulate what
`integrated_pools`, `co2_flux`, `radial_profiles` and `population_statistics`
already compute. Moving them would put DataFrames **and CSV** into
`Project.toml` — which today lists only Distributions and Plots — to buy
packaging, not physics. The scripts currently work only because `Pkg.activate`
leaves the default environment on `LOAD_PATH`; package code cannot rely on
that, so the dependency would have to be real.

**Residual.** `run_diameter_sweep` still recomputes `P_0 = (4/3)πr₀³ρ_POM`
rather than taking `pom_population`'s `P_0_per_particle`. Same formula, same
ρ_POM, so they agree — but it is a second definition of one number.

### D3. Stale claims in tracked documentation — **resolved 2026-07-28**

| File | Claim | Actual | State |
|---|---|---|---|
| `README.md:4,123` | "524/524 tests passing" | 1363 pass, 1 deliberate fail | ✅ |
| `README.md` project tree | missing `postprocessing/`, `tessellation.jl`, `initial_conditions.jl`, `de_gryze/`; listed a `scripts/diagnostics_30day.jl` that does not exist | — | ✅ rewritten |
| `paper/de_gryze/README.md:14` | 23 POM classes | 5 | ✅ |
| `paper/de_gryze/README.md:7,16` | 28-day / 60-day | 45-day run; paper says 21 | ✅ |
| `paper/de_gryze/README.md:8` | "sandy loam → clay" | sandy loam → clay **loam** | ✅ |
| `paper/de_gryze/README.md:3` | wrong paper title and author list | conflated the 2005 SBB companion | ✅ |
| `docs/REFERENCE.md:4` | Authoritative physics: `manuscript_updated.tex` | file does not exist | ✅ (fixed earlier this session) |
| `manuscript_changes.md:222` | 1289 pass / 6 fail / 3 error | superseded | ✅ marked superseded, not silently edited |

---

## E. Dead code left this session *(mine)*

| Location | Issue |
|---|---|
| `postprocess_dataframe.jl:539` | `D_ret = [D_agg[i,t] for i in 1:n_sizes]` — a pointless copy after the `max()` was removed |
| `postprocess_dataframe.jl:636` | `MWD_fixed_weight_mm` computed, consumed by nothing |
| `postprocess_dataframe.jl:613,642` | `WAS_*` columns now redundant with the `pct_*` class columns; a comment says so, both were kept |
| `postprocess_dataframe.jl` | `f_agg_vol`, `POM_mass_frac` emitted, consumed by nothing |

---

## F. Spec drift created this session *(mine)* — **resolved 2026-07-28**

`paper/de_gryze/degryze_EJSS_2006_spec.md` §6 cited line numbers that had all
shifted (L78, L79, L111, L198–204, L269), and rows 1 (temperature), 2 (sieve
classes) and 7 (`f_C_POM`) were resolved but still marked 🔴. §0 item 2 warned
about the 20 °C setting that had since been corrected.

§6 is re-checked against the file and re-scored, and a new §6.0 lists what is
actually still open. §0 items 1–3 now carry the fix next to the warning rather
than being deleted, so the reason survives with the result. The §6 header now
states that line numbers move on every edit and the symbol is the reference.

---

## G. Why A–C keep recurring

The single-aggregate model and the population upscaling are entangled across the
`src/` boundary:

- **ω is derived in a paper script** but changes initial conditions inside
  `src/physics/initial_conditions.jl` (all pools divided by ω, POM deliberately
  not). The decision lives in `src/`; the number that drives it lives in a script.
- **"What is an aggregate" is defined in two places.** `r_agg` in
  `src/postprocessing/aggregate_radius.jl`; the aggregate's mass, its POM core,
  and its sieve class in `paper/de_gryze/postprocess_dataframe.jl`.
- **The population construction** (`N_POM`, `P_0_per_particle`, packing) is
  script-level, but nothing about it is experiment-specific.

Because the nearest file containing similar logic was already in the wrong place,
each local edit landed there too. The fix is not to move the edits — it is to move
the boundary.

---

## Remediation order

1. **`docs/REFERENCE.md` §3.3 + §9** — internal contradiction, documentation only,
   no compile risk. Do first.
2. **Delete `_backup_20260727/`, root logs, `scripts/diag_dfm.jl`** — removes a
   trap that can silently revert fixes.
3. **`soils.jl`** — drop the struct, return the existing model types.
4. **`src/physics/tessellation.jl`** — one home for ω, packing, `N_POM`; both
   scripts call it.
5. **`postprocess_dataframe.jl` → `src/`** — ✅ done 2026-07-28, narrowed: the
   physics moved, the DataFrame packaging stayed. See D2 above for why.
6. **Stale README / spec line references** — ✅ done 2026-07-28, after 3–5.

Items 1–2 carry no compile risk. Items 3–5 cannot be verified without running
Julia.

---

## Status — 2026-07-28

| Item | State |
|---|---|
| 1. REFERENCE.md §3.3 + §9 | ✅ |
| 2. `_backup_20260727/`, root logs, `diag_dfm.jl` | moved to `_to_delete/`; **`rm -rf _to_delete` still needs running on the Mac** |
| 3. `soils.jl` → `degryze_soils.jl` | ✅ |
| 4. `src/physics/tessellation.jl` | ✅ |
| 5. population physics → `src/postprocessing/population.jl` | ✅ |
| 6. stale README / spec line references | ✅ |

Also still open, unchanged by this pass:

- **D1** — `ρ_POM` is 200.0 in `src/parameters.jl` and the de_gryze scripts but
  19098.6 in `paper/simulations/common.jl`. Same name, same units, 95× apart.
- **C** — `optimize_soil3.jl` has not been migrated to `degryze_soils.jl`: it
  still hard-codes `ρ_bulk = 1300`, `SOC = 0.0221`, `ψ = −29`,
  `f_clay_silt = 0.74` for soil 3, where `degryze_soil(3)` gives 1370 and
  `degryze_ic` derives ψ per soil. Its fitted outputs are invalid for the
  separate reason recorded above.
- **E** — `f_agg_vol` and `POM_mass_frac` are emitted and consumed by nothing.
  Kept deliberately: `f_agg_vol` is the >0.6 interpenetration check and
  `MWD_fixed_weight_mm` is the eq. (1) comparison quantity.
