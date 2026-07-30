# Architecture and housekeeping audit

**Living document. Update items in place; do not open a second audit.**
Started 2026-07-29. Dated audits in `docs/` (`STRUCTURAL_AUDIT_2026-07-27.md`,
`TEST_AUDIT_2026-07-28.md`) are frozen historical records — this one is not.

> **This file overlaps `BACKLOG.md` and one should absorb the other.** Both are
> lists of open work with status. Keeping two is the same proliferation this
> audit criticises in §E. Not resolved unilaterally — see *Open decisions*.

## Status log

| Date | Item | Outcome |
|---|---|---|
| 07-29 | **A2 `M_max`** | **CLOSED.** Field deleted; `maoc_capacity(soil)` is the only definition. `k_ma` redefined dimensionless and anchored to Georgiou et al. (2022): 0.086 HM / 0.048 LM. The old 0.48 was 48 mg/g with the mg→g conversion dropped — 10× high. De Gryze soils are Belgian loess-derived, so all low-activity. REFERENCE §26 erratum 12; BACKLOG item 3 closed. |
| 07-29 | **A1 O₂/Henry in `paper/simulations/`** | **DISPOSITIONED, not repaired.** Owner: the folder is to be rebuilt from scratch once De Gryze yields a trusted parameter set — nothing produced before then is trustable, and the folder structure is duplicative. Repairing the block would be work thrown away. **Do not reuse any output from `paper/simulations/` in the meantime.** |
| 07-29 | **`s_M` / initial condition** | **CLOSED, one parameter removed.** `partition_CM` placed MAOC at `s_M·M_max·f(βC)` while the solver drives toward `M_max·f(βC)` — two isotherms, so the run opened 40 % below equilibrium and sorbed to close a gap nothing created. `s_M` also meant `M/M_eq_full` while documented as `M/M_max`, so 0.6 delivered 0.485, and the `@info` reported the former and looked right. Deleted from `partition_CM`, `InitialConditions` and all call sites; saturation is now an output bounded by `SOC_residual/M_max`, i.e. predicted from measured SOC and texture. `k_L` 1000 → 25000, set by pore-water DOC landing in the observed 10–100 mg/L (48 now, 573 before) — the first constraint ever placed on it. `R` → `SOC_residual`. REFERENCE §5d + §26 erratum 14; manuscript §sec:initial_condition. |
| 07-29 | **ω overlap convention** (found tracing A2 downstream) | **OPEN — see the standing item below.** Analysis committed (REFERENCE §26 erratum 13), repair not. Code unchanged: default is still `f_domain_min = 10`, `ω ≈ 30`. `domain_tessellation` records the geometry with an `@info`, not a warning — overlap is the intended design, so a standing warning on it would be wallpaper. |
| 07-29 | **A1 `paper/simulations/`** | **ARCHIVED** to `paper/_archive/simulations_20260729/`. Likely deletion later. Its O₂/Henry and λ defects go with it unrepaired. |
| 07-29 | **A3 POM dissolution** | **CLOSED — my finding was wrong.** The manuscript already had the code's additive ½(a+b) form; I checked three internal documents and not the physics spec. REFERENCE §12, ARCHITECTURE §4.2 and GUIDE §8 corrected, with the mechanistic argument for additive stated. |
| 07-29 | **A6 stale quantisation warning** | **CLOSED.** Deleted from `fitting/sensitivity.jl`. It told the reader not to use gradient-based optimisation because `r_agg` snaps to grid nodes; the sub-grid interpolation removed that and `test_postprocessing.jl` asserts the opposite. |
| 07-29 | **binder vs `G_c`** | **CLOSED as documentation, no code change.** Not a missing `1/ω`: the binder sums a diluted background contribution and an undiluted residue-derived one, so no single factor applies, and the residue-derived term dominates once decomposition starts. Recorded as a property of the closure — REFERENCE §4.4, `aggregate_radius.jl`, manuscript `r_agg` fifth property. `κ_b` carries the residual and is fitted. |
| 07-30 | **B pass 3 — conservation weight and system carbon** | **CLOSED.** `conservation_weight`/`conservation_weights` in `types.jl` own `4πr²h`; `GridInfo`, `reaction_step.jl` and `api.jl` all call it. `compute_total_carbon` is one loop with `W` / `GridInfo` / `(r_grid, h)` entry points, so the four production call sites pass `grid` and allocate nothing. `total_system_carbon(pools, k)` owns the pools-level total; `carbon_balance_table` gained an `IntegratedPools` method and `result_to_dataframe` now calls it instead of re-deriving. Bitwise result-preserving. New tests pin the stencil property `W_i/(r_i²h²) = 4π/h`, the state-level total against the pools-level one, and the three `compute_total_carbon` entry points against each other. |
| 07-30 | **B pass 2 — DOC partition and gas literals** | **CLOSED, with one non-bitwise change disclosed.** `sorption_capacity`/`C_aqueous` in `maoc.jl` own `θ + ρ_b·k_d`; six `src/` sites routed through them. `R_GAS`, `P_ATM`, `M_O2` are now the only definitions in `src/`. **`P_ATM = 101325` replaced a bare `101000.0`, moving `O2_CONST` +0.32 % (0.273808 → 0.274689)** — the old value was a rounding, the new one is the definition of 1 atm; effect on rates is fourth-decimal because ambient O₂ sits far above `L_B`/`L_F`. Runs before 07-30 differ in the last digits. |
| 07-29 | **A5 `optimize_soil3.jl`**, **A4 `diagnostics/`** | Folded into the same rebuild decision. Neither is worth un-forking if the tree is being restructured; both must be gone or rewired before any of their numbers is quoted again. |

Read-only when written; edits since are logged above and made against the tree, not this file's original text.

## Method and what to trust

Snapshot taken by `tar` on the device and extracted here, so it is not a `device_stage_files` copy — the stale-cache failure of CLAUDE.md §7 does not apply. `src/`, `test/`, `docs/`, `paper/`, `dev_notes/`, `Project.toml`; 130 files.

Five parallel readers covered API surface, duplication, comment hygiene, tests, and documentation consistency. **I independently verified the findings in Section A myself** by opening the files and comparing the code to the documents — those are the ones that would change a number. Sections B–F are the readers' work, spot-checked but not exhaustively re-derived by me. Where a claim is unconfirmed it says so.

Prior audits `docs/STRUCTURAL_AUDIT_2026-07-27.md` and `docs/TEST_AUDIT_2026-07-28.md` were read first; this report says which of their findings are still open rather than re-reporting them as new.

---

## A. Producing wrong numbers now

Six items. These are not tidiness.

### A1. The oxygen/Henry bug has recurred in the manuscript figure script — **DISPOSITIONED 07-29** (folder to be rebuilt; do not quote its output)

`paper/simulations/single_aggregate_physics/run_simulations.jl:143` and `:202`:

```julia
K_H_O = 0.032   # "Approximate Henry's constant at 20°C"
O_aq = O * θ / (θ + K_H_O * θ_a)
```

The solver's own convention is `src/solver/reactions.jl:89`: **`O_aq = O`** — `state.O` is already aqueous. This is the same class of defect as erratum 10, the one where `derived.jl` re-partitioned oxygen and every reported respiration ran at roughly a tenth of the oxygen the solver used. `derived.jl` now carries comments forbidding it and `test_postprocessing.jl:339-377` pins it; this file was never brought along.

Two errors compound. `0.032` is not the Henry constant in this codebase's convention either — `K_H_O2(293.15) ≈ 29.1` (`src/temperature/henry.jl`). `0.032 ≈ 1/31.25`, the reciprocal. With 29.1 the factor is ≈0.065; with 0.032 it is ≈0.984. **The two candidate wrong answers differ by 15×.** Everything downstream in that block — `R_B`, `R_Bb`, `Y_B`, `Γ_B`, `dB/dt`, `R_F`, `Γ_F`, the printed F_m budget — is built on it.

Same file, `:207-210`:

```julia
(F_i_0 + STANDARD_BIO.λ * F_n_0)      # λ on F_n
```

against `src/biology/fungi.jl:82`:

```julia
uptake_biomass = λ * F_i + F_n         # λ on F_i
```

λ is on the wrong pool. With `λ = 0.05` and an initialisation that loads F_i, these differ by up to 20×. The script's own printed label at `:260` repeats the error.

The same block also hand-rolls water content without `alpha_effective` (so it reports θ at a water content the solver never used), uses a constant `Y_F` where the solver uses the space-limited `Y_F_func`, and uses bare `D_Fm0` without the network factor. Verified: `run_simulations.jl:137-138, 193-194, 213, 235`.

**This script generates manuscript figures.** I did not check which figures or whether they are in the current draft.

### A2. ~~`M_max` is two different numbers~~ — **CLOSED 07-29**

Closed by erratum 12. Kept here because the *diagnosis* is the reusable part: it was filed for months as a parameter question and it was not one — it was one quantity with two implementations.

- `src/physics/initial_conditions.jl:382` — `M_max = soil.k_ma * soil.f_clay_silt * soil.ρ_b`
- `src/solver/reactions.jl:164`, `src/api.jl:189`, `src/postprocessing/derived.jl:113` — `soil.M_max`, a field defaulting to `10.0`

For De Gryze soil 3 that is **368 at initialisation against 10 in the solver**. The system starts above its own isotherm ceiling and desorbs from t=0. `degryze_config.jl` already reports the symptom and lowers `κ_d_ref` tenfold to suppress it — so a parameter is currently compensating for a duplication defect. The guard at `initial_conditions.jl:418` compares against the texture value, so it never fires.

`src/biology/maoc.jl:44` documents the texture form as if it were the definition, which makes the docstring the third statement.

### A3. ~~POM dissolution — code and all three documents disagree~~ — **CLOSED 07-29. My finding was wrong.**

`src/carbon/pom_dissolution.jl:36-38`:

```julia
microbial_factor = 0.5 * (B_0/(K_B_P+B_0) + F_n_0/(K_F_P+F_n_0))
```

`docs/ARCHITECTURE.md:116`, `docs/REFERENCE.md` §12 and `docs/GUIDE.md` §8 wrote the two Monod terms **multiplied**. **But `manuscript-4-5.tex` Eq.~(eq:R_P), l.364, has the code's form exactly** — ½(a+b). I checked the three internal documents and not the physics spec, and reported "all three documents" without qualifying which. The manuscript and the code never disagreed; the three internal copies were stale. Fixed in all three, with the mechanistic justification: depolymerisation is extracellular and substrate-side, so either community's enzymes suffice — a product form would assert POM cannot be depolymerised without both, which white-rot fungi and cellulolytic bacteria each falsify.

Not cosmetic. With both factors at 0.5, multiplicative gives 0.25 and the code gives 0.5. As either factor → 0 the multiplicative form goes to zero while the code floors at half the other — **in the code, bacteria alone sustain POM dissolution with zero fungi.** That is a qualitative behavioural difference in the term that drives the entire carbon input.

I cannot tell which is intended. The manuscript decides this, not me.

### A4. `converge_doc.jl` is the cited evidence for solver agreement, and it runs different parameters — **OPEN**, folded into the rebuild

`src/solver/timestepper.jl:19-21` cites the two-solver agreement as justification for keeping both. The measurement lives in `paper/de_gryze/diagnostics/converge_doc.jl`, which does not include `degryze_config.jl` and carries `κ_d_ref = 0.001`, `R_P_max = 2.5`, `k_d_eq = 0.05` (default), `w_E = 0.5` (default), `p_Gc = 0` (default).

Production runs `1e-4`, `1.5`, `0.005`, `0.05`, `1.0`. **The agreement was demonstrated for a model that is no longer the model.** All four scripts in `diagnostics/` have the same problem; `diagnose_speed.jl:42-43` states the requirement explicitly — "These MUST match run_degryze.jl exactly" — and no longer meets it.

### A5. `optimize_soil3.jl` is a complete forked configuration — **OPEN**, folded into the rebuild

It never includes `degryze_config.jl`. Twenty parameters differ, including `ρ_b` 1300 vs 1370, `f_clay_silt` 0.74 vs 0.56, `d_32` never set (placeholder 0.01 vs 0.0064), `k_d_eq` 10× off, `w_E` 10× off, `p_Gc` 0 vs 1, `f_bact` 100× off, `κ_d_ref` 50× off, and the un-renormalised 5-bin POM distribution. Its own comments claim parity ("identical to run_degryze.jl") that has not held for some time.

The file is already known to hold an invalid fit — but that warning lives in `degryze_config.jl` and BACKLOG, **not in `optimize_soil3.jl` itself**, which is where someone about to reuse it would look.

### A6. ~~A stale runtime warning is blocking the calibration path~~ — **CLOSED 07-29**, deleted

`paper/de_gryze/fitting/sensitivity.jl:200-206` prints:

> `r_agg is found by scanning grid nodes ... MWD is therefore a STEP function of every parameter and its elasticity is quantised ... do not put MWD into a gradient-based optimiser until that is fixed.`

`compute_r_agg` has not done that since the sub-grid interpolation went in. `test_postprocessing.jl:166-175` asserts the opposite (`length(unique(vals)) > 15` over a 25-point sweep). This is a false blocker printed at runtime on the exact path you are about to take.

---

## B. Will produce wrong numbers soon — duplication

Jeppson's drift arm. These agree today.

**Configuration.** Consolidation was applied to the two files in view at the time. Five of eleven consumers still hold their own copy: `optimize_soil3.jl` and the four `diagnostics/` scripts (A4, A5). `degryze_config.jl:7-9` records the previous instance of exactly this — the fix was applied to the instance, not the class.

**Van Genuchten — four implementations.** Canonical `water_retention.jl:85-93`; inline copies at `initial_conditions.jl:362-365`, `api.jl:181-184`, and `run_simulations.jl`. The two inline copies use `(-α·ψ)` rather than `α·abs(ψ)` and drop the `ψ ≥ 0 → θ_s` guard. **That path is reachable**: `van_genuchten_inverse` returns exactly `0.0` at saturation, so a soil at WFPS = 1.0 raises a negative base to the power 1.47 and throws `DomainError`. Unconfirmed by execution.

**~~Conservation weight `4πr²h` — three implementations~~ — CLOSED 07-30.** `conservation_weight(r, h)` and `conservation_weights(r_grid, h)` in `types.jl` are the sole definition; `GridInfo` calls the vector form, `reaction_step.jl` and `api.jl` call the scalar form. Bitwise result-preserving. `test_result_struct.jl` pins `grid.W == conservation_weights(...)` and adds the stencil property `W_i/(r_i²h²) = 4π/h` as the tripwire — that constancy is the reason the formula is what it is, and no earlier test stated it.

**~~Total system carbon — three implementations plus a fourth partial~~ — CLOSED 07-30.** Two entry points remain because there are two inputs: `compute_total_carbon(state, W)` for a raw state, `total_system_carbon(pools, k)` for `IntegratedPools`. Each is one implementation, `carbon_balance_table` and `result_to_dataframe` both route through the second, and `test_postprocessing.jl` asserts the two agree — the pin that replaces the collapse that cannot be made. `compute_total_carbon` gained `GridInfo` and weight-vector methods, so the four production call sites now pass `grid` and allocate nothing. `include_co2 = false` covers the normalisation denominator, which is **not** `k = 1` of the full sum unless the first output record is at `t = 0`. Bitwise result-preserving. Adding a tenth pool is now one edit in each of the two.

**~~Arrhenius — three~~ — CLOSED.** `D_fungal_translocation` calls `arrhenius(Ea, T, T_ref)`; the docstring that promised it now tells the truth.

**Langmuir–Freundlich — three.** **~~Sorption retardation `θ + ρ_b·k_d` — eight sites with no owning function~~ — CLOSED.** `sorption_capacity(θ, ρ_b, k_d)` and `C_aqueous(C, θ, soil)` in `maoc.jl` own it; all six `src/` sites route through them. `test_biology.jl` asserts the defining identity and that the relation is linear in `C`, so the ω dilution passes through it.

**Repeated literals.** `src/` is clean: `R_GAS`, `P_ATM`, `M_O2` are the only definitions and `degryze_config.jl` builds `O2_CONST` from them. **Adopting `P_ATM = 101325` was not bitwise result-preserving** — `O2_CONST` moved 0.273808 → 0.274689, +0.32 % — because `101000.0` was a rounding. Every remaining bare `8.314`, `101000.0`, `0.032` sits in `optimize_soil3.jl` or the four `diagnostics/` scripts, which are **now 0.32 % off the production config** as a direct result; that is a second reason A5 should not wait. `0.443` is at 5 sites: `degryze_config.jl:204` (`F_C_POM`, the intended owner), plus independent defaults in `degryze_soils.jl:67` and `postprocess_dataframe.jl:444` and two in `optimize_soil3.jl`. `src/postprocessing/population.jl` correctly requires it with no default. `O2_saturation`, the one function in `src/` whose job ambient O₂ is, is still called by nothing outside its own test.

**Cost arm.** Loop-invariant work in the file already profiled as the bottleneck: `crank_nicolson.jl` rebuilds `1/(r_i²h²)` and both face radii per node per species per half-step — ≈1.6×10⁹ redundant divisions per aggregate-run at the measured step count. `θ_s^(2/3)` is recomputed ≈2.4×10⁸ times, all returning the same number. `reactions.jl:93-127` forms ten rate×Arrhenius products per *node* when the temperature cache is per *step*. `workspace_updates.jl:51` computes `cache.D_Fm` every step and **nothing reads it** — `test_physics.jl:382` says so explicitly.

**Closed 07-29:** `M_max` (two implementations, above) and the `k_ma` unit ambiguity. The `κ_d_ref` override that was masking it is still in place — half its justification is now void and it should be re-examined on the next run.

**Sanctioned and verified.** `mol.jl`'s Laplacian deliberately reproduces `crank_nicolson.jl` term for term; a reader checked all three branches and they agree, including the exact `D₁` cancellation at the flux boundary. **No test pins them together** — the sanction is a comment, not an invariant.

---

## C. Tests

**The default solver has no test.** `grep -rn "run_aggregate_stiff" test/` returns nothing. Untested on the production path: the CO₂ recovery from carbon balance (which *is* the respiration output, and is derived rather than integrated — a sign error there gives a plausible, monotone, entirely wrong curve), the output schedule, `abstol_scale`, `sparse_jac`, the O₂ Dirichlet pin, `t_delay`, the `co2_monotonic` diagnostic, and the `NaN` mass-balance contract.

Two behavioural divergences that no test would catch:

- `n_grid` defaults to **50** in `run_aggregate` and **200** in `run_aggregate_stiff` — so the two solvers solve different discretisations by default, which is precisely the comparison the two-solver design exists to enable.
- With no `output_times`, `run_aggregate` emits on a 1-day interval and `run_aggregate_stiff` emits **exactly two records**. Swapping solvers without passing `output_times` silently loses the time series.
- Tripwire: `test_output_times.jl:63` asserts `abs(mass_balance_error) < 1e-12`, which fails for every stiff record by design. Whoever adds stiff coverage hits this first.

**Still open verbatim from the 2026-07-28 audit:** the 217-assertion `total >= aggregate` block (`test_postprocessing.jl:28-37`) — a subset sum compared to its superset, which passes for any `compute_r_agg` including a broken one, while `pools.r_agg` itself is never asserted anywhere; `test_crank_nicolson.jl:100` with 1400× slack; the "Spatially varying D" testset that uses a uniform field, so `D` never enters the arithmetic, still carrying the comment `# No crashes = success!`; `parameters.jl` has zero validation (`K_B=0`, `n_vg=1` construct happily).

**What has improved, and it is real:** `test_closure.jl` is the strongest file in the suite and includes a self-falsification check. `test_postprocessing.jl:132-180` — the new `p_Gc` coverage — uses a constructed linear profile with closed-form crossings and is ahead of the code rather than behind it. `test_tessellation.jl` and `test_particle_size.jl` test defining properties rather than formulas.

**Mechanics are clean.** `runtests.jl` includes all 20 files; `[targets] test` is complete; nothing would error at runtime.

---

## D. Structure

**Your api.jl observation.** Reading it as "the file called `api` does not contain `run_aggregate_stiff` — it points at `mol_solve.jl`": correct, and it is the symptom of the real thing. `api.jl:2` says "User-facing API"; it holds one of two entry points and the default lives under `solver/`. Worse, **the `solver = :stiff | :split` switch does not exist in `src/` at all** — it lives in `paper/de_gryze/postprocess_dataframe.jl:311-322`. The package exports two sibling functions and leaves solver selection to be reimplemented by every consumer. That is the clearest sign the stiff path is bolted on rather than integrated.

Also: ~25 lines of setup and teardown are duplicated verbatim between `api.jl` and `mol_solve.jl` (grid, environment, state branch, `OutputRecord` packing, `SimulationResult` construction). `timestepper.jl:24` warns against exactly this.

`api.jl` docstring defects, verified: the signature block says `ic::Union{Nothing,InitialConditions}=nothing` where the real signature is `ic::InitialConditions=InitialConditions()`; the "State initialization priority" list has three levels where the code has two; `ω` is documented nowhere. `create_initial_state_legacy` (42 lines) is dead — its own docstring says it is "kept for backward compatibility with existing scripts" and no such script exists.

**Layering inversion.** `solver/mol_solve.jl` is included before `api.jl` and calls `compute_total_carbon`, which is defined in `api.jl`. Julia resolves it at call time, but `solver/` depending on `api.jl` is backwards; those carbon utilities belong below the solver.

**Missing exports that documented text tells you to use.** `EnvironmentalDrivers` is the type parameter of the exported `SimulationResult` and is named in a public docstring, but is not exported. `create_initial_state` is the documented state-construction path and is not exported. `TemperatureCache` / `update_temperature_cache!` have no public alternative and `run_degryze.jl` reaches for them qualified.

**Debris:** `src/physics/initial_conditions.jl_backup` (8.4 kB, February) sits inside the package tree. `Plots` and `Distributions` are declared package dependencies and appear nowhere in `src/` — every user of the library loads Plots.

---

## E. Documentation

**REFERENCE.md contradicts itself in four places**, which matters more than it contradicting the code: §3.2 gives `F_i_min = 1e-4` where §5a and the code say `1e-6` (§3.2 has swapped in `F_n_min`'s value); §3.4 gives `K_E` 2× the code and §5a; §7's bacterial allocation uses `R_B` where the code and §17a use `R_B − R_Bb` — including a sign error in the starvation branch; §10 says `Resp_B = 0` at starvation where the code and §17a give `R_Bb > 0`. §4.2 lists two fields, `clay_fraction` and `silt_fraction`, that do not exist anywhere in `src/`.

**ARCHITECTURE.md needs a decision, not a patch.** It documents `AggregateSolver` and `run_aggregate!` — **neither exists**. Its `SoilProperties` listing still carries `k_F`, `χ` and `a_p`, removed by erratum 11. Its file tree omits the entire `postprocessing/` directory (12 of 34 exports) and four `solver/` files including both halves of the default solver. Its §2 still presents Strang splitting as *the* time integration, with the FBDF upgrade written as a hypothetical contingency.

**GUIDE.md §9 is missing the O₂ capacity divisor** — erratum 8 recurring in the one file that sweep missed. §5 inverts λ (says F_i is the primary uptaker; it contributes at 5%). Same inversion at `parameters.jl:66` and `ARCHITECTURE.md:460`.

**`src/biology/maoc.jl:15` and `:80-81`** still assert `S_C` needs the factor `(θ + ρ_b·k_d)/k_d`, with "NOT just -J_M! Missing factor breaks carbon conservation." Erratum 2 removed that factor in February and the code uses `-J_M`. The docstring survives inside the very file the erratum is about.

**Citations.** Of 34 sampled, 13 broken. ~~`ARCHITECTURE_CLAUDE_CODE.md` cited from `parameters.jl`, `types.jl`, `REFERENCE.md:7`, `GUIDE.md:5`~~ — **FIXED 07-29**, all four now point at `docs/ARCHITECTURE.md` and say it is stale; REFERENCE's header also now names the two live manuscript files and states that nothing is senior. Remaining: Every citation into `reactions.jl` is shifted by about +12 lines (the file grew; six consumers still point at the old positions). BACKLOG item 3 cites `reactions.jl:152`; the line is 164. BACKLOG's "2nd-order claim at timestepper.jl:60" points at a docstring line about `dt_min`.

**`dev_notes/` is gitignored and four `src/` files plus six REFERENCE citations point into it.** On a fresh clone all ten dangle. It is the provenance chain for the Falconer α/β convention — the subject of errata 6, 7 and 8. BACKLOG has this below the hygiene line; it belongs above.

~~**`dev_notes/REFERENCE_PART_V.md`**~~ — **REMOVED 07-29** to `_to_delete/`. Orphaned earlier draft of REFERENCE Part IV carrying a claim REFERENCE has since retracted in as many words. Nothing cited it.

**Point-in-time reports in `docs/`.** `respiration_crash_diagnosis.md` **REMOVED 07-29** to `_to_delete/` — February diagnosis, every `M_max` figure 10× high after erratum 12, compared against a `soil.M_max` field that no longer exists, cited by nothing. `julia_falconer_deviations.md` KEPT: `fungi.jl:17` and `:315` cite it, so it is load-bearing, not tangential; the two wrong-path citations in REFERENCE were corrected. The two dated audits stay as frozen records.

---

## F. Comment hygiene — and my part in it

Measured: `src/` is **68% prose by line** (2,081 code, 4,352 prose). `test/` is 22% and healthy.

`degryze_config.jl` measured as you suspected: **233 lines — 54 code, 161 comment**, three comment lines per line of code. One block runs 66 consecutive comment lines to document five keyword arguments. It also contains a visible botched merge at `:116-119` where `To isolate:` dangles and the instruction is then given again in a different voice. That is mine, from this afternoon.

Three genres, and the third is the one that matters:

**Stale claims.** Nine of eleven confirmed cases date from today's session — the pattern is that a change was made, its rationale was written carefully at the change site, and the four or five places describing the old behaviour were never revisited. `N_POM_BINS` went 5 → 10 and **five files still say five**, including `scan_transport.jl`, whose header states a *conclusion* ("a finer POM size distribution cannot produce a 21-day rise") whose premise was subsequently changed — that is the argument that led to doubling the bins. `scan.jl:16` cites `R_P_max = 2.5` as "the current value", three revisions stale. `run_degryze.jl:194-201` describes assumption A1 and asserts it "contains no reference to any measured MWD" — the exact property A1′ does not have, and the property `setup.jl` depends on knowing.

**Duplication of documents that already exist.** `population.jl` has a 103-line docstring whose `# Physics` section is REFERENCE §5c heading-for-heading. `aggregate_radius.jl` has 211 prose lines to 38 code lines, with the 1/d₃₂ derivation written three times (twice in that file, once in REFERENCE §4.4). Every duplicated pair is a future divergence, and §E shows it has already happened three times.

**Narrative.** Roughly 100 of `degryze_config.jl`'s 161 comment lines belong in a changelog, the spec, or BACKLOG. `sensitivity.jl` tells the "these were labelled controls / that was wrong / the run proved it" story twice, once in the header and once at runtime.

Estimated **~570 lines net removable, ~250 of them from two files**, with no information loss — because the canonical text already exists in REFERENCE.

Two things worth keeping exactly as they are: `setup.jl:99-102`, which prevents fitting to a datum the model was constructed from and has nowhere else to live; and `src/parameters.jl:293-299`, the κ_b provenance, which exists to stop someone restoring `0.0143869` on the belief that it was measured.

**Open questions embedded in code that should be BACKLOG items.** `reactions.jl:124-127` — "ζ is a dimensionless splitting fraction, not a rate, so scaling it by an Arrhenius factor is dimensionally wrong and can drive it above 1 (which would flip the sign of `immobil_n`). Clamped until this is resolved." An unlisted open bug with a live workaround in the code, and ARCHITECTURE §5.1 presents the same behaviour as design. Also `degryze_config.jl:154-156` (nothing sustains background biomass in the far field) — not in BACKLOG.

---

## Proposed order

*(revised 07-29 after the rebuild decision — `paper/simulations/`, `optimize_soil3.jl` and `diagnostics/` are to be rebuilt rather than repaired, so items 2 and 3 became one deferred item.)*

Nothing here is on your critical path for the De Gryze fit except A6, which takes a minute.

1. **A6** — delete the false quantisation warning. One edit, unblocks gradient-based work.
2. ~~**A2**~~ — done 07-29.
3. **Deferred to the rebuild** — A1, A4, A5. Until then nothing under `paper/simulations/`, `optimize_soil3.jl` or `diagnostics/` may be quoted. `converge_doc.jl` is the load-bearing one: the two-solver agreement claim in `timestepper.jl` rests on it and is currently unsupported.
4. **One test calling `run_aggregate_stiff` end-to-end**, plus one pinning `mol_laplacian` against `crank_nicolson_step!`. The default solver and the sanctioned duplication are both currently unguarded.
5. **Re-examine `κ_d_ref`** — it was lowered tenfold partly to mask A2, which is gone.
6. **Stale claims sweep** — the `N_POM_BINS`, `R_P_max` and A1/A1′ cascades. Cheap, and they are actively misleading the work in progress.
7. **Structure** — move the `solver=:stiff|:split` switch into `src/`, extract the shared setup/teardown, fix the two default divergences (`n_grid`, output schedule), export `EnvironmentalDrivers` and `create_initial_state`, delete `create_initial_state_legacy` and `initial_conditions.jl_backup`.
8. **Prose** — `degryze_config.jl` and `aggregate_radius.jl` first; delete in-code copies of text REFERENCE already carries.
9. **Documents** — REFERENCE's four internal contradictions, then the ARCHITECTURE decision, then the citation sweep.

One process rule would have caught seven of the eleven stale claims: **when a value changes, grep the tree for the old value before committing.**

---

## WHERE WE STOPPED — 2026-07-29, resume here

**The ω overlap convention. Decided in principle, not implemented.**

*Fullest treatment: `SI_tessellation.tex` (Overleaf).* It derives ω, proves
exactness for diffusion and linear reactions, identifies the Monod nonlinearity,
bounds the error by `1/ω`, and gives an implementation table the code matches.
**It is a working draft, not a specification** — as is the main text, as is
REFERENCE. Where they disagree, none is presumed right. The audit reached the
Monod result independently, which corroborates it rather than duplicating it.
Read S5–S8 before re-opening anything here.

*Four open questions, none settled* (REFERENCE §26 erratum 13 has the detail):

1. **Is O₂ diluted?** The SI's error analysis dilutes it; the code does not.
   The sink is then ~45× weaker than physical and near-POM anoxia does not form.
   A physics call: O₂ is supplied from the headspace, which argues for leaving it
   undiluted — but then source and sink sit on different scales.
2. **Does the validity argument survive the configuration?** SI S7 rests on
   `C_bg ≈ 0` ("native POM was removed"). `degryze_ic` initialises the full
   measured SOC, 2.14 % → 29.3 µg/mm³. The `1/ω` bound is unaffected; the
   argument that it does not matter is. Revise the argument or the config.
3. **Not yet analysed** (drafting work, not a defect). Not the isotherm (`1/k_L`, `M_max`), the
   sigmoid thresholds (`B_min`, `E_min`, `F_i_min`), the space-limited yields
   (`B_S`, `F_S`), `K_B_P`/`K_F_P`, `K_Fm_net`, `ε_F` in `Π`, or the binder
   against an undiluted `G_c` — which is MWD.
4. **Main text and SI describe different constructions.** `manuscript-4-5.tex`
   l.165–181 uses one domain factor `f_d` and claims exact non-overlap; that
   holds only for `f_d = f_p`, and the SI's subject is `f_m > f_p`. Must be
   reconciled before submission whatever else is decided.

*What was proposed and rejected.* Shrinking to `f_domain = f_pack` (ω = 1). It
was committed on 07-29 and reverted the same day. **Owner's objection, which
stands:** the outer boundary must stay clear of the near-POM gradients, and at
`f_pack` the zero-flux surface sits 1.5 mm from the POM while the day-21 DOC
front travels 3.2 mm. Zero flux at the Voronoi radius is exact only for
identical neighbours; with ten POM size classes it is not, and the error grows
as the cell shrinks. Enlarging the domain does not reduce that error either.

One correction to the record, for whoever picks this up: the **oversized domain**
is documented as deliberate in `dev_notes/manuscript_stability_alignment.md`
§2.4 ("coupled to the overlap correction ω — not a free numerical convenience"),
`manuscript_population_upscaling.md` §X.1, and REFERENCE §5b. But all three
describe the correction as **dilution of concentration at initialisation**.
Correction at the summation stage is nowhere in the documents. The domain choice
is on record; this repair is not.

*Owner's position, 2026-07-29.* **The overlap is required and ω cannot be 1** —
the model will be run on annual root input and field scenarios, not only De
Gryze, and the model volume necessarily exceeds the physical volume. Each domain
starts diluted and then runs as if alone; the dilution stands in for neighbours
sharing the volume. They do not actually share that way, but the alternative is
explicit coupling, which needs packing detail the model does not carry and may
not earn at this level of abstraction. The requirement is **consistency, not
accuracy**: the summation to bulk must account for the overlap.

*The open question is therefore where the `1/ω` is applied*, not whether ω
exists. Candidate: keep the dilution as the physical closure and make the
summation consistent with it — versus initialising physical and correcting only
at summation. Not decided.

What that will require, none of it started:

- Remove the dilution block from `create_initial_state` (Step 8).
- Find every place a per-domain background quantity becomes a per-sample total
  and insert `1/ω` there — `population_statistics`, `integrated_pools`,
  `carbon_balance_table`, and the CO₂ sums. **POM and CO₂ must not be divided**;
  they already are per-particle totals.
- Decide what the far field means when each domain draws on 30× the carbon it
  owns. This is the honest cost of the choice and it has not been quantified.
- Re-run everything. **Every parameter fitted before this is suspect** — `κ_b`
  (20×), `w_E` (10×), `f_bact` (100×), `f_eps` (10×), `κ_d_ref` (10×) were all
  moved to compensate distortions that were ω's.

*Also unresolved and downstream of it:* `κ_s_ref` (0.1 → 0.01, undocumented) and
`κ_d_ref` (0.01 → 1e-4). The ratio went 10 → 100, i.e. 10× stronger hysteresis,
changed silently. The elasticity that justified `κ_d_ref` was measured on the
`M_max` bug. Do not touch these until ω is settled.

---

## Open decisions

- **A3 — POM dissolution.** Product or mean of the two Monod terms? The code and all three documents disagree, and it changes whether bacteria alone can drive dissolution. The manuscript settles this.
- **This file vs `BACKLOG.md`.** Two lists of open work with status. One should absorb the other — my suggestion is that BACKLOG keeps the *decisions and blockers* and this file keeps the *findings with evidence*, or that this file is deleted once every item is closed or moved. Either way, not both indefinitely.
- ~~**A2 — `M_max`**~~ — settled 07-29: derived, `k_ma` dimensionless, LM for De Gryze.
- **ARCHITECTURE.md** — rewrite against the current code, or mark historical and move the still-valid material (§3 discretisation, §5 temperature framework) into REFERENCE?
- **`dev_notes/`** — track it, or stop citing it from `src/`?
- **`optimize_soil3.jl`** — repair or delete? Its fit is already declared unusable.
- **The `ζ` Arrhenius clamp** — `reactions.jl:124-127` says the scaling is dimensionally wrong and clamps it. ARCHITECTURE §5.1 lists it as design. Which?
