# Backlog

**The only list of open work.** Updated as items move; **not** a journal —
history lives in `docs/REFERENCE.md` §26 and the archived audits.

One rule, because breaking it cost a day: **an item is closed by pointing at code
or at the manuscript, never at another note.** Three times on 2026-07-30 a closed
item was re-reported as open because the note had not been updated when the source
was, and the note was read as evidence. If this file and the code disagree, this
file is what is wrong.

Absorbed the housekeeping audit on 2026-07-30 (items 8–14 below). It is archived
at `docs/_archive/AUDIT_20260729.md` and holds the evidence trail for everything it
closed.

---

## Dependency order

**Speed is DONE** (item 1). Duplication is DONE (archived audit §B, four passes).
**Verification and structure are DONE** (items 8, 9, 15). Remaining: **anchors →
decisions → fit.** The architectural work is finished; what is left is
calibration, plus documentation items 12–14.

Items 8, 9 and 15 are DONE (2026-07-30). The solvers agree, the production solver
is under test, the two Laplacians are pinned to each other, and the solver switch
and both diverging defaults are in `src/`.

Fitting free parameters while anchored ones sit off their cited values absorbs
the Group B error into the free ones and makes the deviation invisible. Items 3
and 4 are decisions that change model behaviour and should be settled before,
not during, a fit.

---

## 1. Speed — DONE (2026-07-29)

**Solved by replacing the integrator, not by optimising the old one.**

| | split (reference) | stiff (default) |
|---|---|---|
| wall, 45 d, n=200 | 38.1 s | **1.59 s** |
| accepted steps | 391,773 | 2,088 |
| RHS evaluations | — | 7,045 (3.4/step) |
| Jacobians | — | 39 |
| steps/day, 15–45 d | 8,706 | **8.9** |
| step trend | flat at floor, required step still falling | falls through the run |

**24x on wall clock, 188x on steps.** The scaling matters more than the factor:
the split solver's cost per simulated day grows (|S_C|/C at the POM-adjacent
nodes rises to 2055 /day by day 45 and was still climbing), while the stiff
solver's falls. Multi-year integration is possible with one and not the other.

Practically: one parameter evaluation across five diameters is now ~8 s instead
of ~200 s, so a 500-evaluation fit is about an hour at full n=200 — no need for
a coarse grid. (Measured 07-30 on the current config: 1.16 s per aggregate at
n=200, ~12 s for all ten classes.)

**Diagnosis, for the record.** Δt was set by DOC at the nodes next to the POM
surface, not by the 1e-6 denominator guard (tested and killed) and not by oxygen
(tested, would permit Δt ~7e-3). Both were hypotheses of mine, both wrong. From
about day 15 the criterion demanded Δt below the 1e-4 floor — 4.9e-5 by day 45 —
so the reference solver was taking steps it had itself judged too large and could
not report it, because `n_rejected` is guarded by `dt > dt_min`.

**Agreement.** r_agg identical to 7 figures, CO₂ 1.8e-5, P 1.4e-4, M 1.3e-4,
E 1.3e-3, F_i 3.6e-3, B 2.1e-3. DOC differs by 9.6% at the POM-adjacent nodes and
moves nothing downstream. The six stiff runs across a 100x tolerance range agree
to 1.5e-4; the split solver cannot be convergence-tested at all, because
lowering `dt_min` below what its 10%-change criterion demands changes nothing
(1e-5 and 1e-6 gave bit-identical output).

### Open, small, not blocking

- **`abstol` default is weak.** `max(1e-10, 1e-8·|u₀|)` penalises any pool that
  starts near zero and grows: `F_m` starts at 3.3e-8 and reaches 3.5e-4, so it is
  controlled four orders tighter than its own scale. A per-species characteristic
  scale, or a plain scalar, is probably better. Measured cost is small.
- **`mol_rhs!` allocates 7 length-n vectors per call** so AD sees the right
  element type. `PreallocationTools.DiffCache` is the remedy. Not applied.
- **Split solver's controller cannot be converged by its documented knobs.**
  `dt_min` is inert once below what the criterion demands, and the criterion is a
  10%-change-per-step stability guard, not an error tolerance. Recorded as a
  property of the reference path, not something to fix.

## 1a. Carbon closure — test added (2026-07-29)

`Σ S_pools = −Resp_total` is the model's central invariant and was untested
through 1703 assertions. It is now proven analytically (REFERENCE.md §17a) and
tested pointwise in `test/test_closure.jl` — exact, no integration, no solver.

The proof turns on `softplus(x,ε) − softplus(−x,ε) = x`, exact for any ε. **Any
replacement for softplus must satisfy that antisymmetry or the carbon budget
breaks.** Tested directly, not only through the closure.

## 2. De Gryze parameter optimization

**Open design question: what is shared across the five soils?**

Per-soil, because measured (already wired through `degryze_soils.jl`):
`ρ_b`, `f_clay_silt`, `θ_s`, `d_32` (Table 1 texture), `SOC` (Table 1),
`T_0`, `ψ_0` (from 60 % WFPS at that soil's porosity).

Candidate shared: everything biological, plus `κ_b` and `w_E`.

The case for sharing: the same microbial community and the same binder chemistry
operate in all five soils — only the physical environment differs. If different
biology per soil is needed to fit, the model is absorbing the texture effect
rather than explaining it, and the cross-soil comparison stops being a test.

The case against, for one parameter: soil 5 has CEC 17.4, about 2.6× the others,
which is a mineralogical difference the texture triple does not capture and
which plausibly affects MAOC. That is also precisely where a per-soil free
parameter would silently absorb the texture signal.

**Targets.** Table 3 formation rates (5 scalars) are the cleanest — they
difference out the day-0 texture baseline. Table 2 MWD at day 1 and 21 (10
values) and the CO₂ curves are secondary; note the CO₂ caveat in spec §5.3
(no unamended control, soil 5 respires 119 % of carbohydrate-C).

**Blocked by:** items 1, 3, 4, and the Group B anchors below.

Also: `optimize_soil3.jl` was archived on 2026-07-30
(`paper/_archive/degryze_tooling_20260730/`) rather than migrated. It forked the
configuration and hard-coded both the soil-3 properties and the measured CO₂
series; its stored fit was invalid regardless, having been made against the
broken POM normalization. Calibration is a separate future project: get the fit
visually close by hand first, then build one reusable fitting routine that every
example can use, tested against De Gryze.

## ~~3. `M_max` has two values~~ — **CLOSED 2026-07-29**

Resolved by REFERENCE.md §26 erratum 12. The `SoilProperties.M_max` field is
deleted; `maoc_capacity(soil) = k_ma·f_clay_silt·ρ_b` in `src/biology/maoc.jl`
is the only definition. `k_ma` is now dimensionless and anchored to Georgiou et
al. (2022) — 0.048 low-activity, 0.086 high-activity — which also corrects a 10×
error in the old 0.48 (it was 48 mg/g with the mg→g conversion dropped).

Remaining, and it is a soil property rather than an open decision: which
mineralogy each soil has. De Gryze's five are one Belgian loess-derived profile,
so all low-activity (`degryze_soils.jl`).

Downstream: `degryze_config.jl`'s `κ_d_ref = 1e-4` override was justified partly
by the desorption pulse this caused. That half of its justification is void; the
value was left unchanged so the M_max fix can be measured on its own.

## 4. `delta`, `beta`, `alpha`, `zeta` — calibration, not architecture

**Not a defect and not blocking.** The architecture question was settled
2026-07-30: ζ is a constant splitting fraction with no temperature dependence,
matching Falconer [2005] p. 1728. The values are a calibration matter and are
deferred to the De Gryze fit.

What De Gryze runs today: **δ = 1.0, β_i = β_n = 0, α_i = 0.1, α_n = 0.15,
ζ = 0.2** — every one inherited from the package default, which is the MATLAB
precursor's set, which is Falconer with θ = 1 and β = 0. `degryze_config.jl` sets
none of them.

Falconer's published values are tabulated in `src/parameters.jl` at the default
block, so whoever calibrates has them at the call site rather than in a note.
Two facts from there worth carrying: **ζ = 0.01 in every simulation of both
papers**, against our 0.2; and **θ = 3.0 is the only published value above 1**,
used for the [2005] cases that produce the ring/aggregation patterns, while
[2008] uses 1.0 throughout.

With β = 0 the exponent has no linear loss to cross over against, so δ only
rescales the transition rate. **δ and β are one decision, not two.**

---

## 5. Group B parameters off their cited anchors

`D_B_rel` 1000×, `ω_E` 1500×, `ζ` 20×, `r_B_max` 3.7×, `μ_B` 3.3×.
`REFERENCE.md` §5a, Group B. Return these to their anchors, or record why not,
**before** any fitting.

## 6. `r_agg` step function

`pct_gt2000um` jumps between day 2 and day 5 (`f_agg` 0.010 → 0.249) and then
sits flat. May be downstream of 1, 3 or 5 rather than its own defect — revisit
after those.

## ~~7. CO₂ overshoot~~ — **CLOSED 2026-07-30, by manual tuning**

Was 3599.7 against a measured 2139 µg-C/g at day 21 for soil 3, about 1.7×.
Closed in the course of the hand-tuning campaign, not by a single identified
cause, so there is no defect to point at — the parameters that moved are the ones
listed as working assumptions below.

---

## ~~8a. DOC diverges between the solvers at 45 days~~ — **CLOSED 2026-07-30 by archiving the split solver**

**The finding is about the reference, not about the stiff solver.** A
disagreement with an instrument that evaluates rate laws on negative arguments and
folds its clipping into CO₂ unreported is a fact about that instrument. The stiff
solver at 45 days — and at 22 — is **unverified, not suspect**: there is no
evidence either way, because the only thing that could have produced evidence was
not fit to produce it.

Chasing it showed the reference implementation was not fit to be a reference:
rate laws evaluated on negative arguments, clipping folded into CO₂ without being
reported, lagged coefficients, and first order where its docstring claimed
second. On every accuracy axis the stiff solver is the better formulation, and
repairing all four would have destroyed the one thing the split solver was kept
for — the cheap, independent CO₂ integration.

Archived to `_archive/split_solver_20260730/` with the full disposition:
`timestepper.jl`, `reaction_step.jl`, `diffusion_step.jl`, `crank_nicolson.jl`,
`tridiagonal.jl`, `finite_volumes.jl` (no caller at all), their four test files,
and `verify_solver_agreement.jl`. `Workspace`, `update_effective_diffusion!`,
`run_aggregate` and the `run_aggregate(solver, ...)` dispatch went with them.

**What this costs, stated plainly: there is no integrated mass-balance
measurement in the package today.** `test_closure.jl` still asserts the closure
identity pointwise, which REFERENCE §17a ranks above the integrated check anyway.

## NEXT — cumulative respiration as an explicit state

A version carrying cumulative respired carbon as a state variable, written as a
**separate test case for the mass-balance check**. It restores the measurement
inside the production solver instead of borrowing it from a second integrator.

Carry it **per node**, not as a global scalar: one CO₂ scalar makes the Jacobian
row dense across every node, and the sparse linear solve is already 33 % of
runtime. A per-node field is purely local, keeps the block-tridiagonal structure,
and takes the block from 8×8 to 9×9.

`mass_balance_error` stays in `OutputRecord` and reports NaN until then — the
tests now assert `isnan` rather than a tolerance.

---

## ~~8. Solver agreement (22 days)~~ — **RUN 2026-07-30.** Result in REFERENCE §17a

The two solvers describe the same model, and the split solver's production floor
is what was hiding it.

At `dt_min = 1e-4` the split solver differs from stiff by up to 4.2e-2 (F_m).
Refining to 1e-5 shrinks **all eight** field gaps, median factor 9 — first order
in Δt, i.e. the documented one-step lag in θ and D, not a defect. The split
answer converges *toward* stiff, which is what licenses trusting stiff.
Carbon balance on the split run: −2.7e-13.

The reported quantities were never at risk: at the production floor `r_agg`
agrees to 3.6e-3, CO₂ to 2.2e-4, POM to 6.4e-5. The 4 % is in the fast interior
pools.

**Consequence:** every earlier split-vs-stiff comparison was made at
`dt_min = 1e-4` and therefore *understates* the agreement.

**Not a validation.** Two discretisations agreeing is consistency, not
verification against an analytic solution.

**Re-run** after any change to the source terms or the discretisation, and after
any parameter change that alters stiffness. The claim in §17a is only as current
as the last run.

*The script had two defects, fixed the same day: no warm-up, so it printed a 1.6×
speedup against the documented 24×; and its verdict tested the size of the
movement on refinement, not the direction, so it would have said "trust stiff"
even if the split answer had moved away. The measured gaps were unaffected.*

---

## ~~9. `run_aggregate_stiff` has no tests~~ — **CLOSED 2026-07-30**

`test/test_stiff.jl`, registered in `runtests.jl`:

- **`mol_laplacian` pinned against `crank_nicolson_step!`.** With θ = 0 the
  Crank–Nicolson left-hand matrix is the identity, so `u_new = u + dt·L[u]`
  exactly and `(u_new − u)/dt` recovers its Laplacian to roundoff. Compared node
  for node over three grids, with the inner flux off, on, and reversed, against a
  non-smooth profile and a spatially varying `D` — a stencil error that cancels
  on a linear profile does not cancel there. This is the invariant the comment in
  `mol.jl` had been standing in for.
- **End to end**: result shape, output times honoured, `retcode`, POM consumed,
  recovered CO₂ monotone and non-negative, all eight fields finite and
  non-negative, and `mass_balance_error` **NaN rather than 0.0** — reporting 0.0
  would read as a passing check on something never checked (§17a).
- **Short-horizon agreement** between the two solvers, threshold 0.1. Deliberately
  loose: `mol.jl` documents three differences by construction, and §17a puts the
  22-day gap at 4.2e-2 at the split solver's production floor. 1e-3 is the lag,
  1e-1 is a bug — this guards the bug, so it cannot go flaky.
- **The dispatcher** and the two shared defaults.

## 10. No parameter has a firm value yet

Not a defect, and this replaces the entry that treated `κ_s_ref` and `κ_d_ref` as
one. **Unless a parameter carries a citation it is a working assumption.** The
current set came from hand-tuning against soil 3 and is the best available
starting point, nothing more. That applies to `κ_s_ref` (0.1 → 0.01) and
`κ_d_ref` (0.01 → 1e-4) as much as to `κ_b`, `w_E`, `p_Gc`, `R_P_max`, `f_bact`
and `f_eps`.

Two consequences worth keeping in view rather than acting on:

- The package defaults are earlier guesses, not reasoned values. The best De
  Gryze set should become the defaults once there is one.
- The elasticities in the archived sensitivity sweep were measured while the
  `M_max` bug was live, so they cannot be quoted.

The parameters that *do* carry an anchor, and the size of the gap to it, are
item 5.

---

## ~~11. The `ζ` Arrhenius clamp~~ — **CLOSED 2026-07-30: ζ is constant**

**Decision (owner): no temperature effect on ζ.** There is no credible reason for
the *share* of settling fungal carbon that lands insulated to rise with
temperature. The Arrhenius factor and the `min(·, 1.0)` clamp are both removed,
from `reactions.jl` and from the `derived.jl` post-processing copy. `ARCHITECTURE.md`
§5.1 no longer lists ζ among the rates carrying `Ea_F`.

**This changes results, and not slightly.** De Gryze incubates at 25 °C against
`T_ref` = 20 °C with `Ea_F` = 55 kJ/mol, so `f_fun` = 1.46 and the effective ζ was
**0.292, not 0.2** — the value in the parameter table was never the value in
force. Removing the scaling drops it by 31 %. Less carbon lands insulated, and
`F_i` is the binder in `F_i + w_E·E ≥ G_c`, so expect a lower `r_agg` and a lower
MWD. **The hand-tuned fit will need revisiting.**

---

## 12. `ARCHITECTURE.md` — rewrite or retire

698 lines. It documents an `AggregateSolver` type that does not exist, omits the
default solver entirely, and its `SoilProperties` listing carries fields removed by
erratum 11. Its still-valid material is §3 (discretisation) and §5 (temperature
framework). Decide: fold those two into REFERENCE and archive the file, or rewrite
it against the current code. `src/parameters.jl`, `src/types.jl`, `REFERENCE.md`
and `GUIDE.md` all point at it with a staleness warning, which is not a stable
resting place.

---

## 13. REFERENCE.md internal contradictions

Four, found 2026-07-29 and not yet fixed: §3.2 on `F_i_min`, §3.4 on `K_E`, §7 on
bacterial allocation, §10 on `Resp_B`. Also §4.2 names `clay_fraction` and
`silt_fraction`, which do not exist as fields, and the `reactions.jl` citation
cluster is off by roughly 12 lines. At 2096 lines this file is the one everything
else points at, so a contradiction in it propagates.

---

## 14. Prose volume in `src/`

`src/` is 68 % prose by line; roughly 570 lines are removable, about 250 of them in
`degryze_config.jl` and `aggregate_radius.jl`. Two things to keep exactly as they
are: `src/parameters.jl` on the `κ_b` provenance, which exists to stop someone
restoring `0.0143869` on the belief that it was measured; and the note preventing a
fit to a datum the model was constructed from. Lowest priority here, but it is the
same failure as the documentation sprawl — narrative kept where an invariant
belongs.

---

## ~~15. Structure~~ — **CLOSED 2026-07-30**

All four, verified against the tree:

- **The `:stiff`/`:split` switch is in `src/`.** `run_aggregate(solver, bio, soil,
  ...)` dispatches; `dt_max`/`dt_min` apply to `:split` and are accepted and
  ignored by `:stiff`, so one call site serves both. `postprocess_dataframe.jl`
  calls it and no longer carries its own branch or its own argument validation.
- **`n_grid` defaults agree**: `N_GRID_DEFAULT = 200` in `constants.jl`, used by
  both. It was 50 on the split side and 200 on the stiff side. Every call in
  `test/` and `paper/` passes `n_grid` explicitly — checked, 24 of 24 — so
  nothing in the suite changes behaviour.
- **Default output schedules agree**: `default_output_times(t_start, t_end)`,
  every `min(1, span/10)` days with the end point included. The stiff solver
  previously recorded exactly two points where the split solver recorded a
  trajectory.
- **`create_initial_state_legacy` deleted** — no caller in `src/`, `test/` or
  `paper/`, and it seeded a fixed background DOC of 1e-4 instead of partitioning
  measured SOC. `create_initial_state` and `EnvironmentalDrivers` are now
  exported, so the supported entry point is public.
- **Shared setup extracted**: `setup_run` builds grid, environment and state once
  for both solvers, which had identical copies of that block.

## ~~Structural deviations from Falconer~~ — **RECORDED 2026-07-30**

All three written up as decisions with their justification: REFERENCE §9a (table)
and the manuscript, §sec:fungal_dynamics (one paragraph). `F_i` turnover, ζ
applying to the reaction term only, and the conversion cost as an explicit flux.

---

## Smaller / hygiene

- `rm -rf _to_delete` — holds `_backup_20260727/` (11 pre-patch source copies
  that would silently revert the 2026-07-28 work if restored), plus the
  2026-07-29 `Project.toml` backup and the stray run logs moved there.
- `Pkg.rm("ADTypes"); Pkg.rm("SparseConnectivityTracer")` — added for automatic
  Jacobian sparsity detection, then made unnecessary when the pattern was built
  structurally instead (`mol_jacobian_prototype`). Still in `Project.toml`.
- `commit_session.sh`, `degryze_run.txt` untracked in the repo root.
- `dev_notes/` is gitignored but `REFERENCE.md` and `CLAUDE.md` name files in it
  as required reading (`falconer_answers.md`, the MATLAB notes). A fresh clone
  will not have them. Decide: track, or stop referencing.
- ~~`run_simulations.jl` repaired but not re-run~~ — the
  whole folder was archived 2026-07-29 to
  `paper/_archive/simulations_20260729/`. Nothing in it is to be quoted or
  re-run; the folder is to be rebuilt from scratch once De Gryze yields a trusted
  parameter set.
- Sauter 1926 *Forschungsheft* number to verify against a library record before
  submission.
- `O2_saturation` — the one function in `src/` whose job is ambient O₂ — is called
  by nothing outside its own test.
- Speed, if it is ever wanted: measured 2026-07-30, one aggregate is 1.16 s and the
  sparse linear solve is 33 % of self-time (`klu_l_refactor` 17.9 %,
  `klu_l_solve` 15.2 %). Arithmetic hoisting in the RHS cannot pay — the three
  candidates measured 3.0 %, 0.5 % and ~0 %. The levers are whether `n_grid = 200`
  resolves anything `n_grid = 100` would not, which is a **modelling** question and
  the cheapest possible win, and whether the block-tridiagonal structure warrants a
  specialised solver instead of general sparse KLU. Neither is urgent at 12 s per
  sweep.
- The respiration rate shows a birth-effect-like feature (owner, 2026-07-29,
  parked). Not diagnosed.
- `finite_volumes.jl` — `compute_cell_volumes` has no test and no caller.
  `is_diagonally_dominant` is tested but called from nowhere in `src/`.
- Test audit items not done (`docs/_archive/TEST_AUDIT_2026-07-28.md` §6): the
  217-assertion `total >= aggregate` block, `compute_r_agg` boundary tests
  (G_c = 0.0194 is derivable, so this is cheap and high value), and the
  order-of-accuracy refinement study — which is expected to show first order,
  contradicting the "2nd-order" claim at `timestepper.jl:60`.

---

## Deferred by decision

- ~~**Sieve-confinement `G_c(r)`**~~ — **IMPLEMENTED 2026-07-29** as
  `G_c(r) = G_c(δ_s)·(r/δ_s)^p_Gc` in `critical_binding`. `p_Gc = 0` is the
  package default, so nothing predating it changed. What is NOT settled: the
  exponent is fitted, not derived, and the δ_s curvature limit it stands in for
  is still a stated limitation rather than a solved closure. See REFERENCE.md
  §4.4a. Open question moved here: **what value of `p_Gc`**, and whether the De
  Gryze MWD series can identify it separately from `κ_b`.
- **Two particle-size scalings** (one for EPS, one for hyphae). Rejected: adds a
  parameter, and the single `d_32` scaling covers the observed direction.
