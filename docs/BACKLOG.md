# Backlog

**The only list of open work.** Updated as items move; **not** a journal —
history lives in `docs/REFERENCE.md` §26 and the archived audits.

Two rules.

**An item is closed by pointing at code or at the manuscript, never at another
note.** Breaking it cost a day.

**Do not chase or enforce parameter values** (owner, 2026-07-30). Everything is up
for grabs. Where a reference exists it records what someone else used or measured
— these are not physical constants, and a deviation from one is not a defect. The
first full De Gryze calibration produces the first cut for anything uncited.
Model-architecture consistency is the standing priority; parameter items are
recorded for the calibration, not raised as problems. Three times on 2026-07-30 a closed
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

## 5. Group B parameters differ from their published values — **NOT a defect**

`D_B_rel` 1000×, `ω_E` 1500×, `ζ` 20×, `r_B_max` 3.7×, `μ_B` 3.3×, the Falconer
`α`/`β` set 4–5× low with mobilisation off. `REFERENCE.md` §5a, Group B.

**Reclassified 2026-07-30.** The earlier text said "return these to their anchors
before any fitting". That treated published values as constraints, which they are
not — they record what other authors used or measured. Carried as context for the
calibration, not as work.

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

## ~~17. The outer O₂ boundary is frozen at its `t_start` value~~ — **FIXED 2026-07-30**

`mol_rhs!` computes `O2_gas = p.O2_func(t)` and `O_amb = O2_gas / Tc.K_H_O` every
evaluation (`mol.jl:183,186`) and **discards both**. The outer node's value is set
once in `mol_solve.jl:96` and `du[mol_sid(n, MOL_O)] = 0.0` holds it there.

**Requirement, and the whole of it: the outer boundary reads `O2_func(t)` at time
`t`.** How that function evolves — step, ramp, interpolated table — is the
supplier's business. The model needs no concept of it, and no assumption about
smoothness.

**Fixed by a ghost-node Dirichlet.** `mol_outer_dirichlet(u_n, D_n, r_n, h, value)`
in `mol.jl`: a ghost at `r_n + h` carrying `2·value − u_n`, so the face at
`r_n + h/2` carries `value` exactly. Node `n` is now an ordinary state whose outer
face exchanges with the ambient. Evaluated pointwise — no derivative, no mass
matrix, no new parameter, state layout and Jacobian sparsity unchanged.

**A third thing the pin was hiding:** `du[outer] = 0.0` also discarded `S_O` at
node `n`, so respiration there consumed no oxygen. It does now.

Tested with a **step** driver — the case a derivative-based implementation would
get wrong — plus a constant-driver check that the outer node holds the ambient
rather than drifting, and that the profile falls inward, so the face is a boundary
and not a wall.

Two consequences, both wanted: `state.O[n]` becomes a state that tracks the
ambient rather than a frozen initial value, and the flux across that face is
available as the domain's O₂ sink — to be used, if ever, by post-processing.
**Nothing about that belongs in this model.**

*Rejected: finite-differencing `O_amb(t)` to drive `du[n]`. It makes the model
depend on the supplier being smooth, which is exactly the coupling to avoid.*

Not a parameter question. Checked at the same time and **correct**: `T` and `ψ`
are honoured per evaluation — `update_temperature_cache!` is called with
`T_func(t)` inside the RHS, not cached from `t_start`.

---

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

## 13. REFERENCE.md internal contradictions — **VERIFIED 2026-07-30, all four real**

Checked against the code, not restated. They are **not four independent errors**:
each is an early table contradicting a later section, and in every case the later
section is right. §3.2, §3.4, §7 and §10 predate the softplus regularisation and
the parameter settling; §5a and §17a are current.

| § | says | code says | |
|---|---|---|---|
| 3.2 | `F_i_min = 1.0×10⁻⁴` | `parameters.jl:138` → **1.0e-6** | REFERENCE's own line 567 says 1e-6 too, so it contradicts itself *and* the code, by 100× |
| 3.4 | `K_E = 100·K_B/5` → 2e-3 | `parameters.jl:181` → **1.0e-3**, derived as `50·C_B` | two different derivations; line 560 has the right value with the wrong provenance |
| 7 | hard branch: `R_B > R_Bb` → `Γ_B = Y_B·R_B·(1−γ)`; else `Γ_B = R_B, Γ_E = 0` | softplus of the difference, and `Γ_B` carries a `−σ(−x)` term §7 omits | §7 describes the discontinuous form the regularisation replaced. §17a has the true one |
| 10 | `Resp_B = (1−Y_B)·R_B`, "Starvation: 0" | `Resp_B = R_Bb + σ(x)(1−Y_B)` | §10 drops maintenance and claims zero respiration under starvation, when the code keeps `R_Bb` — the physically meaningful part |

**So the fix is not four patches.** The early tables should be regenerated from
`parameters.jl` and from §17a's proof, which is the section that was written
carefully. Patching line by line reproduces the stratum.

Also §4.2 names `clay_fraction` and `silt_fraction` as fields; neither exists. The
`reactions.jl` citation cluster is off by roughly 12 lines.

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

## 16. Documentation survey — **DONE 2026-07-30. The plan was wrong.**

I proposed archiving `dev_notes/` and six orphans. Reading them first — which is
what the item was blocked on — changed the answer.

### Six untracked files are load-bearing

`dev_notes/` is gitignored, and six files in it are cited from `src/`, from
`README.md` and from `REFERENCE.md`. **A fresh clone gets the citations and not
the sources.**

| file | cited from |
|---|---|
| `falconer_answers.md` (847 lines) | `src/parameters.jl` ×2, `src/biology/fungi.jl`, REFERENCE ×4, GUIDE, CLAUDE |
| `MATLAB_aggregation_logic.md` (195) | `src/postprocessing/aggregate_radius.jl`, REFERENCE, CLAUDE |
| `MATLAB_parameters.md` (385) | REFERENCE, CLAUDE |
| `falconer_questions.md` (202) | **`README.md`** — the front page |
| `claude_code_instructions_aggregate_stability.md` (210) | REFERENCE |
| `archive/REPO_REORGANIZATION.md` (262) | REFERENCE, one rename note |

The move is **promote, not archive**: track these, or delete the citations.
Archiving as originally planned would have broken seven references from `src/`
and the front page.

### Two files would have taken live information with them

- **`patches/O2_state_variable_change.md`** — the only record of why `state.O` is
  aqueous rather than total soil O₂. REFERENCE mentions the consequence twice in
  passing and never the decision. Uncited, and on the original plan I would have
  archived it.
- **`docs/julia_falconer_deviations.md`** — cited from REFERENCE ×2, documents
  live behaviour (`D_eff_fungi_mobile`'s network-dependent diffusivity and its
  Falconer/MATLAB provenance). Internally stale: names `src/physics/fungi.jl` and
  `src/physics/reactions.jl`, neither of which has ever existed. Fold into
  REFERENCE §9a, then retire.

### Two folders nobody has mentioned

- **`scripts/`** — 7 Julia files, 867 lines, February, plus a README. **Every one
  references archived API**; two are split-solver by name. None can run.
  `example_1year_simulation.jl` (226 lines) is the only one whose intent may be
  worth reviving.
- **`patches/`** — 3 patch fragments, 80 lines, February, plus the O₂ document.

### Uncited and superseded — safe to archive (13 files, ~3,900 lines)

`DEVELOPMENT_ROADMAP.md`, `PHASE1_SPEC_1.1-1.3.md`, `PHASE2_SPEC.md` (plans,
executed); the four remaining `dev_notes/archive/` files; `paper_framing.md`
(explicitly "not meant to be direct instruction"); `manuscript_population_upscaling.md`
and `manuscript_stability_alignment.md` (merged — their target sections exist in
`manuscript-4-5.tex`); `O2_debugging_handover.md` and `CLAUDE_CODE_INSTRUCTIONS.md`
(superseded by CLAUDE.md); plus `scripts/` and `patches/`.

`manuscript_changes.md` (root, cited once from `test_postprocessing.jl:337`)
overlaps REFERENCE §26 errata — merge there rather than archive.

### Revised target

Not "36 → 4". **Promote 6, fold 3, archive 23.** Nothing is destructive until you
say so.

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
