# Backlog

Working state. Updated as items move; **not** a journal — history lives in
`docs/REFERENCE.md` §26 and the audit documents.

Last commit: `54efc51` (2026-07-28). This file covers the 2026-07-29 solver work.

---

## Dependency order

**Speed is DONE** (item 1). Remaining: **anchors → decisions → fit.**

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
the n=80 coarse grid `optimize_soil3.jl` uses.

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

Also: `optimize_soil3.jl` has not been migrated to `degryze_soils.jl` — it still
hard-codes `ρ_bulk = 1300`, `SOC = 0.0221`, `ψ = −29`, `f_clay_silt = 0.74` for
soil 3. Its stored fit is invalid regardless, having been made against the
broken POM normalization.

## 3. `M_max` has two values, 29× apart — needs a decision

`initial_conditions.jl:382` computes `M_max = k_ma·f_clay_silt·ρ_b = 288` and
partitions the initial MAOC pool against it. `reactions.jl:152` evolves that
pool against `soil.M_max = 10`. Every run therefore opens with MAOC far above
the isotherm ceiling and produces a spurious desorption pulse. The guard at
`initial_conditions.jl:418` compares against 288 and never fires.

`maoc.jl:44` states `M_max = k_ma·f_clay_silt·ρ_b` as definitional, which would
make 288 intended — but the units do not close as written (`k_ma` is per **gram**
of mineral, `ρ_b` is µg/mm³). Needs the intended definition and its units.
`REFERENCE.md` §5a.

## 4. `delta` — and it is coupled to `beta`

Default is 1.0 (MATLAB linear case); Falconer runs θ = 1.0 as a control and
θ = 3.0 as the nonlinear case. Three values have been in play historically
(2.0 in the docs until 2026-07-27, 1.0 in code, 3.0 in Falconer) and nobody
recorded why 1.0 was chosen.

**δ barely does anything while β = 0.** The exponent sits on the gain, `α·Π^δ`,
and with mobilization disabled there is no competing linear loss for it to cross
over against. Raising δ alone rescales the transition rate rather than producing
the threshold behaviour Falconer describes. **Decide δ and β together.**
`REFERENCE.md` §5a.

## 5. Group B parameters off their cited anchors

`D_B_rel` 1000×, `ω_E` 1500×, `ζ` 20×, `r_B_max` 3.7×, `μ_B` 3.3×.
`REFERENCE.md` §5a, Group B. Return these to their anchors, or record why not,
**before** any fitting.

## 6. `r_agg` step function

`pct_gt2000um` jumps between day 2 and day 5 (`f_agg` 0.010 → 0.249) and then
sits flat. May be downstream of 1, 3 or 5 rather than its own defect — revisit
after those.

## 7. CO₂ overshoot

3599.7 vs measured 2139 µg-C/g at day 21 for soil 3, about 1.7×. Both sides are
total respiration (spec §0a A3), so the discrepancy is real rather than a
partitioning artefact.

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
- `paper/simulations/single_aggregate_physics/run_simulations.jl` was repaired
  but **not re-run**; its `data/` and `.log` predate the break.
- Sauter 1926 *Forschungsheft* number to verify against a library record before
  submission.
- `finite_volumes.jl` — `compute_cell_volumes` has no test and no caller.
  `is_diagonally_dominant` is tested but called from nowhere in `src/`.
- Test audit items not done (`docs/TEST_AUDIT_2026-07-28.md` §6): the
  217-assertion `total >= aggregate` block, `compute_r_agg` boundary tests
  (G_c = 0.0194 is derivable, so this is cheap and high value), and the
  order-of-accuracy refinement study — which is expected to show first order,
  contradicting the "2nd-order" claim at `timestepper.jl:60`.

---

## Deferred by decision

- **Sieve-confinement `G_c(r)`** — the radial scale dependence. Partly addressed:
  texture now enters through `d_32`. What remains is the δ_s curvature limit,
  recorded as a stated limitation in the manuscript rather than modelled. The
  MATLAB precursor had `strength ./ x` commented out; no record of why.
- **Two particle-size scalings** (one for EPS, one for hyphae). Rejected: adds a
  parameter, and the single `d_32` scaling covers the observed direction.
