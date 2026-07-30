# Test suite audit — 2026-07-28

**Question asked:** are the tests meaningful? Not "do they pass" — for every
assertion, *what specific wrong code would make this fail?* An assertion that
cannot answer that is not a test.

**Method.** Six independent auditors, one per group of test files, each reading
the test *and* the source it covers. The numerical-kernel auditor additionally
re-implemented the Crank–Nicolson and Thomas routines in Python and ran mutation
experiments. Findings marked **(verified)** were reproduced independently before
being recorded here.

**Provenance warning.** The first run of this audit was **invalid** and was
discarded. Two causes, both worth recording because they will recur:

1. The file-staging bridge returned **cached copies** of the repository from
   earlier in the session while reporting fresh timestamps. Every staged file's
   checksum differed from disk. Detected only because the auditors' "blockers"
   contradicted a passing test run. **Always checksum a staged snapshot against
   the source before trusting it.**
2. `test_output.txt`, a stale run log at the repository root showing the old
   1219/62/11 result, was mined by several auditors and its five-month-old
   failures reported as current. That file was flagged for deletion in
   `STRUCTURAL_AUDIT_2026-07-27.md` item A3 and had not been removed. It is now
   in `_to_delete/`. **Stale logs in a repository are not clutter; they are a
   source of false findings.**

The second run used a freshly-named tarball snapshot, checksum-verified against
disk, with no logs included.

---

## 1. Headline

The suite reported **1364 assertions**. That number is inflated by loops, and the
loops mostly replay assertions that cannot fail.

| Group | Assertions | Structure |
|---|---|---|
| Post-processing | 432 | **307 (71 %) come from five written lines.** One of them is 217 assertions — 50 % of the file — and cannot fail |
| User API | 256 | 82 % from 16 predicates; the non-negativity family is 120 (47 %), of which 48 are strictly unfalsifiable |
| Temperature | 135 | 90 (67 %) are one `> 0` / `isfinite` loop over nine temperatures |
| Parameters + Types + Environment | 122 | 41 (34 %) load-bearing; 32 cannot fail at all |
| Biology + Reactions | 107 | ~35 vacuous or tautological |
| Solver kernel (tridiag + CN + POM) | 58 | see §3 — the most serious result in the audit |

A count of assertions is not a measure of coverage. Where a loop asserts one
predicate at every grid node, it produces one bit of information, not thirty.

## 2. The three genres that account for nearly all the dead weight

**Asserting what the code guarantees.** `@test all(rec.state.F_i .>= 0.0)` where
`reaction_step.jl:118-129` explicitly clips those pools and `diffusion_step!`
never touches them. `@test all(M_eq .<= M_max)` where `M_eq = M_max·x/(1+x)`.
`@test soil.θ_r <= θ <= soil.θ_s` where `van_genuchten` returns
`θ_r + (θ_s−θ_r)/(≥1)` for any α and n. These pass for every possible
implementation of the physics.

**Restating the implementation.** `@test aq.O_aq[i] ≈ rec.state.O[i]` pins an
assignment statement, thirty times — while the same convention was violated 150
lines below in the same file, untested (§4). `@test α_eff ≈ soil.α_vg *
exp(soil.ω_E*E + soil.ω_F*F)` is `alpha_effective`'s single line retyped. A test
that recomputes the code's formula confirms the arithmetic and says nothing
about whether the formula is right.

**Tolerances with orders of slack.** `@test state.P > P_initial * 0.9` where the
real loss is 8.6e-6 — an **11 500×** margin. `@test abs(J_M_eq) < 0.1` against a
derivable exact value of 3.466e-4 — **289×**. `@test grad_outer < 0.01` against
an actual 7.2e-6 — **1400×**. Each was written as a sanity bound and functions as
an unfalsifiable one.

## 3. The most serious finding: the diffusion solver is barely constrained

`test_crank_nicolson.jl` had 11 assertions. The auditor re-implemented the
scheme and mutated it. **(verified)** Ten of twelve semantic mutations pass all
eleven assertions:

| Mutation | Suite verdict |
|---|---|
| face `D` → harmonic mean | survives |
| face `D` → upwind (left node) | survives |
| face `D` → doubled | survives |
| face radius → left node | survives |
| face radius → geometric mean | survives |
| geometry `r²` → `r³` | survives |
| geometry `r²` → `r¹` | survives |
| **geometry `r²` → `r⁰` (spherical geometry removed entirely)** | **survives** |
| **boundary flux × 2, × 0.5, × 1e-6, × 100** | **survives** |
| θ 0.5 → 1.0 | survives |
| outer face radius mismatched (breaks telescoping) | caught |
| θ 0.5 → 0.0 (unconditionally unstable) | caught |

Why the mass-conservation assertions are nearly free: with `W_i = 4π r_i² h` the
sum `Σ W_i L_i` telescopes for **any** face values whatsoever. It therefore
catches exactly one defect class — the same face evaluated differently in the
two cells that share it — and nothing about the coefficients.

Why the varying-`D` test is empty: it initialises to `fill(1.0, n)`, which is an
exact fixed point of the assembly (`a_i + b_i + c_i = 1`), so `D` never enters.
**(verified)** the state moves by 2e-15 over ten steps. The testset's comment
says it tests face-averaging. It cannot.

The boundary flux — the model's only carbon inlet — was constrained only by
`mass_final > mass_initial`. Note that `crank_nicolson.jl:64` documents the ghost
node as `u[0] + 2h·J/D₀` while the code (`:110`, `:202`) uses `h·J/D[1]`. A factor
of two, in the source's own docstring, with no test able to arbitrate.

`test_tridiagonal.jl:10-17` quotes an architecture review listing three required
analytic tests, with the closed-form solutions given verbatim. Only one was ever
implemented, in the vacuous form above.

## 4. Two real source defects the audit surfaced

Both were found by asking what the tests *fail* to catch, and both were verified
against the source before being acted on.

### 4a. Post-processing disagreed with the solver — **FIXED**

`reactions.jl:77` is `O_aq = O`: since the O_total→C_aq switch the state variable
*is* the aqueous concentration, and the gas–water capacity `θ + K_H·θ_a` appears
on the O₂ **source term** (`reactions.jl:166`), not on the concentration.
`derived.jl:59` documented this correctly. But `derived.jl:205` and `:309` —
inside `respiration_rates` and `carbon_use_efficiency` — re-applied
`O·θ/(θ + K_H·θ_a)`. With `K_H ≈ 34` and θ ≈ 0.4 that is roughly **a tenth** of
the oxygen the solver used, pushing the O₂ Monod factor from ~0.9 to ~0.5.

**Every respiration and CUE figure read out of post-processing was therefore not
the quantity the simulation ran on.** Two further divergences in the same
functions: `C_aq` used `(1 + k_d·ρ_b)` where the solver uses `(θ + ρ_b·k_d)`, and
`ζ_T` was not clamped to 1 after the Arrhenius factor.

The oracle was already in the file and discarded: `derived.jl:188` computed
`src = compute_source_terms(...)` and threw the result away, with
`# Let me check what compute_source_terms actually computes...` left in shipped
source at line 193.

All three divergences are fixed, the dead call and the stray comment removed, and
`test_postprocessing.jl` now asserts `respiration_rates` against
`compute_source_terms` node by node at `rtol=1e-12`.

### 4b. `M_max` has two values, 29× apart — **NOT fixed, needs a decision**

`initial_conditions.jl:382` computes

    M_max = soil.k_ma · soil.f_clay_silt · soil.ρ_b = 0.48 × 0.40 × 1500 = 288.0

and partitions the initial MAOC pool against it. `reactions.jl:152` evolves that
same pool against `soil.M_max`, whose default (`parameters.jl:275`) is **10.0**.
So every run opens with MAOC far above the isotherm ceiling the solver enforces,
and the first timesteps produce a spurious desorption pulse. The `@warn` guard at
`initial_conditions.jl:418` compares `M_0` against 288, so it never fires.
`maoc.jl:44` states `M_max = k_ma · f_clay_silt · ρ_b` as definitional.

**This is not fixable without a decision, and I have not made one.** Neither value
is obviously right: the unit string on `k_ma` is µg-C per *gram* of mineral while
`ρ_b` is µg/mm³, so the product as written is not dimensionally clean either.
Recorded as an open decision in `REFERENCE.md` §5a. Deciding it by picking the
value that looks better would be exactly the failure mode this project has
already had once.

## 5. What was added

All expected values were computed independently before being written into a test.

| Test | What it pins | Verification |
|---|---|---|
| `test/test_tessellation.jl` (new) | The module had **zero** coverage despite being exported and despite ω scaling every background pool. Central assertion is the distribution-free identity `total_POM_C/(V·ρ_b) == I_input`, which holds exactly for any diameters and fractions and simultaneously pins the `f_pack` exponent, both radius/diameter halvings and both `(4/3)π` factors | recomputed for three unrelated distributions, exact to 1e-10 |
| CN: exact discrete flux balance | `ΔM = 4π r₀² J Σdt_half` — replaces `mass_final > mass_initial`, which passed for a flux scaled by 1e-6 or 100 | independent reimplementation, agreement 6.8e-13 |
| CN: analytic steady state | `u(r) = u_out + (J r₀²/D)(1/r − 1/r_max)` under flux-in / Dirichlet-out. This is requirement (2) from the architecture review quoted in `test_tridiagonal.jl`'s own header, specified and never implemented. Also the first test to exercise the `theta = 1.0` path | errors 1.256e-2 (n=50), 3.190e-3 (n=100), 7.970e-4 (n=200); ratios 3.94, 4.00 — second order |
| Thomas vs dense solve | Every prior random-matrix test compared the solver to itself: `thomas` is a copy-wrapper around `thomas!`, and `thomas_factorize!`+`thomas_solve!` run the identical recurrence, so a common-mode error was invisible | `LinearAlgebra.\\` is an independent implementation |
| `van_genuchten_inverse` round-trip | Had **zero** tests despite being the entry point for every WFPS-specified experiment. Round-trip needs no expected values — it is its own reference | — |
| Post-processing oracle | `respiration_rates` vs `compute_source_terms`, node by node. Would have caught §4a on the day it was introduced | fails before the §4a fix, passes after |

## 6. Not done, and why

- **The 217-assertion `total >= aggregate` block** was left in place. Replacing it
  properly means asserting `C_agg` against a hand-summed node set and pinning
  `r_agg`, which is a larger change than an audit should make in passing.
- **`compute_r_agg` boundary tests.** `G_c = τ_w/k_F` is derivable from cited
  constants (`τ_w = √2·μ·v_s/δ_s = 0.0437 Pa`, `k_F = 2.25` → `G_c = 0.0194
  µg/mm³`), and the right test builds a synthetic `OutputRecord` with binding
  exactly at `G_c`, just below, just above, and a non-monotonic profile. Highest
  remaining value; not attempted here because it needs the struct constructed by
  hand and I could not run Julia to check it.
- **Order-of-accuracy study.** The auditor's reading is that the composite scheme
  **cannot** be second order as `timestepper.jl:60` claims: the reaction sub-step
  is Forward Euler, the O₂ diffusion sub-step is Backward Euler, and the POM flux
  is frozen at the start of the step. A dt-refinement test would show slope ≈ 1
  and falsify the docstring. That is a real finding but changing the claim is a
  decision about the method, not a test fix.
- **Deleting weak assertions.** Nothing was deleted. Strong tests were added
  alongside. Removing 300 assertions is a separate, reviewable change.

## 7. Open items worth recording

- `paper/simulations/single_aggregate_physics/run_simulations.jl` was repaired on
  2026-07-28 (it read `trans.trans_i` / `trans.insulation`, gone since February)
  but **has not been re-run**; its `data/` and `.log` files predate the break.
- `src/solver/finite_volumes.jl` — `compute_cell_volumes` has no test and no
  caller anywhere in `src/`. Its own header concedes the volumes are inconsistent
  with the stencil's conservation weights.
- `is_diagonally_dominant` is tested but called from nowhere in `src/`.
- `Y_F_const` and `Y_F_uptake_dependent` are tested; `Y_F_func`, the one the
  solver actually calls, is not.
- `src/parameters.jl` performs **no validation at all** — no `throw`, no
  `@assert`. Every bound in `test_parameters.jl` is asserted against the one
  instance that cannot vary (the default), while user-supplied values, the only
  ones that can be wrong at runtime, are unchecked. Documented keyword overrides
  that construct happily and produce `NaN` downstream include `K_B=0`, `ε_F=0`,
  `ε_Y=0`, `B_min=0`, and `n_vg=1`.
- Four derived-parameter relations are documented in `parameters.jl`
  (`C_B = K_B/5`, `K_Y = 10·r_B_max·C_B/(K_B+C_B)·B_min`, `ε_Y = K_Y/100`,
  `K_E = 50·C_B`). All hold for the defaults; none is tested; a keyword override
  of `K_B` silently breaks all four.
- The mobilization regime (`net < 0`) is never entered through the integrated
  reaction path, because the defaults set `β_i = β_n = 0`. The `abs()` in
  `Resp_F_conv` that `fungi.jl` calls CRITICAL is exercised by exactly one direct
  unit call.
- Every Arrhenius factor is 1.0 in every reaction-layer test, so the temperature
  plumbing in `compute_source_terms` — seven multiplications and the `ζ_T` clamp —
  is untested.
