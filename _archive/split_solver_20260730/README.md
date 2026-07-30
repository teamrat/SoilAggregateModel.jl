# The split solver, archived 2026-07-30

Strang splitting + Crank–Nicolson (`run_aggregate`), the original integrator and,
until today, the reference implementation kept to validate `run_aggregate_stiff`.

## Why it is gone

It was kept for one reason: it accumulates CO₂ independently of the pools, so its
carbon balance is a *measurement* where the stiff solver's is an identity. On
2026-07-30 the 45-day agreement run showed that reason does not survive contact
with the thing it was supposed to validate.

Against the stiff solver at 45 days, refining its `dt_min` by 10× moved seven of
eight fields toward the stiff answer — B by 28×, E by 53× — while **DOC moved
away, 3.2e-2 → 9.2e-2**. Chasing that exposed four defects, none of which is a
step-size effect and none of which refinement removes.

That divergence is evidence about *this* implementation, not about the stiff
solver. It is why the file below was archived; it is not a finding about DOC.

1. **Rate laws evaluated on negative arguments.** `reaction_step.jl` passed the
   raw state to `compute_source_terms`. A negative C gives a negative Monod term
   — bacteria gaining carbon from nothing — and the state is only corrected
   afterwards. `mol.jl` clamps the arguments with `max(u, 0)` instead.
2. **Clipping folded into CO₂.** Any negative pool is zeroed and the carbon
   credited to respiration, so the mass balance closes whether or not the
   clipping was physically right. `clip_carbon` was computed and never reported,
   so the size of the artefact was unknowable.
3. **Lagged coefficients.** θ and the diffusivities came from the state at the
   start of the step and were held across both diffusion halves and the reaction.
4. **First order, not second.** Strang splitting is second order only if each
   sub-step is; `reaction_step!` is Forward Euler. The "2nd-order" claim at
   `timestepper.jl:60` was never tested.

On every accuracy axis the stiff solver is the better formulation. Repairing all
four would have produced a second well-converged discretisation whose
distinguishing feature — the cheap, auditable CO₂ integration — the repairs
themselves destroy.

## What replaces it

Nothing yet, deliberately. **A version carrying cumulative respiration as an
explicit state will be written as a separate test case for the mass-balance
check.** Carried per node rather than as a global scalar: a single CO₂ scalar
makes the Jacobian row dense across every node, and the sparse linear solve is
already 33 % of runtime, while a per-node field keeps the block-tridiagonal
structure and takes the block from 8×8 to 9×9.

Until that exists the integrated mass-balance check does not exist.
`test/test_closure.jl` still asserts the closure identity pointwise, which is the
stronger test of the two — exact, no integration, and it localises a failure to a
node and a state (REFERENCE §17a).

## Contents

`src/` — `timestepper.jl`, `reaction_step.jl`, `diffusion_step.jl`,
`crank_nicolson.jl`, `tridiagonal.jl`, and `finite_volumes.jl` (which had no
caller in `src/` at all).

`test/` — `test_tridiagonal.jl`, `test_crank_nicolson.jl`, `test_reactions.jl`,
`test_timestepper.jl`.

`verify_solver_agreement.jl` — the two-solver comparison. Its 22-day and 45-day
results are recorded in REFERENCE §17a; with one solver there is nothing left to
compare.

Nothing here should be run or quoted.
