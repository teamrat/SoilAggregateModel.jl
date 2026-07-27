# O₂ Debugging Handover
# =====================
# Date: 2026-02-17
# Context: SoilAggregateModel.jl — O₂ state variable switched to C_aq, but O₂ still
# accumulates near POM surface exceeding atmospheric boundary value.


## What was done this session

### Bug 1 (FIXED): POM dissolution size-dependence reversed
- `timestepper.jl` used `bio.P_0` (global default) instead of per-aggregate `state.P_0`
- Added `P_0::Float64` field to `AggregateState`, set during initialization
- Verified: smallest POM now depletes fastest ✓

### Bug 2 (FIXED): O₂ IC/BC mismatch
- `O2_func(t)` returns gas-phase O₂ density (~0.278 μg/mm³)
- IC was set to total O₂ (~0.059), BC was gas-phase (~0.278)
- Fixed: BC now converts gas to appropriate units matching state variable

### Bug 3 (ATTEMPTED, PARTIALLY FIXED): O₂ accumulation near POM
This is the active bug. Three approaches were tried:

**Approach 1 — IRA correction (abandoned):**
Adjust O_total to preserve C_aq when θ changes between timesteps.
Result: O₂ still accumulated near POM. The correction addresses the temporal
artifact (θ changing in time) but not the spatial artifact (diffusing O_total
with spatially varying θ introduces spurious advection).

**Approach 2 — Switch state variable to C_aq (implemented, current state):**
Changed O₂ state variable from O_total to aqueous concentration C_aq.
Changes made in 3 files:

1. `timestepper.jl`:
   - Removed IRA correction
   - BC: `O2_amb = O2_gas / K_H` (was `O2_gas * (θ+K_H*θ_a) / K_H`)
   - POM dissolution: `O_aq_0 = state.O[1]` (was `state.O[1] * θ / (θ+K_H*θ_a)`)

2. `reactions.jl` (compute_source_terms):
   - `O_aq = O` (was `O * θ / (θ + K_H * θ_a)`)
   - `S_O = -α_O * Resp / (θ + K_H * θ_a)` (was `-α_O * Resp`)

3. `initial_conditions.jl`:
   - `O_aq = O2_gas / K_H` (was `O_total = O2_gas * (θ+K_H*θ_a) / K_H`)

**Result:** O₂ values are now in correct range (~0.0097 vs old ~0.059) but
O₂ STILL accumulates near POM over time:
- Day 14: O uniform at ~0.00967 (correct, ≈ O2_gas/K_H)
- Day 30: O at node 1 = 0.0143, O at node 200 = 0.00968 (48% above BC!)

**The O₂ accumulation near POM with only sinks and a Dirichlet BC is physically
impossible. There must be a remaining bug.**


## Current state of the code

### Files modified (relative to /mnt/project/ versions):
- `src/solver/timestepper.jl` — C_aq state variable, no IRA
- `src/solver/reactions.jl` — O_aq = O, S_O divided by capacity
- `src/solver/reaction_step.jl` — unchanged (applies S_O from reactions.jl)
- `src/solver/initial_conditions.jl` — O_aq = O2_gas / K_H
- `src/solver/diffusion_step.jl` — unchanged
- `src/solver/crank_nicolson.jl` — unchanged
- `src/types.jl` — P_0 field added to AggregateState
- `src/api.jl` — P_0 initialization, may have workspace θ init lines (check)

### Key verified facts:
- `compute_source_terms` returns S_O = -0.059 (negative, correct sign)
- `Resp_total` = 0.087 (positive, correct)
- θ is essentially uniform: 0.29899 everywhere (EPS effect on α is negligible
  at current E values ~0.08, by design)
- capacity = θ + K_H * θ_a ≈ 6.08 everywhere
- D_eff_oxygen already computes D̂_O / capacity (correct for C_aq diffusion)
- K_H ≈ 28.77
- O2_amb = O2_gas / K_H ≈ 0.00966


## Where to investigate

### Hypothesis 1: CO₂/O₂ coupling bug
The snapshot shows CO₂ column with large values (711 at day 30). CO₂ is
cumulative and not diffused. But check: is there any code path where CO₂
or respiration somehow feeds back into O?

### Hypothesis 2: Strang splitting asymmetry
The Strang splitting does: diffusion(dt/2) → reaction(dt) → diffusion(dt/2).
The reaction step depletes O₂. The diffusion step redistributes. If the
diffusion is somehow adding O₂ (through a BC inconsistency or coefficient
error), it could accumulate over many steps.

Test: run with O₂ diffusion disabled (set D_O to zero). O should only decrease
from reactions. If it still increases, the bug is in reactions. If it decreases,
the bug is in diffusion.

### Hypothesis 3: Dirichlet BC implementation
Check `crank_nicolson.jl` lines 211-218. For Dirichlet outer BC:
```julia
lower[n-1] = 0.0
diag[n] = 1.0
rhs[n] = value_outer
```
This forces u[n] = value_outer. But the RHS was already set on line 134:
`rhs[n] = u[n] + 0.5 * dt_half * L_explicit` with L_explicit = 0.0.
Then line 217 overwrites it. This should be fine, but verify.

Also: the explicit half of Crank-Nicolson (line 134) uses L_explicit = 0
for the Dirichlet node. But the adjacent node (n-1) uses u[n] in its stencil.
If u[n] was modified by reactions between the two diffusion half-steps, the
explicit stencil sees a different u[n] than the Dirichlet value. This could
create a small source/sink at node n-1.

### Hypothesis 4: Neumann inner BC (zero flux)
The inner BC for O₂ is zero flux: ∂C_aq/∂r = 0 at r₀. This means no O₂
enters through the POM surface. Verify the ghost node implementation doesn't
inadvertently create a flux.

### Hypothesis 5: Numerical artifact from spherical geometry
The spherical Laplacian (1/r²)∂/∂r[r² D ∂u/∂r] has a 1/r² prefactor.
At node 1 (r = 0.625 mm), this is 1/0.39 ≈ 2.56. At node 200 (r = 6.25),
it's 1/39 ≈ 0.026. The inner nodes have ~100× amplification of diffusive
fluxes. A small systematic error in the stencil could accumulate much faster
at inner nodes.


## Diagnostic approach

1. **Ablation test**: disable O₂ diffusion (set workspace.D_O .= 0 before
   diffusion_step!) and check if O₂ decreases monotonically from reactions
   alone. This isolates the problem to diffusion vs reactions.

2. **Conservation test**: compute total O₂ in domain before and after each
   diffusion half-step: `sum(state.O[i] * 4π r[i]² h for i in 1:n)`.
   Compare with expected boundary flux. Any discrepancy points to the CN solver.

3. **Steady-state test**: initialize with uniform C_aq = O2_gas/K_H, set all
   reactions to zero (or use tiny bio rates), run for 100 steps. O should
   remain uniform. If it develops a profile, the diffusion has a bug.

4. **Single CN step trace**: manually call one `crank_nicolson_step!` on the
   O₂ vector with a known flat profile and known uniform D. Verify the
   output is still flat (within machine precision). Then try with a small
   gradient and verify flux direction.


## Parameter context

Current test uses de Gryze silt loam with deliberately weak EPS effect:
```julia
bio = BiologicalProperties(
    κ_s_ref=0.001, κ_d_ref=0.0001,
    F_i_min=1e-6, F_n_min=2e-4, F_m_min=1e-6,
    D_Fn0=0.001, D_Fm0=0.01,
    r_B_max=0.5, r_F_max=0.5, R_P_max=0.1,
    C_B=5.0e-5, μ_B=0.02, μ_F=0.02
)
soil = SoilProperties(ρ_b=ρ_bulk, f_clay_silt=0.74)
```
ψ = -29 kPa, T = 293.15 K, O2_gas from ideal gas law (~0.278 μg/mm³).
5 POM bins: 0.65, 0.95, 1.25, 1.55, 1.85 mm diameter.
200 grid nodes per aggregate. ω = 19.6 dilution factor.


## Key files to read

- `src/solver/crank_nicolson.jl` — the CN solver, most likely location of bug
- `src/solver/diffusion_step.jl` — calls CN for all 5 species
- `src/solver/reactions.jl` — compute_source_terms, verify S_O
- `src/solver/timestepper.jl` — the main loop
- `src/physics/effective_diffusion.jl` — D_eff_oxygen formula
- `O2_state_variable_change.md` — derivation of C_aq approach (in outputs)
