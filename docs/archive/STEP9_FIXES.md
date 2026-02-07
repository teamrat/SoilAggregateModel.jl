# Step 9 Fixes — Instructions for Claude Code

These fixes address three issues found during audit. Apply them in order.
All 523 existing tests must continue to pass after each fix.

---

## Fix 1: Machine-precision carbon conservation (CRITICAL)

### Problem

The 30-day carbon conservation test shows ~2e-5 relative error. The architecture (§16) requires machine precision (~1e-14). The carbon leak comes from non-negativity clipping in `reaction_step!`: when `max(0.0, state.X[i])` truncates a negative value to zero, the clipped carbon vanishes — it is not added to CO₂ or any other pool.

### What to change

**File: `src/solver/reaction_step.jl`**

Replace the current clip-and-forget pattern (lines 76–84) with clip-and-track. For each carbon species (C, B, F_n, F_m, F_i, E, M — NOT oxygen), if the value is negative after Forward Euler, record the clipped amount × volume element and add it to `CO2_cumulative` as "lost carbon". This makes the carbon balance exact by construction.

Specifically, change the non-negativity block to something like:

```julia
# Enforce non-negativity — redirect clipped carbon to CO₂
volume_i = 4.0 * π * r_grid[i]^2 * h  # (already computed below for CO₂)

# For each CARBON species: clip and track
for (val_ref, field) in ((state.C, i), ...)  # pseudocode — expand manually
    if state.X[i] < 0.0
        state.CO2_cumulative += abs(state.X[i]) * volume_i
        state.X[i] = 0.0
    end
end
```

Concretely, for each of the 7 carbon species (C, B, F_n, F_m, F_i, E, M):
```julia
if state.C[i] < 0.0
    clip_total += abs(state.C[i]) * volume_i
    state.C[i] = 0.0
end
```
(Repeat for B, F_n, F_m, F_i, E, M.)

Then after the loop over species:
```julia
state.CO2_cumulative += clip_total
```

For oxygen (O): clip to zero but do NOT add to CO₂ (oxygen is not carbon).

Move the `volume_i` computation ABOVE the clipping block (it's currently below it on line 87).

### Expected result

The 30-day carbon conservation test should now pass at rtol=1e-12 or better. Update the test tolerance on line 196 of `test_timestepper.jl`:
```julia
@test C_total_final ≈ C_total_initial rtol=1e-12
```

Also add a `clip_cumulative` field to `AggregateState` (or to the diagnostics dict) so we can monitor how much carbon was clipped. This is recommended by Architecture §2 ("Track cumulative clipping magnitude as a diagnostic").

---

## Fix 2: Eliminate redundant source term computation (PERFORMANCE)

### Problem

`adapt_timestep` (timestepper.jl lines 194–221) recomputes `compute_source_terms` at every node — the exact same computation `reaction_step!` just performed. This doubles the per-step cost. For 400 aggregates × 10 years, this matters.

### What to change

**Option A (simpler, recommended):** Compute `max_rel_change` inside `reaction_step!` and return it. Then `adapt_timestep` just applies the halve/double logic without recomputing sources.

Modify `reaction_step!` to:
1. After computing `sources` at each node (line 63), compute the per-node relative change:
```julia
threshold = 1e-6
rel_C = abs(sources.S_C * dt / max(C_i, threshold))
rel_B = abs(sources.S_B * dt / max(B_i, threshold))
# ... (all 8 species)
node_max = max(rel_C, rel_B, rel_Fn, rel_Fm, rel_Fi, rel_E, rel_M, rel_O)
max_rel_change = max(max_rel_change, node_max)
```
2. Change the return type: `return max_rel_change` instead of `nothing`.
3. Initialize `max_rel_change = 0.0` at the top of the function.

Modify `adapt_timestep` to accept `max_rel_change::Real` instead of recomputing it:
```julia
function adapt_timestep(max_rel_change::Real, dt::Real, dt_min::Real, dt_max::Real)
    dt_new = dt
    if max_rel_change > 0.10
        dt_new = dt / 2.0
    elseif max_rel_change < 0.01
        dt_new = dt * 2.0
    end
    return max(dt_min, min(dt_max, dt_new))
end
```

Update the call in `run_simulation` (line 121 and 131):
```julia
max_rel = reaction_step!(state, workspace, dt, r_grid, h, bio, soil, ψ)
# ...
dt_new = adapt_timestep(max_rel, dt, dt_min, dt_max)
```

### What to verify

- `reaction_step!` still updates state correctly (existing tests pass).
- `adapt_timestep` still returns sensible values (Test 3 in test_timestepper.jl).
- Update Test 3 to call the new `adapt_timestep` signature. Since the direct-call test now needs a `max_rel_change` value, either:
  (a) Call `reaction_step!` first and use its return value, or
  (b) Pass a known `max_rel_change` value and check the dt logic.

---

## Fix 3: Hardcoded constants in workspace_updates.jl (CLEANUP)

### Problem

`workspace_updates.jl` has three hardcoded values (lines 32, 44, 51) that should come from `constants.jl` or from the parameter structs:
- `T_ref = 298.15` (line 32)
- `D_DOC_ref = 5.0e-6` (line 44) 
- `D_O2_a_ref = 0.2` (line 51)

### What to change

**Option A:** If these constants are already in `constants.jl`, import them. Replace:
```julia
T_ref = 298.15
```
with a reference to the constant (e.g., `T_REF` from constants.jl).

**Option B:** If they're not in `constants.jl`, add them to `BiologicalProperties` or `SoilProperties` as fields, or add them to `constants.jl`. Then reference them here.

**Option C (minimal):** Add a comment explaining these are physical constants that don't change:
```julia
T_ref = 298.15  # K — reference temperature (universal, matches constants.jl)
```

Option A or B is preferred. Option C is acceptable if the others require too many file changes.

For `D_DOC_ref` and `D_O2_a_ref` specifically: check whether `BiologicalProperties` or `SoilProperties` already has fields for these. If so, use `bio.D_DOC_ref` and `bio.D_O2_a_ref` (or whatever the field names are). If not, add them.

---

## Testing checklist

After all three fixes:

1. `julia --project test/runtests.jl` — all tests pass (should be ≥523)
2. 30-day carbon conservation test passes at `rtol=1e-12`
3. `adapt_timestep` test updated for new signature
4. No new NaN or Inf values
5. Verify the hardcoded constants match whatever `constants.jl` defines

---

## Do NOT change

- The Strang splitting order (diffusion–reaction–diffusion) — correct
- The dt/2 passed to `diffusion_step!` — correct  
- The J_P_val frozen for both half-steps — correct (Strang-consistent)
- The `deepcopy` for output snapshots — correct
- The post-hoc adaptive timestep strategy (accept-and-reduce) — acceptable design choice
- The output timing logic — works for dt_max=0.1 with output_interval≥1.0
