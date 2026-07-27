# O₂ State Variable: Switching from O_total to C_aq
# ==================================================
#
# This document describes the change from total soil O₂ (O_total) to
# aqueous O₂ concentration (C_aq) as the state variable, and the
# rationale behind it.


## Problem

The model state variable O was originally total O₂ per unit bulk soil volume:

    O_total = C_aq · θ + C_gas · θ_a

where C_gas = K_H · C_aq (Henry's law, K_H ≈ 30 for O₂).

This formulation is correct only when water content θ is uniform in space
and constant in time. When EPS and F_i modify water retention, θ develops
spatial gradients that create two artifacts:

1. **Spurious advection**: diffusing O_total with spatially varying θ
   introduces a cross-term C_aq · ∇(θ + K_H θ_a) that drives O₂ toward
   high-θ regions (near the POM), producing unphysical O₂ concentrations
   exceeding the atmospheric boundary value.

2. **Storage inconsistency**: when θ changes between timesteps, O_total
   stays fixed but C_aq shifts implicitly, misrepresenting the true
   aqueous concentration available for microbial uptake.


## Solution: C_aq as state variable

Switch to aqueous O₂ concentration as the state variable. The governing
equation is:

    (θ + K_H θ_a) ∂C_aq/∂t = (1/r²) ∂/∂r [r² D̂_O ∂C_aq/∂r] − α_O · Resp

where D̂_O = D_w τ_aq θ + D_a τ_gas K_H θ_a is the total flux-based
diffusion coefficient (preserving dual-phase gas+liquid transport).

Dividing both sides by the capacity (θ + K_H θ_a):

    ∂C_aq/∂t = (1/r²) ∂/∂r [r² D_O^eff ∂C_aq/∂r] − α_O · Resp / (θ + K_H θ_a)

where:

    D_O^eff = D̂_O / (θ + K_H θ_a)
            = D_w τ_aq θ/(θ + K_H θ_a) + D_a τ_gas K_H θ_a/(θ + K_H θ_a)

This is EXACTLY what D_eff_oxygen() already computes. The existing diffusion
coefficient was implicitly designed for C_aq, not O_total.


## Why gas-phase diffusion is preserved

The effective diffusion coefficient D_O^eff contains both pathways:

    D_O^eff = D_w · f_aq · τ_aq  +  D_a · f_gas · τ_gas

where f_aq = θ/(θ + K_H θ_a) and f_gas = K_H θ_a/(θ + K_H θ_a).

In typical unsaturated soil: D_a ≈ 10⁴ × D_w, so the gas term dominates
despite the tortuosity penalty. Switching to C_aq does not lose the gas
pathway — it is embedded in D_O^eff through the f_gas weighting.

Physically: a gradient in C_aq implies a proportional gradient in C_gas
(via Henry's law). Gas diffuses along its own gradient in air-filled pores.
The effective D captures both transport mechanisms acting on the single
driving force ∇C_aq.


## Changes required

### 1. reactions.jl — compute_source_terms()

O_aq computation simplifies (O is now C_aq directly):

    BEFORE:  O_aq = O * θ / (θ + K_H * θ_a)
    AFTER:   O_aq = O

S_O must be divided by capacity:

    BEFORE:  S_O = -α_O * Resp_total
    AFTER:   S_O = -α_O * Resp_total / (θ + K_H * θ_a)

### 2. initial_conditions.jl — create_initial_state()

Initialize O to aqueous concentration:

    BEFORE:  O_total = O2_gas * (θ + K_H * θ_a) / K_H
    AFTER:   O_aq = O2_gas / K_H
    state.O .= O_aq

Note: O₂ is NOT diluted by ω (unchanged from before).

### 3. timestepper.jl — run_simulation()

O₂ boundary condition simplifies:

    BEFORE:  O2_amb = O2_gas * (workspace.θ[end] + K_H * workspace.θ_a[end]) / K_H
    AFTER:   O2_amb = O2_gas / K_H

Remove the IRA correction entirely — not needed when C_aq is the state
variable (C_aq doesn't change when θ changes; only the total does).

POM dissolution O_aq computation simplifies:

    BEFORE:  O_aq_0 = state.O[1] * θ_0 / (θ_0 + K_H * θ_a_0)
    AFTER:   O_aq_0 = state.O[1]

### 4. diffusion_step.jl — NO CHANGE

D_eff_oxygen already returns the correct coefficient for C_aq.

### 5. effective_diffusion.jl — NO CHANGE

D_eff_oxygen formula is already D̂_O / (θ + K_H θ_a).

### 6. Postprocessing

Any code that reads state.O and expects total O₂ needs updating.
In practice, most diagnostics want C_aq (for Monod terms, anoxic zone
delineation), so this simplifies postprocessing.

The spatial snapshot column "O" now represents C_aq [μg/mm³].
To recover O_total for mass balance purposes:

    O_total = state.O * (θ + K_H * θ_a)


## Consequences for anoxic zone formation

With C_aq as state variable, the anoxic zone emerges correctly from two
mechanisms:

1. **Reduced diffusivity**: where EPS increases θ, θ_a decreases, shutting
   down the fast gas pathway in D_O^eff. O₂ resupply slows.

2. **Increased capacity sink**: the reaction sink α_O·Resp/(θ + K_H θ_a)
   is LARGER where θ is high (capacity is lower), meaning each unit of
   respiration depletes C_aq faster in wet regions.

Both effects concentrate O₂ depletion near the POM, creating the expanding
anoxic zone that drives the bacterial-to-fungal transition.


## What this does NOT change

- Carbon balance (O₂ is not part of carbon accounting)
- D_eff_oxygen formula
- Diffusion solver (crank_nicolson.jl)
- Diffusion step structure (diffusion_step.jl)
- Boundary condition types (still Neumann inner, Dirichlet outer)
