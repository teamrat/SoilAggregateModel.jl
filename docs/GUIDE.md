# Soil Aggregate Model — Guide

**Last Updated**: 2026-02-06  
**Audience**: Users and developers of the Julia v2 implementation  
**Companion documents**: REFERENCE.md (quick lookup), ARCHITECTURE_CLAUDE_CODE.md (solver design)

---

## Table of Contents

### Part I: Model Theory
1. [Physical Setting](#1-physical-setting)
2. [State Variables](#2-state-variables)
3. [Carbon Partitioning](#3-carbon-partitioning)
4. [Microbial Dynamics](#4-microbial-dynamics)
5. [Fungal Architecture](#5-fungal-architecture)
6. [EPS](#6-eps)
7. [MAOC Formation](#7-maoc-formation)
8. [POM Dissolution](#8-pom-dissolution)
9. [Oxygen](#9-oxygen)
10. [Temperature Framework](#10-temperature-framework)
11. [Water Retention and Diffusion](#11-water-retention-and-diffusion)
12. [Boundary Conditions](#12-boundary-conditions)
13. [Conservation Law](#13-conservation-law)

### Part II: Numerical Method
14. [Strang Splitting](#14-strang-splitting)
15. [Crank–Nicolson Diffusion](#15-cranknnicolson-diffusion)
16. [Reaction Step](#16-reaction-step)
17. [Adaptive Time-Stepping](#17-adaptive-time-stepping)
18. [Non-Negativity Protection](#18-non-negativity-protection)

### Part III: Code Organization
19. [File Structure](#19-file-structure)
20. [Parameter Structs](#20-parameter-structs)
21. [Workspace and Memory](#21-workspace-and-memory)
22. [Testing Strategy](#22-testing-strategy)

### Part IV: Developer Notes
23. [Design Principles](#23-design-principles)
24. [Known Corrections from Manuscript](#24-known-corrections-from-manuscript)
25. [Future Work](#25-future-work)

---

# Part I: Model Theory

## 1. Physical Setting

The model describes a single soil aggregate: a spherical domain with a POM (particulate organic matter) core at the center. The POM core has fixed radius r₀ (it does not shrink as POM decomposes; voids form within the core). The soil matrix extends from r₀ to r_max.

Computational grid: uniform radial spacing r_i = r₀ + i·h, i = 0, ..., n−1, h = (r_max − r₀)/(n−1). All PDEs are discretized in spherical coordinates with this grid.

The aggregate is an open system: oxygen diffuses in from the boundary, CO₂ from respiration diffuses out. Carbon that is respired leaves the system permanently — this is not a mass balance error, it is the intended behavior.

---

## 2. State Variables

Nine state variables evolve in time. Five are transported by diffusion and solved via Crank–Nicolson tridiagonal systems. Three are immobile and advanced pointwise. One (POM) is a scalar ODE.

**Diffusing (5)**:
- **C** — total mobile carbon. This is the *bulk* concentration (aqueous + fast-sorbed). See §3 for how microbes and MAOC see different derived concentrations.
- **B** — bacterial biomass. Small diffusion coefficient (chemotaxis, proportional to D_C).
- **F_n** — non-insulated fungi. Hyphal tip extension, diffuses through pore water.
- **F_m** — mobile fungi. Translocation network. Diffuses with **no tortuosity** (transport is internal to hyphae, not through pore space).
- **O** — oxygen. Dual-phase: aqueous + gas, with Henry's law partitioning.

**Immobile (3)**:
- **F_i** — insulated fungi. Coated by new growth, bound to mineral surfaces. Source of aggregate stability.
- **E** — EPS. Produced by bacteria during growth.
- **M** — MAOC. Mineral-associated organic carbon, formed by slow stabilization of the fast-sorbed pool.

**Scalar (1)**:
- **P** — POM. Well-mixed within the fixed core. Dissolves enzymatically at the surface.

---

## 3. Carbon Partitioning

This is subtle and critical. The state variable C represents total mobile carbon per mm³ bulk soil. But different processes see different fractions of it.

The fast equilibrium partitioning gives:
```
C_aq = C / (θ + ρ_b · k_d)           aqueous concentration [μg/mm³ water]
C_eq = k_d · C_aq                     sorbed concentration  [μg/g solid]
```

This partitioning is instantaneous (local equilibrium assumption). The key decisions:

- **Microbial uptake** (R_B, R_F) uses **C_aq** — microbes consume what is dissolved.
- **EPS degradation** (R_rec_E) uses **C_aq** — scavenging responds to dissolved substrate scarcity.
- **MAOC equilibrium** (M_eq) uses **C_eq** — the Langmuir-Freundlich isotherm describes mineral surface reactions.
- **Diffusion PDE** uses **C** — the total mobile pool diffuses, retarded by the factor θ/(θ+ρ_b·k_d).
- **MAOC coupling in S_C** is simply **−J_M** (no (θ+ρ_b·k_d)/k_d factor). This is because J_M is defined as dM/dt per mm³ bulk soil, and carbon lost from C equals carbon gained by M, one-to-one. See REFERENCE.md §19, correction #2.

---

## 4. Microbial Dynamics

### Bacteria

Bacterial uptake follows dual-Monod kinetics on C_aq and O_aq, with water stress:

$$R_B = r_{B,\max}(T) \cdot \frac{C_{aq}}{K_B + C_{aq}} \cdot \frac{O_{aq}}{L_B + O_{aq}} \cdot B \cdot \exp(\nu_B \, \psi) \cdot h_B$$

The model distinguishes growth and starvation. A "basal" uptake R_Bb is computed at the maintenance concentration C_B. When R_B > R_Bb (growth), the uptake is partitioned into biomass growth Γ_B = Y_B·R_B·(1−γ) and EPS production Γ_E = Y_B·R_B·γ. When R_B ≤ R_Bb (starvation), all uptake goes to maintenance: Γ_B = R_B, Γ_E = 0.

Water stress uses **positive ν_B** with **negative ψ**: exp(ν_B · ψ) < 1 in unsaturated soil. This is a deliberate sign convention.

### Fungi

Fungal uptake uses the same dual-Monod form but only (F_i + λ·F_n) contributes to uptake (not F_m):

$$R_F = r_{F,\max}(T) \cdot \frac{C_{aq}}{K_F + C_{aq}} \cdot \frac{O_{aq}}{L_F + O_{aq}} \cdot (F_i + \lambda \, F_n) \cdot \exp(\nu_F \, \psi)$$

All assimilated carbon Γ_F = Y_F·R_F enters F_m (the translocation network), then redistributes to F_i and F_n via transitions. Only F_i dies; F_n converts to F_i via insulation. This three-compartment architecture is detailed in §5.

---

## 5. Fungal Architecture

Three compartments capture the biology of fungal hyphae:

- **F_m** (mobile) — the translocation network. Receives all growth. Diffuses freely (no tortuosity — transport is internal to hyphae). Distributes resources to immobile compartments.
- **F_n** (non-insulated) — newly settled hyphae. λ fraction participates in uptake. Converts to F_i via insulation (rate ζ). Can be mobilized back to F_m.
- **F_i** (insulated) — old hyphae coated by new growth and EPS. Primary contributor to uptake. Dies at rate μ_F. Binds mineral particles (aggregate stability).

### Transitions

The mobile-to-immobile ratio Π = F_m/(F_i + F_n + ε_F) controls direction:

```
net_i = η · (β_i·Π − α_i·Π^δ) · F_i
net_n = η · (β_n·Π − α_n·Π^δ) · F_n
```

When Π is small (low mobile fraction), net terms are positive (immobilization — F_m → F_i or F_n). When Π is large (high mobile fraction), the α·Π^δ term dominates and net terms are negative (mobilization — F_i or F_n → F_m). The exponent δ > 1 creates a nonlinear threshold.

The conversion efficiency η < 1 means (1−η) of transferred biomass is lost to respiration. This respiration (Resp_F^conv) uses **abs()** on the net transfer sum to ensure it is always positive, regardless of transfer direction.

One-way insulation: trans_ni = ζ·F_n (F_n → F_i), not reversible.

---

## 6. EPS

EPS is produced by bacteria during growth: Γ_E = Y_B · R_B · γ.

EPS degrades when dissolved substrate is scarce:
$$R_{rec,E} = \mu_E^{\max}(T) \cdot \frac{K_E}{K_E + C_{aq}} \cdot E \cdot h_E$$

The inhibition term K_E/(K_E + C_aq) means degradation is suppressed when C_aq is high (plenty of substrate available) and maximal when C_aq → 0. This uses **C_aq** (not the total C), because the sorbed fraction is unavailable for microbial use.

EPS modifies water retention by decreasing the van Genuchten α parameter (more EPS → soil holds more water at a given potential).

---

## 7. MAOC Formation

The two-stage MAOC model separates fast reversible partitioning from slow structural stabilization.

**Stage 1** (instantaneous): Carbon partitions between aqueous (C_aq) and sorbed (C_eq = k_d·C_aq) phases. This is fast equilibrium, not explicitly solved — just algebraic.

**Stage 2** (kinetic): The sorbed pool C_eq drives slow formation of structurally stabilized MAOC (M) via a Langmuir-Freundlich isotherm:

$$M_{eq} = M_{\max} \cdot \frac{(k_L \cdot C_{eq})^{n_{LF}}}{1 + (k_L \cdot C_{eq})^{n_{LF}}}$$

The rate law uses softplus regularization (smooth approximation to max(0, x)):

$$J_M = \kappa_s(T) \cdot \varphi_\varepsilon(M_{eq} - M) \;-\; \kappa_d(T) \cdot \varphi_\varepsilon(M - M_{eq})$$

**Hysteresis**: κ_s > κ_d means sorption is faster than desorption. A wetting-drying cycle concentrates DOC during drying (θ↓ → C_aq↑ → C_eq↑ → enhanced MAOC formation), then dilutes during rewetting — but M responds slowly, retaining gains. This is a mechanistic explanation for the Birch effect.

**Moisture coupling** is fully mechanistic: θ appears in C_aq = C/(θ + ρ_b·k_d), so drying directly concentrates the aqueous phase and elevates C_eq without any empirical scalar.

**Temperature coupling**: κ_s and κ_d have separate activation energies. Since ℰ_{a,d} > ℰ_{a,s}, warming narrows the hysteresis gap — a testable prediction.

---

## 8. POM Dissolution

POM dissolution is enzymatic, driven by microbial activity at the POM surface. The flux density at r = r₀:

$$J_P = R_P^{\max}(T) \cdot \frac{P}{P_0} \cdot \frac{B_0}{K_{B,P}+B_0} \cdot \frac{F_{n,0}}{K_{F,P}+F_{n,0}} \cdot \frac{\theta_0}{\theta_P+\theta_0} \cdot \frac{O_{aq,0}}{L_P+O_{aq,0}}$$

J_P [μg-C/mm²/day] is the surface-specific rate. The total rate R_P = 4πr₀²·J_P [μg-C/day] is what appears in dP/dt = −R_P.

J_P enters the carbon PDE as a Neumann boundary condition at the inner boundary: −D_C ∂C/∂r|_{r₀} = J_P. It does **not** appear as a volumetric source term — dissolved carbon enters the domain through the POM surface as a flux.

---

## 9. Oxygen

Oxygen exists in both aqueous and gas phases. The effective diffusion combines both pathways:

$$D_O = D_{O_2,w}(T) \cdot \tau_{aq} \cdot \frac{\theta}{\theta + K_H \cdot \theta_a} + D_{O_2,a}(T) \cdot \tau_{gas} \cdot \frac{K_H \cdot \theta_a}{\theta + K_H \cdot \theta_a}$$

Henry's law constant K_H(T) controls the partitioning. K_H increases with warming (less O₂ dissolves in warm water), which is physically correct and important for the anoxic zone predictions.

Oxygen is consumed by all respiration: S_O = −α_O·(Resp_B + Resp_F + Resp_F^conv).

The outer boundary has Dirichlet: O(r_max) = O_amb(t). The inner boundary is no-flux.

---

## 10. Temperature Framework

All temperature dependence uses Arrhenius scaling with six distinct activation energies, not Q10. This is a deliberate upgrade from the MATLAB version.

$$k(T) = k_\text{ref} \cdot \exp\!\left[\frac{\mathcal{E}_a}{R}\left(\frac{1}{T_\text{ref}} - \frac{1}{T}\right)\right]$$

| ℰ_a | Default (J/mol) | Applied to |
|-----|---------|------------|
| ℰ_{a,B} | 60,000 | μ_max_B, m_B, r_B_max |
| ℰ_{a,F} | 55,000 | All fungal rates (growth, death, transitions, D_Fn, D_Fm) |
| ℰ_{a,E} | 50,000 | μ_E_max |
| ℰ_{a,s} | 25,000 | κ_s |
| ℰ_{a,d} | 40,000 | κ_d |
| ℰ_{a,P} | 60,000 | R_P_max |

All fungal biological rates share a single ℰ_{a,F} because these processes are mediated by common cellular machinery and separate energies would be unconstrained by data.

Diffusion coefficients use mechanistic temperature laws: Stokes-Einstein+VFT for DOC in water, Han-Bartels for O₂ in water, Chapman-Enskog for O₂ in air. Henry's law uses van't Hoff.

All temperature-dependent quantities are computed **once per timestep** and stored in a `TemperatureCache` struct. This avoids redundant exponential evaluations.

---

## 11. Water Retention and Diffusion

Water content θ(ψ) follows van Genuchten, modified by EPS and insulated fungi:

$$\alpha_\text{eff}(r) = \alpha_0 \cdot \exp(\omega_E \cdot E(r) + \omega_F \cdot F_i(r))$$

where ω_E, ω_F < 0 (EPS and F_i *decrease* α, meaning the soil retains more water at a given potential).

θ and derived diffusion coefficients are updated **once per timestep** (not per RHS evaluation), because E and F_i evolve on biological timescales (~days).

Tortuosity uses two different forms:
- Aqueous: τ_aq = θ²/θ_s^(2/3) — a deliberate simplification of Millington-Quirk. Not a bug.
- Gas: τ_gas = θ_a^(10/3)/θ_s² — standard Millington-Quirk.

---

## 12. Boundary Conditions

| Species | r = r₀ (inner) | r = r_max (outer) |
|---------|-----------------|---------------------|
| C | Neumann: −D_C ∂C/∂r = J_P (POM flux) | Neumann: ∂C/∂r = 0 |
| B | Neumann: ∂B/∂r = 0 | Neumann: ∂B/∂r = 0 |
| F_n | Neumann: ∂F_n/∂r = 0 | Neumann: ∂F_n/∂r = 0 |
| F_m | Neumann: ∂F_m/∂r = 0 | Neumann: ∂F_m/∂r = 0 |
| O | Neumann: ∂O/∂r = 0 | Dirichlet: O = O_amb(t) |

Ghost nodes handle Neumann; row elimination handles Dirichlet.

Note: C outer boundary is Neumann (no-flux), not Dirichlet. Carbon does not enter or leave through the outer boundary — the aggregate is closed to carbon at r_max.

---

## 13. Conservation Law

Summing all source/sink terms (excluding diffusion, which redistributes but conserves):

$$S_P + S_C + S_B + S_{F_i} + S_{F_n} + S_{F_m} + S_E + S_M = -(Resp_B + Resp_F + Resp_F^{conv})$$

The only carbon that leaves the system is CO₂ from respiration (which consumes O₂ and escapes through the boundary). This identity is verified in test_biology.jl at rtol = 1e-10.

For the full domain integral:
```
P(0) = P(t) + ∫(C+B+F_i+F_n+F_m+E+M)·4πr²dr + CO₂_cumulative
```

---

# Part II: Numerical Method

## 14. Strang Splitting

The operator splitting is:

```
Each Δt:
  1. Diffusion half-step (Δt/2): Crank–Nicolson for C, B, F_n, F_m, O
  2. Reaction full-step  (Δt)  : Pointwise source/sink at each node + POM scalar
  3. Diffusion half-step (Δt/2): Same as step 1
```

This gives second-order accuracy O(Δt²). The splitting decouples transport (tridiagonal, O(n)) from reactions (pointwise, O(n)), avoiding a monolithic 8n+1 implicit system.

---

## 15. Crank–Nicolson Diffusion

Each diffusing species satisfies:
$$u^{k+1} - u^k = \frac{\Delta t}{4}\left(L[u^k] + L[u^{k+1}]\right)$$

where L is the spherical Laplacian with spatially-varying diffusion. This rearranges to a tridiagonal system solved by the Thomas algorithm in O(n). Five tridiagonal solves per half-step.

The spherical Laplacian at node i:
$$L_i = \frac{1}{r_i^2 h^2}\left[r_{i+1/2}^2 D_{i+1/2}(u_{i+1} - u_i) - r_{i-1/2}^2 D_{i-1/2}(u_i - u_{i-1})\right]$$

Face-averaged diffusion: D_{i+1/2} = (D_i + D_{i+1})/2.

---

## 16. Reaction Step

At each grid node independently:
1. Compute C_aq, C_eq from C and θ
2. Compute all uptake rates (R_B, R_F), maintenance (R_Bb)
3. Determine growth/starvation regime
4. Compute assimilation (Γ_B, Γ_E, Γ_F), respiration (Resp_B, Resp_F)
5. Compute fungal transitions (net_i, net_n, trans_ni) and Resp_F_conv
6. Compute recycling (R_rec_B, R_rec_F, R_rec_E)
7. Compute MAOC dynamics (M_eq, J_M)
8. Assemble source terms (S_C, S_B, S_Fi, S_Fn, S_Fm, S_E, S_M, S_O)
9. Advance with Forward Euler: u^{k+1} = u^k + Δt · S

Then advance POM scalar: P^{k+1} = P^k − Δt · R_P.

Accumulate CO₂: CO₂ += Δt · Σ_i (Resp_B + Resp_F + Resp_F_conv)_i · V_i.

---

## 17. Adaptive Time-Stepping

```
max_change = max_i |S_i · Δt / u_i|    over all nodes and species

if max_change > 0.10 → halve Δt
if max_change < 0.01 → double Δt
bracket: Δt_min = 1e-4 days, Δt_max = 0.1 days
```

---

## 18. Non-Negativity Protection

Three sigmoid functions prevent state variables from going negative:

$$h(x) = \frac{\exp(\beta \, x)}{\exp(\beta \, x) + \exp(\beta \, x_{\min})}$$

with β = 50/x_min. Applied to:
- h_B: multiplies R_Bb (maintenance) and R_rec_B (death)
- h_Fi: multiplies R_rec_F (fungal death)
- h_E: multiplies R_rec_E (EPS degradation)

After each reaction step, non-negativity is enforced by clipping. Cumulative clipping magnitude is tracked as a diagnostic.

The Π ratio uses ε_F protection: Π = F_m/(F_i + F_n + ε_F) to avoid division by zero.

---

# Part III: Code Organization

## 19. File Structure

See README.md for full tree. Key design:

- **temperature/** — pure functions, no state. Compute scaling factors and physical constants.
- **physics/** — water retention and effective diffusion. Depend on θ, soil properties, temperature cache.
- **biology/** — all reaction terms. One file per organism group. Pure functions operating on scalar inputs (per-node).
- **solver/** — tridiagonal, Crank–Nicolson, reaction step, Strang loop. This is the numerical engine.
- **api.jl** — user entry point. Constructs solver, runs simulation, returns output.

---

## 20. Parameter Structs

Two main parameter structs (immutable):

**BiologicalProperties**: all rate constants, half-saturations, yields, activation energies at T_ref. Contains bacterial, fungal, EPS, MAOC, POM, and oxygen parameters. See REFERENCE.md §3 for complete listing.

**SoilProperties**: van Genuchten hydraulic parameters, bulk density, clay/silt fractions, MAOC sorption parameters, reference diffusion coefficients, aggregate stability parameters. See REFERENCE.md §4.

**EnvironmentalDrivers{FT,Fψ,FO}**: parametric type for T(t), ψ(t), O₂(t). Constants compile to literals; functions specialize at call site. No performance penalty for time-varying drivers.

---

## 21. Workspace and Memory

Pre-allocated workspace avoids heap allocations in the hot loop:

```julia
struct Workspace
    # Tridiagonal system (reused for each of 5 species)
    lower, diag, upper, rhs :: Vector{Float64}

    # Spatially varying (updated once per timestep)
    θ, θ_a, D_C, D_B, D_Fn, D_O :: Vector{Float64}

    # Temperature cache
    f_T :: TemperatureCache
end
```

The state is struct-of-arrays (`AggregateState` with 8 vectors + 2 scalars), optimal for tridiagonal solves where each species is a contiguous vector.

Target: zero allocations in the time-stepping hot loop. Memory < 1 MB per aggregate at n = 1500.

---

## 22. Testing Strategy

Bottom-up, one layer at a time:

| Step | Module | Tests |
|------|--------|-------|
| 1 | types, parameters, environment | Struct construction, defaults, type stability |
| 2 | temperature | Arrhenius returns 1.0 at T_ref; Han-Bartels units; Henry sign; VFT viscosity |
| 3 | tridiagonal | Thomas solver vs known solutions |
| 4 | physics | Water retention curve, effective diffusion at known θ |
| 5 | biology | All rate functions with hand-computed inputs; conservation identity |
| 6 | POM dissolution | Boundary flux with known biomass, θ, O₂ |
| 7 | Crank–Nicolson | Analytical steady-state diffusion in sphere |
| 8 | Reaction step | Zero-diffusion batch reactor |
| 9 | Strang loop | 30-day integration test |
| 10 | API | User-facing interface, output, carbon balance check |

Each step has its own test file. Tests run with `julia runtests.jl`.

Conservation test (test_biology.jl): verifies S_C + S_B + S_Fi + S_Fn + S_Fm + S_E + S_M = −(Resp_B + Resp_F + Resp_F_conv) at rtol = 1e-10, with realistic parameters (ρ_b = 1200).

---

# Part IV: Developer Notes

## 23. Design Principles

1. **Manuscript is authoritative for physics**. If code and manuscript disagree on an equation, the manuscript wins (with exceptions documented in REFERENCE.md §19).

2. **No unit conversions in kernels**. Everything is μg/mm³, mm, days, kPa, K, J/mol. Conversions happen only at the user API boundary.

3. **No default parameters in kernels**. Every parameter is passed explicitly. Defaults exist only in constructor functions.

4. **Temperature-dependent quantities computed once per timestep**. Stored in TemperatureCache. Never recomputed during Strang sub-steps.

5. **θ updated once per timestep**. E and F_i evolve on biological timescales; updating θ during sub-steps is wasted work.

6. **Zero heap allocations in hot loop**. All workspace is pre-allocated. Functions are mutation-based (!) where needed.

7. **Struct-of-arrays** for state. Each species is a contiguous vector. Optimal for tridiagonal solves.

---

## 24. Known Corrections from Manuscript

These issues were found during the 2026-02-05/06 code audit. The code has the corrected versions. The manuscript may still need updating.

### 24.1 MAOC Coupling in S_C

**Manuscript**: S_C = ... − J_M · (θ + ρ_b·k_d)/k_d  
**Correct**: S_C = ... − J_M

The factor was derived under the assumption that J_M represents a rate per mm³ of solid phase, but J_M is actually defined as dM/dt per mm³ bulk soil. When J_M μg/mm³/day of carbon transfers from the mobile pool C to the mineral pool M, the total mobile pool C decreases by exactly J_M — no scaling needed. With typical parameters (θ=0.35, ρ_b=1200, k_d=0.1), the erroneous factor ≈ 1204, which would over-drain C by three orders of magnitude.

### 24.2 Henry's Law Sign

**Manuscript**: K_H(T) = K_H_ref · exp[**−**ΔH_sol/R · (1/T − 1/T_ref)]  
**Correct**: K_H(T) = K_H_ref · exp[ΔH_sol/R · (1/T − 1/T_ref)]   (no leading negative)

With ΔH_sol = −12,000 J/mol (exothermic dissolution), the correct formula gives K_H increasing with warming (less O₂ dissolves in warm water). The erroneous sign would make warm soils *more* oxygenated at the aggregate core — opposite to reality.

### 24.3 EPS Degradation Substrate

**Manuscript (old)**: K_E/(K_E + C)  
**Correct**: K_E/(K_E + C_aq)

Scavenging responds to dissolved concentration, not total mobile carbon. Using C would underestimate scavenging pressure because the sorbed fraction is unavailable.

### 24.4 Resp_F^conv Absolute Value

**Manuscript (old)**: (1−η)·[net_i/η + net_n/η]  
**Correct**: (1−η)·|net_i/η + net_n/η|

During mobilization, the bracketed expression goes negative. Without abs(), respiration becomes negative — nonphysical.

### 24.5 POM Dissolution Notation

The manuscript now distinguishes J_P [μg/mm²/day] (flux density) from R_P = 4πr₀²·J_P [μg/day] (total rate). The Neumann BC uses J_P; the POM ODE uses R_P.

---

## 25. Future Work

### Solver Assembly (current priority)
Steps 7–10 of the implementation order: Crank–Nicolson, reaction step, Strang loop, API.

### Cohort Management
Multiple aggregates formed at different times from different POM input events. Requires efficient data structures for 400+ concurrent aggregates and scheduling of POM additions.

### Profile Coupling
Two-way coupling with a soil profile model (Richards equation + heat transport + gas transport). One-way (profile → aggregate: T, ψ) is straightforward. Two-way (aggregate → profile: O₂ consumption, CO₂ production) requires careful time-stepping coordination.

### Performance Validation
Target: single aggregate < 30 seconds for 10-year run at n = 1500. Verify zero allocations in hot loop with `@allocated`. Benchmark tridiagonal solve cost.

---

**End of Guide**
