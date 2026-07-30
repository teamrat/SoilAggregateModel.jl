# Soil Aggregate Model v2 — Implementation Architecture

**Status**: Final — authoritative guide for implementation  
**Companion**: `manuscript.tex` is authoritative for all physics equations  
**Rule**: If this document and the manuscript disagree on physics, the manuscript wins.  
**Units**: μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol throughout.

---

## 0. Performance Targets

- 400 aggregates × n=1500 × 10 years
- Single aggregate: < 30 seconds for 10-year run at n=1500
- Memory: < 1 MB per aggregate (enables 1000+ concurrent)
- Zero heap allocations in the time-stepping hot loop

---

## 1. Problem Structure

9 state variables. 5 diffuse, 3 are immobile local ODEs, 1 is a scalar ODE:

| Variable | Symbol | Transport | Solver role |
|----------|--------|-----------|-------------|
| Dissolved organic carbon | C | Diffusion (retarded by sorption) | Tridiagonal PDE |
| Bacteria | B | Diffusion (chemotaxis, D_B ∝ D_C) | Tridiagonal PDE |
| Non-insulated fungi | F_n | Diffusion (hyphal tip extension) | Tridiagonal PDE |
| Mobile fungi | F_m | Diffusion (internal translocation, constant D) | Tridiagonal PDE |
| Oxygen | O | Diffusion (aqueous + gas phase) | Tridiagonal PDE |
| Insulated fungi | F_i | Immobile | Local ODE at each node |
| EPS | E | Immobile | Local ODE at each node |
| MAOC | M | Immobile | Local ODE at each node |
| POM | P | Spatially uniform | Single scalar ODE |

**Key insight**: 5 independent tridiagonal solves (O(n) each) + pointwise reactions replaces one monolithic 8n+1 implicit system.

Note: B and F_n have small diffusion coefficients relative to C and O. If profiling shows the extra tridiagonal solves are not worth the cost, these can be demoted to local ODEs. Start with all 5.

---

## 2. Time Integration: Strang Splitting

```
For each timestep Δt:
    1. Diffusion half-step (Δt/2): Solve C, B, F_n, F_m, O via Crank-Nicolson (5 tridiagonal)
    2. Reaction full-step (Δt):    Advance all 8 spatial variables + P via source/sink terms
    3. Diffusion half-step (Δt/2): Solve C, B, F_n, F_m, O via Crank-Nicolson (5 tridiagonal)
```

Strang splitting gives second-order accuracy O(Δt²).

### Adaptive Timestep

- If max(|rate × Δt / state|) > 0.10 at any node → halve Δt
- If max(|rate × Δt / state|) < 0.01 everywhere → double Δt
- Bracket: Δt_min = 1e-4 days, Δt_max = 0.1 days

### Non-negativity

After each reaction step, enforce non-negativity on all state variables. Track cumulative clipping magnitude as a diagnostic. If clipping exceeds 0.1% of total carbon in any step, reduce Δt_max. The h_B, h_{F_i}, h_E sigmoids (Section 4.4) prevent most negativity, but numerical overshoot is still possible with explicit integration.

### Fallback

If Strang splitting proves insufficiently accurate for stiff MAOC/microbial coupling, upgrade to DifferentialEquations.jl `KenCarp4` (IMEX-ARK). The reaction/diffusion separation is the same; only the time integrator changes.

---

## 3. Spatial Discretization

### 3.1 Grid

Uniform radial grid: r_i = r_POM + i·h, i = 0, ..., n-1, h = (r_max - r_POM) / (n-1).

(Future option: logarithmic grid near POM surface. Tridiagonal structure preserved; only coefficients change.)

### 3.2 Spherical Laplacian

For species u with spatially varying diffusion D(r):

$$L_i = \frac{1}{r_i^2 h^2}\left[r_{i+1/2}^2 D_{i+1/2}(u_{i+1} - u_i) - r_{i-1/2}^2 D_{i-1/2}(u_i - u_{i-1})\right]$$

Face-averaged coefficients: D_{i+1/2} = (D_i + D_{i+1}) / 2.

### 3.3 Crank-Nicolson Diffusion Half-Step

$$u^{n+1} - u^n = \frac{\Delta t}{4}\left(L[u^n] + L[u^{n+1}]\right)$$

Rearranges to a tridiagonal system Au^{n+1} = b, solved by Thomas algorithm in O(n). No iteration needed. Five tridiagonal solves per half-step = O(5n) total.

### 3.4 Boundary Conditions

| Species | r = r_POM (inner) | r = r_max (outer) |
|---------|-------------------|-------------------|
| C | Neumann: −D_C ∂C/∂r = R_P/(4πr₀²) (POM flux per unit area) | Neumann: ∂C/∂r = 0 |
| B | Neumann: ∂B/∂r = 0 | Neumann: ∂B/∂r = 0 |
| F_n | Neumann: ∂F_n/∂r = 0 | Neumann: ∂F_n/∂r = 0 |
| F_m | Neumann: ∂F_m/∂r = 0 | Neumann: ∂F_m/∂r = 0 |
| O | Neumann: ∂O/∂r = 0 | Dirichlet: O = O_amb(t) |

The POM dissolution flux enters as a proper Neumann BC at the first grid node, NOT as a volumetric source.

Implementation: ghost node for Neumann; row elimination for Dirichlet.

---

## 4. Reaction Step

### 4.1 Overview

At each grid point i independently, compute source/sink terms and advance. ALL 8 spatial variables get source/sink updates during the reaction step. The 3 immobile species (F_i, E, M) are advanced entirely within the reaction step. The 5 diffusing species get reactions here and transport in the diffusion step.

### 4.2 POM Dissolution (Scalar, depends on state at r = r₀)

POM dissolution is enzymatic, driven by microbial activity at the POM surface. The dissolution rate depends on local biomass, water content, and oxygen:

$$R_P = 4\pi r_0^2 \; R_P^{\max} \; \frac{P}{P_0} \; \tfrac{1}{2}\left(\frac{B_0}{K_{B,P} + B_0} + \frac{F_{n,0}}{K_{F,P} + F_{n,0}}\right) \; \frac{\theta_0}{\theta_P + \theta_0} \; \frac{O_{aq,0}}{L_P + O_{aq,0}}

**Additive, not multiplicative.** The two microbial terms are averaged, not
multiplied. Depolymerisation is extracellular and acts on the substrate, so
enzymes from either community suffice: at saturating biomass of one community
alone the rate is half maximal, and both at saturation recover the full rate. A
product form would assert that POM cannot be depolymerised without both
populations present, which is false — white-rot fungi mineralise lignin without
bacteria, and cellulolytic bacteria act without fungi. The equal weighting is a
stated default and is calibratable; the additive *form* is not. Matches
`manuscript-4-5.tex` Eq.~(\ref{eq:R_P}) and
`src/carbon/pom_dissolution.jl:36-38`.


$$

where subscript 0 denotes values at the first grid node (r = r_POM), and K_{B,P}, K_{F,P} are half-saturation constants for bacterial and fungal contribution to POM dissolution.

R_P scales with temperature via ε_{a,P} (POM activation energy, see Section 5).

**Note**: R_P is a scalar (total mass rate, μg-C/day). The Neumann BC flux is R_P / (4πr₀²) per unit area.

### 4.3 Source/Sink Terms (per node)

**`docs/REFERENCE.md` §14 is the single source for these.** The block below is a
summary; if it disagrees with §14, §14 wins, because §14 is transcribed from
`src/solver/reactions.jl` and re-verified against it.

```
S_P  = −R_P
S_C  = −R_B − R_F + R_rec − J_M
S_B  = Γ_B − R_rec,B
S_Fi = immobil_i − R_rec,F
S_Fn = immobil_n
S_Fm = Γ_F − immobil_i − immobil_n − Resp_F_conv
S_E  = Γ_E − R_rec,E
S_M  = J_M
S_O  = −α_O·(Resp_B + Resp_F + Resp_F_conv) / (θ + K_H(T)·θ_a)
```

where `immobil_i`, `immobil_n` and `Resp_F_conv` come from
`fungal_transitions` (`src/biology/fungi.jl`) — see REFERENCE §14 for their
definitions. **There is no separate `±ζ·F_n` insulation term**: ζ is a
dimensionless splitting fraction applied inside `fungal_transitions`, not a
one-way drain rate.

*Corrected 2026-07-28. This block previously carried the inverted α/β
convention, the deleted `ζ·F_n` drain, the pre-audit
`−J_M·(θ+ρ_b·k_d)/k_d` coupling (REFERENCE §26 erratum 2, fixed in the code in
February), and an `S_O` missing its O₂ capacity divisor. Four stale claims in
eight lines, in the document `CLAUDE.md` names as authoritative.*

Cumulative CO₂ (diagnostic, not a state variable):
```
dCO2/dt = Resp_B + Resp_F + Resp_F_conv   (summed over all nodes × volume elements)
```

### 4.4 Non-Negativity Sigmoids

Three smooth threshold functions prevent state variables from going negative:

**h_B** (bacteria):
$$h_B = \frac{\exp(\beta_B \, B)}{\exp(\beta_B \, B) + \exp(\beta_B \, B_{\min})}$$
with β_B = 50/B_min. Applied to: R_{B,b} (maintenance) and R_{rec,B} (death).

**h_{F_i}** (insulated fungi):
$$h_{F_i} = \frac{\exp(\beta_F \, F_i)}{\exp(\beta_F \, F_i) + \exp(\beta_F \, F_{i,\min})}$$
with β_F = 50/F_{i,min}. Applied to: R_{rec,F} (fungal death).

**h_E** (EPS):
$$h_E = \frac{\exp(\beta_E \, E)}{\exp(\beta_E \, E) + \exp(\beta_E \, E_{\min})}$$
with β_E = 50/E_min. Applied to: R_{rec,E} (EPS degradation).

### 4.5 Π Protection

The mobile-to-immobile ratio Π = F_m / (F_i + F_n) needs protection against division by zero:

$$\Pi = \frac{F_m}{F_i + F_n + \varepsilon_F}$$

with ε_F a small constant (e.g., 1e-10 μg/mm³).

### 4.6 MAOC: Smooth Switching

Replace the non-differentiable max(0, ·) with softplus regularization:

$$J_M = \kappa_s(T) \cdot \phi_\varepsilon(M_{eq} - M) - \kappa_d(T) \cdot \phi_\varepsilon(M - M_{eq})$$

$$\phi_\varepsilon(x) = \varepsilon \cdot \ln(1 + e^{x/\varepsilon})$$

Use ε_maoc = 0.01 μg/mm³. Indistinguishable from max(0,x) for |x| > 0.1 but C∞ smooth.

### 4.7 Local Reaction Integrator

Start with Forward Euler. If instability arises (MAOC switching, low O₂), upgrade to implicit midpoint or SDIRK2 per node. The local system is 8×8, so direct Jacobian solve is trivial.

---

## 5. Temperature Framework

Three distinct mechanisms — do NOT conflate them.

### 5.1 Arrhenius for Biological Rates

$$k(T) = k_{ref} \cdot \exp\left[\frac{\mathcal{E}_a}{R}\left(\frac{1}{T_{ref}} - \frac{1}{T}\right)\right]$$

**Symbol**: Use 𝓔_a (mathcal E) for activation energy to avoid clash with E (EPS state variable).

Six distinct activation energies:

| Process | Symbol | Default (J/mol) | Applied to |
|---------|--------|-----------------|------------|
| Bacterial metabolism | 𝓔_{a,B} | 60,000 | μ_max_B, m_B, r_{B,max} |
| Fungal metabolism | 𝓔_{a,F} | 55,000 | μ_max_F, m_F, α_i, α_n, β_i, β_n, D_Fn0, D_Fm0 — **not ζ**, which is a splitting fraction and is not temperature-scaled (2026-07-30) |
| EPS degradation | 𝓔_{a,E} | 50,000 | μ_E_max (k_EPS) |
| MAOC sorption | 𝓔_{a,s} | 25,000 | κ_s |
| MAOC desorption | 𝓔_{a,d} | 40,000 | κ_d |
| POM dissolution | 𝓔_{a,P} | 60,000 | R_P_max |

**Key decision**: All fungal biological rates (growth, death, transitions, translocation, hyphal extension) share a single 𝓔_{a,F}. Rationale: these processes are mediated by common cellular machinery; separate activation energies would be unconstrained by available data.

**Key prediction**: 𝓔_{a,d} > 𝓔_{a,s} → warming narrows MAOC hysteresis.

### 5.2 Diffusion in Water

**DOC**: Stokes-Einstein via water viscosity ratio:

$$D_{C0}(T) = D_{C0,ref} \cdot \frac{T}{T_{ref}} \cdot \frac{\eta(T_{ref})}{\eta(T)}$$

Water viscosity (Vogel-Fulcher-Tammann):

$$\ln[\eta(T) / \text{mPa·s}] = -3.7188 + \frac{578.919}{T - 137.546}$$

Valid 273–373 K.

**O₂ in water**: Han & Bartels (1996) empirical:

$$\log_{10}[D_{O_2,w} / \text{cm}^2\text{s}^{-1}] = -4.410 + \frac{773.8}{T} - \left(\frac{506.4}{T}\right)^2$$

Convert to mm²/day: multiply cm²/s by 10² × 86400 = 8.64 × 10⁶.

### 5.3 Diffusion in Air

Chapman-Enskog:

$$D_a(T) = D_{a,ref} \cdot (T / T_{ref})^{1.75}$$

### 5.4 Henry's Law Constants

Van't Hoff:

$$K_H(T) = K_{H,ref} \cdot \exp\left[-\frac{\Delta H_{sol}}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right]$$

O₂: K_{H,ref} = 31.25 (dimensionless at 298K), ΔH_sol = −12,000 J/mol.

### 5.5 Implementation: TemperatureCache

All temperature-dependent quantities are computed **once per timestep** and stored:

```julia
mutable struct TemperatureCache
    # Arrhenius factors (dimensionless multipliers on reference rates)
    f_bac::Float64      # bacteria: μ_max_B, m_B, r_B_max
    f_fun::Float64      # fungi: μ_max_F, m_F, transitions, D_Fn0, D_Fm0
    f_eps::Float64      # EPS: μ_E_max
    f_maoc_s::Float64   # MAOC sorption: κ_s
    f_maoc_d::Float64   # MAOC desorption: κ_d
    f_pom::Float64      # POM: R_P_max
    # Pure-phase diffusion coefficients [mm²/day]
    D_O2_w::Float64
    D_DOC_w::Float64
    D_O2_a::Float64
    D_Fm::Float64       # = D_Fm0_ref × f_fun (spatially uniform)
    # Henry's law
    K_H_O::Float64
end
```

The functions that compute these (Arrhenius, Stokes-Einstein, etc.) are **module-level functions**, NOT stored in structs. This avoids abstract `Function` type performance penalty. Users who want custom temperature functions can override at the API level.

---

## 6. Environmental Drivers

Temperature, matric potential, and ambient O₂ are externally imposed functions of time. They may be constant, time-varying, or interpolated from data.

```julia
struct EnvironmentalDrivers{FT, Fψ, FO}
    T::FT       # T(t) → temperature [K]
    ψ::Fψ       # ψ(t) → matric potential [kPa]
    O2::FO      # O2(t) → boundary O₂ concentration [μg/mm³]
end
```

Parametric types ensure Julia specializes: constant functions are as fast as literal constants.

Convenience constructors:
```julia
# All constant
EnvironmentalDrivers(293.15, -33.0, 0.21)
# → internally wraps scalars as t -> value

# Mixed
EnvironmentalDrivers(t -> 293.15 + 5*sin(2π*t/365), -33.0, 0.21)

# From data
using Interpolations
EnvironmentalDrivers(
    linear_interpolation(t_data, T_data),
    linear_interpolation(t_data, ψ_data),
    0.21
)
```

---

## 7. Water Content

θ(r) depends on ψ(t), E(r), F_i(r) via modified van Genuchten (manuscript Section 2.4.8):

$$\alpha_{eff}(r) = \alpha_0 \cdot \exp(\omega_E \cdot E(r) + \omega_F \cdot F_i(r))$$

Updated **once per timestep**, NOT per RHS evaluation. E and F_i evolve on biological timescales (~days). θ is held fixed during Strang sub-steps.

---

## 8. Effective Diffusion Coefficients

Computed once per timestep from θ(r), θ_a(r), and temperature-dependent pure-phase values:

```
D_C[i]  = D_DOC_w(T) · τ(θ[i]) · θ[i]/(θ[i] + ρ_b·k_d)
D_B[i]  = D_B_rel · D_C[i]
D_Fn[i] = D_Fn0_ref · f_fun · τ(θ[i])
D_Fm    = D_Fm0_ref · f_fun                   (scalar, no tortuosity)
D_O[i]  = D_O2_w(T)·θ[i]/(θ[i]+K_H·θ_a[i])·θ[i]²/θ_s^(2/3)
          + D_O2_a(T)·K_H·θ_a[i]/(θ[i]+K_H·θ_a[i])·θ_a[i]^(10/3)/θ_s²
```

where τ(θ) = θ²/θ_s^(2/3) (Millington-Quirk). Tortuosity is temperature-independent.

Note: D_Fm has no tortuosity — internal translocation occurs within the hyphal network, not through pore space. This means mobile fungi can transport through dry regions. The manuscript should state this explicitly.

---

## 9. State Layout and Memory

### 9.1 State: Struct-of-Arrays

```julia
mutable struct AggregateState
    C::Vector{Float64}      # n — dissolved organic carbon
    B::Vector{Float64}      # n — bacteria
    F_n::Vector{Float64}    # n — non-insulated fungi
    F_m::Vector{Float64}    # n — mobile fungi
    O::Vector{Float64}      # n — oxygen
    F_i::Vector{Float64}    # n — insulated fungi
    E::Vector{Float64}      # n — EPS
    M::Vector{Float64}      # n — MAOC
    P::Float64              # scalar — POM
    CO2_cumulative::Float64 # diagnostic — total CO₂ respired
end
```

Struct-of-arrays: each species is a contiguous vector. Optimal for the tridiagonal solve (dominant cost). Reaction computation at node i requires 8 strided reads, but L1-cacheable at n=1500 (8 × 8 bytes × 1500 ≈ 96 KB).

### 9.2 Workspace (Pre-Allocated)

```julia
struct Workspace
    # Tridiagonal system (reused for each of 5 diffusing species)
    lower::Vector{Float64}   # n-1
    diag::Vector{Float64}    # n
    upper::Vector{Float64}   # n-1
    rhs::Vector{Float64}     # n

    # Spatially varying quantities (updated once per timestep)
    θ::Vector{Float64}       # n — water content
    θ_a::Vector{Float64}     # n — air-filled porosity
    D_C::Vector{Float64}     # n — effective C diffusion
    D_B::Vector{Float64}     # n — effective B diffusion
    D_Fn::Vector{Float64}    # n — effective F_n diffusion
    D_O::Vector{Float64}     # n — effective O diffusion
    # D_Fm is scalar — stored in TemperatureCache

    # Temperature cache
    f_T::TemperatureCache
end
```

**Zero allocations in the hot loop.**

### 9.3 Solver Struct

```julia
struct AggregateSolver{FT, Fψ, FO}
    # Grid
    n::Int
    r::Vector{Float64}      # node positions [mm]
    h::Float64              # grid spacing [mm]
    r_POM::Float64          # inner boundary [mm]
    r_max::Float64          # outer boundary [mm]

    # State (mutable)
    state::AggregateState

    # Workspace (pre-allocated)
    workspace::Workspace

    # Parameters (immutable)
    biology::BiologicalProperties
    soil::SoilProperties

    # Environment (callable)
    env::EnvironmentalDrivers{FT, Fψ, FO}
end
```

---

## 10. Parameter Structs

### 10.1 BiologicalProperties

```julia
struct BiologicalProperties
    # --- Bacterial ---
    r_B_max::Float64        # Max specific uptake rate at T_ref [1/day]
    K_B::Float64            # Half-saturation for DOC [μg/mm³]
    L_B::Float64            # Half-saturation for O₂ [μg/mm³]
    ν_B::Float64            # Water potential sensitivity [1/kPa]
    Y_B_max::Float64        # Maximum growth yield [-]
    K_Y::Float64            # Half-saturation for yield [-]
    γ::Float64              # EPS allocation fraction [-]
    C_B::Float64            # Basal carbon requirement [μg/mm³]
    μ_B::Float64            # Mortality rate at T_ref [1/day]
    B_min::Float64          # Minimum viable biomass [μg/mm³]
    Ea_B::Float64           # Activation energy [J/mol]

    # --- Fungal ---
    r_F_max::Float64        # Max specific uptake rate at T_ref [1/day]
    K_F::Float64            # Half-saturation for DOC [μg/mm³]
    L_F::Float64            # Half-saturation for O₂ [μg/mm³]
    ν_F::Float64            # Water potential sensitivity [1/kPa]
    Y_F::Float64            # Growth yield [-] (or Y_F_max + K_YF if uptake-dependent)
    μ_F::Float64            # Mortality rate at T_ref [1/day]
    F_i_min::Float64        # Minimum viable insulated biomass [μg/mm³]
    Ea_F::Float64           # Activation energy [J/mol] — shared by ALL fungal rates

    # --- Fungal transitions (Falconer naming; see REFERENCE §3.3, §5a) ---
    α_i::Float64            # IMMOBILIZATION rate → insulated [1/day] (gain; carries δ)
    α_n::Float64            # IMMOBILIZATION rate → non-insulated [1/day] (gain; carries δ)
    β_i::Float64            # MOBILIZATION rate from insulated [1/day] (loss; linear in Π)
    β_n::Float64            # MOBILIZATION rate from non-insulated [1/day] (loss; linear in Π)
    delta::Float64          # Immobilization exponent — Falconer's θ [-]
    η_conv::Float64         # Conversion efficiency [-]
    ζ::Float64              # Insulation SPLITTING FRACTION on the F_n tendency [-] — NOT a rate
    λ::Float64              # Fraction of F_n at uptake surfaces [-]
    D_Fn0::Float64          # Hyphal extension diffusivity at T_ref [mm²/day]
    D_Fm0::Float64          # Internal translocation rate at T_ref [mm²/day]
    ε_F::Float64            # Π denominator protection [μg/mm³]

    # --- EPS ---
    μ_E_max::Float64        # Max EPS degradation rate at T_ref [1/day]
    K_E::Float64            # Substrate inhibition concentration [μg/mm³]
    E_min::Float64          # Minimum EPS for h_E sigmoid [μg/mm³]
    Ea_EPS::Float64         # Activation energy [J/mol]

    # --- MAOC ---
    κ_s_ref::Float64        # Sorption rate at T_ref [1/day]
    κ_d_ref::Float64        # Desorption rate at T_ref [1/day]
    Ea_MAOC_sorb::Float64   # Activation energy, sorption [J/mol]
    Ea_MAOC_desorb::Float64 # Activation energy, desorption [J/mol]
    ε_maoc::Float64         # Softplus smoothing width [μg/mm³]

    # --- POM ---
    R_P_max::Float64        # Max dissolution rate at T_ref [μg-C/mm²/day]
    P_0::Float64            # Initial POM mass [μg-C]
    r_0::Float64            # POM radius [mm]
    θ_P::Float64            # Half-saturation water content for dissolution [-]
    L_P::Float64            # Half-saturation O₂ for dissolution [μg/mm³]
    K_B_P::Float64          # Half-saturation bacteria for dissolution [μg/mm³]
    K_F_P::Float64          # Half-saturation fungi for dissolution [μg/mm³]
    ρ_POM::Float64          # POM carbon density [μg-C/mm³]
    Ea_POM::Float64         # Activation energy [J/mol]

    # --- Oxygen ---
    α_O::Float64            # Respiratory quotient [μg-O₂/μg-C]

    # --- Reference ---
    T_ref::Float64          # Reference temperature [K]
end
```

### 10.2 SoilProperties

```julia
struct SoilProperties
    # Van Genuchten
    θ_r::Float64            # Residual water content [-]
    θ_s::Float64            # Saturated water content [-]
    α_vg::Float64           # van Genuchten α [1/kPa]
    n_vg::Float64           # van Genuchten n [-]

    # EPS/fungi modification of water retention
    ω_E::Float64            # EPS effect on α (negative) [mm³/μg]
    ω_F::Float64            # Fungi effect on α (negative) [mm³/μg]

    # Equilibrium sorption
    k_d_eq::Float64         # Linear partition coefficient [mm³/μg]
    ρ_b::Float64            # Bulk density [μg/mm³]

    # MAOC capacity (Langmuir-Freundlich)
    k_L::Float64            # Langmuir affinity [mm³/μg]
    n_LF::Float64           # Freundlich exponent [-]
    k_ma::Float64           # MAOC capacity per g clay+silt [g-C/g-mineral, dimensionless]
                            # M_max = maoc_capacity(soil) = k_ma·f_clay_silt·ρ_b
    f_clay_silt::Float64    # Clay+silt mass fraction [-]

    # Reference diffusion at T_ref [mm²/day]
    D_C0_ref::Float64       # DOC in water
    D_O2_w_ref::Float64     # O₂ in water
    D_O2_a_ref::Float64     # O₂ in air
    D_B_rel::Float64        # Bacterial motility relative to D_C [-]

    # Aggregate stability
    k_F::Float64            # Specific binding strength [kPa/(μg/mm³)]
    χ::Float64              # Particle adhesion length scale [mm]
    a_p::Float64            # Particle radius [mm]
end
```

---

## 11. Initialization

```julia
function initialize_state(n, biology, soil) -> AggregateState
    # Uniform profiles at background values
    C  = fill(C_background, n)       # small, near-zero
    B  = fill(biology.B_min, n)      # minimum viable
    F_n = fill(F_n_background, n)    # small
    F_m = fill(F_m_background, n)    # small
    O  = fill(O_amb, n)              # ambient oxygen everywhere
    F_i = fill(F_i_background, n)    # small
    E  = fill(0.0, n)                # no EPS initially
    M  = fill(M_eq_background, n)    # equilibrium MAOC from background C
    P  = biology.P_0                 # full POM mass
    CO2 = 0.0
    return AggregateState(C, B, F_n, F_m, O, F_i, E, M, P, CO2)
end
```

M_eq_background computed from Langmuir-Freundlich isotherm at background C concentration.

---

## 12. Output

Record at user-specified output times (first and last always included):

```julia
struct OutputRecord
    t::Float64                      # time [days]
    state::AggregateState           # full state (deep copy)
    mass_balance_error::Float64     # diagnostic
end
```

Output contains state variables at all grid points. Post-processing (aggregate radius, pool partitioning inside/outside r_agg, etc.) is done on demand from these snapshots — NOT computed during the simulation.

---

## 13. File Structure

```
src/
├── SoilAggregateModel.jl        # Module definition, exports
│
├── types.jl                      # AggregateState, Workspace, TemperatureCache, OutputRecord
├── parameters.jl                 # BiologicalProperties, SoilProperties constructors + defaults
├── environment.jl                # EnvironmentalDrivers{FT,Fψ,FO} + convenience constructors
│
├── temperature/
│   ├── arrhenius.jl              # arrhenius(Ea, T, T_ref) → factor
│   ├── diffusion_pure.jl        # D_O2_water(T), D_DOC_water(T, D_ref, T_ref), D_O2_air(T, D_ref, T_ref)
│   │                            # Stokes-Einstein+VFT, Han-Bartels, Chapman-Enskog
│   └── henry.jl                 # K_H_O2(T) via van't Hoff
│
├── physics/
│   ├── water_retention.jl        # θ(ψ, E, F_i) — modified van Genuchten
│   ├── effective_diffusion.jl    # All D_eff computations (fills workspace arrays)
│   └── aggregate_stability.jl    # r_agg(t) diagnostic (post-processing only)
│
├── biology/
│   ├── bacteria.jl               # R_B, R_Bb, h_B, Γ_B, Γ_E, Resp_B — per node
│   ├── fungi.jl                  # R_F, Π, transitions, Resp_F, Resp_F_conv, h_Fi — per node
│   ├── eps.jl                    # R_rec_E, h_E — per node
│   └── maoc.jl                   # J_M with softplus, M_eq from Langmuir-Freundlich — per node
│
├── carbon/
│   └── pom_dissolution.jl        # R_P (scalar, reads state at node 0)
│
├── solver/
│   ├── tridiagonal.jl            # Thomas algorithm (in-place, overwrites rhs with solution)
│   ├── crank_nicolson.jl         # Assemble tridiagonal + solve for one species, one half-step
│   ├── diffusion_step.jl         # Call crank_nicolson for all 5 diffusing species
│   ├── reactions.jl              # Compute all source/sink terms at one node (no allocation)
│   ├── reaction_step.jl          # Loop over nodes + POM scalar + CO₂ accumulation
│   └── timestepper.jl            # Strang splitting main loop + adaptive Δt
│
└── api.jl                        # run_aggregate() — user entry point, output collection
```

---

## 14. Main Loop

```julia
function run_aggregate!(solver::AggregateSolver, tspan;
                        output_times=Float64[], dt_initial=0.01)
    t, t_end = tspan
    dt = dt_initial
    ws = solver.workspace
    outputs = OutputRecord[]

    while t < t_end
        dt = min(dt, t_end - t)

        # Current environment
        T     = solver.env.T(t)
        ψ     = solver.env.ψ(t)
        O2_bc = solver.env.O2(t)

        # Temperature cache (once per step)
        update_temperature_cache!(ws.f_T, T, solver.biology, solver.soil)

        # Water content (once per step)
        update_water_content!(ws.θ, ws.θ_a, ψ, solver.state, solver.soil)

        # Effective diffusion (once per step)
        update_effective_diffusion!(ws, solver.soil, ws.f_T)

        # === Strang splitting ===
        diffusion_step!(solver, dt/2, O2_bc)     # 5 tridiagonal solves
        reaction_step!(solver, dt)                # pointwise + POM scalar + CO₂
        diffusion_step!(solver, dt/2, O2_bc)     # 5 tridiagonal solves

        t += dt
        dt = adapt_timestep(solver, dt)
        maybe_record!(outputs, solver, t, output_times)
    end

    return outputs
end
```

---

## 15. Implementation Order

Build bottom-up, testing each layer before moving to the next:

1. **types.jl, parameters.jl, environment.jl** — structs with default constructors
2. **temperature/*.jl** — Arrhenius, Stokes-Einstein+VFT, Han-Bartels, Chapman-Enskog, van't Hoff. Unit test each.
3. **solver/tridiagonal.jl** — Thomas algorithm. Test against known solutions.
4. **physics/*.jl** — water retention, effective diffusion. Test with known θ, T.
5. **biology/*.jl** — all reaction terms, one node at a time. Test against manuscript equations with known inputs.
6. **carbon/pom_dissolution.jl** — test with known biomass, θ, O₂.
7. **solver/crank_nicolson.jl + diffusion_step.jl** — test with analytical steady-state diffusion in sphere.
8. **solver/reactions.jl + reaction_step.jl** — test with zero-diffusion scenario (pure batch reactor).
9. **solver/timestepper.jl** — Strang splitting, adaptive Δt. Integration test: 30-day run.
10. **api.jl** — user-facing interface, output, carbon balance check.

Each step should have its own test file before proceeding.

---

## 16. Validation

- **Carbon conservation**: P + ∫(C+B+F_i+F_n+F_m+E+M)·4πr²dr + CO₂_cumulative = P₀ (to machine precision)
- **Regression**: Match v1 code on 30-day benchmark at T = T_ref, constant ψ, constant O₂
- **Steady-state diffusion**: Analytical solution for constant-source sphere
- **Performance**: @allocated == 0 in hot loop; @time < 30s for n=1500, 10-year
