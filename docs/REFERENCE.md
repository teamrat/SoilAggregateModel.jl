# Soil Aggregate Model — Reference Manual

**Last Updated**: 2026-02-09
**Authoritative source for physics**: `manuscript_updated.tex`
**Authoritative source for solver design**: `ARCHITECTURE_CLAUDE_CODE.md`
**Units throughout**: μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol

---

## Table of Contents

### Part I: Variables and Parameters
1. [State Variables](#1-state-variables)
2. [Physical / Diffusion Constants](#2-physical--diffusion-constants)
3. [Biological Parameters](#3-biological-parameters)
4. [Soil Properties](#4-soil-properties)
5. [Environmental Drivers](#5-environmental-drivers)

### Part II: Rates, Fluxes, and Source Terms
6. [Concentration Partitioning (C / C_aq / C_eq)](#6-concentration-partitioning)
7. [Uptake and Assimilation](#7-uptake-and-assimilation)
8. [Recycling and Death](#8-recycling-and-death)
9. [Fungal Transitions](#9-fungal-transitions)
10. [Respiration](#10-respiration)
11. [MAOC Dynamics](#11-maoc-dynamics)
12. [POM Dissolution](#12-pom-dissolution)
13. [Diffusion Coefficients](#13-diffusion-coefficients)
14. [Source/Sink Summary (S terms)](#14-sourcesink-summary)

### Part III: Function Catalog
15. [By Module](#15-function-catalog-by-module)

### Part IV: Computational Methods
16. [Softplus Regularization](#16-softplus-regularization)
17. [Non-Negativity Clipping and CO₂ Correction](#17-non-negativity-clipping-and-co2-correction)
18. [Conservation Weights and Spherical Laplacian](#18-conservation-weights-and-spherical-laplacian)
19. [Strang Splitting and POM-Diffusion Coupling](#19-strang-splitting-and-pom-diffusion-coupling)
20. [Adaptive Timestepping](#20-adaptive-timestepping)
21. [Sigmoid Threshold Functions](#21-sigmoid-threshold-functions)
22. [Space-Limited Yield](#22-space-limited-yield)

### Part V: Conventions
23. [Units](#23-units)
24. [Sign Conventions](#24-sign-conventions)
25. [Naming Conventions](#25-naming-conventions)
26. [Errata and Corrections](#26-errata-and-corrections)

---

# Part I: Variables and Parameters

## 1. State Variables

| Symbol | Code | Description | Units | Domain | BC inner (r = r₀) | BC outer (r = r_max) |
|--------|------|-------------|-------|--------|--------------------|----------------------|
| C | `C` | Total mobile carbon (aqueous + fast-sorbed) | μg-C/mm³ | Spatial PDE | Neumann: −D_C ∂C/∂r = J_P | Neumann: ∂C/∂r = 0 |
| B | `B` | Bacterial biomass | μg-C/mm³ | Spatial PDE | Neumann: ∂B/∂r = 0 | Neumann: ∂B/∂r = 0 |
| F_n | `F_n` | Non-insulated fungi | μg-C/mm³ | Spatial PDE | Neumann: ∂F_n/∂r = 0 | Neumann: ∂F_n/∂r = 0 |
| F_m | `F_m` | Mobile fungi (translocation) | μg-C/mm³ | Spatial PDE | Neumann: ∂F_m/∂r = 0 | Neumann: ∂F_m/∂r = 0 |
| O | `O` | Oxygen (total: aqueous + gas) | μg-O₂/mm³ | Spatial PDE | Neumann: ∂O/∂r = 0 | Dirichlet: O = O_amb(t) |
| F_i | `F_i` | Insulated fungi | μg-C/mm³ | Local ODE | — | — |
| E | `E` | EPS | μg-C/mm³ | Local ODE | — | — |
| M | `M` | MAOC | μg-C/mm³ | Local ODE | — | — |
| P | `P` | POM (well-mixed core) | μg-C | Scalar ODE | — | — |

State vector length: 8n + 1.  CO₂_cumulative is a diagnostic, not a state variable.

**Key**: C is the *total* mobile pool. Microbes see C_aq; MAOC equilibrium sees C_eq. See §6.

---

## 2. Physical / Diffusion Constants

Reference temperature T_ref = 293.15 K (20°C) unless otherwise noted.

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| D_C0_ref | `D_C0_ref` | 86.4 | mm²/day | Glucose in water at T_ref |
| D_O2_w_ref | `D_O2_w_ref` | 173.7 | mm²/day | O₂ in water at T_ref (Han-Bartels) |
| D_O2_a_ref | `D_O2_a_ref` | 1.52×10⁶ | mm²/day | O₂ in air at T_ref (Chapman-Enskog) |
| K_H_O2_ref | — | 31.25 | — | Henry constant at 298.15 K (C_gas/C_aq) |
| ΔH_sol_O2 | — | −12,000 | J/mol | Enthalpy of O₂ dissolution |
| R | `R_GAS` | 8.314 | J/(mol·K) | Gas constant |

Temperature dependence:
- **D_C0(T)**: Stokes-Einstein via VFT viscosity ratio (valid 273–373 K)
- **D_O2_w(T)**: Han & Bartels (1996) empirical formula
- **D_O2_a(T)**: Chapman-Enskog, ∝ (T/T_ref)^1.75
- **K_H(T)**: van't Hoff: `K_H_ref · exp[ΔH_sol/R · (1/T − 1/T_ref)]` — **no leading negative sign**

---

## 3. Biological Parameters

### 3.1 Bacterial

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| r_B_max | `r_B_max` | 2.0 | day⁻¹ | Max specific uptake rate at T_ref |
| K_B | `K_B` | 1.0×10⁻⁴ | μg/mm³ | Half-sat for DOC (uses C_aq) |
| L_B | `L_B` | 1.29×10⁻³ | μg-O₂/mm³ | Half-sat for O₂ |
| ν_B | `ν_B` | 5.8×10⁻⁴ | kPa⁻¹ | Water stress sensitivity (**positive**; exp(ν·ψ) < 1 when ψ < 0) |
| Y_B_max | `Y_B_max` | 0.7 | — | Maximum growth yield |
| K_Y | `K_Y` | — | μg/mm³ | Half-sat for yield (if uptake-dependent) |
| γ | `γ` | 0.2 | — | EPS allocation fraction |
| C_B | `C_B` | — | μg/mm³ | Basal carbon requirement |
| μ_B | `μ_B` | 0.012 | day⁻¹ | Mortality rate at T_ref |
| B_min | `B_min` | 1.0×10⁻⁴ | μg/mm³ | Minimum viable biomass |
| ℰ_a,B | `Ea_B` | 60,000 | J/mol | Activation energy |

### 3.2 Fungal

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| r_F_max | `r_F_max` | 2.0 | day⁻¹ | Max specific uptake rate at T_ref |
| K_F | `K_F` | 1.0×10⁻⁴ | μg/mm³ | Half-sat for DOC (uses C_aq) |
| L_F | `L_F` | 1.29×10⁻³ | μg-O₂/mm³ | Half-sat for O₂ |
| ν_F | `ν_F` | 7.58×10⁻⁵ | kPa⁻¹ | Water stress sensitivity (**positive**; ν_F < ν_B → fungi more drought-tolerant) |
| Y_F | `Y_F` | 0.6 | — | Growth yield |
| μ_F | `μ_F` | 0.012 | day⁻¹ | Mortality rate at T_ref (**only F_i dies**) |
| F_i_min | `F_i_min` | 1.0×10⁻⁴ | μg/mm³ | Minimum viable insulated biomass |
| ℰ_a,F | `Ea_F` | 55,000 | J/mol | Activation energy (shared by ALL fungal rates) |

### 3.3 Fungal Transitions

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| α_i | `α_i` | 0.1 | day⁻¹ | Mobilization rate, insulated |
| α_n | `α_n` | 0.15 | day⁻¹ | Mobilization rate, non-insulated |
| β_i | `β_i` | 0.0 | day⁻¹ | Immobilization rate, insulated |
| β_n | `β_n` | 0.0 | day⁻¹ | Immobilization rate, non-insulated |
| δ | `delta` | 2.0 | — | Mobilization exponent (must be > 1) |
| η | `η_conv` | 0.8 | — | Conversion efficiency; (1−η) lost to respiration |
| ζ | `ζ` | 0.2 | day⁻¹ | Insulation rate (F_n → F_i) |
| λ | `λ` | 0.05 | — | Reduced uptake fraction for insulated hyphae (λ ≪ 1) |
| D_Fn0 | `D_Fn0` | — | mm²/day | Hyphal extension diffusivity at T_ref |
| D_Fm0 | `D_Fm0` | 1.0 | mm²/day | Translocation rate at T_ref (no tortuosity) |
| ε_F | `ε_F` | 1e-10 | μg/mm³ | Π denominator protection |

Fungal uptake: only (λ·F_i + F_n) contributes. All Γ_F enters F_m. Only F_i dies.

### 3.4 EPS

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| μ_E_max | `μ_E_max` | 0.002 | day⁻¹ | Max EPS degradation rate at T_ref |
| K_E | `K_E` | 100·K_B/5 | μg/mm³ | Substrate inhibition threshold |
| E_min | `E_min` | — | μg/mm³ | Minimum EPS for h_E sigmoid |
| ℰ_a,E | `Ea_EPS` | 50,000 | J/mol | Activation energy |

### 3.5 MAOC

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| κ_s_ref | `κ_s_ref` | 0.1 | day⁻¹ | Sorption rate at T_ref |
| κ_d_ref | `κ_d_ref` | 0.01 | day⁻¹ | Desorption rate at T_ref (κ_s > κ_d → hysteresis) |
| ℰ_a,s | `Ea_MAOC_sorb` | 25,000 | J/mol | Activation energy, sorption |
| ℰ_a,d | `Ea_MAOC_desorb` | 40,000 | J/mol | Activation energy, desorption (ℰ_a,d > ℰ_a,s → warming narrows hysteresis) |
| ε_maoc | `ε_maoc` | 0.01 | μg/mm³ | Softplus smoothing width |

### 3.6 POM

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| R_P^max | `R_P_max` | 1.0 | μg-C/mm²/day | Surface-specific max dissolution rate at T_ref |
| θ_P | `θ_P` | 0.1 | — | Half-sat water content |
| L_P | `L_P` | 1.29×10⁻³ | μg-O₂/mm³ | Half-sat O₂ |
| K_B,P | `K_B_P` | — | μg/mm³ | Half-sat bacteria for dissolution |
| K_F,P | `K_F_P` | — | μg/mm³ | Half-sat fungi for dissolution |
| ρ_POM | `ρ_POM` | 200 | μg-C/mm³ | POM carbon density |
| ℰ_a,P | `Ea_POM` | 60,000 | J/mol | Activation energy |

### 3.7 Temperature Response (Arrhenius)

All biological rates scale via:

$$k(T) = k_\text{ref} \cdot \exp\!\left[\frac{\mathcal{E}_a}{R}\left(\frac{1}{T_\text{ref}} - \frac{1}{T}\right)\right]$$

Six distinct activation energies: ℰ_a,B (bacteria), ℰ_a,F (fungi — shared by all fungal rates), ℰ_a,E (EPS), ℰ_a,s (MAOC sorption), ℰ_a,d (MAOC desorption), ℰ_a,P (POM dissolution).

### 3.8 Oxygen

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| α_O | `α_O` | 2.2 | μg-O₂/μg-C | Respiratory quotient |

---

## 4. Soil Properties

### 4.1 Hydraulic (van Genuchten)

| Symbol | Code | Sandy loam | Clay loam | Loam | Units |
|--------|------|------------|-----------|------|-------|
| θ_s | `θ_s` | 0.5 | 0.43 | 0.43 | — |
| θ_r | `θ_r` | 0.06 | 0.095 | 0.078 | — |
| α_VG | `α_vg` | 0.1133 | 0.01899 | 0.03598 | kPa⁻¹ |
| n_VG | `n_vg` | 1.47 | 1.31 | 1.56 | — |

Water content: θ(ψ) = θ_r + (θ_s − θ_r)·[1 + (α·|ψ|)^n]^{−m}, m = 1 − 1/n.

Modified by EPS and insulated fungi: α_eff = α₀ · exp(ω_E·E + ω_F·F_i), where ω_E, ω_F < 0 (EPS/fungi increase water retention by decreasing α).

### 4.2 Physical

| Symbol | Code | Sandy loam | Units |
|--------|------|------------|-------|
| ρ_b | `ρ_b` | 1500 | μg/mm³ (= kg/m³) |
| f_clay | `clay_fraction` | 0.18 | — |
| f_silt | `silt_fraction` | 0.22 | — |

### 4.3 MAOC Sorption

| Symbol | Code | Sandy loam | Units | Notes |
|--------|------|------------|-------|-------|
| k_d | `k_d_eq` | 0.05 | mm³/μg | Linear partition coefficient |
| M_max | `M_max` | 10.0 | μg-C/mm³ | Max MAOC capacity |
| k_L | `k_L` | 10.0 | mm³/μg | Langmuir affinity |
| n_LF | `n_LF` | 0.7 | — | Freundlich exponent (< 1 → heterogeneous sites) |
| k_ma | `k_ma` | 0.48 | μg-C/g-mineral | Mineral activity coefficient |
| f_cs | `f_clay_silt` | 0.40 | — | Clay + silt mass fraction |

Clay-dependent capacity: M_max = k_ma · ρ_b · f_clay_silt.

### 4.4 Aggregate Stability

| Symbol | Code | Units | Notes |
|--------|------|-------|-------|
| k_F | `k_F` | kPa/(μg/mm³) | Specific binding strength |
| χ | `χ` | mm | Particle adhesion length scale |
| a_p | `a_p` | mm | Particle radius |

---

## 5. Environmental Drivers

```julia
struct EnvironmentalDrivers{FT, Fψ, FO}
    T::FT       # T(t) → temperature [K]
    ψ::Fψ       # ψ(t) → matric potential [kPa]  (negative in unsaturated soil)
    O2::FO      # O2(t) → boundary O₂ concentration [μg/mm³]
end
```

Parametric types: constants are as fast as literals.

---

# Part II: Rates, Fluxes, and Source Terms

## 6. Concentration Partitioning

The single state variable C represents total mobile carbon. Three derived concentrations:

```
C      = state variable (total mobile carbon per mm³ bulk soil)     [μg/mm³]
C_aq   = C / (θ + ρ_b · k_d)                                      [μg/mm³ water]
C_eq   = k_d · C_aq = k_d · C / (θ + ρ_b · k_d)                  [μg/g solid]
```

| Where used | Concentration |
|-----------|---------------|
| State variable / diffusion PDE | C |
| Bacterial uptake R_B (Monod) | C_aq |
| Fungal uptake R_F (Monod) | C_aq |
| EPS degradation inhibition | C_aq |
| Bacterial maintenance R_Bb | C_B (constant parameter) |
| MAOC equilibrium M_eq | C_eq |
| S_C coupling to MAOC | −J_M (no factor) |

---

## 7. Uptake and Assimilation

**Bacterial uptake** (per node):
$$R_B = r_{B,\max}(T) \cdot \frac{C_{aq}}{K_B + C_{aq}} \cdot \frac{O_{aq}}{L_B + O_{aq}} \cdot B \cdot \exp(\nu_B \, \psi) \cdot h_B$$

**Bacterial maintenance**:
$$R_{B,b} = r_{B,\max}(T) \cdot \frac{C_B}{K_B + C_B} \cdot \frac{O_{aq}}{L_B + O_{aq}} \cdot B \cdot h_B$$

Growth vs starvation: if R_B > R_Bb (growth regime), then Γ_B = Y_B · R_B · (1−γ), Γ_E = Y_B · R_B · γ. If R_B ≤ R_Bb (starvation), Γ_B = R_B, Γ_E = 0.

**Fungal uptake** (per node):
$$R_F = r_{F,\max}(T) \cdot \frac{C_{aq}}{K_F + C_{aq}} \cdot \frac{O_{aq}}{L_F + O_{aq}} \cdot (\lambda \, F_i + F_n) \cdot \exp(\nu_F \, \psi)$$

All Γ_F = Y_F · R_F enters F_m.

---

## 8. Recycling and Death

| Rate | Formula | Notes |
|------|---------|-------|
| R_rec,B | μ_B(T) · B · h_B | Bacterial mortality |
| R_rec,F | μ_F(T) · F_i · h_Fi | **Only F_i dies** |
| R_rec,E | μ_E(T) · K_E/(K_E + C_aq) · E · h_E | Uses **C_aq**, not C |

R_rec = R_rec,B + R_rec,F + R_rec,E (total recycling back to C pool).

---

## 9. Fungal Transitions

Relative mobile fraction: Π = F_m / (F_i + F_n + ε_F)

Net transfer functions (already contain η):
```
net_i = η · (β_i · Π − α_i · Π^δ) · F_i     (positive = immobilization)
net_n = η · (β_n · Π − α_n · Π^δ) · F_n     (positive = immobilization)
```

Insulation: trans_ni = ζ · F_n (one-way, F_n → F_i)

Source terms from transitions:
- S_Fi gets: +net_i + trans_ni
- S_Fn gets: +net_n − trans_ni
- S_Fm gets: −net_i − net_n

---

## 10. Respiration

| Rate | Formula | Notes |
|------|---------|-------|
| Resp_B | (1 − Y_B) · R_B | Growth: uptake minus yield. Starvation: 0. |
| Resp_F | (1 − Y_F) · R_F | Fungal uptake respiration |
| Resp_F^conv | (1−η) · \|net_i/η + net_n/η\| | Uses **abs()** to prevent negative respiration |

Total O₂ consumption: S_O = −α_O · (Resp_B + Resp_F + Resp_F^conv)

**Carbon conservation identity** (per node, excluding diffusion):
$$S_C + S_B + S_{F_i} + S_{F_n} + S_{F_m} + S_E + S_M = -(Resp_B + Resp_F + Resp_F^{conv})$$

---

## 11. MAOC Dynamics

**Two-stage model**:

Stage 1 (fast equilibrium): C ↔ C_eq via k_d (instantaneous)

Stage 2 (slow kinetic): C_eq → M via Langmuir-Freundlich isotherm

$$M_{eq} = M_{\max} \cdot \frac{(k_L \cdot C_{eq})^{n_{LF}}}{1 + (k_L \cdot C_{eq})^{n_{LF}}}$$

Rate law with softplus regularization:
$$J_M = \kappa_s(T) \cdot \varphi_\varepsilon(M_{eq} - M) \;-\; \kappa_d(T) \cdot \varphi_\varepsilon(M - M_{eq})$$

where φ_ε(x) = ε · ln(1 + e^{x/ε}) ≈ max(0, x) smoothly.

**Coupling to S_C**: `S_C = ... − J_M` (just J_M, no (θ+ρ_b·k_d)/k_d factor).

---

## 12. POM Dissolution

**Flux density** at POM surface (r = r₀):
$$J_P = R_P^{\max}(T) \cdot \frac{P}{P_0} \cdot \frac{B_0}{K_{B,P}+B_0} \cdot \frac{F_{n,0}}{K_{F,P}+F_{n,0}} \cdot \frac{\theta_0}{\theta_P+\theta_0} \cdot \frac{O_{aq,0}}{L_P+O_{aq,0}}$$

Subscript 0 = values at first grid node. J_P has units [μg-C/mm²/day].

**Total dissolution rate**: R_P = 4πr₀² · J_P [μg-C/day]

**Where each is used**:
- POM ODE: dP/dt = −R_P
- Inner Neumann BC for C: −D_C ∂C/∂r|_{r₀} = J_P
- CO₂ diagnostic: R_P is total carbon entering domain per day

---

## 13. Diffusion Coefficients

Computed once per timestep from θ(r) and temperature-dependent pure-phase values.

| Species | D_eff formula | Notes |
|---------|--------------|-------|
| C | D_C0(T) · τ_aq · θ/(θ + ρ_b·k_d) | Retarded by sorption |
| B | D_B_rel · D_C | Chemotaxis, proportional to D_C |
| F_n | D_Fn0 · f_fun(T) · τ_aq | Hyphal tip extension |
| F_m | D_Fm0 · f_fun(T) | **No tortuosity** — internal translocation |
| O | D_O2_w(T)·θ·τ_aq/(θ+K_H·θ_a) + D_O2_a(T)·K_H·θ_a·τ_gas/(θ+K_H·θ_a) | Aqueous + gas phase |

Tortuosity:
- Aqueous: τ_aq = θ²/θ_s^(2/3)  (deliberate simplification, not standard M-Q)
- Gas: τ_gas = θ_a^(10/3)/θ_s²  (standard Millington-Quirk)

---

## 14. Source/Sink Summary

```
S_P  = −R_P

S_C  = −R_B − R_F + R_rec − J_M     (+ diffusion + POM Neumann BC)

S_B  = Γ_B − R_rec,B

S_Fi = ζ·F_n + η·(β_i·Π − α_i·Π^δ)·F_i − R_rec,F

S_Fn = η·(β_n·Π − α_n·Π^δ)·F_n − ζ·F_n

S_Fm = Γ_F − η·(β_i·Π − α_i·Π^δ)·F_i − η·(β_n·Π − α_n·Π^δ)·F_n − Resp_F^conv   (+ diffusion)

S_E  = Γ_E − R_rec,E

S_M  = J_M

S_O  = −α_O·(Resp_B + Resp_F + Resp_F^conv)     (+ diffusion)
```

**Conservation**: S_P + ∫(S_C + S_B + S_Fi + S_Fn + S_Fm + S_E + S_M)·4πr²dr = −∫(Resp_B + Resp_F + Resp_F^conv)·4πr²dr

---

# Part III: Function Catalog

## 15. Function Catalog by Module

### temperature/arrhenius.jl
| Function | Signature | Returns |
|----------|-----------|---------|
| `arrhenius` | `(Ea, T, T_ref)` | Dimensionless scaling factor |

### temperature/diffusion_pure.jl
| Function | Returns | Method |
|----------|---------|--------|
| `D_DOC_water(T, D_ref, T_ref)` | mm²/day | Stokes-Einstein + VFT viscosity |
| `D_O2_water(T)` | mm²/day | Han & Bartels (1996) |
| `D_O2_air(T, D_ref, T_ref)` | mm²/day | Chapman-Enskog T^1.75 |
| `water_viscosity_VFT(T)` | ln(η/mPa·s) | Vogel-Fulcher-Tammann |
| `D_fungal_translocation(T, D_ref, T_ref, Ea)` | mm²/day | Arrhenius |

### temperature/henry.jl
| Function | Returns | Notes |
|----------|---------|-------|
| `henry_vant_hoff(K_H_ref, ΔH_sol, T, T_ref)` | Dimensionless K_H | **No leading negative sign** |
| `K_H_O2(T)` | Dimensionless K_H | Convenience for O₂ |
| `O2_saturation(T, P_atm=0.21)` | μg-O₂/mm³ | Air-saturated aqueous concentration |

### physics/water_retention.jl
| Function | Returns |
|----------|---------|
| `θ(ψ, E, F_i, soil)` | Water content (van Genuchten, modified α) |

### physics/effective_diffusion.jl
| Function | Returns |
|----------|---------|
| `D_eff_DOC(D_w, θ, θ_s, ρ_b, k_d)` | Effective DOC diffusion |
| `D_eff_O2(D_w, D_a, θ, θ_s, K_H)` | Effective O₂ diffusion (aq + gas) |
| `D_eff_bacteria(D_C, D_B_rel)` | Effective bacterial motility |
| `D_eff_fungi_noninsulated(D_Fn0, f_T, θ, θ_s)` | Effective F_n diffusion |
| `D_eff_fungi_mobile(D_Fm0, f_T)` | Effective F_m diffusion (no tortuosity) |

### biology/bacteria.jl
| Function | Computes |
|----------|----------|
| `R_B(C_aq, O_aq, B, ψ, ...)` | Bacterial uptake rate |
| `R_Bb(C_B, O_aq, B, ...)` | Bacterial maintenance rate |
| `h_B(B, B_min)` | Non-negativity sigmoid |
| `Y_B_func(C_aq, ...)` | Yield (constant or uptake-dependent) |
| `Gamma_B(R_B, R_Bb, Y_B, γ)` | Growth allocation |
| `Gamma_E(R_B, R_Bb, Y_B, γ)` | EPS allocation |
| `Resp_B(R_B, R_Bb, Y_B)` | Bacterial respiration |
| `R_rec_B(B, μ_B, h_B)` | Bacterial death/recycling |

### biology/fungi.jl
| Function | Computes |
|----------|----------|
| `Pi_protected(F_i, F_n, F_m, ε_F)` | Π = F_m/(F_i+F_n+ε_F) |
| `R_F(C_aq, O_aq, F_i, F_n, ψ, λ, ...)` | Fungal uptake rate |
| `Y_F(...)` | Fungal yield |
| `Gamma_F(R_F, Y_F)` | Fungal growth (→ F_m) |
| `Resp_F(R_F, Y_F)` | Fungal uptake respiration |
| `h_Fi(F_i, F_i_min)` | Non-negativity sigmoid |
| `R_rec_F(F_i, μ_F, h_Fi)` | Fungal death (F_i only) |
| `fungal_transitions(F_i, F_n, F_m, Π, ...)` | Returns (net_i, net_n, trans_ni) |
| `Resp_F_conv(net_i, net_n, η)` | Conversion respiration, uses **abs()** |

### biology/eps.jl
| Function | Computes |
|----------|----------|
| `h_E(E, E_min)` | Non-negativity sigmoid |
| `R_rec_E(E, C_aq, K_E, μ_E, h_E)` | EPS recycling (uses **C_aq**) |

### biology/maoc.jl
| Function | Computes |
|----------|----------|
| `softplus(x, ε)` | ε·ln(1+exp(x/ε)), numerically stable |
| `M_eq_langmuir_freundlich(C_eq, M_max, k_L, n_LF)` | MAOC equilibrium capacity |
| `J_M(M, M_eq, κ_s, κ_d, ε_maoc)` | Net MAOC formation rate |

### solver/tridiagonal.jl
| Function | Notes |
|----------|-------|
| `thomas_solve!(a, b, c, d)` | In-place Thomas algorithm, overwrites d with solution |

---

# Part IV: Computational Methods

This section documents numerical techniques used in the implementation that depart from, extend, or supplement the manuscript equations. These are not approximations introduced for convenience — each addresses a specific mathematical or physical requirement (smoothness, positivity, conservation) and has been verified to preserve the model's conservation identity to machine precision.

---

## 16. Softplus Regularization

**Problem.** The manuscript uses piecewise-linear operations — $\max(0, x)$ and $\min(0, x)$ — to switch between growth and starvation regimes in the bacterial allocation equations. These have discontinuous derivatives at $x = 0$, which introduces stiffness: the ODE system's Jacobian jumps, forcing the adaptive timestepper to take unnecessarily small steps near the transition.

**Solution.** Replace $\max(0, x)$ with the softplus function and $\min(0, x)$ with its complement:

$$
\varphi_\varepsilon(x) = \varepsilon \ln\!\bigl(1 + e^{x/\varepsilon}\bigr)
$$

implemented in numerically stable form:

$$
\varphi_\varepsilon(x) = \begin{cases}
x + \varepsilon \ln(1 + e^{-x/\varepsilon}) & x > 0 \\[4pt]
\varepsilon \ln(1 + e^{x/\varepsilon}) & x \leq 0
\end{cases}
$$

The key substitutions are:

| Manuscript | Implementation |
|------------|----------------|
| $\max(0, R_{\text{diff}})$ | $\varphi_\varepsilon(R_{\text{diff}})$ |
| $\min(0, R_{\text{diff}})$ | $-\varphi_\varepsilon(-R_{\text{diff}})$ |

where $R_{\text{diff}} = R_B - R_{B,b}$ (uptake minus maintenance).

**Conservation identity.** The softplus satisfies the exact algebraic identity:

$$
\varphi_\varepsilon(x) - \varphi_\varepsilon(-x) = x
$$

for all $x$ and $\varepsilon > 0$. This is not a numerical approximation — it holds to machine precision because:

$$
\varepsilon\ln(1+e^{x/\varepsilon}) - \varepsilon\ln(1+e^{-x/\varepsilon}) = \varepsilon\ln\frac{1+e^{x/\varepsilon}}{1+e^{-x/\varepsilon}} = \varepsilon \cdot \frac{x}{\varepsilon} = x
$$

As a consequence, the bacterial carbon balance

$$
\Gamma_B + \Gamma_E + \text{Resp}_B = R_B
$$

is preserved exactly under softplus substitution, because all allocation terms decompose into $\varphi_\varepsilon(R_{\text{diff}})$ and $\varphi_\varepsilon(-R_{\text{diff}})$ whose sum telescopes to $R_{\text{diff}}$.

**Where used.** Bacterial allocation ($\Gamma_B$, $\Gamma_E$, $\text{Resp}_B$), yield function ($Y_B$), and MAOC sorption/desorption switching ($J_M$). The smoothing width $\varepsilon_Y = K_Y / 100$ is small enough that the softplus is indistinguishable from $\max(0, \cdot)$ outside a narrow transition zone of width $\sim 3\varepsilon$.

---

## 17. Non-Negativity Clipping and CO₂ Correction

**Problem.** Forward Euler can overshoot, driving a pool negative when the consumption rate exceeds the available substrate within one timestep. Negative concentrations are unphysical and, if left uncorrected, contaminate subsequent source-term evaluations.

**Solution.** After each Forward Euler update, negative pools are clipped to zero. The critical subtlety is in the CO₂ bookkeeping.

**The accounting logic.** Consider a node where the Forward Euler update produces $C^{n+1} = C^n + \Delta t \cdot S_C < 0$. The source terms satisfy the conservation identity:

$$
\sum_X S_X + \text{Resp}_{\text{total}} = 0
$$

so the CO₂ accumulation line $\text{CO}_2 \mathrel{+}= \Delta t \cdot \text{Resp}_{\text{total}} \cdot V_i$ already accounts for all carbon leaving the pools. The conservation identity holds exactly after the Euler update — total carbon (pools + CO₂) is unchanged.

Now clipping sets $C$ from a negative value to zero, increasing the pool's carbon content by $|C^{n+1}| \cdot V_i$. To maintain conservation, CO₂ must decrease by the same amount:

```julia
if state.C[i] < 0.0
    clip_carbon += state.C[i] * volume_i   # negative value → reduces CO₂
    state.C[i] = 0.0
end
...
state.CO2_cumulative += dt * Resp_total * volume_i + clip_carbon
```

Note that `clip_carbon` carries the **signed** (negative) value of the clipped pool times volume. This subtracts from CO₂, exactly compensating the carbon "created" by zeroing the negative pool.

**Physical interpretation.** Forward Euler over-estimated respiration because it assumed constant rates throughout the step. In reality, as the substrate approached zero, consumption would have slowed. The CO₂ correction removes the excess respiration that could not have physically occurred — it is not a fudge but the correct mass balance adjustment for an explicit integrator that overshoots.

**The original bug (diagnosed February 2026).** The initial implementation used `abs(state.C[i])` instead of `state.C[i]`, which added a positive value to CO₂ rather than subtracting. This doubled the clipped carbon: the pool gained $|C^{n+1}| \cdot V_i$ (from negative to zero) and CO₂ also gained $|C^{n+1}| \cdot V_i$, creating $2|C^{n+1}| \cdot V_i$ of carbon from nothing. The bug produced a conservation error of approximately 0.014%/year under adaptive timestepping (0.7%/year under fixed $\Delta t = 0.001$ d), was invisible at low biomass (no clipping events), and was masked by the adaptive stepper's tendency to take small steps that minimized overshoot. It was isolated through systematic elimination of all other conservation pathways (diffusion, source-term identity, POM coupling) and confirmed by tracing per-step carbon budgets at high biomass density.

---

## 18. Conservation Weights and Spherical Laplacian

**Problem.** In spherical geometry, the natural cell volume for node $i$ is the shell $(4\pi/3)(r_{i+1/2}^3 - r_{i-1/2}^3)$. However, using geometrically exact shell volumes for discrete integration does not produce exact conservation with the Crank–Nicolson stencil, because the ratio $V_i / (r_i^2 h^2)$ varies across nodes.

**Solution.** The conservation weights are chosen to match the discrete Laplacian stencil:

$$
W_i = 4\pi \, r_i^2 \, h
$$

These are **not** the true shell volumes — they are the weights that make the discrete conservation identity hold exactly. The spherical Laplacian stencil at node $i$ is:

$$
L_i = \frac{1}{r_i^2 h^2}\Bigl[r_{i+1/2}^2 D_{i+1/2}(u_{i+1} - u_i) - r_{i-1/2}^2 D_{i-1/2}(u_i - u_{i-1})\Bigr]
$$

Multiplying by $W_i = 4\pi r_i^2 h$ gives:

$$
W_i \cdot L_i = \frac{4\pi}{h}\Bigl[r_{i+1/2}^2 D_{i+1/2}(u_{i+1} - u_i) - r_{i-1/2}^2 D_{i-1/2}(u_i - u_{i-1})\Bigr]
$$

Summing over all nodes, each interior flux appears once with $+$ and once with $-$, producing exact telescoping:

$$
\sum_i W_i \cdot L_i = \frac{4\pi}{h}\bigl[\text{outer boundary flux} - \text{inner boundary flux}\bigr]
$$

This guarantees that diffusion conserves mass exactly (to machine precision) when integrated with weights $W_i$. The same weights must be used everywhere: in `compute_total_carbon`, in the CO₂ accumulation (`volume_i`), and in post-processing integrals.

**Implementation note.** The weights are precomputed in `GridInfo` and stored as `grid.W`. A separate `compute_cell_volumes` function exists for geometrically exact shell volumes (used only for visualization and post-processing, never for conservation accounting).

---

## 19. Strang Splitting and POM-Diffusion Coupling

**Architecture.** Each timestep applies second-order Strang splitting:

1. Diffusion half-step ($\Delta t / 2$): Crank–Nicolson for $C$, $B$, $F_n$, $F_m$, $O$
2. Reaction full-step ($\Delta t$): Forward Euler for all pools + POM scalar + CO₂
3. Diffusion half-step ($\Delta t / 2$): identical to step 1

This decouples transport (tridiagonal, $O(n)$ per species) from reactions (pointwise, $O(n)$), avoiding a monolithic $8n+1$ implicit system.

**POM dissolution coupling.** POM dissolution creates a flux $J_P$ at the inner boundary that enters the dissolved carbon pool $C$ via the Neumann boundary condition. The same dissolution rate decreases the POM scalar in the reaction step. To ensure exact mass balance:

- $J_P$ is computed **once** at the beginning of each timestep from the current state
- The **same** $J_P$ is used for both diffusion half-steps (as a Neumann BC) and for the POM decrease ($P \mathrel{-}= \Delta t \cdot R_P$ where $R_P = J_P \cdot 4\pi r_0^2$)
- Carbon entering $C$ through the boundary = $J_P \cdot 4\pi r_0^2 \cdot \Delta t$ (summed over both half-steps)
- Carbon leaving $P$ = $R_P \cdot \Delta t = J_P \cdot 4\pi r_0^2 \cdot \Delta t$
- These are identical by construction, regardless of whether $J_P$ is "stale" relative to the post-reaction state

Freezing $J_P$ does introduce a splitting error (the solution trajectory differs from the unsplit system), but it does **not** introduce a conservation error. The splitting error is $O(\Delta t^2)$, consistent with Strang splitting.

---

## 20. Adaptive Timestepping

**Criterion.** After each reaction step, the maximum relative change across all nodes and species is computed:

$$
\rho = \max_{i, X} \frac{|S_X \cdot \Delta t|}{\max(u_X^{(i)},\, \tau)}
$$

where $\tau = 10^{-6}$ prevents division by near-zero values.

**Adjustment.**

| Condition | Action |
|-----------|--------|
| $\rho > 0.10$ | Halve $\Delta t$ |
| $\rho < 0.01$ | Double $\Delta t$ |
| otherwise | Keep $\Delta t$ |

Bounds: $\Delta t_{\min} = 10^{-4}$ d, $\Delta t_{\max} = 0.1$ d.

**Note on conservation.** The adaptive stepper does **not** reject and re-do steps. When $\rho$ exceeds the threshold, the current step is accepted and the next step uses the reduced $\Delta t$. This is a "predict then correct" strategy: the occasional large-$\rho$ step introduces $O(\Delta t)$ trajectory error but does not break conservation (the per-step conservation identity holds for any $\Delta t$). In practice, the stepper maintains $\rho \approx 0.01\text{–}0.10$, corresponding to 1–10% relative change per step.

**Interaction with clipping.** Smaller $\Delta t$ means less overshoot, fewer clipping events, and better trajectory accuracy. This is why the adaptive stepper produced 50$\times$ better conservation (0.014%/year) than fixed $\Delta t = 0.001$ (0.7%/year) before the clipping bug was fixed — the adaptive stepper was taking smaller steps in regions of stiff dynamics, reducing the number of clipping corrections needed.

---

## 21. Sigmoid Threshold Functions

**Problem.** Biological rates (mortality, maintenance, insulation) should vanish smoothly as biomass approaches a minimum viable threshold, not discontinuously.

**Solution.** Three sigmoid functions enforce soft lower bounds:

$$
h(x) = \frac{\exp(\beta\, x)}{\exp(\beta\, x) + \exp(\beta\, x_{\min})}
$$

with $\beta = 50 / x_{\min}$, giving a smooth transition from 0 to 1 over a width of approximately $4 x_{\min} / 50$. Applied to:

| Function | Variable | Threshold | Effect |
|----------|----------|-----------|--------|
| $h_B(B)$ | Bacteria | $B_{\min}$ | Shuts off maintenance $R_{B,b}$ and mortality near extinction |
| $h_{F_i}(F_i)$ | Insulated fungi | $F_{i,\min}$ | Shuts off fungal mortality near extinction |
| $h_E(E)$ | EPS | $E_{\min}$ | Shuts off EPS recycling near depletion |

These prevent numerical extinction artifacts where a pool oscillates around zero due to competing production and consumption terms. The sigmoid ensures a smooth, monotonic approach to zero rather than oscillatory overshoot.

---

## 22. Space-Limited Yield

**Problem.** Without a carrying capacity mechanism, bacterial biomass grows exponentially until substrate is exhausted, producing unrealistic densities ($> 6$ $\mu$g/mm$^3$, equivalent to $> 6$ kg/m$^3$ bacterial carbon).

**Solution.** Growth yield decreases with local biomass density via a Monod-form space limitation factor:

$$
Y_B = Y_{B,\max} \cdot \frac{\varphi_\varepsilon(R_{\text{diff}})}{\varphi_\varepsilon(R_{\text{diff}}) + K_Y} \cdot \frac{B_S}{B + B_S}
$$

At $B \gg B_S$, yield approaches zero and all uptake is respired as maintenance. This produces an emergent carrying capacity that depends on local substrate and oxygen availability, rather than a hard cap. The analogous fungal form uses total fungal biomass $F_i + F_n + F_m$ with half-saturation $F_S$.

**Conservation.** Because reduced yield routes carbon to respiration rather than growth, the conservation identity $\Gamma_B + \Gamma_E + \text{Resp}_B = R_B$ continues to hold exactly. No carbon is created or destroyed by space limitation — it only changes the partition between biomass and CO₂.

---

# Part V: Conventions

## 23. Units

| Quantity | Unit | Equivalence |
|----------|------|-------------|
| Concentration | μg/mm³ | = kg/m³ = g/L |
| Length | mm | |
| Time | day | |
| Temperature | K | 293.15 = 20°C |
| Water potential | kPa | negative in unsaturated soil |
| Activation energy | J/mol | |
| Diffusion | mm²/day | |
| Rate constant | day⁻¹ | |

**No unit conversions inside kernels.**

---

## 24. Sign Conventions

- **State variables**: always ≥ 0
- **ψ (water potential)**: negative in unsaturated soil (e.g., −30 kPa at field capacity)
- **ν_B, ν_F**: positive values; stress factor = exp(ν·ψ) < 1
- **ω_E, ω_F**: negative values; EPS/fungi decrease α (increase retention)
- **R terms** (R_B, R_F, R_rec, R_P): always positive (magnitude of process)
- **J_M**: positive when sorption dominates, negative when desorption dominates
- **S terms**: sign indicates net production (+) or consumption (−) of that pool
- **dP/dt**: always ≤ 0 (POM only decreases)

---

## 25. Naming Conventions

| Pattern | Meaning |
|---------|---------|
| `R_*` | Rate (always positive magnitude) |
| `S_*` | Source/sink (signed, in PDE RHS) |
| `Gamma_*` / `Γ_*` | Assimilation / growth allocation |
| `Resp_*` | Respiration (always positive) |
| `R_rec_*` | Recycling / death rate |
| `D_*` | Diffusion coefficient |
| `h_*` | Non-negativity sigmoid |
| `f_*` | Arrhenius temperature factor |
| `_ref` suffix | Value at reference temperature |

---

## 26. Errata and Corrections

Issues found during 2026-02-05/06 audit. The **code** has the corrected versions; the manuscript may still need updating.

| # | Issue | Manuscript | Code (correct) |
|---|-------|-----------|----------------|
| 1 | Henry's law sign | Extra negative in van't Hoff exponent | `exp[ΔH_sol/R · (1/T − 1/T_ref)]`, no leading negative |
| 2 | MAOC coupling in S_C | −J_M · (θ+ρ_b·k_d)/k_d | −J_M (no factor; J_M is per mm³ bulk soil) |
| 3 | EPS degradation substrate | K_E/(K_E + C) | K_E/(K_E + C_aq) |
| 4 | Resp_F^conv sign | No abs() | abs() wraps net transfer sum |
| 5 | POM notation | Ambiguous R_P units | J_P [μg/mm²/day] vs R_P = 4πr₀²·J_P [μg/day] |

---

**End of Reference**
