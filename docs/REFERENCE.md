# Soil Aggregate Model — Reference Manual

**Last Updated**: 2026-02-06  
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

### Part IV: Conventions
16. [Units](#16-units)
17. [Sign Conventions](#17-sign-conventions)
18. [Naming Conventions](#18-naming-conventions)
19. [Errata and Corrections](#19-errata-and-corrections)

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
| λ | `λ` | 0.05 | — | Fraction of F_n at uptake surfaces (λ ≪ 1) |
| D_Fn0 | `D_Fn0` | — | mm²/day | Hyphal extension diffusivity at T_ref |
| D_Fm0 | `D_Fm0` | 1.0 | mm²/day | Translocation rate at T_ref (no tortuosity) |
| ε_F | `ε_F` | 1e-10 | μg/mm³ | Π denominator protection |

Fungal uptake: only (F_i + λ·F_n) contributes. All Γ_F enters F_m. Only F_i dies.

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
$$R_F = r_{F,\max}(T) \cdot \frac{C_{aq}}{K_F + C_{aq}} \cdot \frac{O_{aq}}{L_F + O_{aq}} \cdot (F_i + \lambda \, F_n) \cdot \exp(\nu_F \, \psi)$$

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

# Part IV: Conventions

## 16. Units

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

## 17. Sign Conventions

- **State variables**: always ≥ 0
- **ψ (water potential)**: negative in unsaturated soil (e.g., −30 kPa at field capacity)
- **ν_B, ν_F**: positive values; stress factor = exp(ν·ψ) < 1
- **ω_E, ω_F**: negative values; EPS/fungi decrease α (increase retention)
- **R terms** (R_B, R_F, R_rec, R_P): always positive (magnitude of process)
- **J_M**: positive when sorption dominates, negative when desorption dominates
- **S terms**: sign indicates net production (+) or consumption (−) of that pool
- **dP/dt**: always ≤ 0 (POM only decreases)

---

## 18. Naming Conventions

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

## 19. Errata and Corrections

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
