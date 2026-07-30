# Soil Aggregate Model — Reference Manual

**Last Updated**: 2026-07-29

**Physics** — the manuscript on Overleaf, `Dropbox/Apps/Overleaf/Aggregation
2026 rev3/`: `manuscript-4-5.tex` (main text, the current revision — 4-2, 4-3,
4-4 and `manuscript.tex` are older) and `SI_tessellation.tex` (domain
tessellation and the ω overlap correction). `manuscript_updated.tex` does not
exist and never did.

**Everything is working draft.** The manuscript, this file and the code are
three descriptions of one model, none of them senior. Where they disagree,
record it — §26 — rather than assuming the other one is right.

**Solver design** — `docs/ARCHITECTURE.md`, formerly `ARCHITECTURE_CLAUDE_CODE.md`
(renamed 2026-02, see `dev_notes/archive/REPO_REORGANIZATION.md`). **Stale:** it
documents an `AggregateSolver` type that does not exist, omits the default
solver entirely, and its `SoilProperties` listing still carries fields removed
by erratum 11. `docs/AUDIT.md` tracks the decision on whether to rewrite or
retire it. Do not treat it as current.
**Units throughout**: μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol

---

## Table of Contents

### Part I: Variables and Parameters
1. [State Variables](#1-state-variables)
2. [Physical / Diffusion Constants](#2-physical--diffusion-constants)
3. [Biological Parameters](#3-biological-parameters)
4. [Soil Properties](#4-soil-properties)
5. [Environmental Drivers](#5-environmental-drivers)
5a. [Parameter Provenance and Anchors](#5a-parameter-provenance-and-anchors)
5b. [Domain Tessellation and Population Upscaling](#5b-domain-tessellation-and-population-upscaling)
5c. [Aggregate Mass, Sieve Classes and Sample Statistics](#5c-aggregate-mass-sieve-classes-and-sample-statistics)

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
17a. [Carbon Closure, and What Actually Verifies It](#17a-carbon-closure-and-what-actually-verifies-it)
18. [Conservation Weights and Spherical Laplacian](#18-conservation-weights-and-spherical-laplacian)
19. [Strang Splitting and POM-Diffusion Coupling](#19-strang-splitting-and-pom-diffusion-coupling)
20. [Adaptive Timestepping](#20-adaptive-timestepping)
20a. [Method of Lines and the Stiff Solver](#20a-method-of-lines-and-the-stiff-solver)
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
- **K_H(T)**: van't Hoff: `K_H_ref · exp[ΔH_sol/R · (1/T − 1/T_ref)]` — **no leading negative sign**. This is `arrhenius_ratio(ΔH_sol, T, T_ref)`, and `henry_vant_hoff` calls it: van't Hoff and Arrhenius are one exponential with ΔH_sol in place of Ea, so the package has one and not three.

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
Naming follows Falconer (2005) eq. 2.4 / (2008) Box 1: **α is immobilization**
(mobile → sessile, the gain, and it carries the exponent); **β is mobilization**
(sessile → mobile, the loss, linear in Π). These were inverted before 2026-07.
See §5a for provenance.

| Symbol | Code | Value | Units | Notes |
|--------|------|-------|-------|-------|
| α_i | `α_i` | 0.1 | day⁻¹ | **Immobilization** rate → insulated (gain; carries δ) |
| α_n | `α_n` | 0.15 | day⁻¹ | **Immobilization** rate → non-insulated (gain; carries δ) |
| β_i | `β_i` | 0.0 | day⁻¹ | **Mobilization** rate from insulated (loss; linear in Π) |
| β_n | `β_n` | 0.0 | day⁻¹ | **Mobilization** rate from non-insulated (loss; linear in Π) |
| δ | `delta` | 1.0 | — | Immobilization exponent — Falconer's θ. 1.0 = MATLAB linear case; Falconer's nonlinear runs use 3.0 |
| η | `η_conv` | 0.8 | — | Conversion efficiency (Falconer's γ); (1−η) lost to respiration |
| ζ | `ζ` | 0.2 | — | Insulation **splitting fraction** on the F_n tendency. Dimensionless, NOT a rate. Clamped to ≤ 1 after the Arrhenius factor (`reactions.jl:115`) |
| λ | `λ` | 0.05 | — | Reduced uptake fraction for insulated hyphae (λ ≪ 1) |
| D_Fn0 | `D_Fn0` | 0.01 | mm²/day | Hyphal extension diffusivity at T_ref |
| D_Fm0 | `D_Fm0` | 1.0 | mm²/day | Translocation at T_ref. No tortuosity, but **network-dependent** — scaled by (F_n+F_i)/(F_n+F_i+K_Fm_net) |
| K_Fm_net | `K_Fm_net` | 2e-3 | μg/mm³ | Half-saturation for the D_Fm network factor |
| ε_F | `ε_F` | 1e-4 | μg/mm³ | Π denominator protection |

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
| κ_d_ref | `κ_d_ref` | 0.01 | day⁻¹ | Desorption rate at T_ref. κ_s/κ_d is the hysteresis strength |
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
| k_L | `k_L` | 10.0 | mm³/μg | Langmuir affinity |
| n_LF | `n_LF` | 0.7 | — | Freundlich exponent (< 1 → heterogeneous sites) |
| k_ma | `k_ma` | 0.048 | — | MAOC capacity per unit clay+silt, g-C/g-mineral |
| f_cs | `f_clay_silt` | 0.40 | — | Clay + silt mass fraction |

**There is no `M_max` field.** The capacity is computed in one place —
`maoc_capacity(soil)` in `src/biology/maoc.jl`:

    M_max = k_ma · f_clay_silt · ρ_b        [-]·[-]·[µg/mm³] = [µg-C/mm³]

`k_ma` is **dimensionless** (g-C per g of clay+silt), which is what makes the
identity close. `K_MA_HIGH_ACTIVITY = 0.086` and `K_MA_LOW_ACTIVITY = 0.048`
(`src/constants.jl`) are Georgiou et al. (2022), 86 ± 9 and 48 ± 6 mg-C
g⁻¹-mineral for high- and low-activity minerals, converted mg/g → g/g. The
default is the low-activity value. See §26 erratum 12.

### 4.4 Aggregate Stability

*Rewritten 2026-07-28. The previous entry listed `k_F` in kPa (the code used Pa),
plus two fields — `χ` and `a_p` — that no source file read. See §26 erratum 11.*

| Symbol | Code | Units | Notes |
|--------|------|-------|-------|
| κ_b | `κ_b` | Pa·mm/(μg/mm³) | Specific binding strength per unit binder carbon. **Fitted** |
| w_E | `w_E` | — | Binding weight of EPS relative to insulated hyphae, per unit carbon. **Fitted** |
| d₃₂ | `d_32` | mm | Sauter mean particle diameter. **Measured** — from texture via `sauter_from_texture` |
| p_Gc | `p_Gc` | — | Size dependence of the threshold, `G_c ∝ (r/δ_s)^p_Gc`. **Fitted**; default 0 = no size dependence |

**The criterion.** Local cohesive strength against the wall shear stress of the
oscillatory Stokes layer:

    τ_c(r) = (κ_b / d₃₂) · [F_i(r) + w_E·E(r)]  ≥  τ_w(r)

which reduces to a concentration threshold:

    F_i(r) + w_E·E(r)  ≥  G_c(r) = (τ_w · d₃₂ / κ_b) · (r/δ_s)^p_Gc

With `p_Gc = 0` — the package default — the threshold does not depend on radial
position and this is the pre-2026-07-29 form `G_c = τ_w·d₃₂/κ_b` exactly. See
§4.4a.

`τ_w = √2·μ·v_s/δ_s` with `δ_s = √(2ν_w/Ω)`; for the standard sieving protocol
(`L_s = 13` mm, `f_s = 34` min⁻¹) this gives **δ_s = 0.7510 mm** and
**τ_w = 0.04367 Pa**. That chain is unchanged and independently verified.

#### Why strength scales as 1/d₃₂

Failure separates particles at their bonded contacts. Particles of diameter `d`
tile a failure surface, so the number of bonded contacts crossing unit area
scales as `1/d²`. Strength is that count times the force per bond. Where the
bond force scales linearly with particle size — the geometry of a binder bridge
held at a contact, as for a capillary bridge — strength `∝ (1/d²)·d = 1/d`.

At fixed binder concentration a **finer soil is stronger**: more contacts in the
same failure area. Equivalently a coarser soil needs more binder to hold, so
`G_c ∝ d₃₂`.

This is the form derived and tested against three sand sizes by **Albalasmeh and
Ghezzehei (2014)**, whose model makes aggregate stability depend on soil texture
together with the strength, density and mass fraction of the binding organic
matter. It is a deduction plus that prior result, **not** a fit to the De Gryze
soils — those are the test, and the comparison is reported in
`paper/de_gryze/degryze_EJSS_2006_spec.md`.

**The caveat.** The 1/d result requires the bond force to scale linearly with
particle size. If the bridge instead stays geometrically similar so that force
`∝ d²`, the size dependence cancels entirely. The linear case is the capillary
analogue; EPS is a hydrogel forming bridges at contacts, so the analogy is
reasonable but it is an analogy.

#### Why the Sauter mean

The bond count follows interfacial area, so the correct average is the one that
preserves the volume-to-surface ratio — the Sauter mean `d₃₂ = 1/Σ(fᵢ/dᵢ)`
(**Sauter, 1926**). A geometric or arithmetic mean is dominated by the coarse
fraction, which supplies the fewest contacts per unit mass, and gives a
measurably worse ordering across the five De Gryze soils.

`d₃₂` is dominated by the finest class, so the representative diameter assumed
for clay is the most consequential input (`TEXTURE_CLASS_DIAMETERS`, clay =
1 µm). **Ordering between soils is robust to that choice; the spread is not.**

#### 4.4a The size dependence of the threshold

*Added 2026-07-29. Implemented in `critical_binding`; `p_Gc = 0` is the default,
so every result predating this stands unchanged.*

    G_c(r) = G_c(δ_s) · (r/δ_s)^p_Gc

`p_Gc > 0` makes a larger aggregate need more binder per unit volume — it is
weaker at fixed binder concentration.

**The exponent is fitted.** Two mechanisms act in that direction and neither
fixes a value:

* *Flaw statistics.* A larger body samples more of the flaw-size distribution,
  so its strength falls with size (Weibull). For soil aggregates specifically,
  measured tensile strength falls with aggregate diameter — **Dexter and Watts**
  report roughly `Y ∝ d^-0.5` to `d^-1`, i.e. `p_Gc` between 0.5 and 1.
* *Boundary-layer protrusion.* `δ_s = 0.7510 mm` sits inside the size range the
  sieve series resolves. An aggregate smaller than δ_s lies inside the Stokes
  layer and sees a velocity difference growing with its own radius; a larger one
  protrudes into the free stream. The flat-wall closure that gives `τ_w` does not
  hold across that transition. This is the limitation already recorded under
  `compute_r_agg`; `p_Gc` is an empirical stand-in for it, **not** a derivation
  of the curvature correction.

A third possibility — that the disruptive stress genuinely is size-independent,
which is what a Stokes-drag balance on a sphere in uniform shear gives, since
the force scales as `r²` and so does the failure area — is `p_Gc = 0`, inside
the same family. The parameter spans a structural question the model does not
settle; it is not a correction bolted onto a settled one.

**`p_Gc < 0` is a form already rejected.** Before the 2026-07 rewrite the code
tested `[F_i + w_E·E]·r ≥ G_c`, which is `p_Gc = -1`: larger aggregates
inherently more stable. It was removed because its derivation distributed Stokes
drag over the aggregate surface — the wrong failure model — and because it made
the criterion easier to satisfy the further out it was tested, so the aggregate
had no reason to stop growing (`dev_notes/claude_code_instructions_aggregate_stability.md`).
`p_Gc > 0` is the opposite sign and is the one the aggregate-strength
measurements support. Setting `p_Gc < 0` reinstates the rejected form; the
parameter permits it because the family is continuous, not because it is open.

**Why δ_s is the pivot.** `G_c(δ_s) = τ_w·d₃₂/κ_b` for every `p_Gc`, so changing
the exponent rotates the threshold about a fixed point instead of also moving
its level. δ_s is measured from the sieving protocol, so no new fitted length
enters. `κ_b` sets the level and `p_Gc` the shape, and they are separable for
that reason — but a run that changes both at once cannot attribute anything to
either.

**What motivated it.** Measured 2026-07-29 at `p_Gc = 0` on De Gryze soil 3:
`r_agg` sits at `r_0` through day 2, jumps on day 3, then is flat for 42 days.
The flat threshold cuts the binder profile out in the far field where the
profile is flat, while the binder accumulates 6-fold near the POM — nothing
grades growth. A threshold rising with `r` moves the crossing inward into the
region still accumulating, which is the only mechanism in this model that can
give a radius growing over 21 days. The MATLAB precursor had the same idea as a
commented-out `strength ./ x`.

#### What is fitted and what is not

`κ_b` and `w_E` are **fitted**. The argument above fixes the *form* — strength
proportional to binder concentration and inversely to particle size — and says
nothing about the coefficients.

`κ_b = 0.16` since 2026-07-29 (via 0.012, 0.008, 0.024, 0.072, 0.1, 0.15, 0.18). The previous `0.0143869` was **not measured and
not fitted to data**: it was `2.25 × d₃₂(soil 3)`, chosen to reproduce the legacy
`G_c = τ_w/2.25` exactly, and that `k_F = 2.25` had a dimensionally inconsistent
derivation (§26 erratum 11). It preserved a number rather than a result, and
matching it constrained nothing. 0.16 is a fitted value like any other and
carries no more authority than the evidence behind it.

`w_E = 0.5` was hard-coded in `aggregate_radius.jl` before the 2026-07 revision
and is now a named parameter so it can be fitted rather than hidden. See §5a.

#### Two scalings in the binder when domains overlap

`F_i + w_E·E` sums contributions at two scales whenever `ω > 1`:

| contribution | scaling | why |
|---|---|---|
| pre-existing background binding agents | reduced by `ω` | apportioned across the domains that share that soil (§5d) |
| residue-derived binding agents | not scaled by `ω` | produced from an undiluted surface flux; near the residue, set by release against diffusion rather than by domain size |

Both are correct. Aggregation may be present at `t = 0`, and newly produced EPS
adds to what the soil already held, so the pre-existing fraction must be
distributed across domains rather than counted once per domain.

`G_c` is a single physical threshold, so the criterion is exact only where one
contribution dominates. Once residue decomposition is underway that is the
residue-derived one — `F_i` sits at `F_i_min` and the crossing is set by EPS from
residue carbon. The residual discrepancy is carried by `κ_b`, which is fitted.

**No single factor corrects this.** Multiplying the binder by `ω` would inflate
the dominant term, which is not diluted by `ω`.

`p_Gc` is **fitted** as well, and its default of 0 is a default, not a finding —
it preserves the pre-2026-07-29 model, nothing more. §4.4a says what is known
about its value.

The package default `d_32 = 0.01` mm is a **nominal placeholder**. Every real
problem must override it.

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

## 5a. Parameter Provenance and Anchors


> **Every value in this package without a citation is a working assumption,
> including all "defaults".** They are initial guesses, not reasoned settings,
> and carry no more authority than a value set in a problem folder. Do not
> describe one as reasonable, sensible, physical or calibrated unless a source
> is named here. The intention is that a good De Gryze fit eventually becomes
> the package starting point.

**Purpose.** Records, for every parameter that can move a model output, whether a
literature or precursor value exists and whether the code currently uses it.
Sections 2–4 give values; this section gives their standing.

**What counts as an anchor.** A published source, or a value in the MATLAB
precursor that itself carries a citation. A value inherited from the precursor
*without* a citation is **not** an anchor — it is one prior choice among many.
Derived quantities (`K_Y`, `ε_Y`, `C_B`, `K_E`) inherit the standing of what they
are derived from.

**Source documents.** `dev_notes/MATLAB_parameters.md` and
`dev_notes/MATLAB_aggregation_logic.md` (parse of the 2018 precursor, with its
citations); `dev_notes/falconer_answers.md` (extraction from Falconer 2005/2008).

> ⚠ `dev_notes/` is listed in `.gitignore`. Those three files are **untracked**
> and are the only record of most citations below. This section exists partly so
> the provenance survives independently of them.

---

### Group A — anchored, and the anchor is in use

| Parameter | Value | Anchor |
|---|---|---|
| `K_B` | 1e-4 | Maier 2015 Box 3.1 — range 1e-5–3e-3 |
| `L_B`, `L_F` | 1.29e-3 | Maier 2015 — range 1e-5–3e-3 |
| `α_O` | 2.2 | Maier 2015 Eq 3.19 (2.667) / Ex 3.5 (2.167) |
| `ν_F` | 7.58e-5 | Tresner & Hayes 1971, via Ghezzehei et al. 2018 |
| `θ_r` | 0.06 | UNSODA, average silt loam |
| `D_C0_ref` | 86.4 | Flury & Gimmi, in Dane & Topp (Eds.) 2002 |
| `λ` | 0.05 | MATLAB; Falconer λ₂ = 0.01 (2005), 0.1 (2008) |

### Group B — anchored, and the code deviates

The actionable list. Correcting these requires no new judgement.

| Parameter | Anchor | `parameters.jl` | `run_degryze.jl` | Deviation |
|---|---|---|---|---|
| `D_B_rel` | Wu et al. 2006, AEM 72:4987–4994 → **0.01** | 0.001 | **1e-5** | **1000×** |
| `ζ` | MATLAB **0.01**; Falconer **0.01** | **0.2** | — | **20×** |
| `r_B_max` | Maier 2015 — **0.48–2.16** | 2.0 ✓ | **8.0** | 3.7× above max |
| `μ_B` | Wang et al. 2014, Biogeosciences 11:1817 → **0.012** | 0.012 ✓ | **0.0036** | 3.3× below |
| `ν_B` | Ghezzehei et al. 2018 → **2e-4** | **5.8e-4** | — | 2.9× above |
| `α_i, α_n, β_i, β_n` | Falconer 2005 via MATLAB — **0.6, 0.8, 0.1, 0.2** | 0.1, 0.15, 0, 0 | — | gains 4–5× low; mobilization disabled |
| `Y_B_max` | Maier 2015 Fig 3.8 → 0.5 | 0.7 | **0.4** | no range published |
| `θ_s` | UNSODA silt loam → 0.43 | **0.5** | — | 16% |
| `n_vg` | UNSODA silt loam → 1.39 | **1.47** | — | 6% |
| `delta` | Falconer θ: 1.0 (linear control) or 3.0 (nonlinear runs), `falconer_answers.md` §B2; MATLAB 1 | 1.0 | — | **OPEN DECISION — see below** |

Note that `run_degryze.jl` moves four of these off their anchors simultaneously
(`r_B_max`, `μ_B`, `Y_B_max`, `D_B_rel`), and `R_P_max` was subsequently tuned
against that configuration.

### Open decision: `delta`

**Status: undecided. The default is 1.0 and nothing has settled that it should
stay there.**

Falconer runs θ = 1.0 as a *linear control* and θ = 3.0 as the nonlinear case,
and calls θ > 1 the regime of interest (`dev_notes/falconer_answers.md` §B2).
Our default of 1.0 reproduces the MATLAB reference bit-for-bit — MATLAB is
Falconer with θ = 1 and β = 0, a strict special case — so it is a
porting-verification choice, not a scientific one.

At least three values have been in play historically: **2.0** (documented in
§3.3 until 2026-07-27), **1.0** (the code), **3.0** (Falconer). Whoever chose
1.0 did not record why, and `docs/julia_falconer_deviations.md` — the
document whose whole subject is where Julia departs from Falconer — does not
mention `delta` at all.

Note that δ only bites when β > 0. With mobilization disabled (β_i = β_n = 0,
the current default) the gain is `α·Π^δ` with nothing to compete against it, so
raising δ rescales the transition rate rather than creating the threshold
behaviour Falconer describes. **δ and β should be decided together, not
separately.**

This is recorded here and nowhere else. It was previously encoded as a failing
assertion in `test/test_parameters.jl` — see §26.

### Resolved: `M_max` *(2026-07-29 — see §26 erratum 12)*

Was an open decision with two values in the same run — 288 from
`k_ma·f_clay_silt·ρ_b` at initialisation against `soil.M_max = 10.0` in the
solver. Both were wrong.

There is now **one** definition, `maoc_capacity(soil)`, and the `M_max` field is
gone. `k_ma` is dimensionless and anchored to Georgiou et al. (2022). For the
De Gryze soils at low-activity mineralogy, `M_max` = 30.9, 44.1, 36.8, 38.2 and
48.5 µg-C/mm³ for soils 1–5.

The residual open question is **mineralogy, not units**: whether a given soil is
high- or low-activity. That is a per-soil property with two candidate values, not
a free parameter.

### Group C — no anchor, and outcome-relevant

Fitting belongs here. Roughly fifteen parameters against two observables.

| Parameter | Value | Controls | Precursor |
|---|---|---|---|
| **`κ_b`** | 0.16 | sets the LEVEL of `G_c` (inversely — `G_c = τ_w·d₃₂/κ_b`) | none. Was 0.0143869, reverse-engineered to reproduce the previous `G_c = 0.0194` at soil 3, itself from `k_F = 2.25` whose stated derivation `σ_h·d_p/(4ρ_h)` was dimensionally inconsistent (§26 erratum 11). Changed 2026-07-29 — matching the legacy number constrained nothing |
| **`w_E`** | 0.5 (run: 0.05) | relative value of EPS vs hyphal carbon | none; inherited as a hard-coded literal from the MATLAB precursor (`strength = Fi + 0.5*E`) |
| **`p_Gc`** | 0.0 (run: 1.0) | size dependence of `G_c`; sets whether `r_agg` can grow at all | MATLAB `strength ./ x`, commented out, no record of why. Dexter & Watts put tensile strength at `d^-0.5` to `d^-1`, i.e. 0.5-1 |
| **`γ`** | 0.2 | EPS allocation → the dominant binding agent | MATLAB 0.05, uncited |
| `D_Fn0` | 0.01 (run: 1e-5) | hyphal front speed | MATLAB 1e-11, uncited |
| `D_Fm0` | 1.0 (run: 1e-3) | translocation | MATLAB 1e-5, uncited |
| `μ_E_max`, `K_E` | 0.002, 0.001 | EPS turnover | MATLAB `mu_E` = 0.01, uncited |
| `R_P_max` | 1.0 (run: 1.5) | substrate supply rate | MATLAB split soluble/enzymatic, "not proper partitioning" |
| `κ_s_ref`, `κ_d_ref` | 0.1, 0.01 | MAOC buffering | no MAOC pool in MATLAB |
| `k_L`, `n_LF` | 10.0, 0.7 | MAOC isotherm shape | no MAOC pool in MATLAB |
| `k_ma` | 0.048 | MAOC capacity via `maoc_capacity` | **Group A** — Georgiou et al. (2022), 48 ± 6 mg-C/g-mineral (low-activity); 86 ± 9 for high-activity |
| `B_S`, `F_S` | 0.5, 0.2 | space-limited yield | Julia-era; MATLAB `logistic_B` commented out |
| `ε_F` | 1e-4 | caps Π | Julia-era; MATLAB `PI = Fm/(Fi+Fn)`, unregularized |
| `F_i_min`, `F_n_min`, `F_m_min` | 1e-6, 1e-4, 1e-6 | fungal floors — set the starting point | MATLAB IC: `Fi = 0`, `E = 0` |
| `f_insulated` | 0.5 | initial F_i / F_n split | Falconer E1: **not stated**; MATLAB starts F_i = 0 |
| `ω_E`, `ω_F` | −0.001, −0.0005 | EPS / fungal effect on water retention | MATLAB `WRC_change` = −1.5, **cited** to Ghezzehei & Albalasmeh 2015 — different parameterization; the citation did not survive the port |

### Group D — not parameters, but they move the result

Hard-coded in `src/postprocessing/aggregate_radius.jl`, unreachable from any
configuration. They fix `τ_w = 0.0437 Pa` for every soil and every run:

| Quantity | Value | Note |
|---|---|---|
| `L_s` | 13.0 mm | sieving stroke length — assay property |
| `f_s` | 34/60 Hz | sieving frequency — assay property |
| `μ` | 1.002e-3 Pa·s | water viscosity at 20 °C |
| `ν_w` | 1.004 mm²/s | water kinematic viscosity at 20 °C |

`L_s` and `f_s` describe one apparatus. A different sieving protocol changes
`τ_w` and therefore `G_c`, with no way to express that through the parameter
structs.

Also hard-coded in `paper/de_gryze/run_degryze.jl`: `POM_mean = 1.25`,
`POM_sigma = 0.23` (the paper reports only "0.5–2 mm"; the normal distribution
is an invention), `f_domain_min = 10.0`, `ρ_POM = 200`.

---

### Implication for calibration

Group B should be corrected **before** any fitting exercise, not during it. A fit
performed over Group C while Group B sits off its anchors will absorb the Group B
errors into unanchored parameters, and the deviation becomes invisible.


## 5b. Domain Tessellation and Population Upscaling

`src/physics/tessellation.jl`. The solver integrates one domain around one POM
particle; comparing against a bulk-soil measurement requires the geometry that
links the two. Both quantities below follow from measured inputs, so neither
belongs in an experiment script.

### Domain geometry — `domain_tessellation`

| Quantity | Formula | Meaning |
|---|---|---|
| `φ_POM` | `I_input · ρ_b / ρ_POM` | POM volume fraction of bulk soil [-] |
| `f_pack` | `(1/φ_POM)^(1/3)` | packing-cell diameter / POM diameter [-] |
| `f_domain` | `max(f_pack, f_domain_min)` | model domain radius / POM radius [-] |
| `ω` | `(f_domain/f_pack)³` | overlap correction [-] |

The fullest treatment is `SI_tessellation.tex` (Overleaf) — proofs, the error
bound and an implementation table. **It is a working draft, not a specification**;
where it and the code disagree, neither is presumed right. §26 erratum 13 lists
the disagreements.

`f_domain_min` (default 10) exists so the radial grid resolves the near-POM
gradients. Where it binds, neighbouring domains overlap and each holds more than
its share of surrounding soil — `ω` is that excess.

### The ω convention

`create_initial_state` divides **background carbon** by `ω`: C, B, F_n, F_m,
F_i, E, M. **POM is not divided**, being a lumped scalar at the domain centre
rather than carbon spread through the soil.

The same asymmetry holds downstream. `Σ nᵢ·CO₂ᵢ` and `Σ nᵢ·P₀ᵢ` are already
per-unit-soil totals and must not be divided by `ω` again; aggregate volume
fractions are built from physical diameters and physical counts and carry no
`ω` either. Dividing twice understates the amendment by a factor of `ω`.

### Population — `pom_population`

    V_pack[i] = (4/3)π(dᵢ·f_pack/2)³        soil owned by one particle
    N_POM[i]  = f_POM[i] · V_soil / V_pack[i]
    P_0[i]    = (4/3)π(dᵢ/2)³ · ρ_POM

**The cells tile the soil:** `Σ Nᵢ·V_pack,ᵢ = V_soil` exactly, provided
`Σ f_POM = 1`. This identity is why the population bookkeeping in §5c can
subtract totals rather than looping per class — the two are the same arithmetic.
It also means `f_POM` **must be normalised**: a distribution truncated at the bin
edges recovers an `I_input` low by exactly the shortfall.

**Verification:** `total_POM_C / (V_soil · ρ_b · 1e-6)` must reproduce the
experiment's stated amendment rate in µg-C per g soil.

This is an **exact, distribution-free identity**, not a sanity check:

    Σ Nᵢ·P₀,ᵢ = Σ [fᵢ·V / ((4/3)π(dᵢ f_pack/2)³)] · (4/3)π(dᵢ/2)³·ρ_POM
              = (V·ρ_POM / f_pack³) · Σ fᵢ
              = V·ρ_POM·φ_POM                    (f_pack³ = 1/φ_POM)
              = V·I_input·ρ_b

so it holds for any diameters and any normalised fractions, and it pins the
`f_pack` exponent, both radius/diameter halvings and both `(4/3)π` factors at
once. Asserted in `test/test_tessellation.jl` (added 2026-07-28; this module
previously had no tests at all).

### `log_interpolate_fraction(lo, mid, hi)`

`ln(mid/lo) / ln(hi/lo)` — the fraction of a size class below `mid` assuming
equal mass per logarithmic interval. Used where a measured size fraction spans
two reporting classes and the split is unreported. It is an assumption about
particle-size distributions and must be recorded as such wherever applied.

### `van_genuchten_inverse(θ, α, n, θ_r, θ_s)`

`src/physics/water_retention.jl`. Returns the ψ [kPa] giving water content θ:

    ψ = -(S_e^(-1/m) - 1)^(1/n) / α,    S_e = (θ-θ_r)/(θ_s-θ_r),  m = 1-1/n

For experiments specifying water content or WFPS rather than a potential. The
result depends on α, n, θ_r and θ_s, so it must be evaluated with the same soil
the run will use.


## 5c. Aggregate Mass, Sieve Classes and Sample Statistics

`src/postprocessing/population.jl`. §5b builds the population; this section
turns it into something a wet-sieving measurement can be compared with.

### `ρ_b` is the amended sample's bulk density

Throughout §5c, `ρ_b` is the bulk density of the **amended** sample — soil plus
the residue mixed into it — so that `ρ_b · V` is the total sample mass. This is
a deliberate simplification: it removes any need to track the amendment
separately when forming `f_agg` and the class shares. The amendment is ~1 % of
soil mass, so the amended and bare-soil readings differ by ~1 %.

The same `ρ_b` also sets the aggregate shell density (`aggregate_mass`), the POM
volume fraction in the tessellation (§5b), and porosity via `θ_s = 1 − ρ_b/ρ_s`.
The convention therefore propagates, at that ~1 % level, into all three. Where a
measured bulk density is available for the amended soil it should be used
directly; where only the bare-soil value is reported, using it here is the
approximation being accepted.

### Mass of one aggregate — `aggregate_mass`

    m = ρ_b·(4/3)π·max((D_agg/2)³ − r_0³, 0)   mineral shell
      + P_C / f_C_POM                           remaining residue

Mass, not volume. `r_agg` can sit still for days while the residue inside is
respired away, so volume weighting discards the only time dependence retained
mass has, and simultaneously treats a POM-rich aggregate as if it were solid
mineral. `f_C_POM` is the residue carbon fraction and is experiment-specific —
the function has **no default** for it or for `ρ_b`.

### Sieve classes — `sieve_class`

`k` apertures define `k+1` classes, ascending. The convention is `D < aperture`
passes, matching wet sieving. Class 1 is below the first aperture; class `k+1`
is at or above the last.

### Sample statistics — `population_statistics`

Takes `n_sizes × n_times` matrices of `D_agg`, `POM`, `CO₂` and `CO₂ flux`,
plus `r_0` and the particle counts `n_dist` from §5b.

| Output | Definition |
|---|---|
| `MWD_agg_only` | `Σᵢ nᵢmᵢDᵢ / Σᵢ nᵢmᵢ` — over aggregates **only** |
| `f_agg` | aggregate share of sample mass `ρ_b·V_soil` |
| `f_agg_vol` | aggregate share of sample volume — diagnostic |
| `POM_mass_frac` | residue share of retained coarse mass |
| `CO2_total`, `CO2_flux_total` | `Σᵢ nᵢ·CO₂ᵢ`, `Σᵢ nᵢ·flux ᵢ` |
| `class_pct` | `(k+1) × n_times`, % of sample mass per class |
| `MWD_fixed_weight` | `Σₖ nominalₖ·class_pctₖ / 100` |

### The well-mixed closure, and where it stops being defensible

The classes are defined by **the assay's sieve series**, not by sand/silt/clay.
The algorithm is universal; only the cutoffs are problem-specific. Mapping a
measured texture onto a particular sieve series belongs to the experiment.

Given no information about which particles a growing aggregate takes up, the
least-assuming closure is that the matrix is well mixed and material is
incorporated **in proportion to its relative abundance**. The remainder's
composition is then unchanged and only its total falls, which is what
`cls[k] += (1 − f_agg)·mineral_class_fractions[k]` encodes. This is a statement
about what the model does **not** track, not a claim about uptake physics.

**The continuum/discrete boundary.** The model is a continuum — a homogeneous
shell of bulk-density matrix around a POM core. The class accounting is discrete
— it speaks of particles of named sizes. The two are compatible only while the
shell is much thicker than the particles it is said to contain. A 0.6 mm
aggregate cannot engulf a 0.5 mm grain, and no bookkeeping rule makes it able to.

At the De Gryze operating point this is not a marginal concern. Day 21, soil 3:

| Quantity | Value |
|---|---|
| mass-weighted shell thickness | **485–711 µm** (`shell_thickness_mm`) |
| coarse mineral class | 250–2000 µm |
| share of that class coarser than the shell | ~50–68 % |
| **share of the whole sample coarser than its own shell** | **~14 %** |
| `f_agg` from day 5 | **0.25** |

So roughly a seventh of the sample is mineral that the closure assigns to
aggregate interiors it could not physically fit inside, while a quarter of the
sample is claimed as aggregate. The **total** mass balance is unaffected — the
columns still sum to 100 — but the *composition* of the remainder is not
defensible in detail.

`shell_thickness_mm` is returned from `population_statistics` for exactly this
check: compare it against the coarse end of `sieve_sizes` in every run.

**No size-restricted uptake rule is applied, deliberately.** Restricting
absorption to particles finer than the aggregate would splice a discrete
particle model into a continuum one. That is the incoherence restated, not a fix
for it. Resolving it properly means choosing one representation.

### Per-class occupancy — what the bulk formula hides

Subtracting total shell mass from total sample mass is exactly equivalent to
subtracting per particle and summing over classes, by the tiling identity in
§5b. The equivalence is convenient, but it means **one size class overrunning
its own share of soil is silently absorbed by the others** — the total stays
consistent while a class has eaten more than it owns.

`population_statistics` therefore accepts `cell_volume_mm3` (from
`pom_population`'s `V_pack`) and returns `cell_occupancy[i,t]`, the aggregate
volume of class `i` divided by the cell it owns. **A value ≥ 1 is a failure of
the independent-domain picture**, not a large number. `NaN` when the cell
volumes are not supplied.

De Gryze soil 3 at day 21, for scale: 0.47 for the 0.65 mm class down to 0.17
for the 1.85 mm class. The smallest POM class is the one to watch, because it
owns the least soil.

### `f_agg` is not clamped

`f_agg > 1` means the model produced more aggregate than there is sample. It is
reported raw and warns, because a clamp would leave `f_agg` reading exactly 1.0
while the class columns summed past 100. Only the mineral top-up is floored, so
the class shares stay non-negative and the disagreement stays visible.

### What the model does not predict

The model predicts aggregates. It carries **no particle-size distribution for
the unaggregated matrix**. Reporting class shares as fractions of the whole
sample therefore requires that distribution to be supplied
(`mineral_class_fractions`, `k+1` fractions summing to 1); the unaggregated
remainder `(1 − f_agg)` is then spread across the classes with it and the
columns sum to 100 %. Supply nothing and the class columns carry aggregate mass
alone, summing to `f_agg`.

`class_nominal_mm` without `mineral_class_fractions` is rejected: a
fixed-weight sum over shares that are not whole-sample shares is not comparable
with a published MWD.

`MWD_agg_only` and a published MWD are **different quantities**. The first is a
mean over POM-centred aggregates; the second is a fixed-weight sum over sieve
classes covering the whole sample. A fixed-weight MWD is a convex combination
of its nominals, so it saturates at the largest and is bounded below by the
smallest whatever the aggregates do.

### ω does not appear here

Every sum in this file is built from physical particle counts and physical
diameters and is already a per-sample total. `ω` corrects background-carbon
*concentration* inside an oversized domain and is applied once, at
initialization (§5b). Applying it again understates the amendment by that
factor.

### Placement

`paper/de_gryze/postprocess_dataframe.jl` holds a `population_outputs` wrapper
that extracts columns from a sweep DataFrame, calls `population_statistics`,
and names the class columns. It defines nothing. It sits outside `src/` only
because the package does not depend on DataFrames or CSV.

## 5d. Initial condition — the SOC partition

**Undocumented until 2026-07-29** — it sets two of nine state variables and
appeared in neither the manuscript nor this file. Now in both: `manuscript-4-5.tex`
§`sec:initial_condition`.

Background SOC is apportioned at `t = 0` by `create_initial_state`
(`src/physics/initial_conditions.jl`):

```
SOC_vol       = SOC · ρ_b                                  measured × measured
biotic pools  = f_bact, f_fungi, f_eps  × SOC_vol           three fractions
SOC_residual  = SOC_vol − (biotic pools)
```

`SOC_residual` then splits between DOC and MAOC by `partition_CM`, which solves
mass balance against the sorption isotherm simultaneously:

```
C + M = SOC_residual                    and        M = M_max · f_LF(β·C)
β = k_L·k_d / (θ + ρ_b·k_d)                        M_max = maoc_capacity(soil)
```

One equation, one unknown after substitution; solved by fixed-point iteration in
the DOC-limited regime and bisection in the mineral-limited one.

**`J_M` at the partition point is not exactly zero.** The softplus regularisation
has `φ_ε(0) = ε·ln2`, so a smoothed switch carries a permanent sorption bias

    J_M(M_eq, M_eq) = (κ_s − κ_d)·ε_maoc·ln2

= 6.24e-4 µg mm⁻³ day⁻¹ at the package rates. Over 21 days that moves `M` by
0.05 %, against 184× that much if the state is placed at 60 % of `M_eq` — so the
floor is the right target and the distinction matters only when reading a
diagnostic. `test_api.jl` asserts its value rather than asserting zero.

### Saturation is an output

`M/M_max` is not specified. It is whatever sorption equilibrium with the
available carbon produces, and it is bounded:

> as `k_L → ∞`, `C → 0` and `M/M_max → SOC_residual/M_max`

— a ceiling fixed by **measured SOC and measured texture alone**. So the model
predicts each soil's saturation deficit and the prediction is testable against
reported values (Georgiou et al. 2022). For the De Gryze soils at low-activity
mineralogy the ceilings are 0.79, 0.53, 0.80, 0.87, 0.82 for soils 1–5 — a
spread of 0.34 across five soils from one meadow, generated by texture.

### Sorption equilibrium is not SOC steady state

`M_eq(C)` is fast, local and physicochemical — `κ_s ≈ 0.1 day⁻¹` against a
21-day incubation, so mineral surfaces equilibrate with the pore water almost
immediately. Starting there is not a claim that inputs balance outputs. It is
only a claim that the fast process has run.

### What the partition cannot do

It has two pools to hold background carbon and no third. Any SOC the isotherm
will not take is declared dissolved, because nothing else can hold it. Soils
whose native POM has *not* been removed therefore need either a particulate
background pool or an explicitly unrepresented fraction — see §26 erratum 14.
`k_L` is the parameter that decides how much the minerals take, and the only
constraint on it is that pore-water DOC land in the observed 10–100 mg/L.

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

Net **signed** tendencies. α (immobilization, gain) carries the exponent;
β (mobilization, loss) is linear. Positive = net immobilization:
```
net_i = (α_i · Π^δ − β_i · Π) · F_i
net_n = (α_n · Π^δ − β_n · Π) · F_n
```

Conversion efficiency applied:
```
trans_i = η · net_i
trans_n = η · net_n
```

ζ is a **splitting fraction** on the F_n tendency, not an independent drain. A
fraction ζ of trans_n is redirected to F_i; the remainder stays with F_n. Because
both are signed, a mobilizing F_n also draws a ζ share off F_i:
```
immobil_i = trans_i + ζ · trans_n
immobil_n = (1 − ζ) · trans_n
```

Conversion respiration — `abs()` is required so mobilization also carries a cost
and carbon still closes:
```
Resp_F_conv = (1 − η) · |net_i + net_n|
```

Source terms (`reactions.jl:165–167`):
```
S_Fi = immobil_i − R_rec_F
S_Fn = immobil_n
S_Fm = Γ_F − immobil_i − immobil_n − Resp_F_conv
```

Mobilization is currently disabled (β_i = β_n = 0), so `S_Fn ≥ 0` and F_n has no
sink of any kind — Falconer gives b_n no death term either (falconer_answers.md
§C1). With β > 0 the bracket becomes signed and F_n can lose biomass.

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
$$J_P = R_P^{\max}(T) \cdot \frac{P}{P_0} \cdot \tfrac{1}{2}\left(\frac{B_0}{K_{B,P}+B_0} + \frac{F_{n,0}}{K_{F,P}+F_{n,0}}\right) \cdot \frac{\theta_0}{\theta_P+\theta_0} \cdot \frac{O_{aq,0}}{L_P+O_{aq,0}}$$

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

Transcribed from `src/solver/reactions.jl:160–172` and
`src/biology/fungi.jl:324–350`, verified 2026-07-28.

The fungal terms are written through the intermediates the code actually
computes, not expanded. Expanding them is what let this section drift out of
step with §9 for five months — see §26.

```
Fungal intermediates (fungi.jl, fungal_transitions):

  net_i     = (α_i·Π^δ − β_i·Π)·F_i         signed: + immobilizing, − mobilizing
  net_n     = (α_n·Π^δ − β_n·Π)·F_n
  immobil_i = η·net_i + ζ·η·net_n            ζ splits the F_n tendency
  immobil_n = (1 − ζ)·η·net_n
  Resp_F^conv = (1 − η)·|net_i + net_n|      abs(): mobilization also costs C


S_P  = −R_P

S_C  = −R_B − R_F + R_rec − J_M                      (+ diffusion + POM Neumann BC)

S_B  = Γ_B − R_rec,B

S_Fi = immobil_i − R_rec,F

S_Fn = immobil_n                                     (+ diffusion)

S_Fm = Γ_F − immobil_i − immobil_n − Resp_F^conv     (+ diffusion)

S_E  = Γ_E − R_rec,E

S_M  = J_M

S_O  = −α_O·(Resp_B + Resp_F + Resp_F^conv) / (θ + K_H(T)·θ_a)   (+ diffusion)
```

**There is no separate insulation term.** `S_Fn` has no `−ζ·F_n` and `S_Fi` has
no `+ζ·F_n`. ζ is a dimensionless splitting fraction applied *inside*
`fungal_transitions` to the F_n tendency, not a one-way drain rate. The
independent-drain form was removed in February 2026; if you find `trans_ni`
anywhere, it is stale.

**F_n has no death term.** Only F_i recycles (`R_rec,F` uses `F_i` alone). That
matches Falconer, who gives neither sessile pool an independent loss
(`dev_notes/falconer_answers.md` §C1) — though the F_i death term itself is
ours, not Falconer's.

**The O₂ divisor is not optional.** `S_O` is divided by the O₂ capacity
`θ + K_H(T)·θ_a`, because O is stored as a bulk concentration while respiration
consumes from the aqueous phase (`reactions.jl:170`).

**Conservation**:

    S_P + ∫(S_C + S_B + S_Fi + S_Fn + S_Fm + S_E + S_M)·4πr²dr
        = −∫(Resp_B + Resp_F + Resp_F^conv)·4πr²dr

A useful sub-check: `immobil_i + immobil_n = η·(net_i + net_n)`, so ζ cancels
out of `S_Fi + S_Fn + S_Fm = Γ_F − Resp_F^conv − R_rec,F` exactly. **If a change
to ζ alters total fungal carbon, ζ has stopped being a splitting fraction** —
which is the failure mode the Arrhenius clamp at `reactions.jl:116` guards
against.

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

### physics/tessellation.jl
| Function | Computes |
|----------|----------|
| `domain_tessellation(; ρ_POM, I_input, ρ_b, f_domain_min)` | φ_POM, f_pack, f_domain, ω |
| `pom_population(diam, f_POM, tess; ρ_POM, soil_volume_mm3)` | N_POM, P_0 per particle, total POM C |
| `log_interpolate_fraction(lo, mid, hi)` | equal-mass-per-log-interval split |

### postprocessing/population.jl
| Function | Computes |
|----------|----------|
| `sieve_class(D, sieve_sizes)` | class index, ascending, `D < aperture` passes |
| `aggregate_mass(D_agg, r_0, P_C; ρ_b, f_C_POM)` | shell + core dry mass [µg] |
| `population_statistics(D_agg, POM, CO2, flux, r_0, n_dist; ...)` | MWD, f_agg, class shares, population CO₂ |

### physics/particle_size.jl
| Function | Computes |
|----------|----------|
| `sauter_mean_diameter(f, d)` | `d₃₂ = 1/Σ(fᵢ/dᵢ)` — preserves volume/surface ratio |
| `sauter_from_texture(sand, silt, clay; ...)` | `d₃₂` [mm] from a USDA texture triple |
| `TEXTURE_CLASS_DIAMETERS` | assumed class midpoints (clay 1 µm, silt 10 µm, sand 316 µm) |

### physics/water_retention.jl
| Function | Returns |
|----------|---------|
| `θ(ψ, E, F_i, soil)` | Water content (van Genuchten, modified α) |
| `van_genuchten_inverse(θ, α, n, θ_r, θ_s)` | ψ [kPa] from water content — inverse of the above |

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
| `fungal_transitions(F_i, F_n, F_m, Π, ...)` | Returns **(immobil_i, immobil_n, Resp_F_conv)** |
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

## 17a. Carbon Closure, and What Actually Verifies It

### The identity

The seven pool source terms must sum to minus total respiration, at every node,
for any state:

$$
S_C + S_B + S_{F_n} + S_{F_m} + S_{F_i} + S_E + S_M = -\,\text{Resp}_\text{total}
$$

This is the model's carbon closure. It is a property of how the rate laws
compose in `compute_source_terms`, **not** something any solver enforces. If it
were broken — by a changed yield expression, a dropped recycling term, a sign
error in a transition — every carbon number the model produces would be wrong,
and nothing would raise an error.

### Proof

Summing the seven terms as written in `reactions.jl`:

- `J_M` cancels between `S_C` and `S_M`.
- `immobil_i` cancels between `S_{F_m}` and `S_{F_i}`; `immobil_n` between
  `S_{F_m}` and `S_{F_n}`. Transitions move carbon between fungal pools and
  create none.
- `R_rec_total = R_rec_B + R_rec_F + R_rec_E` in `S_C` cancels the three
  recycling losses.

leaving $-R_B - R_F + \Gamma_B + \Gamma_F + \Gamma_E - \text{Resp}_{F,\text{conv}}$,
so with $\text{Resp}_\text{total} = \text{Resp}_B + \text{Resp}_F + \text{Resp}_{F,\text{conv}}$
the requirement splits in two.

**Fungal**, exact by inspection:

$$\Gamma_F - R_F = Y_F R_F - R_F = -(1-Y_F)R_F = -\text{Resp}_F$$

**Bacterial.** With $x = R_B - R_{Bb}$ and $\sigma(\cdot) = \text{softplus}(\cdot,\varepsilon_Y)$,

$$\Gamma_B + \Gamma_E = Y_B\sigma(x)(1-\gamma) - \sigma(-x) + Y_B\sigma(x)\gamma = Y_B\sigma(x) - \sigma(-x)$$

and $\text{Resp}_B = R_{Bb} + \sigma(x)(1-Y_B)$. Requiring
$\Gamma_B + \Gamma_E - R_B = -\text{Resp}_B$ reduces, after $Y_B\sigma(x)$
cancels from both sides, to

$$\sigma(x) - \sigma(-x) = x$$

which holds **exactly, for any ε**:

$$\varepsilon\ln\frac{1+e^{x/\varepsilon}}{1+e^{-x/\varepsilon}} = \varepsilon\ln e^{x/\varepsilon} = x$$

The softplus regularisation (§16) therefore preserves the closure identically.
That is not luck — it is the property $\max(0,x) - \max(0,-x) = x$ surviving the
smoothing, and it is the reason ε can be chosen freely without breaking the
carbon budget. **Any future replacement for softplus must satisfy the same
antisymmetry relation.**

The closure is exact analytically. In floating point it holds to roundoff.

### What each solver's carbon balance actually tells you

`compute_carbon_balance_error` reports
$(P + \sum_i W_i \text{pools}_i + \text{CO}_2 - C_0)/C_0$. What that number means
depends entirely on where CO₂ came from, and the two solvers differ.

**Split solver — a real, if diluted, test.** It advances the pools by
$\Delta t\,S_k$ and accumulates CO₂ separately as
$\sum_i \Delta t\,\text{Resp}_i W_i$. These are two independent computations.
The per-step change in the total is
$\Delta t \sum_i W_i(\sum_k S_{k,i} + \text{Resp}_i)$, which vanishes **only if
the closure identity holds**. Its measured $-5.1\times10^{-13}$ over 391,773
steps is therefore evidence that the identity holds along the whole trajectory —
roundoff on a quantity that would otherwise drift.

Two parts of that number are vacuous and should not be credited: the spherical
stencil telescopes (§18) so transport conserves by construction, and the POM→DOC
transfer cancels by construction (§20a). A third part is worse than vacuous —
the clipping correction (§17) adds clipped carbon back to CO₂, which restores
the balance whether or not the clipping was physically right, and so masks
closure violations wherever clipping is active.

**Stiff solver — no test at all.** It recovers CO₂ as $C_0$ minus current total
carbon (§20a). The balance is then identically zero by definition, which is why
`mass_balance_error` is reported as `NaN` rather than `0.0`. Removing CO₂ from
the state removed the independent computation that made the comparison mean
something.

### The three available tests, ranked

1. **Pointwise identity.** Evaluate `compute_source_terms` at a node and assert
   $|\sum_k S_k + \text{Resp}_\text{total}| \le 10^{-12}\times\text{scale}$.
   Exact, no integration, no solver, microseconds, and it localises a failure to
   a node and a state. This is the primary test and **it does not yet exist** —
   none of the current assertions check it.

2. **Split-solver balance error.** An integrated backstop. Weaker (diluted over
   the trajectory, and blinded by clipping) but it samples the states the
   trajectory actually visits, which a hand-chosen unit test may not. Keep it.

3. **Comparing the two solvers' CO₂.** Confounded: the two run different
   trajectories with different integrators, so a disagreement cannot be
   attributed to the closure rather than to the integration. It is a useful
   cross-check between two independent implementations, but it is not a closure
   test, and item 1 subsumes what it would tell you.

### Why the split solver is kept

Not for accuracy, not for speed, and not because it reproduced the MATLAB
precursor. It is kept because it computes CO₂ independently and the stiff solver
structurally cannot, so it carries test 2 — and because two independent
implementations over shared physics are the only convergence evidence this model
has. Once test 1 exists, test 2 becomes a redundancy rather than the only line
of defence.

### Implementation note

`R_diff = R_B - R_Bb` is computed three times per node — once in
`compute_source_terms`, and again inside `Gamma_B` and `Gamma_E`. Same inputs, so
no drift is possible today, but it is the same expression in three places
(CLAUDE.md §8).

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

**Implementation note.** `conservation_weight(r, h)` in `src/types.jl` is the single definition of $W = 4\pi r^2 h$; `conservation_weights(r_grid, h)` is its vector form. `GridInfo` calls it once and stores the result as `grid.W`, which is what every caller with a grid in hand should use. `compute_total_carbon` takes either the weight vector, a `GridInfo`, or `(r_grid, h)` and runs one loop. A separate `compute_cell_volumes` function exists for geometrically exact shell volumes $(4\pi/3)(r_+^3 - r_-^3)$ (used only for visualization and post-processing, never for conservation accounting) — it is a different quantity and is not interchangeable with $W$.

The system total appears at two levels because it has two inputs: `compute_total_carbon(state, W)` sums a raw `AggregateState` against the weights, and `total_system_carbon(pools, k)` in `src/postprocessing/integration.jl` sums the already-integrated `IntegratedPools`. Each is a single implementation, `carbon_balance_table` and `result_to_dataframe` both go through the second, and `test_postprocessing.jl` asserts the two agree.

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

**Note on conservation.** The adaptive stepper does **not** reject and re-do steps. When $\rho$ exceeds the threshold, the current step is accepted and the next step uses the reduced $\Delta t$. This is a "predict then correct" strategy: the occasional large-$\rho$ step introduces $O(\Delta t)$ trajectory error but does not break conservation (the per-step conservation identity holds for any $\Delta t$). 

**Correction (2026-07-28).** This section previously claimed that "in practice, the stepper maintains $\rho \approx 0.01\text{–}0.10$, corresponding to 1–10% relative change per step." That was measured and is false for most of a run. For De Gryze soil 3, $\rho$ at $\Delta t = 10^{-4}$ d is 9e-4 at day 0.25 and 0.150 by day 21, rising to 0.206 by day 45. From roughly day 15 onward the criterion demands $\Delta t < \Delta t_{\min}$ — 6.7e-5 d at day 21, 4.9e-5 d at day 45 — so the solver takes steps it has itself judged too large. It cannot report this: `n_rejected` is guarded by `dt > dt_min`, which is false at the floor, and a 45-day run reports `n_rejected = 9` out of 391,773 steps. The limiter is DOC at the nodes adjacent to the POM surface, where $|S_C|/C$ reaches 2055 /day. See §20a.

**Interaction with clipping.** Smaller $\Delta t$ means less overshoot, fewer clipping events, and better trajectory accuracy. This is why the adaptive stepper produced 50$\times$ better conservation (0.014%/year) than fixed $\Delta t = 0.001$ (0.7%/year) before the clipping bug was fixed — the adaptive stepper was taking smaller steps in regions of stiff dynamics, reducing the number of clipping corrections needed.

---

## 20a. Method of Lines and the Stiff Solver

`run_aggregate_stiff` (`src/solver/mol.jl`, `src/solver/mol_solve.jl`) integrates
the same model as `run_aggregate` with an implicit stiff solver instead of Strang
splitting. Both paths call the same `compute_source_terms` and discretise space
identically; there is one implementation of the biology.

### Why a second solver exists

The criterion in §20 is not an accuracy estimate. It bounds the relative change
the **explicit** reaction step is allowed to make, which makes it a stability
guard: the step is set by the fastest pool anywhere in the domain, regardless of
whether that pool matters to the answer.

Measured (De Gryze soil 3, $D$ = 1.25 mm, $n$ = 200):

| $t$ [d] | limiter | $\lvert S\rvert/u$ [1/d] | node | $\Delta t$ permitted [d] |
|---|---|---|---|---|
| 0 | $F_m$ | 222 | 1 | 4.5e-4 |
| 1 | $O$ | 15.3 | 1 | 6.5e-3 |
| 5 | $C$ | 333 | 15 | 3.0e-4 |
| 21 | $C$ | 1501 | 2 | 6.7e-5 |
| 45 | $C$ | 2055 | 1 | 4.9e-5 |

$S_C$ is roughly flat over the last stretch (−6.63 → −6.39 μg-C/mm³/d) while the
pool it is divided by thins. The ratio is still climbing at day 45. **Cost per
simulated day therefore grows without bound**, which no constant-factor
optimisation addresses and which rules the scheme out for multi-year runs.

An implicit method sizes its step from accuracy rather than stability, so the
step grows once the fast pools reach quasi-steady state and cost tracks activity
rather than elapsed time.

### Formulation

State vector, node-major with species fastest, so the Jacobian is
block-tridiagonal with 8×8 blocks:

$$
u_{8(i-1)+k}, \quad k \in \{C, B, F_n, F_m, O, F_i, E, M\}, \quad i = 1 \ldots n
$$

with $P$ at index $8n+1$; $8n+1$ states in total. Species-major ordering would
give five tridiagonal blocks glued by dense reaction coupling, which fills in
badly under LU.

The spatial operator reproduces `crank_nicolson.jl` term for term — same
node-centred spherical stencil, same **arithmetic** face average of $D$, same
ghost node at the flux boundary. A finite-volume flux form with harmonic face
averages is arguably better and is deliberately **not** used, because then a
disagreement between the two solvers could not be attributed to time
integration. Change the space discretisation separately, and test it separately.

POM coupling is exact rather than merely consistent: $J_P$ is computed once and
enters $\mathrm{d}P/\mathrm{d}t = -4\pi r_0^2 J_P$ and the node-1 boundary term
with opposite signs. With $W_i = 4\pi r_i^2 h$ the node-1 gain is
$W_1 \cdot (1/r_1^2h^2) \cdot r_0^2 h J_P = 4\pi r_0^2 J_P$, identical.

Oxygen at the outer node is pinned ($\mathrm{d}u/\mathrm{d}t = 0$, initial value
$O_{amb}$). This is the same net condition as the split solver, which overwrites
that node on every diffusion half-step — not an approximation of it. Its
respiration still counts toward CO₂, matching `reaction_step.jl`.

### Three differences that are not cosmetic

1. **Coefficient lagging.** $\theta$ and the effective diffusivities are
   evaluated at the current state here. `timestepper.jl` computes them once per
   step from the state at the *start* of the step and holds them fixed across
   both diffusion halves and the reaction, so the split solver lags its own
   coefficients by one step. The two agree as $\Delta t \to 0$.
2. **No negativity clipping.** §17 clips each pool at zero and credits the
   clipped carbon to CO₂. That exists to catch Forward Euler overshoot. Any CO₂
   the split solver reports that originated in clipping is a numerical artefact,
   and the difference between the two runs measures it.
3. **Cumulative respired carbon is not integrated.** It is recovered at output
   times from the carbon balance — initial total carbon minus current total
   carbon. It is a whole-domain quantity wanted a few dozen times per run, so
   carrying it through the solver as a per-node quadrature builds machinery for
   a resolution nothing uses. It also has a dense Jacobian row, and a dense row
   forces the column colouring that sparse forward-mode AD uses to fall back to
   one RHS evaluation per column — measured at 91 % of a 45-day run's cost when
   it was carried as a state.

   The consequence is that `mass_balance_error` is `NaN` on this path. The
   recovery makes the balance identically zero, so reporting a number would
   assert a check that was never performed. What remains available at output
   resolution is monotonicity — respiration cannot be negative, so the recovered
   quantity cannot fall — and that is reported as
   `diagnostics["co2_monotonic"]`.

Because of these, agreement between the solvers should be judged against the
split solver run at a $\Delta t_{\min}$ small enough to be converged, not against
its production settings.

### Solver configuration

Default is `FBDF` with a detected sparse Jacobian and a KLU factorisation. BDF is
the documented choice for systems of this size, but is also documented to
tolerate less stiffness than the ESDIRK family; `KenCarp47` is the fallback if
Newton iterations fail.

Sparsity is detected on a state with every pool forced strictly positive. The
rate laws are guarded with $\max(0, u)$; on a state where a pool sits at zero the
guard short-circuits, the tracer never records the coupling, and the pattern
comes out missing entries — which degrades the Newton solve silently.

`abstol` defaults to a **vector**, $\max(10^{-10}, 10^{-8}|u_{0,i}|)$ per state.
The pools span roughly eight orders of magnitude, so a scalar tolerance either
over-resolves MAOC or under-resolves DOC. Note the weakness: basing it on the
INITIAL value penalises any pool that starts near zero and grows. `F_m` starts
at $3.3\times10^{-8}$ (`F_m_min`/ω) and reaches $3.5\times10^{-4}$, so it is
controlled four orders tighter than its own scale. Not yet resolved.

### Known inefficiency

`mol_rhs!` allocates seven length-$n$ vectors per evaluation ($\theta$,
$\theta_a$, and five diffusivities) so that forward-mode AD sees them at the
right element type. `PreallocationTools.DiffCache` is the standard remedy and is
not yet applied. Recorded here rather than left silent.

### Dependencies

`OrdinaryDiffEqBDF`, `OrdinaryDiffEqSDIRK`, `LinearSolve`, `SciMLBase`,
`SparseConnectivityTracer`, `ADTypes`, `SparseArrays`. These are hard
dependencies of the package rather than of a script, because the integrator is
core model machinery — placing it in `paper/` would require every driver to
carry its own copy (CLAUDE.md §4). If per-session load time proves costly, the
principled remedy is a package extension (weak dependency), not relocation.

`SourceTerms` was made parametric in its element type at the same time. A
hard-coded `Float64` there would silently force the stiff solver onto finite
differences.

---

## 21. Sigmoid Threshold Functions

**Problem.** Biological rates (mortality, maintenance, insulation) should vanish smoothly as biomass approaches a minimum viable threshold, not discontinuously.

**Solution.** Three sigmoid functions enforce soft lower bounds:

$$
h(x) = \frac{\exp(\beta\, x)}{\exp(\beta\, x) + \exp(\beta\, x_{\min})} = \frac{1}{1 + \exp[-\beta(x - x_{\min})]}
$$

with $\beta = 50 / x_{\min}$, giving a smooth transition from 0 to 1 over a width of approximately $4 x_{\min} / 50$. Applied to:

| Function | Variable | Threshold | Effect |
|----------|----------|-----------|--------|
| $h_B(B)$ | Bacteria | $B_{\min}$ | Shuts off maintenance $R_{B,b}$ and mortality near extinction |
| $h_{F_i}(F_i)$ | Insulated fungi | $F_{i,\min}$ | Shuts off fungal mortality near extinction |
| $h_E(E)$ | EPS | $E_{\min}$ | Shuts off EPS recycling near depletion |

These prevent numerical extinction artifacts where a pool oscillates around zero due to competing production and consumption terms. The sigmoid ensures a smooth, monotonic approach to zero rather than oscillatory overshoot.

**Implementation.** The two forms above are algebraically identical; the code uses the right-hand one, because the left-hand one overflows for $x \gg x_{\min}$. `sigmoid_threshold(x, x_min, steepness)` in `src/math_utils.jl` is the single definition, with `steepness` $= \beta x_{\min}$ so it is scale-free in $x$. `h_B`, `h_E` and `h_Fi` call it with `SIGMOID_STEEPNESS = 50`, and the POM activation delay `pom_delay_factor(t, t_delay)` calls it with `POM_DELAY_STEEPNESS = 10` — the delay is the same switch centred on $t_{delay}$ with a width of $t_{delay}/10$, which was not previously visible from either write-out. Neither steepness has a citation; both are working assumptions carried over from the MATLAB code.

`pom_delay_factor` returns exactly 1 when `t_delay ≤ 0`, which is every configuration in `paper/`. Both `timestepper.jl` and `mol.jl` call it, so the two integrators cannot disagree on when POM turns on.

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

### Errata in THIS document — 2026-07

These are cases where the reference manual itself was wrong, not the
manuscript. Recorded because the same failure recurred twice.

| # | Section | Was | Now | Found |
|---|---|---|---|---|
| 6 | §3.3, §9 | α labelled **mobilization**, β **immobilization**; exponent on the loss term | α = immobilization (gain, carries δ); β = mobilization (loss, linear) — matching Falconer (2005) eq. 2.4 / (2008) Box 1 and `parameters.jl` | 2026-07-27 |
| 7 | §9 | `trans_ni = ζ·F_n`, a one-way drain | ζ is a splitting fraction applied inside `fungal_transitions`; the drain form was deleted in February 2026 | 2026-07-27 |
| 8 | §14 | Carried **both** errata 6 and 7 — inverted α/β **and** `±ζ·F_n` — and omitted the O₂ capacity divisor from `S_O` | Transcribed from `reactions.jl` through the code's own intermediates | 2026-07-28 |

**Why erratum 8 happened.** §3.3 and §9 were fixed on 2026-07-27 and the item
was reported closed without §14 being checked. §14 restates the same source
terms in expanded form, so it holds a *second copy* of every convention in §9 —
and a fix applied to one copy does not reach the other.

Two rules follow, and they are the reason this errata table exists:

1. **Search the whole document for the symbol, not the section you remember.**
   `grep -n "α_i\|β_i\|ζ" docs/REFERENCE.md` would have caught it in seconds.
2. **§14 now defers to the intermediates** (`immobil_i`, `immobil_n`,
   `Resp_F^conv`) instead of expanding them. One statement of a convention
   cannot contradict another statement of it if there is only one.

**The repo-wide sweep that followed** (2026-07-28) found the same two
conventions stale in six more places. All fixed:

| File | Was |
|---|---|
| `docs/ARCHITECTURE.md` §4.3 | inverted α/β, `±ζ·F_n`, pre-audit `−J_M·(θ+ρ_b·k_d)/k_d` (erratum 2, fixed in code in February), `S_O` missing its capacity divisor |
| `docs/ARCHITECTURE.md` struct listing | α labelled "Mobilization rate", β "Immobilization rate", ζ "Insulation rate [1/day]" |
| `docs/GUIDE.md` §5 | inverted α/β, η misplaced inside `net`, `trans_ni = ζ·F_n` described as one-way insulation |
| `docs/GUIDE.md` reaction-step listing | `(net_i, net_n, trans_ni)` |
| `docs/REFERENCE.md` §15 catalog | `fungal_transitions` documented as returning `(net_i, net_n, trans_ni)`; it returns `(immobil_i, immobil_n, Resp_F_conv)` |
| `src/physics/initial_conditions.jl:344` | comment gave `S_Fn` with α/β swapped (the arithmetic beneath it used α_n's value under β_n's name, so the conclusion held) |
| `paper/simulations/single_aggregate_physics/run_simulations.jl` | **live code**, read `trans.trans_i` / `trans.trans_n` / `trans.insulation` — fields that have not existed since February. Would throw on first use |

### Erratum 9 — a modelling choice encoded as a test *(2026-07-28)*

`test/test_parameters.jl:31` asserted `bio.delta > 1.0` inside the default-constructor
smoke test, against a default of 1.0. It had been failing for an unknown period.

It was never a valid test. Its neighbours in that testset assert **invariants**
(`Y_B_max <= 1.0`, `Ea_B > 0.0`) — things whose violation means the parameter set
is wrong. `delta > 1.0` asserts a **choice**: δ = 1 is Falconer's own linear
control, not a broken value. A default constant also cannot regress; it changes
only when someone edits `parameters.jl`, which git already records.

Replaced with `@test bio.delta > 0.0`, which *is* an invariant (δ ≤ 0 makes the
gain term constant or divergent). The open question moved to §5a, where open
questions belong.

**The cost of having left it red:** with one known failure, `1 failed` reads as
"the known one" and a genuine regression arriving as `2 failed` goes unnoticed.
A permanently red suite does not preserve a question — it spends the signal that
would tell you about the next real defect.

---

### Erratum 10 — post-processing disagreed with the solver *(2026-07-28)*

`respiration_rates` and `carbon_use_efficiency` (`src/postprocessing/derived.jl`)
re-derive the uptake and transition terms from the stored state. That second
implementation had drifted from `compute_source_terms` in three places:

| Quantity | Was (derived.jl) | Solver (reactions.jl) |
|---|---|---|
| `O_aq` | `O·θ/(θ + K_H·θ_a)` (:205, :309) | `O` (:77) — state.O is already aqueous |
| `C_aq` | `C/(1 + k_d·ρ_b)` (:54) | `C/(θ + ρ_b·k_d_eq)` (:75) |
| `ζ_T` | `ζ·f_fun`, unclamped (:227) | `min(ζ·f_fun, 1.0)` (:115) |

The O₂ one is the consequential one: with `K_H ≈ 34` and θ ≈ 0.4 the factor is
≈ 0.11, so reported respiration and CUE were evaluated at roughly a tenth of the
oxygen the solver used. **Every respiration figure taken from post-processing
before this date is not the quantity the simulation ran on.**

All three are corrected. The oracle that catches this class of drift was already
present and discarded — `derived.jl:188` computed `compute_source_terms(...)` and
threw the result away — so `test_postprocessing.jl` now asserts
`respiration_rates` against `compute_source_terms` node by node at `rtol=1e-12`.

**The general lesson.** Any quantity computed twice will diverge. Where a second
implementation is unavoidable, the test must compare the two, not check that each
is internally consistent — 432 assertions in that file checked structure,
non-negativity and self-sums, and none of them could see a 10× error in the
driving concentration.

---

### Erratum 12 — uninitialised state and workspace arrays *(2026-07-28)*

`AggregateState(n)` and `Workspace(n)` allocated every array with
`Vector{Float64}(undef, n)`. Any pool a caller forgot to set therefore read
whatever happened to be on the heap: wrong, and **non-deterministic**.

Found through a test in `test_physics.jl` that set `E` and `F_i` but not `F_n`,
then read `F_n` through `D_eff_fungi_mobile`. It passed or failed depending on
the run, with no code change between — two consecutive runs of the same suite
disagreed.

Both constructors now fill with **`NaN`**. Of the three options:

| fill | behaviour on a forgotten pool |
|---|---|
| `undef` | wrong answer, non-deterministic, silent |
| `0.0` | wrong answer, deterministic, silent — "no fungi" is perfectly plausible |
| **`NaN`** | **propagates and fails at the first assertion that touches it** |

This matches the convention `TemperatureCache` already used for the Arrhenius
factors, and `test_types.jl` now asserts the contract for all eight state pools
and the workspace scratch arrays. Zero would be the worst choice for the
workspace in particular, where a zero diffusivity is a physically meaningful
value that would quietly stop transport.

`create_initial_state` assigns all eight pools plus `P`, `P_0` and
`CO2_cumulative`, so no production path depends on the fill. Both constructors
run once per simulation, never in the hot loop.

### Erratum 11 — the binding-strength derivation *(2026-07-28)*

`k_F = σ_h·d_p/(4ρ_h)` did not balance dimensionally: `σ_h d_p/ρ_h` has units
Pa·m⁴·kg⁻¹ where Pa·m³·kg⁻¹ is required — one extra length. The stated value
2.25 was recovered only by expressing `d_p` in metres while `ρ_h` was in
µg mm⁻³; in the model's own length unit the identical formula gives 2250.

The cause was not a slipped exponent. The prose described **tensile rupture**
(hyphae snapping across a surface, hence σ_h) while the algebra multiplied by a
**pull-out** length — two different failure mechanisms spliced together. Working
the stated argument through (isotropic strands, ⟨|cos θ|⟩ = ½ twice) gives
`k_F = σ_h/(4ρ_h)` with no `d_p`, which is dimensionally clean and reproduces
the accompanying claim that the result is independent of hyphal diameter.

Two further inconsistencies in the same quantity: `docs/REFERENCE.md` §4.4
listed `k_F` in **kPa** while the code used **Pa**, and the fields `χ` and `a_p`
were carried in `SoilProperties` with no reader anywhere in `src/`.

Replaced by the bond-counting form in §4.4, with the texture dependence carried
by the Sauter mean diameter and the two coefficients declared as fitted.

---

### Erratum 12 — `M_max` was two quantities, and both were wrong *(2026-07-29)*

`SoilProperties` carried an `M_max` field defaulting to **10.0** µg-C/mm³, read
by `reactions.jl`, `api.jl` and `derived.jl`. `initial_conditions.jl` ignored it
and computed `k_ma · f_clay_silt · ρ_b` instead. For De Gryze soil 3 that was
**368 against 10**.

Consequence: the state was initialised against one ceiling and evolved against
another an order of magnitude lower, so `M(0) ≫ M_eq` at every node and every
run opened with a desorption flux that had no equilibrium to reach. The guard in
`create_initial_state` compared `M_0` against the texture value, not the
solver's, so it never fired. In `paper/de_gryze/degryze_config.jl` the symptom
was being suppressed by lowering `κ_d_ref` tenfold — a parameter compensating
for a duplication defect.

The 368 was independently wrong. `k_ma = 0.48` with the unit "µg-C per gram of
mineral" is Georgiou et al. (2022)'s 48 mg-C g⁻¹ with the mg → g conversion
dropped. The formula needs a **mass fraction**; read literally the stated unit
puts the product short by 10⁶, and read as the code used it, 10× high.

Fixed by deleting the field, making `maoc_capacity(soil)` the only definition,
and redefining `k_ma` as dimensionless with the two Georgiou values in
`src/constants.jl`. `test_parameters.jl` now asserts the field's **absence** and
bounds `k_ma`, so neither arm can return silently.

This is Jeppson's drift arm in both directions at once: two implementations of
one quantity, and a third statement in a docstring (`maoc.jl:44` presented the
texture formula as definitional while the code it documented used the field).

---

### Erratum 13 — ω: four open questions between the drafts *(2026-07-29)*

Not an erratum in the sense of the others — nothing here is settled and nothing
is being corrected. It records where two working drafts and the code disagree
about the overlap correction, so the disagreements are not rediscovered.

`SI_tessellation.tex` (Overleaf) has the fullest treatment: it derives `f_p`,
`f_m`, `ω = (f_m/f_p)³` and the `C_bg/ω` dilution, proves exactness for
diffusion and linear reactions, identifies the Monod nonlinearity, bounds the
error by `1/ω`, and gives an implementation table. **The code matches that
table.** The two proofs stand on their own — they do not depend on which
document is senior.

An independent trace on 2026-07-29 reached the same Monod result by another
route, which is corroboration, and then raised four things the SI does not
settle.

**1. Is oxygen diluted?** The SI's error analysis writes
`(O/ω)/(L_B + O/ω)` — it dilutes O₂. The code does not
(`create_initial_state`, Step 8: "O₂ is NOT diluted — it is a boundary
condition"). So `S_O = -α_O·Resp/capacity` divides a respiration built from
diluted carbon by an undiluted oxygen capacity. Measured on soil 3, the O₂ sink
is ~45× weaker than the physical one and `monod_O` never leaves 0.87, so
near-POM anoxia does not form. **A physics decision, not a bug report:** O₂ is a
boundary condition supplied from the headspace, which is an argument for leaving
it undiluted — but then the sink term and the source term are on different
scales. Unresolved.

**2. Does the validity argument survive the configuration?** SI S7 rests the
case on background carbon being negligible: "the native POM was removed, so
`C_bg ≈ 0`". `degryze_ic` initialises the full measured SOC — 2.14 % for soil 3,
`SOC_vol = 29.3` µg/mm³, partitioned across DOC, biomass, EPS and MAOC. That is
not trace seeding. The `1/ω` bound is unaffected; the argument that the bound
does not matter is. Either the argument needs revising or the configuration
does.

**3. Terms not yet analysed.** The SI treats the Monod form only. Also nonlinear
in a diluted state against an undiluted constant: the Langmuir–Freundlich
isotherm (`1/k_L`, `M_max`), the sigmoid thresholds (`B_min`, `E_min`,
`F_i_min` — widths of `X_min/50`, effectively switches), the space-limited
yields (`B_S`, `F_S`), the POM-dissolution Monods (`K_B_P`, `K_F_P`), the
translocation network factor (`K_Fm_net`) and `ε_F` inside `Π`. Drafting work,
not a defect.

The binder against `G_c` is a separate case and is **not** a missing `1/ω`: the
binder is a sum of a diluted background contribution and an undiluted
residue-derived one, so no single factor applies. Documented as a property of the
closure in §4.4.

**4. Main text and SI describe different constructions.** `manuscript-4-5.tex`
§2 (l.165–181) uses a single domain factor `f_d`, sets `r_max,j = f_d·r_0,j`,
and states the tessellation "satisfies the non-overlapping constraint exactly".
That identity holds only for `f_d = f_p`. The SI's entire subject is
`f_m > f_p`, where domains do overlap. This one has to be reconciled before
submission whatever else is decided.

Also: the SI's worked example tabulates `ρ_POM = 100` µg-C/mm³ (`f_p = 2.55`,
`ω = 60.2`) with 200 in a footnote. The code uses 200 (`f_p = 3.21`,
`ω = 30.1`).

**The overlap itself is not in question.** It is required — the model will be
run on annual root input and field scenarios where model volume necessarily
exceeds physical volume — and explicit domain coupling would need packing detail
the model does not carry. The open question is where the `1/ω` is applied. See
`docs/AUDIT.md`.

---

### Erratum 14 — `s_M`: two isotherms for one quantity *(2026-07-29)*

`partition_CM` placed the initial MAOC at `M = s_M · M_max · f_LF(β·C)` with
`s_M = 0.6` for De Gryze. `reactions.jl` drives `M` toward
`M_max · f_LF(β·C)` — **no `s_M`**. So the state opened at 60 % of what the
solver calls equilibrium and sorbed from `t = 0` to close a gap nothing in the
experiment had created. Two implementations of one quantity, at any `ω`.

Two further defects in the same parameter:

- **The name did not match the code.** `s_M` was documented as "fraction of
  mineral capacity occupied", i.e. `M/M_max` — the convention Georgiou et al.
  (2022) report saturation deficits in. The code made it `M/M_eq_full`, the
  fraction of *local equilibrium*. Since `f_LF < 1` always, the achieved
  capacity-saturation was always below the number requested: `s_M = 0.6` on soil
  3 delivered **0.485**.
- **The diagnostic hid it.** The `@info` line reported `M/M_eq_full`, so it
  printed 60.0 % and looked correct.

`s_M` also had no entry in §5a, §2–§4, the manuscript or the SI. It was a free
parameter that partitioned background carbon and was written down nowhere.

**Fixed** by deleting it — from `partition_CM`, from `InitialConditions`, and
from every call site. The partition now starts on the isotherm, so `J_M(0) = 0`,
and saturation is an output with a measured ceiling (§5d). One fewer degree of
freedom.

**Consequence, and it needed a second change.** Removing `s_M` moved soil 3 from
`C = 11.14, M = 17.86` to `C = 4.09, M = 24.91` — better, but pore-water DOC was
still 573 mg/L against an observed 10–100. Cause: the partition has two pools and
no third, so carbon the isotherm will not take is declared dissolved by default.
De Gryze destroyed native POM, so background SOC there *is* mineral-associated —
the isotherm was simply refusing it, at `k_L = 1000`. `k_L` raised to 25000,
which puts DOC at 48 mg/L with the minerals holding 28.7 of 29.0. **`k_L` is
Group C with no anchor; that DOC range is the first constraint ever placed on
it**, and it is a physical one rather than a fit to a compared observable.

`k_d_eq` was also implicated and is *not* changed: it moves solution DOC strongly
(50× range → 1398 down to 29 mg/L) and MAOC not at all, because
`C_eq = k_d·C/(θ+ρ_b·k_d) → C/ρ_b` saturates. So it cannot substitute for `k_L`,
and its 0.05 → 0.005 cut on 2026-07-29 stands on the transport argument alone.

**Renamed** `partition_CM`'s first argument `R` → `SOC_residual`. The manuscript
uses `\mathcal{R}` for total respiration and `R_*` is a rate throughout `src/`;
a carbon stock called `R` collided with both.

**Documented**, for the first time, in three places at once: §5d here,
`manuscript-4-5.tex` §`sec:initial_condition` (which also states the structural
limitation for intact soils), and `test_api.jl`, which asserts mass balance,
`J_M(0) = 0` against the solver's own isotherm, monotonicity of saturation in
`k_L`, convergence to the `SOC_residual/M_max` ceiling, the DOC band that set
`k_L`, and the absence of the `s_M` field.

---

Historical documents that *quote* the old form were left alone on purpose:
`docs/STRUCTURAL_AUDIT_2026-07-27.md`, `dev_notes/falconer_questions.md`
(the brief describes what we had when we asked), `docs/julia_falconer_deviations.md`
and everything under `dev_notes/archive/`.

---

**End of Reference**
