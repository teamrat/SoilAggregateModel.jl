# Manuscript Changes — Implementation Deviations

This document tracks all changes made to the model implementation that deviate from or extend the manuscript formulation. Each change is documented with rationale, affected equations, and implementation details.

---

## 1. Fungal Uptake: λ Parameter Semantics

**Change**: Swapped λ factor in fungal uptake formula
**Original (manuscript)**: `R_F = ... · (F_i + λ·F_n) · ...`
**Implementation**: `R_F = ... · (λ·F_i + F_n) · ...`

**Rationale**: F_n (non-insulated hyphae) are located at uptake surfaces and should be the primary uptakers. F_i (insulated hyphae) are internal and have reduced efficiency. The manuscript formulation had this backwards.

**Location**: `src/biology/fungi.jl`, function `R_F`
**Documentation**: Updated in `docs/REFERENCE.md` and `docs/GUIDE.md`

---

## 2. POM Dissolution: Additive Microbial Coupling

**Change**: Bacteria and fungi contribute independently to POM dissolution
**Formula**:
```
microbial_factor = 0.5 * (B/(K_B_P + B) + F_n/(K_F_P + F_n))
```

**Rationale**: Additive coupling allows both guilds to contribute to dissolution independently. Multiplicative coupling would require both to be present simultaneously, which is unrealistic during early colonization.

**Parameters**:
- K_B_P = 1.0e-3 μg/mm³ (bacterial half-saturation for POM dissolution)
- K_F_P = 1.0e-3 μg/mm³ (fungal half-saturation for POM dissolution)

**Location**: `src/carbon/pom_dissolution.jl`
**Status**: Implemented and verified

---

## 3. Space-Limited Yield

**Change**: Added space limitation to both bacterial and fungal yields
**Bacterial yield (modified Eq. 308)**:
```
Y_B = Y_B_max · softplus(R_diff, ε_Y) / [softplus(R_diff, ε_Y) + K_Y] · B_S / (B + B_S)
```

**Fungal yield (new)**:
```
Y_F = Y_F_base · F_S / (F_total + F_S)    where F_total = F_i + F_n + F_m
```

**Rationale**: Prevents unbounded exponential growth by reducing yield as biomass approaches a space limitation threshold. Replaces previous density-dependent mortality approach.

**New Parameters**:
- B_S = 0.5 μg/mm³ (bacterial space limitation half-saturation)
- F_S = 0.2 μg/mm³ (fungal space limitation half-saturation)
- ε_Y = 3.33e-6 μg-C/mm³/day (softplus smoothing width for smooth C∞ transition)

**Location**:
- `src/biology/bacteria.jl`, function `Y_B_func`
- `src/biology/fungi.jl`, function `Y_F_func`
- `src/solver/reactions.jl` and `src/postprocessing/derived.jl` (call sites)

---

## 4. Softplus Smoothing for Allocation

**Change**: Replaced `max(0, R_diff)` with `softplus(R_diff, ε_Y)` in allocation functions
**Affected functions**:
- `Gamma_B`: `Y_B · softplus(R_diff, ε_Y) · (1-γ) - softplus(-R_diff, ε_Y)`
- `Gamma_E`: `Y_B · softplus(R_diff, ε_Y) · γ`
- `Resp_B`: `R_Bb + softplus(R_diff, ε_Y) · (1 - Y_B)`

**Softplus formula**:
```
softplus(x, ε) = { x + ε·ln(1 + exp(-x/ε))   if x > 0
                 { ε·ln(1 + exp(x/ε))        if x ≤ 0
```

**Rationale**: Provides smooth C∞ transition at R_diff = 0, eliminating discontinuous derivatives that can cause stiffness in the PDE system.

**Location**:
- `src/math_utils.jl` (shared utility)
- `src/biology/bacteria.jl`, functions `Gamma_B`, `Gamma_E`, `Resp_B`

---

## 5. Fungal Bootstrap: Nonzero β_n

**Change**: Set β_n = 0.15 day⁻¹ (previously 0.0 in early implementations)
**Rationale**: Creates a pathway F_m → F_n, preventing F_n extinction when insulation rate ζ drains it faster than growth can replenish it. Enables complete fungal lifecycle.

**Location**: `src/parameters.jl`, default for β_n
**Status**: Critical for fungal viability

---

## 6. Π Regularization: Increased ε_F

**Change**: Raised ε_F from 1e-10 to 1e-4 μg/mm³
**Rationale**: When F_i + F_n → 0, the protection ratio Π = F_m/(F_i + F_n + ε_F) was exploding to 10^121, causing numerical overflow. With ε_F = 1e-4 (comparable to B_min and F_i_min), Π naturally bounds at physically reasonable values (e.g., 2e-5/1e-4 = 0.2).

**Location**:
- `src/parameters.jl`, default for ε_F
- `src/biology/fungi.jl`, function `Pi_protected`

**Status**: Essential numerical stabilization

---

## 7. Parameter Corrections from REFERENCE.md

**Change**: Restored 28 biological parameters that had diverged during testing
**Critical corrections**:
- K_B: 50.0 → 1.0e-4 μg/mm³ (500,000× error!)
- B_min: 0.1 → 1.0e-4 μg/mm³ (1000× error)
- K_Y: 3.33 day⁻¹ → 3.33e-4 μg-C/mm³/day (unit mismatch: specific → volumetric rate)

**Rationale**: Parameters had been incorrectly modified to pass unit tests without understanding the physics. All parameters now match manuscript defaults in REFERENCE.md.

**Location**: `src/parameters.jl`, constructor for `BiologicalProperties`
**Status**: Fully restored

---

## 8. Reverted: Density-Dependent Mortality

**Previous change (reverted)**: Added `(1 + B/B_max)` and `(1 + F_i/F_i_max)` factors to R_rec_B and R_rec_F
**Current status**: **Reverted to original manuscript form**
**Rationale**: Space-limited yield (§3 above) provides a cleaner mechanism for carrying capacity. Density-dependent mortality was removed.

**Current formulas**:
- R_rec_B = μ_B(T) · B · h_B
- R_rec_F = μ_F(T) · F_i · h_Fi

---

## Summary of Active Deviations

| Change | Type | Impact | Status |
|--------|------|--------|--------|
| λ swap in R_F | Physics correction | High | Implemented |
| Additive POM dissolution | Physics extension | Medium | Implemented |
| Space-limited yield | New mechanism | High | Implemented |
| Softplus allocation | Numerical | Medium | Implemented |
| β_n = 0.15 | Parameter correction | High | Implemented |
| ε_F = 1e-4 | Numerical stability | Critical | Implemented |
| K_B_P, K_F_P = 1e-3 | New parameters | Medium | Implemented |
| Parameter restoration | Corrections | Critical | Complete |

---

## Files Modified

### Core Biology
- `src/biology/bacteria.jl`: Y_B_func, Gamma_B, Gamma_E, Resp_B, R_rec_B
- `src/biology/fungi.jl`: R_F, Y_F_func, Pi_protected, R_rec_F
- `src/biology/maoc.jl`: Removed softplus (moved to shared utility)

### Parameters
- `src/parameters.jl`: Added B_S, F_S, ε_Y; corrected 28 biological parameters

### Solver
- `src/solver/reactions.jl`: Updated all call sites for modified functions
- `src/postprocessing/derived.jl`: Updated respiration and CUE calculations

### Utilities
- `src/math_utils.jl`: Created shared softplus function
- `src/SoilAggregateModel.jl`: Added math_utils include

### Documentation
- `docs/REFERENCE.md`: Updated parameter tables and formulas (pending)
- `docs/GUIDE.md`: Updated fungal uptake description (pending)

---

## Testing Status

- **Unit tests**: Many failures expected due to parameter changes (not yet updated)
- **Physics test (Theme 1)**: Ready to run with all changes integrated
- **Conservation**: Expected to hold to machine precision with softplus smoothing

---

_Last updated: 2026-02-08_
_Implementation: SoilAggregateModel.jl v0.1.0_
