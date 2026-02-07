# parameters.jl — Bug Report and Correction Guide

**Date**: 2026-02-06  
**Scope**: All default parameter values in `parameters.jl`, cross-referenced against `manuscript_updated.tex` (authoritative) and `ARCHITECTURE_CLAUDE_CODE.md`.  
**Rule**: Manuscript wins on physics. Architecture wins on implementation detail. If both agree and code disagrees, code is wrong.

---

## HOW TO USE THIS DOCUMENT

For each bug: read the **Problem**, **Evidence**, **Fix**, and **Cascade** sections. Apply fixes in the order listed (Section 1 first — ρ_b cascades everywhere). After fixing `parameters.jl`, grep for each affected quantity in all `.jl` files to verify no hardcoded duplicates exist.

Run the full test suite after each fix. If tests fail, the failure reveals secondary bugs masked by the original error.

---

## UNIT SYSTEM REMINDER

| Quantity | Unit | Equivalence |
|----------|------|-------------|
| Concentration | μg/mm³ | = kg/m³ = g/L |
| Length | mm | |
| Time | day | |
| Pressure | kPa | |
| Temperature | K | |
| Energy | J/mol | |
| Diffusion | mm²/day | 1 cm²/s = 8.64×10⁶ mm²/day |
| Density | μg/mm³ | 1 g/cm³ = 10³ μg/mm³ (NOT 10⁶) |

**Critical conversion**: 1 g/cm³ = 1 g/(10 mm)³ = 1 g / 10³ mm³ = 10⁻³ g/mm³ = 10³ μg/mm³.  
NOT 10⁶. The factor-of-10³ error in ρ_b likely arose from confusing μg/mm³ with kg/m³ (they ARE equal) but then separately computing 1.3 g/cm³ × 10⁶ μg/g = 1.3×10⁶, forgetting the cm³→mm³ conversion (÷10³).

---

## 1. ρ_b — CRITICAL (cascades everywhere)

### Problem
```julia
ρ_b = 1.3e6,        # [μg/mm³] = 1.3 g/cm³    ← WRONG
```

### Evidence
1.3 g/cm³ = 1.3 × 10³ μg/mm³ = 1300.0 μg/mm³.

Proof: 1 g = 10⁶ μg. 1 cm³ = 10³ mm³.  
So 1 g/cm³ = 10⁶ μg / 10³ mm³ = 10³ μg/mm³.

Alternatively: μg/mm³ = (10⁻⁹ kg)/(10⁻⁹ m³) = kg/m³. And 1.3 g/cm³ = 1300 kg/m³ = 1300 μg/mm³. ✓

### Fix
```julia
ρ_b = 1.3e3,        # [μg/mm³] = 1.3 g/cm³
```

### Cascade
ρ_b appears in:
- `C_aq = C / (θ + ρ_b·k_d)` — aqueous concentration (reactions.jl line 80)
- `C_eq = k_d·C_aq` — equilibrium sorbed (reactions.jl line 81)
- `D_eff_DOC` — retardation factor θ/(θ+ρ_b·k_d) (effective_diffusion.jl line 42)
- `D_eff_bacteria` — inherits from D_C
- `M_max = k_ma·f_clay_silt·ρ_b` — if formula is used (see Bug #7)
- `S_C` MAOC coupling factor (θ+ρ_b·k_d)/k_d (see Bug #3)

With ρ_b = 1.3e6 (wrong): ρ_b·k_d = 650,000 → C_aq ≈ C/650,000 → all microbes starved.  
With ρ_b = 1300 (correct): ρ_b·k_d = 650 → C_aq = C/650 → retardation ~2000×.  

NOTE: After fixing ρ_b, k_d_eq = 0.5 gives retardation factor R = 1 + 650/0.3 ≈ 2167. This is physically high but not impossible for strongly sorbing soils. Revisit after initial runs — see Bug #4.

---

## 2. ρ_POM — CRITICAL (same conversion error)

### Problem
```julia
ρ_POM = 1.0e6,      # POM carbon density [μg-C/mm³]    ← WRONG
```

### Evidence
POM is plant residue. Bulk density ~1.0–1.4 g/cm³, carbon fraction ~0.45.  
POM carbon density ≈ 0.45–0.63 g/cm³ = 450–630 μg/mm³.

Same error pattern as ρ_b: multiplied by 10⁶ instead of 10³.

### Fix
```julia
ρ_POM = 500.0,       # POM carbon density [μg-C/mm³] ≈ 0.5 g-C/cm³
```

Adjust to match your literature value. Plausible range: 300–700 μg/mm³.

### Cascade
ρ_POM is used in POM geometry (radius↔mass conversion):
- `P = (4/3)π r₀³ ρ_POM` or equivalently `r₀ = (3P/(4π ρ_POM))^(1/3)`
- With ρ_POM 1000× too large: for a given P_0, the inferred aggregate domain is 10× too small, or for a given r_0, the inferred mass is 1000× too large.

Check whether P_0 = 1000 μg-C is consistent with r_0 = 0.1 mm at corrected ρ_POM:
- `P = (4/3)π(0.1)³ × 500 = (4/3)π × 0.001 × 500 = 2.09 μg-C`

So P_0 = 1000 μg-C with r_0 = 0.1 mm and ρ_POM = 500 are inconsistent. Either:
- r_0 = (3×1000/(4π×500))^(1/3) = (0.477)^(1/3) = 0.781 mm (POM radius ~0.8 mm if P_0 = 1000), or
- P_0 = 2.09 μg-C if r_0 = 0.1 mm and ρ_POM = 500

**Decision needed**: Which two of {P_0, r_0, ρ_POM} are primary? The third should be computed.

---

## 3. S_C MAOC coupling factor — CODE IS CORRECT, manuscript/architecture need update

### Status: NOT A CODE BUG — document correction needed

### Analysis
The code (reactions.jl line 162) is correct:
```julia
S_C_val = -R_B_val - R_F_val + R_rec_total - J_M_val
```

$J_M \equiv \partial M / \partial t$ is defined in total-volume units [μg-C/mm³/day], same basis as all other terms in $S_C$. When MAOC absorbs carbon at rate $J_M$, the total dissolved carbon pool $C$ (which includes both $C_{aq}$ and the equilibrium-sorbed $C_{eq}$) decreases by exactly $J_M$. No coupling factor is needed.

Proof via carbon conservation: $S_C + S_M = (-J_M) + (J_M) = 0$ for MAOC terms. ✓

If the factor $(θ + ρ_b k_d)/k_d$ were included, the sum would be $-J_M \cdot R_f + J_M \neq 0$, **breaking** conservation.

### Action needed: UPDATE MANUSCRIPT AND ARCHITECTURE
Both manuscript (line 429) and architecture (§4.3, line 131) incorrectly include the factor. The corrected $S_C$ equation should read:
```
S_C = −R_B − R_F + R_rec − J_M
```

Also update `maoc.jl` comments (lines 15–16, 115–116) which reference the incorrect factor — these "CRITICAL" warnings are stale and misleading.

---

## 4. k_d_eq — REVIEW AFTER ρ_b FIX

### Problem
```julia
k_d_eq = 0.5,       # [mm³/μg]
```

### Evidence
With corrected ρ_b = 1300: ρ_b·k_d = 650. Retardation factor at θ = 0.3: R ≈ 2167.

Literature K_d for dissolved organic carbon: typically 1–100 L/kg.

Unit conversion: 1 L/kg = 1000 mL / 1000 g = 1 cm³/g = 1000 mm³ / 10⁶ μg = **10⁻³ mm³/μg**.

So literature range: k_d = 10⁻³ to 10⁻¹ mm³/μg.

k_d_eq = 0.5 is at the high end (500 L/kg equivalent). For strongly sorbing soils with high organic matter, this is physically possible but extreme. Typical mineral soils: k_d ≈ 0.001–0.01 mm³/μg (1–10 L/kg).

### Recommendation
- If k_d = 0.01: ρ_b·k_d = 13, retardation = 44. C_aq = C/44. Moderate sorption.
- If k_d = 0.1: ρ_b·k_d = 130, retardation = 434. Strong sorption.
- If k_d = 0.5: ρ_b·k_d = 650, retardation = 2167. Extreme sorption — nearly all DOC immobile.

**Decision needed**: What retardation factor is physically appropriate for your target soil? This controls effective DOC diffusion and microbial substrate availability. A value of k_d = 0.01–0.1 mm³/μg is more typical.

### Cascade
k_d_eq affects:
- C_aq computation
- Effective DOC diffusion (D_C has θ/(θ+ρ_b·k_d) factor)
- S_C MAOC coupling factor
- K_E effectiveness (see Bug #5)
- All microbial uptake rates (via C_aq)

---

## 5. K_E — needs recalibration for C_aq

### Problem
```julia
K_E = 200.0,        # Substrate inhibition concentration [μg/mm³]
```

### Evidence
Manuscript (line 391) writes the EPS degradation formula with $C$ (total dissolved carbon):
```
R_rec_E = μ_E^max(T) · K_E/(K_E + C) · E · h_E
```

Code (eps.jl line 74) correctly uses C_aq instead of C:
```julia
inhibition = K_E / (K_E + C_aq)
```

This correction is documented as "MANUSCRIPT_CHANGES #3". The issue is that K_E = 200 was likely set assuming comparison against total C, not C_aq. With any meaningful retardation, C_aq ≪ C, so K_E/(K_E + C_aq) ≈ 1.0 always — EPS degradation is never inhibited by dissolved carbon abundance.

### Fix
K_E should be set relative to expected C_aq values. If typical C_aq ≈ 0.1–10 μg/mm³ (depends on k_d fix), then K_E should be in the same range:
```julia
K_E = 5.0,          # Substrate inhibition for EPS [μg/mm³] — in C_aq units
```

**Decision needed**: Set K_E to the C_aq concentration at which you want EPS degradation to be 50% inhibited. Literature guidance sparse — this is a model-specific parameter.

---

## 6. K_Y units documentation — HIGH (misleading)

### Problem
```julia
K_Y::Float64            # Half-saturation for yield [-]       ← WRONG UNITS
...
K_Y = 1.0,
```

### Evidence
In bacteria.jl:
```julia
Y_B_max * R_diff_pos / (R_diff_pos + K_Y)
```
R_diff = R_B - R_Bb has units [μg-C/mm³/day]. For Monod to be dimensionally consistent, K_Y must have same units.

### Fix
```julia
K_Y::Float64            # Half-saturation for yield [μg-C/mm³/day]
```

The default value K_Y = 1.0 μg-C/mm³/day may or may not be appropriate. Check: if R_diff at half-max uptake is ~1 μg-C/mm³/day, does that match expected bacterial growth rates? At r_B_max = 5.0/day, B = 1 μg/mm³, C_aq at half-saturation (K_B), and full O₂: R_B ≈ 5 × 0.5 × 1 × 1 = 2.5 μg/mm³/day. If R_Bb ≈ 0.5, then R_diff ≈ 2.0, and Y_B = 0.5 × 2/(2+1) = 0.33. This seems reasonable — yield at ~2/3 of max. Acceptable default.

---

## 7. M_max inconsistency with k_ma formula — MEDIUM

### Problem
```julia
M_max = 100.0,      # [μg/mm³]               ← hardcoded
k_ma = 0.05,        # [μg-C/g-mineral]
f_clay_silt = 0.5,  # [-]
```

Manuscript (line 266): M_max = k_ma · f_clay_silt · ρ_b

### Evidence
With corrected ρ_b = 1300 μg/mm³ = 1.3 × 10⁻³ g/mm³:

M_max = 0.05 [μg-C/g-mineral] × 0.5 × 1.3e-3 [g/mm³] = 3.25e-5 μg-C/mm³

This is essentially zero. The formula doesn't work because k_ma is in μg-C per g-mineral, but ρ_b is in μg/mm³ not g/mm³.

### Root cause
Units mismatch. The formula M_max = k_ma · f_clay_silt · ρ_b requires consistent units. Options:

**Option A**: Keep M_max as a directly specified parameter (100 μg/mm³ seems high — literature MAOC is typically 10–50 mg-C/g-soil). In model units: 10–50 mg-C/g-soil × ρ_b_in_g/mm³ = 10–50 × 1.3e-3 = 0.013–0.065 μg-C/mm³. Hmm, that's tiny.

Wait — let me redo this more carefully. MAOC is typically 10–50 mg-C/g-soil = 10–50 μg-C/mg-soil. In volumetric: multiply by ρ_b in mg/mm³:

ρ_b = 1300 μg/mm³ = 1.3 mg/mm³

M_max = 10–50 [μg-C/mg-soil] × 1.3 [mg-soil/mm³] = 13–65 μg-C/mm³

So M_max = 100 is high but not unreasonable for the *capacity* (actual stocks are below capacity). Abramoff et al. report soils at 21–42% of capacity.

**Option B**: Fix k_ma units. If k_ma should give M_max in μg/mm³ when multiplied by f_clay_silt and ρ_b in μg/mm³:

k_ma [μg-C/μg-mineral] × f_clay_silt [-] × ρ_b [μg/mm³] → μg-C/mm³

k_ma = M_max / (f_clay_silt × ρ_b) = 100 / (0.5 × 1300) = 0.154 [μg-C/μg-mineral]

But the comment says [μg-C/g-mineral], not [μg-C/μg-mineral]. Converting: 0.154 μg-C/μg-mineral = 154,000 μg-C/g-mineral = 154 mg-C/g-mineral. Too high.

### Fix recommendation
**Simplest**: Remove the formula. Keep M_max as a direct input parameter. Remove or document k_ma and f_clay_silt as informational only. The formula adds confusion without value when M_max can be directly constrained from measurements.

If keeping the formula: redefine k_ma in units that work. With ρ_b in μg/mm³:
```julia
k_ma = 0.154,       # [μg-C/μg-mineral] = 154 mg-C/g-mineral (MAOC capacity per unit mineral)
```
Then M_max = 0.154 × 0.5 × 1300 = 100 μg/mm³. ✓

**Decision needed**: Direct M_max or computed from formula?

---

## 8. D_B_rel — manuscript/code disagreement

### Problem
```julia
D_B_rel = 0.001,              # Bacterial motility << DOC diffusion
```

### Evidence
Manuscript (line 507) states: "D_{B,rel} ≈ 0.5 based on experimental observations of bacterial diffusivity relative to dissolved organic substrates."

Code has 0.001 — factor of 500 discrepancy.

### Analysis
D_B_rel = 0.5 means bacteria diffuse at half the rate of DOC. This is physically reasonable for motile bacteria in saturated conditions (swimming speeds ~30 μm/s, comparable to DOC molecular diffusion).

D_B_rel = 0.001 means bacteria are essentially immobile relative to DOC. This may be the intended behavior for soil bacteria in unsaturated conditions where motility is severely limited.

### Fix
**Decision needed**: Which value reflects your intended physics? If manuscript is authoritative:
```julia
D_B_rel = 0.5,               # Bacterial motility relative to D_C [-] (manuscript)
```

If the code value was a deliberate correction from an unrealistic manuscript value:
- Update manuscript to D_B_rel ≈ 0.001
- Add comment explaining the change

---

## 9. λ — placeholder value

### Problem
```julia
λ = 0.5,            # Fraction of F_n at uptake surfaces [-]
```

### Evidence
Manuscript (line 279): "λ ≪ 1 reflects the small fraction of non-insulated hyphae at active uptake surfaces."

λ = 0.5 is not ≪ 1. With this value, F_n contributes half its biomass to uptake — nearly equivalent to F_i.

### Fix
```julia
λ = 0.05,           # Fraction of F_n at uptake surfaces [-] (λ ≪ 1)
```

Or even smaller (0.01–0.1 range). This determines how much fungal uptake comes from non-insulated vs. insulated hyphae.

---

## 10. R_GAS duplicate definition — LOW (will cause warning/error)

### Problem
Both `arrhenius.jl` (line 10) and `henry.jl` (line 12) define:
```julia
const R_GAS = 8.314
```

Also `diffusion_pure.jl` (line 211) defines a local `R_GAS = 8.314`.

### Fix
Define R_GAS once in a constants file or in the module header. Import everywhere else.

```julia
# constants.jl
const R_GAS = 8.314  # [J/(mol·K)]
```

Remove definitions from arrhenius.jl, henry.jl, and diffusion_pure.jl.

---

## 11. Diffusion reference values — VERIFY

### Current values
```julia
D_C0_ref = 1.0e-5 * 8.64e6,   # DOC: ~1e-5 cm²/s → 86.4 mm²/day
D_O2_w_ref = 2.0e-5 * 8.64e6, # O₂ in water: ~2e-5 cm²/s → 172.8 mm²/day
D_O2_a_ref = 0.2 * 8.64e6,    # O₂ in air: ~0.2 cm²/s → 1.728e6 mm²/day
```

### Evidence
These are stored as reference values in SoilProperties but diffusion_pure.jl implements the full physical formulas (Stokes-Einstein+VFT, Han-Bartels, Chapman-Enskog). The TemperatureCache should use the formula-based functions, not these reference values directly.

Verify that the TemperatureCache actually calls D_DOC_water(), D_O2_water(), D_O2_air() from diffusion_pure.jl, using these _ref values only as the reference point for Stokes-Einstein scaling. If the code just multiplies _ref × Arrhenius factor (wrong for diffusion), the temperature dependence is incorrect.

### Status
Not a bug in parameters.jl per se, but verify the temperature coupling is using the correct physical formulas from diffusion_pure.jl.

---

## SUMMARY: REQUIRED ACTIONS

### Must fix immediately (blocks all simulation):
1. **ρ_b**: 1.3e6 → 1300.0
2. **ρ_POM**: 1.0e6 → ~500.0 (confirm with literature value)

### Must fix in manuscript/architecture (code is correct):
3. **S_C MAOC factor**: Remove (θ + ρ_b·k_d)/k_d from S_C equation in manuscript line 429 and architecture §4.3 line 131. Also clean up stale "CRITICAL" comments in maoc.jl.

### Must fix before calibration:
4. **k_d_eq**: Review after ρ_b fix. Likely reduce from 0.5 to 0.01–0.1
5. **K_E**: Recalibrate for C_aq units. Likely reduce from 200 to 1–10
6. **K_Y**: Fix units documentation from [-] to [μg-C/mm³/day]
7. **M_max/k_ma**: Resolve formula vs. direct value. Fix units.
8. **D_B_rel**: Resolve manuscript (0.5) vs. code (0.001) discrepancy
9. **λ**: Reduce from 0.5 to 0.01–0.1

### Cleanup:
10. **R_GAS**: Consolidate to single definition
11. **Diffusion refs**: Verify TemperatureCache uses physical formulas

---

## VERIFICATION PROTOCOL

After applying fixes:

1. **Unit check**: For every parameter with physical units, verify: value × unit = physically meaningful quantity. Use the conversion 1 g/cm³ = 10³ μg/mm³.

2. **Carbon conservation**: Run 30-day simulation. Verify P₀ = P(t) + ∫(C+B+F_i+F_n+F_m+E+M)·4πr²dr + CO₂(t) to machine precision.

3. **Steady-state sanity**: At constant T, ψ, O₂:
   - C_aq should be ~0.1–10 μg/mm³ (not 10⁻⁴ or 10³)
   - Bacterial biomass should grow initially then stabilize
   - MAOC should accumulate slowly
   - Total respiration should be positive and bounded

4. **Diffusion sanity**: D_C_eff at θ = 0.3 should be ~0.01–0.1 mm²/day (not 10⁻⁶ or 10²)

5. **Retardation check**: Print C_aq/C ratio at several nodes. Should be ~10⁻² to 10⁻¹ (not 10⁻⁶).
