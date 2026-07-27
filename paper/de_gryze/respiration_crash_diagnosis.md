# Respiration Boom-Bust Diagnosis (De Gryze Simulation)

Symptom: cumulative CO₂ flux rises exponentially then crashes abruptly at ~day 12,
regardless of parameter combinations tested.

---

## 1. Initial DOC is Essentially Zero

Trace the initialization chain for the de Gryze parameters
(`run_degryze.jl`, `initial_conditions.jl:188`).

**Physical SOC (before ω dilution):**

```
SOC_vol = 0.022 × 1300               = 28.6 μg-C/mm³
B_0   = 0.01 × 28.6                  = 0.286
F_i_0 = 0.01 × 28.6                  = 0.286
E_0   = 0.005 × 28.6                 = 0.143
biotic total                         ≈ 0.715 μg-C/mm³
Residual R = 28.6 − 0.715            ≈ 27.885 μg-C/mm³
```

**MAOC partition** (`partition_CM`, `initial_conditions.jl:103`):

```
M_max = k_ma × f_clay_silt × ρ_b = 0.48 × 0.74 × 1300  = 461.8 μg-C/mm³
s_M × M_max = 0.4 × 461.8                               = 184.7 μg-C/mm³
```

Since `s_M × M_max = 184.7 >> R = 27.885`, the Langmuir-Freundlich isotherm
always wants to absorb more than the available residual. The iteration hits the cap
`M_new = min(..., R − C_floor)` on the first step and converges to:

```
M_0 ≈ 27.885 μg-C/mm³     (nearly all residual)
C_0 ≈ C_floor = 1e-6 μg-C/mm³     (enforced minimum)
```

**After ω = 28.8 dilution** (`initial_conditions.jl:257`):

```
C_0_model = 1e-6 / 28.8  ≈ 3.5e-8 μg-C/mm³   ← essentially zero
M_0_model = 27.885 / 28.8 ≈ 0.968 μg-C/mm³
B_0_model = 0.286  / 28.8 ≈ 0.010 μg-C/mm³
```

**Effective DOC available to microbes** (retardation factor, `reactions.jl`):

```
C_aq = C / (θ + ρ_b × k_d_eq)
     = 3.5e-8 / (0.292 + 1300 × 0.05)
     = 3.5e-8 / 65.3
     ≈ 5.4e-10 μg-C/mm³
```

This is five orders of magnitude below `K_B = 1e-4`. Bacterial uptake from native
SOC is negligible at t = 0.

---

## 2. POM is the Only Carbon Source — and It is Undiluted

POM mass is **not divided by ω** (`initial_conditions.jl:277`). For d = 1.25 mm:

```
P_0 = (4/3)π(0.625)³ × 200 ≈ 204 μg-C
```

Estimated initial POM dissolution flux (d = 1.25 mm bin):

| Factor | Value |
|---|---|
| moisture: θ/(θ_P + θ) = 0.292/(0.10 + 0.292) | 0.745 |
| O₂: O_aq/(L_P + O_aq) = 0.00967/(0.00129 + 0.00967) | 0.882 |
| bacterial: B/(K_B_P + B) = 0.010/(0.001 + 0.010) | 0.908 |
| fungal: F_n/(K_F_P + F_n) ≈ (tiny F_n_min after dilution) | ~0.007 |
| microbial_factor = 0.5 × (0.908 + 0.007) | 0.457 |

```
J_P ≈ 0.5 × 1.0 × 0.457 × 0.745 × 0.882  ≈ 0.150 μg-C/mm²/day
R_P = J_P × 4π × r_0²  = 0.150 × 4π × 0.625² ≈ 0.74 μg-C/day
```

This concentrated DOC flux into node 1 drives a rapid concentration spike that feeds
exponential microbial growth. Without this source, microbes would not grow at all.

---

## 3. The Boom-Bust Mechanism

### Phase 1 (days 0–~5): Slow buildup

DOC accumulates near POM surface. Uptake is Monod-limited (`C_aq << K_B`).
Respiration grows slowly. Microbes are largely dormant in the outer domain.

### Phase 2 (days ~5–12): Exponential bloom

DOC near POM exceeds `K_B`. Bacteria enter near-saturated Monod kinetics.
Effective net growth rate:

```
g_net ≈ Y_B × r_B_max − μ_B ≈ 0.7 × 1.0 − 0.02 = 0.68 day⁻¹
doubling time ≈ ln(2) / 0.68 ≈ 1.0 day
```

Respiration rises exponentially. This is the observed initial increase.

### Phase 3 (~day 12): Crash

Microbial carbon demand exceeds POM supply. Sustainable biomass at steady state:

```
B_sustainable = R_P / (μ_B / Y_B) = 0.74 / (0.02 / 0.7) ≈ 25.9 μg-C per aggregate
```

Once biomass exceeds this, DOC depletes faster than POM dissolves. Because Monod
kinetics are highly nonlinear, uptake collapses nearly instantaneously when
`C_aq → 0`. Death rate `μ_B × B` continues, so biomass falls rapidly.

### Why the crash is so sharp

Three amplifiers make the transition abrupt rather than gradual:

1. **Fast DOC diffusion** (~86 mm²/day × tortuosity). A DOC deficit at the POM
   surface propagates across the 6.25 mm domain in hours — the entire domain
   crashes simultaneously.

2. **No MAOC buffer**: `κ_d_ref = 0.0001 day⁻¹` (10× below already-low default).
   MAOC holds ~0.97 μg-C/mm³ but releases it on a timescale of 1/κ_d_ref = 10,000
   days. For practical purposes it is inert on the ~12-day experiment horizon.

3. **Nonlinear Monod kinetics**: When `C_aq` drops from `K_B = 1e-4` to `1e-7`,
   uptake falls by 1000×. The transition from growth to starvation is essentially a
   step function.

---

## 4. Root Cause: Structural, Not Parameter-Sensitive

Multiple parameter combinations all crash because the crash is **structural** — a
consequence of the domain architecture, not parameter tuning.

| Quantity | Physical value | Model value | Ratio |
|---|---|---|---|
| DOC | 28.6 μg-C/mm³ | ~3.5e-8 μg-C/mm³ | 8 × 10⁸ |
| MAOC | 27.9 μg-C/mm³ | 0.97 μg-C/mm³ | 28.8 |
| POM | 204 μg-C | **204 μg-C** (undiluted) | 1 |

The ω = 28.8 dilution correctly represents that each model domain contains only 1/ω
of the background soil carbon. But POM is not diluted because it is physically
localized at the domain center. This creates an imbalance: one undiluted POM
particle in an essentially carbon-free matrix. Microbial growth feeds entirely off
POM → boom-bust is guaranteed regardless of `r_B_max`, `μ_B`, or `R_P_max`.

The experimental data likely shows a gradual inflection (not a crash) because:
- Incubation data integrates over many particles and a heterogeneous soil
- MAOC desorption provides DOC buffering over weeks
- Background SOC is never truly zero in practice

---

## 5. Suggested Levers (in order of likely impact)

### a. Increase `κ_d_ref` — MAOC buffering (highest priority)

Current: `0.0001`. Default: `0.01`. Try `0.001`–`0.01`.

Faster MAOC desorption means that when DOC crashes near the POM surface, MAOC
releases DOC to fill the gap. This converts the abrupt crash into a plateau. On a
10,000-day timescale (current setting) MAOC is completely inert; on a 100-day
timescale (κ_d = 0.01) it provides meaningful buffering within the 45-day experiment.

### b. Use `t_delay` — microbial equilibration window (already implemented)

Setting `t_delay = 3` (days) in `run_simulation` suppresses POM dissolution via the
sigmoid switch until microbes equilibrate. This prevents an initial DOC spike from
driving an early overshoot and shifts the peak later.

### c. Lower `B_S` — engage space limitation before the crash

Current: `B_S = 0.5 μg/mm³`. Estimated peak biomass ≈ 0.025 μg/mm³ << B_S, so
space limitation never engages. Try `B_S = 0.01`–`0.05`. This caps growth before
the crash and smooths the transition.

### d. Increase `K_B` — wider Monod half-saturation

`K_B = 1e-4 → 5e-4 μg/mm³`. A higher half-saturation means the collapse of C_aq
translates to a more gradual uptake reduction rather than a near-step-function cliff.

### e. Increase `κ_s_ref` / reduce retardation — allow faster C recycling

Current retardation factor: `θ + ρ_b × k_d_eq = 65.3`. This makes POM-derived DOC
65× less available to microbes than its volumetric concentration suggests. If DOC
sorbs strongly to minerals, it is also less bioavailable. Check whether `k_d_eq =
0.05` is appropriate for a silt loam (literature range: 0.02–0.2 mm³/μg).

---

## 6. Recommended Diagnostic Sequence

1. Run baseline (current parameters) — confirm crash at ~day 12.
2. Increase `κ_d_ref` to `0.005` — check if plateau replaces crash.
3. If yes: tune `r_B_max` and `R_P_max` to match peak timing and magnitude.
4. If no: add `t_delay = 3` — check shift in crash timing.
5. Compare cumulative CO₂ and instantaneous flux with de Gryze data per size class.
