# Soil Aggregate Model — Julia v2

**Last Updated**: 2026-02-06  
**Status**: Kernel modules verified against manuscript; solver assembly in progress

---

## Overview

A mechanistic model of soil aggregate formation and carbon cycling. The model solves 9 coupled equations (5 diffusion PDEs + 3 local ODEs + 1 scalar ODE) on a fixed spherical grid, using Strang-split Crank–Nicolson diffusion and pointwise reactions with adaptive time-stepping.

**Target**: 400 aggregates × n = 1500 × 10 years, with zero heap allocations in the hot loop.

---

## State Variables

| Variable | Symbol | Transport | Description |
|----------|--------|-----------|-------------|
| Dissolved organic carbon | C | Diffusion (retarded) | Total mobile carbon (aqueous + fast-sorbed) |
| Bacteria | B | Diffusion | Bacterial biomass |
| Non-insulated fungi | F_n | Diffusion | Hyphal tips, extends via diffusion |
| Mobile fungi | F_m | Diffusion | Translocation network (no tortuosity) |
| Oxygen | O | Diffusion (aq + gas) | Dissolved + gas-phase oxygen |
| Insulated fungi | F_i | Immobile (local ODE) | EPS-coated, binds particles |
| EPS | E | Immobile (local ODE) | Extracellular polymeric substances |
| MAOC | M | Immobile (local ODE) | Mineral-associated organic carbon |
| POM | P | Scalar ODE | Particulate organic matter, well-mixed core |

Total state: 8n + 1 (+ cumulative CO₂ diagnostic).

---

## Project Structure

```
src/
├── types.jl                  # AggregateState, Workspace, TemperatureCache
├── parameters.jl             # BiologicalProperties, SoilProperties + defaults
├── environment.jl            # EnvironmentalDrivers{FT,Fψ,FO}
│
├── temperature/
│   ├── arrhenius.jl          # arrhenius(Ea, T, T_ref)
│   ├── diffusion_pure.jl     # Stokes-Einstein+VFT, Han-Bartels, Chapman-Enskog
│   └── henry.jl              # K_H(T) via van't Hoff (sign-corrected)
│
├── physics/
│   ├── water_retention.jl    # θ(ψ, E, F_i) — modified van Genuchten
│   └── effective_diffusion.jl # D_C, D_B, D_Fn, D_Fm, D_O
│
├── biology/
│   ├── bacteria.jl           # R_B, R_Bb, h_B, Γ_B, Γ_E, Resp_B
│   ├── fungi.jl              # R_F, Π, transitions, Resp_F, Resp_F_conv, h_Fi
│   ├── eps.jl                # R_rec_E (uses C_aq), h_E
│   └── maoc.jl               # J_M with softplus, M_eq Langmuir-Freundlich
│
├── carbon/
│   └── pom_dissolution.jl    # J_P (flux density) and R_P = 4πr₀²·J_P
│
├── solver/
│   ├── tridiagonal.jl        # Thomas algorithm (in-place)
│   ├── crank_nicolson.jl     # CN half-step for one species
│   ├── diffusion_step.jl     # 5 tridiagonal solves
│   ├── reactions.jl          # All source/sink terms at one node
│   ├── reaction_step.jl      # Loop over nodes + POM + CO₂
│   └── timestepper.jl        # Strang splitting + adaptive Δt
│
└── api.jl                    # run_aggregate() entry point

test/
├── runtests.jl
├── test_types.jl
├── test_parameters.jl
├── test_temperature.jl
├── test_environment.jl
├── test_tridiagonal.jl
├── test_physics.jl
└── test_biology.jl
```

---

## Numerical Method

**Strang splitting** (second-order):

```
Each Δt:
  1. Diffusion half-step (Δt/2) — Crank–Nicolson, 5 tridiagonal solves, O(n) each
  2. Reaction full-step (Δt)   — pointwise at each node + POM scalar + CO₂
  3. Diffusion half-step (Δt/2) — same as step 1
```

Adaptive Δt: halve if max relative change > 10%, double if < 1%.

---

## Key Design Decisions (from audit)

Decisions made during the manuscript→code verification (2026-02-05/06):

1. **MAOC coupling in S_C**: `S_C = ... − J_M` (not `−J_M·(θ+ρ_b·k_d)/k_d`). The manuscript formula was wrong — the factor is unity because J_M is already per mm³ bulk soil.

2. **Henry's law sign**: `K_H(T) = K_H_ref · exp[ΔH_sol/R · (1/T − 1/T_ref)]` — no leading negative sign. The manuscript formula had an erroneous extra negative.

3. **EPS degradation uses C_aq**: `R_rec_E = μ_E · K_E/(K_E + C_aq) · E · h_E`, not raw state variable C. Scavenging responds to dissolved concentration.

4. **Resp_F^conv uses abs()**: Prevents nonphysical negative respiration during mobilization-dominated regimes.

5. **ν_B, ν_F are positive**: Water stress = exp(ν·ψ) where ψ < 0, so positive ν gives suppression < 1.

6. **Temperature**: Arrhenius with 6 distinct activation energies (not Q10). See REFERENCE.md §3.7.

7. **Tortuosity**: τ_aq = θ²/θ_s^(2/3) is a deliberate simplification, not standard Millington-Quirk.

---

## Documentation

- **[GUIDE.md](GUIDE.md)** — Theory, usage, developer guide
- **[REFERENCE.md](REFERENCE.md)** — Variables, parameters, functions (quick lookup)
- **[ARCHITECTURE_CLAUDE_CODE.md](ARCHITECTURE_CLAUDE_CODE.md)** — Implementation architecture (authoritative for solver design)
- **Manuscript** (`manuscript_updated.tex`) — Authoritative for all physics equations

---

## Units

All internal computations use: **μg, mm, days, kPa, K, J/mol**. No unit conversions inside kernels.

| Quantity | Unit | Equivalence |
|----------|------|-------------|
| Concentration | μg/mm³ | = kg/m³ |
| Diffusion | mm²/day | |
| Rates | day⁻¹ | |
| Water potential | kPa | negative in unsaturated soil |

---

## Background

- **2016–2025**: MATLAB implementation (pdepe). Established core physics.
- **2026-01**: Julia v1 (DifferentialEquations.jl/CVODE monolithic solver).
- **2026-02**: Julia v2 (Strang splitting, zero-allocation design). Added MAOC pool, Arrhenius temperature framework, enzymatic POM dissolution.

**Theoretical basis**: Ghezzehei, T.A. et al. (2026), manuscript in preparation.

---

## Contact

**Teamrat A. Ghezzehei**  
University of California, Merced  
taghezzehei@ucmerced.edu

---

## License

[To be determined]
