# SoilAggregateModel.jl

**Last Updated**: 2026-02-06
**Status**: ✅ Complete — 524/524 tests passing, carbon conservation to machine precision

---

## Overview

A mechanistic model of soil aggregate formation and carbon cycling. The model solves 9 coupled equations (5 diffusion PDEs + 3 local ODEs + 1 scalar ODE) on a fixed spherical grid, using Strang-split Crank–Nicolson diffusion and pointwise reactions with adaptive time-stepping.

**Performance target**: 400 aggregates × 1500 grid points × 10 years, with zero heap allocations in the hot loop.

**Key features**:
- Custom Crank-Nicolson solver with Thomas algorithm (O(n) per species)
- Strang splitting (second-order operator splitting)
- Adaptive timestep control
- Exact carbon conservation (< 10⁻¹² relative error)
- Temperature-dependent rates via Arrhenius kinetics

---

## Quick Start

### Installation

```bash
git clone https://github.com/teamrat/aggregate.git
cd aggregate
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run Tests

```bash
julia --project=. test/runtests.jl
```

### Basic Usage

```julia
using SoilAggregateModel

# Set up parameters
bio = BiologicalProperties()
soil = SoilProperties()

# Define environmental drivers (constant or time-varying)
T(t) = 293.15  # Temperature [K]
ψ(t) = -10.0   # Water potential [kPa]
O2(t) = 0.3    # Ambient O₂ [μg/mm³]

# Run 30-day simulation
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))

# Access results
for output in result.outputs
    println("t = ", output.t, " days: POM = ", output.state.P, " μg-C")
end
```

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
aggregate/
├── src/                      # Source code (SoilAggregateModel module)
│   ├── SoilAggregateModel.jl # Main module file
│   ├── types.jl              # AggregateState, Workspace, TemperatureCache
│   ├── parameters.jl         # BiologicalProperties, SoilProperties
│   ├── constants.jl          # Physical constants (R_GAS)
│   ├── environment.jl        # EnvironmentalDrivers{FT,Fψ,FO}
│   │
│   ├── temperature/          # Temperature dependencies
│   │   ├── arrhenius.jl      # Arrhenius factor
│   │   ├── diffusion_pure.jl # Stokes-Einstein, Han-Bartels, Chapman-Enskog
│   │   └── henry.jl          # Henry's law K_H(T)
│   │
│   ├── physics/              # Soil physics
│   │   ├── water_retention.jl    # Modified van Genuchten θ(ψ, E, F_i)
│   │   └── effective_diffusion.jl # Tortuosity-limited diffusion
│   │
│   ├── biology/              # Biological processes
│   │   ├── bacteria.jl       # Bacterial growth, mortality, recycling
│   │   ├── fungi.jl          # Fungal growth, transitions, translocation
│   │   ├── eps.jl            # EPS production and recycling
│   │   └── maoc.jl           # MAOC sorption (Langmuir-Freundlich)
│   │
│   ├── carbon/               # Carbon cycling
│   │   └── pom_dissolution.jl # Enzymatic POM breakdown
│   │
│   ├── solver/               # Numerical solver
│   │   ├── tridiagonal.jl        # Thomas algorithm
│   │   ├── crank_nicolson.jl     # CN diffusion step
│   │   ├── finite_volumes.jl     # Conservation weights
│   │   ├── reactions.jl          # Source/sink terms
│   │   ├── diffusion_step.jl     # 5-species diffusion
│   │   ├── reaction_step.jl      # Local reactions + POM + CO₂
│   │   ├── workspace_updates.jl  # Update caches once per timestep
│   │   └── timestepper.jl        # Strang splitting + adaptive Δt
│   │
│   └── api.jl                # Public API (run_aggregate)
│
├── test/                     # Test suite (524 tests)
│   ├── runtests.jl           # Main test runner
│   ├── test_types.jl         # Data structures
│   ├── test_parameters.jl    # BiologicalProperties, SoilProperties
│   ├── test_environment.jl   # EnvironmentalDrivers
│   ├── test_temperature.jl   # Arrhenius, diffusion, Henry's law
│   ├── test_tridiagonal.jl   # Thomas algorithm
│   ├── test_crank_nicolson.jl # CN solver + BCs
│   ├── test_physics.jl       # Water retention, effective diffusion
│   ├── test_biology.jl       # Bacteria, fungi, EPS, MAOC
│   ├── test_pom.jl           # POM dissolution
│   ├── test_reactions.jl     # Source/sink computation
│   ├── test_timestepper.jl   # Time integration, adaptive Δt
│   └── test_api.jl           # User-facing API
│
├── docs/                     # Documentation
│   ├── ARCHITECTURE.md       # Implementation architecture (authoritative)
│   ├── GUIDE.md              # Theory, usage, developer guide
│   ├── REFERENCE.md          # Variables, parameters, functions
│   ├── PARAMETER_BUG_REPORT.md
│   ├── REPO_REORGANIZATION.md
│   └── archive/              # Historical documents
│
├── scripts/                  # Diagnostic scripts
│   └── diagnostics_30day.jl  # 30-day simulation with detailed diagnostics
│
├── paper/                    # Manuscript materials
│   ├── simulations/          # Figure generation scripts
│   │   ├── common.jl         # Shared parameters
│   │   └── README.md
│   ├── figures/              # Output figures (syncs to Overleaf)
│   └── data/                 # Simulation outputs
│
├── CLAUDE.md                 # Instructions for Claude Code
├── Project.toml              # Julia package manifest
└── README.md                 # This file
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

- **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** — Implementation architecture (authoritative for solver design)
- **[docs/GUIDE.md](docs/GUIDE.md)** — Theory, usage, developer guide
- **[docs/REFERENCE.md](docs/REFERENCE.md)** — Variables, parameters, functions (quick lookup)
- **[docs/PARAMETER_BUG_REPORT.md](docs/PARAMETER_BUG_REPORT.md)** — Bug fixes and parameter corrections
- **Manuscript** — Authoritative for all physics equations (see docs/archive/)

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
- **2026-01**: Julia port attempt (abandoned — poorly implemented).
- **2026-02**: **First working Julia version**. Custom Strang splitting solver, zero-allocation design. Added MAOC pool, Arrhenius temperature framework, enzymatic POM dissolution.

**Development** (2026-02):
- Complete rewrite with custom Crank-Nicolson solver
- Verification against manuscript equations
- Parameter corrections (ρ_b, ρ_POM, λ, K_Y units)
- 524 tests, carbon conservation to machine precision

**Theoretical basis**: Ghezzehei, T.A. et al. (2026), manuscript in preparation.

---

## Contact

**Teamrat A. Ghezzehei**  
University of California, Merced  
taghezzehei@ucmerced.edu

---

## License

[To be determined]
