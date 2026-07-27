# Diagnostic and Validation Scripts

This folder contains standalone scripts for performance benchmarking, conservation validation, diagnostic testing, and example simulations for the SoilAggregateModel.

## Example Simulation

### `example_1year_simulation.jl`
**Purpose**: Complete example showing how to run a 1-year simulation and export results to CSV
**What it does**:
- Runs 365-day simulation with daily outputs
- Extracts time-series for all carbon pools using `integrated_pools()` function
- Calculates CO2 flux using `co2_flux()` function
- Exports 3 CSV files for analysis in R/Python
- Demonstrates correct usage of package postprocessing functions
- Includes R code examples for plotting

**Usage**: `julia scripts/example_1year_simulation.jl`

**Output files** (saved to `output/`):
- `timeseries_pools.csv` - All carbon pools over time (C, B, F_i, F_n, F_m, E, M, POM, CO2)
- `timeseries_fluxes.csv` - CO2 flux [μg-C/day]
- `spatial_profile_day365.csv` - Final spatial distribution at day 365

**Key features**:
- Uses existing package functions (`integrated_pools`, `co2_flux`, `carbon_balance_table`)
- Correct terminology: F_insulated, F_noninsulated, F_mobile
- Correct units: O2 in μg/mm³ (not molar fraction), T in K, ψ in kPa

---

## Production Tools

### `performance_diagnostic.jl`
**Purpose**: Benchmark solver performance for production planning
**What it does**:
- Tests 1-year runs at different dt_max values (0.1, 0.5, 1.0)
- Tests 5-year baseline run for scaling analysis
- Reports step counts, wall time, and conservation errors
- Determines if solver is stiffness-limited or dt_max-limited

**Usage**: `julia --project=. scripts/performance_diagnostic.jl`

**Key output**: Per-aggregate runtime projections for population-scale simulations (400 aggregates × 10 years)

---

## Conservation Validation

### `fixed_dt_365day_test.jl`
**Purpose**: Verify conservation at machine precision over 365 days
**What it does**:
- Runs 365-day simulation with fixed dt=0.001
- Tests conservation without adaptive timestepper effects
- Reports final conservation error (target: < 10^-12)

**Usage**: `julia --project=. scripts/fixed_dt_365day_test.jl`

**Expected**: ~10^-14 relative error (machine precision)

### `adaptive_timestepper_conservation_test.jl`
**Purpose**: Verify conservation with adaptive timestepping
**What it does**:
- Runs 365-day simulation with adaptive dt (dt_max=0.1)
- Tests accept-then-adjust strategy
- Reports conservation error and step count

**Usage**: `julia --project=. scripts/adaptive_timestepper_conservation_test.jl`

**Expected**: ~10^-11 to 10^-14 relative error

### `diffusion_neumann_test.jl`
**Purpose**: Verify Neumann boundary condition at POM surface
**What it does**:
- Tests single diffusion step with Neumann BC
- Compares numerical flux to expected flux from POM dissolution
- Verifies Crank-Nicolson implementation

**Usage**: `julia --project=. scripts/diffusion_neumann_test.jl`

**Expected**: Flux error < 10^-14

---

## Bug Discovery Diagnostics

### `high_biomass_strang_test.jl`
**Purpose**: Test Strang splitting conservation at high biomass
**What it does**:
- Loads day-60 state (B ≈ 4.5 μg/mm³, high clipping events)
- Runs 100 D-R-D Strang steps
- Checks per-step conservation residual
- **This script found the clipping accounting bug (Feb 2026)**

**Usage**: `julia --project=. scripts/high_biomass_strang_test.jl`

**Historical**: Before fix: 2.5×10^-5 μg-C leak per step. After fix: 1.1×10^-12 (machine precision)

### `identity_check_day60.jl`
**Purpose**: Verify conservation identity algebraically
**What it does**:
- Loads day-60 state
- Computes conservation identity: S_C + S_B + S_Fi + S_Fn + S_Fm + S_E + S_M + Resp_total
- Checks if source term algebra is exact
- **Proved the bug was NOT in biology formulas, isolating it to clipping accounting**

**Usage**: `julia --project=. scripts/identity_check_day60.jl`

**Expected**: Identity error < 10^-20 (proves algebra is correct)

---

## Notes

- All scripts assume you're in the project root directory
- All scripts use `using SoilAggregateModel` (package must be instantiated)
- Scripts contain `## Setup` markers for VS Code Julia extension cell execution
- Timings are hardware-dependent (scripts run on M3 Pro/Ultra)

---

Last updated: 2026-02-10
