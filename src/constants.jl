"""
    constants.jl

Shared physical constants used throughout the model.

All constants are in SI-derived units compatible with the model's unit system:
μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol

# Module Structure
This file should be included EXACTLY ONCE at the top of the main module file or test file,
before any other implementation files. Individual implementation files assume these constants
are already defined and should NOT include this file directly.

Correct pattern:
```julia
module SoilAggregateModel
    include("constants.jl")      # Include once here
    include("temperature/arrhenius.jl")  # Uses R_GAS
    include("temperature/henry.jl")      # Uses R_GAS
    # ...
end
```

Incorrect pattern:
```julia
# In arrhenius.jl
include("constants.jl")  # ❌ Don't do this in implementation files
```
"""

# Gas constant [J/(mol·K)]
const R_GAS = 8.314

# Standard atmosphere [Pa] and O₂ molar mass [kg/mol]. Both were bare literals,
# with `101000.0` disagreeing with the `1 atm = 101325` used in henry.jl by
# 0.32 % — which is why ambient gas-phase O₂ existed as five different numbers.
# 101325 is the definition. Adopting it moved `O2_CONST` in degryze_config.jl
# from 0.273808 to 0.274689, so runs before 2026-07-30 differ in the last
# digits. The scripts that still carried the old literals were archived on
# 2026-07-30 (`paper/_archive/degryze_tooling_20260730/`), so this is now the
# only definition anywhere that runs.
const P_ATM = 101325.0
const M_O2  = 0.032

# Mineral MAOC capacity coefficient [g-C per g of clay+silt] — DIMENSIONLESS.
#
# Georgiou et al. (2022): MOC_max = 86 ± 9 mg-C/g-mineral for high-activity
# minerals and 48 ± 6 for low-activity, converted here from mg/g to g/g.
#
# Dimensionless is not cosmetic. `M_max = k_ma · f_clay_silt · ρ_b` closes in
# µg-C/mm³ ONLY if k_ma is a mass fraction; with the "µg-C per gram" unit the
# package carried until 2026-07-29 the identity is short by 1e6. See
# REFERENCE.md §26 erratum 12.
# Radial nodes, default for BOTH entry points. They disagreed until 2026-07-30 —
# `run_aggregate` defaulted to 50 and `run_aggregate_stiff` to 200, so the same
# call resolved to a different grid depending on which solver ran. 200 is the
# production value; every call in `paper/` and `test/` passes `n_grid` explicitly,
# so this default governs casual use only.
const N_GRID_DEFAULT = 200

const K_MA_HIGH_ACTIVITY = 0.086
const K_MA_LOW_ACTIVITY  = 0.048

# Steepness of the sigmoid threshold `sigmoid_threshold(x, x_min, steepness)`,
# expressed as β·x_min so it is scale-free in x.
#
# 50 is used by the three viability thresholds h_B, h_E and h_Fi: the switch
# spans roughly x_min ± 0.1·x_min. It was a bare literal at all three sites.
# 10 is used by the POM activation delay, whose width was written as
# `0.1·t_delay`, i.e. β = 1/(0.1·t_delay) = 10/t_delay. Neither number has a
# citation; both are working assumptions inherited from the MATLAB code.
const SIGMOID_STEEPNESS    = 50.0
const POM_DELAY_STEEPNESS  = 10.0
