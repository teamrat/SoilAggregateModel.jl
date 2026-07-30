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
# digits. `optimize_soil3.jl` and the four `diagnostics/` scripts still carry
# the old literals and so are now 0.32 % off from the production config.
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
const K_MA_HIGH_ACTIVITY = 0.086
const K_MA_LOW_ACTIVITY  = 0.048
