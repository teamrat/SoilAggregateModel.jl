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
