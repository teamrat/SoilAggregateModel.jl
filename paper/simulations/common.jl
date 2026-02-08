# common.jl — Canonical parameter sets and utilities for Phase 2
#
# Usage: include("../common.jl") at the top of each theme's run_simulations.jl
#
# This file defines the EXACT parameter sets used for single-aggregate analyses.
# These values must not be changed after figures are generated without re-running
# all affected themes.

using SoilAggregateModel
using Printf

# ═══════════════════════════════════════════════════════════════
# STANDARD reference: medium-textured temperate soil
# Used for: radial profiles, carbon fate, stability lag, life cycles
# ═══════════════════════════════════════════════════════════════

STANDARD_BIO = BiologicalProperties()    # all defaults
STANDARD_SOIL = SoilProperties()         # all defaults

# Standard environmental conditions
STANDARD_T = 293.15        # 20°C
STANDARD_ψ = -33.0         # field capacity [kPa]
STANDARD_O2 = 0.27         # ambient O₂ [μg/mm³] (≈21% atmospheric)

# Standard POM
STANDARD_r_0 = 0.25        # 0.5 mm diameter POM [mm]
STANDARD_r_max = 5.0       # domain extends to 5 mm [mm]
STANDARD_n_grid = 100      # grid resolution

# ═══════════════════════════════════════════════════════════════
# Output time schedules
# ═══════════════════════════════════════════════════════════════

# 60-month run with monthly output (for carbon fate, MAOC, stability)
TIMES_60mo = Float64.([0; collect(1:60) .* 30.44])

# Dense early + sparse late (for capturing initial dynamics)
TIMES_DENSE = Float64.(vcat(
    0:0.5:7,                        # every 0.5 days for first week
    8:1:30,                         # daily for first month
    collect(2:6) .* 30.44,          # monthly months 2-6
    collect(9:3:60) .* 30.44        # every 3 months after that
))

# Snapshot times for radial profiles (months → days)
SNAPSHOT_MONTHS = [1, 6, 12, 24, 48]
SNAPSHOT_TIMES = Float64.(SNAPSHOT_MONTHS .* 30.44)

# ═══════════════════════════════════════════════════════════════
# Helper: write a NamedTuple or struct of vectors to CSV
# ═══════════════════════════════════════════════════════════════

"""
    write_csv(filename, data::NamedTuple)

Write a NamedTuple of equal-length vectors to CSV.
Compatible with R's readr::read_csv().

# Arguments
- `filename`: Output CSV path
- `data`: NamedTuple with equal-length vector fields

# Example
```julia
data = (t=[0.0, 1.0, 2.0], x=[1.0, 2.0, 3.0])
write_csv("output.csv", data)
```
"""
function write_csv(filename, data)
    names = keys(data)
    n = length(data[first(names)])

    open(filename, "w") do io
        # Header
        println(io, join(string.(names), ","))
        # Data rows
        for i in 1:n
            vals = [data[k][i] for k in names]
            println(io, join([@sprintf("%.10g", v) for v in vals], ","))
        end
    end
    println("  Wrote $filename ($n rows, $(length(names)) columns)")
end
