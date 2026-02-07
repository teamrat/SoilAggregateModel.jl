# common.jl — Shared parameter sets and utilities for paper figures
#
# Usage: include("common.jl") at the top of each figure script
#
# This file defines the EXACT parameter sets used in the published paper.
# Do not modify after publication without documenting the change.

using SoilAggregateModel
using CairoMakie  # or whatever plotting package you choose

# === Canonical parameter set (medium-textured soil, temperate climate) ===
# These override the package defaults where needed for the paper simulations.

const PAPER_BIO = BiologicalProperties(
    # Override defaults as needed for paper
    # Example: r_B_max = 5.0, ...
)

const PAPER_SOIL = SoilProperties(
    # Override defaults as needed for paper
    # Example: k_d_eq = 0.05, ...
)

# === Standard environmental conditions ===
const T_STANDARD = 293.15    # 20°C
const ψ_STANDARD = -33.0     # Field capacity [kPa]

# === Figure output ===
const FIGURE_DIR = joinpath(@__DIR__, "..", "figures")
mkpath(FIGURE_DIR)

"""
    savepaperfig(fig, name)

Save figure to paper/figures/ as both PDF and PNG.
PDF goes to Overleaf via Dropbox sync.
"""
function savepaperfig(fig, name)
    save(joinpath(FIGURE_DIR, name * ".pdf"), fig)
    save(joinpath(FIGURE_DIR, name * ".png"), fig, px_per_unit=3)
    println("Saved: $name.pdf and $name.png")
end
