# finite_volumes.jl
# Geometrically exact finite-volume cell volumes for spherical geometry
#
# NOTE: These volumes are NOT used for discrete conservation in the current implementation.
# The spherical Laplacian stencil requires conservation weights W_i = 4πr_i²h to ensure
# exact telescoping (W_i/(r_i²h²) = constant). The geometrically exact shells computed
# here do not have this property. This file is retained for reference and potential
# future use (e.g., post-processing, visualization).

"""
    compute_cell_volumes(r_grid::Vector{Float64})

Compute geometrically exact finite-volume cell volumes for spherical geometry.

# Arguments
- `r_grid::Vector{Float64}`: Radial grid points [mm] (cell centers)

# Returns
- `volume::Vector{Float64}`: Cell volumes [mm³] such that ∑(u[i] * volume[i])
  equals the exact spherical integral

# Algorithm
Each cell i extends from r_minus to r_plus:
- Node 1 (inner): r_minus = r_grid[1] (POM surface), r_plus = (r_grid[1] + r_grid[2])/2
- Interior nodes: r_minus = (r_grid[i-1] + r_grid[i])/2, r_plus = (r_grid[i] + r_grid[i+1])/2
- Node n (outer): r_minus = (r_grid[n-1] + r_grid[n])/2, r_plus = r_grid[n]

Volume = (4π/3)(r_plus³ - r_minus³)

# Notes
- These volumes ensure discrete conservation when used consistently
- The inner boundary is at r_grid[1], not r_grid[1] - h/2
- Summing with these weights gives the correct total volume

# Manuscript reference
Architecture §3: Finite volume discretization for spherical geometry
"""
function compute_cell_volumes(r_grid::Vector{Float64})
    n = length(r_grid)
    volume = Vector{Float64}(undef, n)

    for i in 1:n
        # Cell boundaries (interfaces)
        r_minus = if i == 1
            r_grid[1]  # Inner boundary at POM surface
        else
            (r_grid[i-1] + r_grid[i]) / 2.0  # Midpoint
        end

        r_plus = if i == n
            r_grid[n]  # Outer boundary at aggregate surface
        else
            (r_grid[i] + r_grid[i+1]) / 2.0  # Midpoint
        end

        # Spherical shell volume: (4π/3)(r_plus³ - r_minus³)
        volume[i] = (4.0 * π / 3.0) * (r_plus^3 - r_minus^3)
    end

    return volume
end
