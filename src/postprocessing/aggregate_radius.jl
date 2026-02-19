# aggregate_radius.jl
# Aggregate stability radius from binding agent profiles

"""
    compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)

Stable aggregate radius [mm] from oscillatory wall shear stress criterion.

The aggregate is the region where biologically produced binding agents
generate sufficient cohesive strength to resist hydrodynamic disaggregation
during wet sieving. The criterion compares the local cohesive strength
τ_c(r) = k_F·F_i(r) + k_E·E(r) against the wall shear stress τ_w
from the oscillatory Stokes boundary layer (Batchelor, 1967, §4.3, §5.14).

Dividing through by k_F with k_E = k_F/2:

    F_i(r) + ½·E(r) ≥ G_c = √2·μ·v_s / (k_F·δ_s)

where:
- μ: dynamic viscosity of water [Pa·s]
- v_s = π·f_s·L_s: maximum sieving velocity [mm/s]
- δ_s = √(2·ν_w/Ω): Stokes boundary layer thickness [mm]
- Ω = 2π·f_s: angular frequency [rad/s]
- k_F: specific binding strength [Pa/(μg/mm³)]

G_c is a concentration threshold [μg/mm³] independent of position r.

# Arguments
- `record::OutputRecord`: State snapshot at a specific time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Biological and soil parameters

# Returns
- `r_agg::Float64`: Stable aggregate radius [mm]
  - Returns 0.0 if criterion is not met at any node
  - Returns r_grid[i] of the outermost node satisfying the criterion

# Notes
- Scans ALL nodes and takes the outermost passing one (no early break),
  because binding agent profiles can be non-monotonic.
- The flat-wall approximation is accurate for aggregates with diameter > ~3 mm
  (δ_s ≈ 0.75 mm). For smaller aggregates, curvature effects cause the model
  to overestimate stability. This bias is partially absorbed into calibrated k_F.

# References
- Batchelor, G. K. (1967). An Introduction to Fluid Dynamics. Cambridge
  University Press. §4.3 (oscillating plate solution), §5.14 (extension to
  oscillating rigid bodies).
"""
function compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)
    # Standard wet sieving parameters (Eijkelkamp apparatus)
    L_s = 13.0          # stroke length [mm]
    f_s = 34.0 / 60.0   # frequency [Hz] (34 oscillations/min)
    μ = 1.002e-3         # dynamic viscosity of water at 20°C [Pa·s]
    ν_w = 1.004          # kinematic viscosity of water at 20°C [mm²/s] = 1.004e-6 m²/s

    # Maximum sieving velocity [mm/s]
    v_s = π * f_s * L_s  # ≈ 23.1 mm/s

    # Angular frequency [rad/s]
    Ω = 2π * f_s  # ≈ 3.56 rad/s

    # Stokes boundary layer thickness [mm]
    δ_s = sqrt(2.0 * ν_w / Ω)  # ≈ 0.75 mm

    # Wall shear stress [Pa] — peak amplitude from exact flat-plate solution
    # τ_w = √2 · μ · v_s / δ_s
    # (the √2 arises from the phase difference between velocity and stress
    #  at the wall in the oscillatory Stokes solution; Batchelor §4.3)
    #
    # NOTE: v_s is in mm/s but μ is in Pa·s = kg/(m·s).
    # Convert v_s to m/s and δ_s to m for consistent Pa units:
    τ_w = sqrt(2.0) * μ * (v_s * 1e-3) / (δ_s * 1e-3)  # [Pa]

    # Specific binding strength [Pa/(μg/mm³)]
    k_F = params.soil.k_F

    # Critical threshold [μg/mm³]
    G_c = τ_w / k_F

    # State variables
    F_i = record.state.F_i
    E = record.state.E

    # Compute binding concentration at all nodes
    binding = F_i .+ 0.5 .* E

    # Find outermost node where binding >= G_c
    r_agg = 0.0
    for i in grid.n:-1:1  # scan from outer to inner
        if binding[i] >= G_c
            r_agg = grid.r_grid[i]
            break  # found the outermost, done
        end
    end

    return r_agg
end
