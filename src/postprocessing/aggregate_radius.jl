# aggregate_radius.jl
# Aggregate stability radius from binding agent profiles

"""
    compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)

Stable aggregate radius [mm] from manuscript Section 2.4.7 (Eq. in line 753).

Scans outward from r_0. The stable radius is the largest r where:
    (F_i(r) + ½·E(r)) · (r + χ·a_p) ≥ G_c

where G_c = 3μv / (2k_F) is the critical binding strength threshold.

# Arguments
- `record::OutputRecord`: State snapshot at a specific time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Biological and soil parameters

# Returns
- `r_agg::Float64`: Stable aggregate radius [mm]
  - Returns 0.0 if criterion is not met at any node (including r_0)
  - Returns r_grid[i] of the outermost node satisfying the criterion

# Physical basis
From Stokes drag force balance during wet sieving:
- Hydrodynamic stress must not exceed binding strength
- Binding strength = (F_i + 0.5·E) × effective radius
- Effective radius r_eff = r + χ·a_p accounts for particle adhesion
- G_c threshold from sieving velocity v = π·f·L

# Units
- μ: dynamic viscosity [Pa·s = kg/(m·s)]
- v: sieving velocity [mm/s]
- k_F: specific binding strength [kPa/(μg/mm³)]
- G_c: [Pa·mm] = [kg/s²] = [kg·mm/s²]
- All consistent with model's μg/mm³/days unit system

# Manuscript reference
Section 2.4.7, Aggregate Stability Mechanics
"""
function compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)
    # Sieve shaker parameters (standard wet sieving)
    L = 30.0        # stroke length [mm]
    f = 50.0 / 60.0  # frequency [Hz] = [1/s] (50 oscillations/min)
    μ = 1.002e-3     # dynamic viscosity of water at 20°C [Pa·s = kg/(m·s)]

    # Maximum sieving velocity [mm/s]
    v = π * f * L  # ≈ 78.5 mm/s

    # Soil parameters
    k_F = params.soil.k_F      # [kPa/(μg/mm³)]
    χ = params.soil.χ          # [mm]
    a_p = params.soil.a_p      # [mm]

    # Critical binding strength G_c
    # Units check:
    # LHS of criterion: (F_i + 0.5·E) × r_eff = [μg/mm³] × [mm] = [μg/mm²]
    # So G_c must have units [μg/mm²]
    #
    # G_c = 3μv/(2k_F)
    # μ: [Pa·s] = [kg/(m·s)]
    # v: [mm/s]
    # k_F: [kPa/(μg/mm³)] = [1000 Pa/(μg/mm³)]
    #
    # G_c = 3 × [kg/(m·s)] × [mm/s] / (2 × 1000 [Pa/(μg/mm³)])
    #     = 3μv / (2000 k_F) × [(kg·mm/m)/s²] / [(kg/(m·s²))/(μg/mm³)]
    #     = 3μv / (2000 k_F) × [μg/mm²]
    #
    # So: G_c = 3μv / (2000 k_F) where k_F is in kPa units
    G_c = (3.0 * μ * v) / (2000.0 * k_F)  # [μg/mm²]

    # State variables
    F_i = record.state.F_i
    E = record.state.E
    r = grid.r_grid

    # Scan outward from r_0 to find outermost node satisfying criterion
    r_agg = 0.0

    for i in 1:grid.n
        # Effective radius
        r_eff = r[i] + χ * a_p

        # Binding strength
        binding_strength = (F_i[i] + 0.5 * E[i]) * r_eff

        # Check criterion
        if binding_strength >= G_c
            r_agg = r[i]
        else
            # Once criterion fails, stop (monotonic scan)
            break
        end
    end

    return r_agg
end
