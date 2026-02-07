"""
    eps.jl

Extracellular polymeric substance (EPS) degradation.

EPS is produced by bacteria (Γ_E) and degraded enzymatically when soluble
carbon is scarce. Degradation is inhibited when C_aq is abundant.

# Dependencies
Requires: BiologicalProperties

# Critical C usage
- R_rec_E uses C_aq (MANUSCRIPT_CHANGES #3), NOT raw C
- Compute C_aq = C / (θ + ρ_b·k_d) ONCE per node, pass to function
"""

"""
    h_E(E::Real, E_min::Real)

Sigmoid threshold function for EPS viability.

# Formula (manuscript Eq. 403)
    h_E = exp(β_E·E) / [exp(β_E·E) + exp(β_E·E_min)]

where β_E = 50/E_min.

# Arguments
- `E`: EPS concentration [μg/mm³]
- `E_min`: Minimum EPS for sigmoid [μg/mm³]

# Returns
- Threshold factor [-] (0 ≤ h_E ≤ 1)

# Notes
- As E → 0: h_E → 0 (shuts off degradation)
- As E → ∞: h_E → 1 (full degradation rate)
- Applied to R_rec_E (EPS degradation)
- Analogous to h_B and h_Fi
"""
function h_E(E::Real, E_min::Real)
    β_E = 50.0 / E_min
    # Numerically stable sigmoid: 1 / (1 + exp(-β_E * (E - E_min)))
    # Avoids overflow when E >> E_min
    1.0 / (1.0 + exp(-β_E * (E - E_min)))
end

"""
    R_rec_E(μ_E_max_T::Real, C_aq::Real, K_E::Real, E::Real, E_min::Real)

EPS degradation (recycling) rate with substrate inhibition.

# Formula (manuscript Eq. 397, CORRECTED)
    R_rec_E = μ_E^max(T) · K_E/(K_E + C_aq) · E · h_E

# Arguments
- `μ_E_max_T`: Max degradation rate at current T [1/day]
- `C_aq`: Aqueous DOC concentration [μg/mm³ water] ← CRITICAL: uses C_aq, not C
- `K_E`: Substrate inhibition concentration [μg/mm³]
- `E`: EPS concentration [μg/mm³]
- `E_min`: Minimum EPS for sigmoid [μg/mm³]

# Returns
- Degradation rate [μg-C/mm³/day]

# Notes
- **CRITICAL**: Uses C_aq (MANUSCRIPT_CHANGES #3), NOT raw C state variable
- Substrate inhibition: degradation increases as C_aq decreases
- As C_aq → 0: R_rec_E → μ_E_max_T · E · h_E (maximum scavenging)
- As C_aq → ∞: R_rec_E → 0 (EPS preserved when DOC abundant)
- μ_E_max_T already includes Arrhenius factor
- Returns carbon to DOC pool (contributes to R_rec in S_C)
"""
function R_rec_E(μ_E_max_T::Real, C_aq::Real, K_E::Real, E::Real, E_min::Real)
    inhibition = K_E / (K_E + C_aq)
    h_E_val = h_E(E, E_min)

    μ_E_max_T * inhibition * E * h_E_val
end
