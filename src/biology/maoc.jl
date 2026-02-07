"""
    maoc.jl

Mineral-associated organic carbon (MAOC) sorption/desorption kinetics.

Implements smooth switching between rate-limited sorption and desorption toward
local equilibrium determined by Langmuir-Freundlich isotherm.

# Dependencies
Requires: SoilProperties

# Critical C usage
- M_eq uses C_eq (equilibrium sorbed concentration), NOT C or C_aq
- Compute C_eq = k_d · C_aq = k_d·C / (θ + ρ_b·k_d) ONCE per node
- S_C coupling: contribution is -J_M · (θ + ρ_b·k_d) / k_d, NOT just -J_M

# Softplus regularization
- Replaces max(0, x) with smooth C∞ approximation
- Prevents discontinuous derivatives at M = M_eq
- Uses numerically stable form to avoid overflow
"""

"""
    softplus(x::Real, ε::Real)

Smooth approximation to max(0, x) with C∞ continuity.

# Formula (manuscript Eq. 243)
    φ_ε(x) = ε · ln(1 + exp(x/ε))

Numerically stable implementation:
    x > 0:  x + ε·ln(1 + exp(-x/ε))
    x ≤ 0:  ε·ln(1 + exp(x/ε))

# Arguments
- `x`: Input value (can be any real)
- `ε`: Smoothing width [same units as x]

# Returns
- Smoothed value (always ≥ 0)

# Notes
- As ε → 0: softplus(x, ε) → max(0, x)
- For |x| >> ε: softplus(x, ε) ≈ max(0, x) to within ε
- Typical ε = 0.01 μg/mm³ for MAOC
- Numerically stable form avoids overflow for large |x|
"""
function softplus(x::Real, ε::Real)
    if x > 0.0
        # Stable for large positive x
        x + ε * log(1.0 + exp(-x / ε))
    else
        # Stable for large negative x
        ε * log(1.0 + exp(x / ε))
    end
end

"""
    M_eq_langmuir_freundlich(C_eq::Real, M_max::Real, k_L::Real, n_LF::Real)

Equilibrium MAOC concentration from Langmuir-Freundlich isotherm.

# Formula (manuscript Eq. 264)
    M_eq = M_max · (k_L·C_eq)^n_LF / [1 + (k_L·C_eq)^n_LF]

# Arguments
- `C_eq`: Equilibrium sorbed DOC concentration [μg/mm³ solid]
- `M_max`: Maximum sorption capacity [μg/mm³]
- `k_L`: Langmuir affinity constant [mm³/μg]
- `n_LF`: Freundlich heterogeneity exponent [-]

# Returns
- Equilibrium MAOC [μg/mm³]

# Notes
- **CRITICAL**: Uses C_eq = k_d·C_aq, NOT C or C_aq directly
- n_LF < 1: heterogeneous sites (typical for soil minerals)
- n_LF = 1: recovers standard Langmuir isotherm
- M_max depends on soil texture: M_max = k_ma · f_clay_silt · ρ_b
- As C_eq → 0: M_eq → 0
- As C_eq → ∞: M_eq → M_max
"""
function M_eq_langmuir_freundlich(C_eq::Real, M_max::Real, k_L::Real, n_LF::Real)
    k_L_C = k_L * C_eq
    k_L_C_n = k_L_C^n_LF

    M_max * k_L_C_n / (1.0 + k_L_C_n)
end

"""
    J_M(M::Real, M_eq::Real, κ_s_T::Real, κ_d_T::Real, ε_maoc::Real)

Net MAOC flux with smooth switching between sorption and desorption.

# Formula (manuscript Eq. 237)
    J_M = κ_s(T)·φ_ε(M_eq - M) - κ_d(T)·φ_ε(M - M_eq)

where φ_ε is softplus regularization.

# Arguments
- `M`: Current MAOC concentration [μg/mm³]
- `M_eq`: Equilibrium MAOC (from M_eq_langmuir_freundlich) [μg/mm³]
- `κ_s_T`: Sorption rate at current T [1/day]
- `κ_d_T`: Desorption rate at current T [1/day]
- `ε_maoc`: Softplus smoothing width [μg/mm³] (default 0.01)

# Returns
- Net MAOC flux [μg-C/mm³/day]

# Notes
- When M < M_eq: sorption dominates (J_M > 0, M increases)
- When M > M_eq: desorption dominates (J_M < 0, M decreases)
- κ_s_T, κ_d_T already include Arrhenius factors
- Softplus ensures C∞ smoothness near M = M_eq
- **CRITICAL S_C coupling**: Contribution to S_C is -J_M·(θ + ρ_b·k_d)/k_d
  (NOT just -J_M! Missing factor breaks carbon conservation.)
"""
function J_M(M::Real, M_eq::Real, κ_s_T::Real, κ_d_T::Real, ε_maoc::Real)
    sorption = κ_s_T * softplus(M_eq - M, ε_maoc)
    desorption = κ_d_T * softplus(M - M_eq, ε_maoc)

    sorption - desorption
end
