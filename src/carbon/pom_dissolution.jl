# pom_dissolution.jl
# POM dissolution: enzymatic breakdown at the POM surface
# Flux density J_P [μg-C/mm²/day] and total rate R_P [μg-C/day]

"""
    J_P(P, P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T, K_B_P, K_F_P, θ_P, L_P)

Compute POM dissolution flux density at the POM surface [μg-C/mm²/day].

# Arguments
- `P::Real`: Current POM mass [μg-C]
- `P_0::Real`: Initial POM mass [μg-C]
- `B_0::Real`: Bacterial biomass at r=r₀ [μg-C/mm³]
- `F_n_0::Real`: Non-insulated fungi at r=r₀ [μg-C/mm³]
- `θ_0::Real`: Water content at r=r₀ [-]
- `O_aq_0::Real`: Aqueous oxygen concentration at r=r₀ [μg/mm³]
- `R_P_max_T::Real`: Max dissolution rate at temperature T [μg-C/mm²/day]
- `K_B_P::Real`: Half-saturation bacteria [μg-C/mm³]
- `K_F_P::Real`: Half-saturation fungi [μg-C/mm³]
- `θ_P::Real`: Half-saturation water content [-]
- `L_P::Real`: Half-saturation O₂ [μg/mm³]

# Returns
- Flux density [μg-C/mm²/day]

# Manuscript reference
Equation for R_P in Section on POM dissolution (Architecture §4.2)
"""
function J_P(P::Real, P_0::Real, B_0::Real, F_n_0::Real, θ_0::Real, O_aq_0::Real,
             R_P_max_T::Real, K_B_P::Real, K_F_P::Real, θ_P::Real, L_P::Real)
    # POM availability (normalized by initial mass)
    pom_factor = P / P_0

    # Microbial enzyme contributions (Monod kinetics)
    # Additive coupling: bacteria and fungi contribute independently
    bacterial_contribution = B_0 / (K_B_P + B_0)
    fungal_contribution = F_n_0 / (K_F_P + F_n_0)
    microbial_factor = 0.5 * (bacterial_contribution + fungal_contribution)

    # Moisture limitation
    moisture_factor = θ_0 / (θ_P + θ_0)

    # Oxygen limitation
    oxygen_factor = O_aq_0 / (L_P + O_aq_0)

    # Total flux density
    R_P_max_T * pom_factor * microbial_factor * moisture_factor * oxygen_factor
end

"""
    R_P(J_P_val, r_0)

Compute total POM dissolution rate [μg-C/day] from flux density.

# Arguments
- `J_P_val::Real`: Flux density at POM surface [μg-C/mm²/day]
- `r_0::Real`: POM radius [mm]

# Returns
- Total dissolution rate [μg-C/day]

# Notes
The total rate is the flux density integrated over the spherical surface:
R_P = 4π·r₀²·J_P
"""
function R_P(J_P_val::Real, r_0::Real)
    4.0 * π * r_0^2 * J_P_val
end

"""
    pom_delay_factor(t::Real, t_delay::Real)

Multiplier on `R_P_max` implementing the POM activation delay [-].

    t_delay ≤ 0  ->  1.0 exactly (no delay)
    t_delay > 0   ->  sigmoid_threshold(t, t_delay, POM_DELAY_STEEPNESS)

The switch is centred on `t = t_delay` (factor 0.5 there) with a transition
width of `t_delay / POM_DELAY_STEEPNESS`, i.e. 10 % of the delay: `t_delay = 10`
reaches 90 % activation at `t ≈ 12.2`. Its purpose is to let the microbial pools
equilibrate before POM input begins.

Both solvers call this. It was written out separately in `timestepper.jl` and
`mol.jl`, which is the one place the two integrators could silently disagree on
when POM turns on — and the audit's open question about them is whether they
agree.
"""
function pom_delay_factor(t::Real, t_delay::Real)
    t_delay > 0.0 || return 1.0
    return sigmoid_threshold(t, t_delay, POM_DELAY_STEEPNESS)
end
