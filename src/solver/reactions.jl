# reactions.jl
# Compute all source/sink terms at a single grid node

"""
    SourceTerms

Source/sink terms for all 8 state variables at a single node [all in μg/mm³/day].

Fields:
- `S_C`: DOC source/sink
- `S_B`: Bacteria source/sink
- `S_Fn`: Non-insulated fungi source/sink
- `S_Fm`: Mobile fungi source/sink
- `S_Fi`: Insulated fungi source/sink
- `S_E`: EPS source/sink
- `S_M`: MAOC source/sink
- `S_O`: Oxygen source/sink
- `Resp_total`: Total respiration (for CO₂ accumulation) [μg-C/mm³/day]

Parametric in the element type so the same function serves both the
operator-split solver (`Float64`) and the method-of-lines right-hand side,
whose Jacobian is built by forward-mode AD and therefore evaluates every rate
law on dual numbers. A hard-coded `Float64` here would silently force the stiff
solver onto finite differences.
"""
struct SourceTerms{T<:Real}
    S_C::T
    S_B::T
    S_Fn::T
    S_Fm::T
    S_Fi::T
    S_E::T
    S_M::T
    S_O::T
    Resp_total::T
end

function SourceTerms(S_C::Real, S_B::Real, S_Fn::Real, S_Fm::Real, S_Fi::Real,
                     S_E::Real, S_M::Real, S_O::Real, Resp_total::Real)
    v = promote(S_C, S_B, S_Fn, S_Fm, S_Fi, S_E, S_M, S_O, Resp_total)
    return SourceTerms{typeof(v[1])}(v...)
end

"""
    compute_source_terms(C, B, F_n, F_m, F_i, E, M, O, θ, θ_a, ψ,
                        bio::BiologicalProperties, soil::SoilProperties,
                        temp_cache::TemperatureCache)

Compute all source/sink terms at a single grid node.

# Arguments
- State variables at the node:
  - `C`: Total DOC [μg-C/mm³]
  - `B`: Bacteria [μg-C/mm³]
  - `F_n`: Non-insulated fungi [μg-C/mm³]
  - `F_m`: Mobile fungi [μg-C/mm³]
  - `F_i`: Insulated fungi [μg-C/mm³]
  - `E`: EPS [μg-C/mm³]
  - `M`: MAOC [μg-C/mm³]
  - `O`: Aqueous oxygen concentration [μg/mm³]
  - `θ`: Water content [-]
  - `θ_a`: Air-filled porosity [-]
  - `ψ`: Water potential [kPa]
- Parameters:
  - `bio`: BiologicalProperties
  - `soil`: SoilProperties
  - `temp_cache`: TemperatureCache with temperature-dependent values

# Returns
- `SourceTerms`: Struct with all 8 source terms + total respiration

# Notes
- All biology functions are imported from biology/*.jl
- Follows exact computation order from test_biology.jl (verified against manuscript)
- Zero allocation (all computations scalar)
- CRITICAL: ζ splitting is handled inside fungal_transitions() — no separate insulation term
  (matches MATLAB / Falconer 2005, 2008)

# Manuscript reference
Architecture §4: Reaction step, source/sink terms
"""
function compute_source_terms(C::Real, B::Real, F_n::Real, F_m::Real, F_i::Real,
                              E::Real, M::Real, O::Real, θ::Real, θ_a::Real, ψ::Real,
                              bio::BiologicalProperties, soil::SoilProperties,
                              temp_cache::TemperatureCache)
    # === STEP 1: Compute C_aq, C_eq, O_aq ONCE ===
    C_aq = C / (θ + soil.ρ_b * soil.k_d_eq)
    C_eq = soil.k_d_eq * C_aq
    O_aq = O    

    # === STEP 2: Bacterial terms ===
    # Temperature-dependent rates
    r_B_max_T = bio.r_B_max * temp_cache.f_bac
    μ_B_T = bio.μ_B * temp_cache.f_bac

    # Uptake and maintenance
    R_B_val = R_B(C_aq, O_aq, B, ψ, r_B_max_T, bio.K_B, bio.L_B, bio.ν_B)
    R_Bb_val = R_Bb(bio.C_B, O_aq, B, r_B_max_T, bio.K_B, bio.L_B, bio.B_min)

    # Yield (depends on uptake - maintenance difference, and space limitation)
    R_diff = R_B_val - R_Bb_val
    Y_B_val = Y_B_func(R_diff, bio.Y_B_max, bio.K_Y, B, bio.B_S, bio.ε_Y)

    # Growth
    Γ_B_val = Gamma_B(R_B_val, R_Bb_val, Y_B_val, bio.γ, bio.ε_Y)
    Γ_E_val = Gamma_E(R_B_val, R_Bb_val, Y_B_val, bio.γ, bio.ε_Y)

    # Respiration and recycling
    Resp_B_val = Resp_B(R_Bb_val, R_diff, Y_B_val, bio.ε_Y)
    R_rec_B_val = R_rec_B(μ_B_T, B, bio.B_min)

    # Total recycled carbon from bacteria
    R_rec_bacteria = R_rec_B_val

    # === STEP 3: Fungal terms ===
    # Temperature-dependent rates
    r_F_max_T = bio.r_F_max * temp_cache.f_fun
    μ_F_T = bio.μ_F * temp_cache.f_fun
    # α = immobilization (gain), β = mobilization (loss) — Falconer naming
    α_i_T = bio.α_i * temp_cache.f_fun
    α_n_T = bio.α_n * temp_cache.f_fun
    β_i_T = bio.β_i * temp_cache.f_fun
    β_n_T = bio.β_n * temp_cache.f_fun
    # NOTE: ζ is a dimensionless SPLITTING FRACTION, not a rate, so scaling it
    # by an Arrhenius factor is dimensionally wrong and can drive it above 1
    # (which would flip the sign of immobil_n). Clamped until this is resolved.
    ζ_T = min(bio.ζ * temp_cache.f_fun, 1.0)

    # Protection ratio
    Π_val = Pi_protected(F_m, F_i, F_n, bio.ε_F)

    # Uptake (NOTE: R_F takes F_i, F_n as separate arguments, NOT F_m)
    R_F_val = R_F(C_aq, O_aq, F_i, F_n, bio.λ, ψ, r_F_max_T, bio.K_F, bio.L_F, bio.ν_F)

    # Yield (space-limited)
    Y_F_val = Y_F_func(bio.Y_F, F_i, F_n, F_m, bio.F_S)

    # Growth
    Γ_F_val = Gamma_F(Y_F_val, R_F_val)

    # Respiration
    Resp_F_val = Resp_F(R_F_val, Y_F_val)

    # Transitions (ζ splitting already applied inside fungal_transitions)
    trans = fungal_transitions(F_i, F_n, F_m, Π_val, α_i_T, α_n_T, β_i_T, β_n_T,
                               ζ_T, bio.delta, bio.η_conv, bio.ε_F)

    # Recycling (death)
    R_rec_F_val = R_rec_F(μ_F_T, F_i, bio.F_i_min)

    # Total recycled carbon from fungi
    R_rec_fungi = R_rec_F_val

    # === STEP 4: EPS terms ===
    μ_E_max_T = bio.μ_E_max * temp_cache.f_eps
    R_rec_E_val = R_rec_E(μ_E_max_T, C_aq, bio.K_E, E, bio.E_min)

    # Total recycled carbon
    R_rec_total = R_rec_bacteria + R_rec_fungi + R_rec_E_val

    # === STEP 5: MAOC terms ===
    κ_s_T = bio.κ_s_ref * temp_cache.f_maoc_s
    κ_d_T = bio.κ_d_ref * temp_cache.f_maoc_d
    M_eq_val = M_eq_langmuir_freundlich(C_eq, soil.M_max, soil.k_L, soil.n_LF)
    J_M_val = J_M(M, M_eq_val, κ_s_T, κ_d_T, bio.ε_maoc)

    # === STEP 6: Total respiration ===
    Resp_total_val = Resp_B_val + Resp_F_val + trans.Resp_F_conv

    # === STEP 7: Compute source terms ===
    # CRITICAL: S_C uses corrected formula (no factor on J_M)
    # CRITICAL: ζ splitting handled in fungal_transitions — immobil_i includes ζ·trans_n,
    #           immobil_n = (1-ζ)·trans_n. No separate insulation term.
    #           (Matches MATLAB lines 389-390, 404-406; Falconer 2005/2008)
    S_C_val = -R_B_val - R_F_val + R_rec_total - J_M_val
    S_B_val = Γ_B_val - R_rec_B_val
    S_Fn_val = trans.immobil_n
    S_Fm_val = Γ_F_val - trans.immobil_i - trans.immobil_n - trans.Resp_F_conv
    S_Fi_val = trans.immobil_i - R_rec_F_val
    S_E_val = Γ_E_val - R_rec_E_val
    S_M_val = J_M_val
    capacity_O = θ + temp_cache.K_H_O * θ_a
    S_O_val = -bio.α_O * Resp_total_val / capacity_O

    SourceTerms(S_C_val, S_B_val, S_Fn_val, S_Fm_val, S_Fi_val, S_E_val, S_M_val,
                S_O_val, Resp_total_val)
end
