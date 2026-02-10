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
"""
struct SourceTerms
    S_C::Float64
    S_B::Float64
    S_Fn::Float64
    S_Fm::Float64
    S_Fi::Float64
    S_E::Float64
    S_M::Float64
    S_O::Float64
    Resp_total::Float64
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
  - `O`: Total oxygen [μg/mm³]
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
- CRITICAL: Uses corrected formulas (no MAOC factor in S_C, trans_n already has η)

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
    O_aq = O * θ / (θ + temp_cache.K_H_O * θ_a)

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
    α_i_T = bio.α_i * temp_cache.f_fun
    α_n_T = bio.α_n * temp_cache.f_fun
    β_i_T = bio.β_i * temp_cache.f_fun
    β_n_T = bio.β_n * temp_cache.f_fun
    ζ_T = bio.ζ * temp_cache.f_fun

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

    # Transitions (insulation, mobilization, translocation, conversion respiration)
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
    # CRITICAL: S_Fn uses trans.trans_n directly (already contains η)
    S_C_val = -R_B_val - R_F_val + R_rec_total - J_M_val
    S_B_val = Γ_B_val - R_rec_B_val
    S_Fn_val = trans.trans_n - trans.insulation
    S_Fm_val = Γ_F_val - trans.trans_i - trans.trans_n - trans.Resp_F_conv
    S_Fi_val = trans.insulation + trans.trans_i - R_rec_F_val
    S_E_val = Γ_E_val - R_rec_E_val
    S_M_val = J_M_val
    S_O_val = -bio.α_O * Resp_total_val

    SourceTerms(S_C_val, S_B_val, S_Fn_val, S_Fm_val, S_Fi_val, S_E_val, S_M_val,
                S_O_val, Resp_total_val)
end
