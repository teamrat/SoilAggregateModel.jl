"""
    fungi.jl

Fungal uptake, allocation, transitions, and turnover.

Implements three-pool fungal dynamics (F_i, F_n, F_m) with mobilization/
immobilization transitions and conversion respiration cost.

# Dependencies
Requires: BiologicalProperties

# Critical notes
- Resp_F_conv uses abs() to handle mobilization correctly (MANUSCRIPT_CHANGES #2)
- Π = F_m / (F_i + F_n + ε_F) has division-by-zero protection
- All transition rates share single Ea_F activation energy
"""

"""
    Pi_protected(F_m::Real, F_i::Real, F_n::Real, ε_F::Real)

Mobile-to-immobile biomass ratio with division-by-zero protection.

# Formula (manuscript Eq. 350)
    Π = F_m / (F_i + F_n + ε_F)

# Arguments
- `F_m`: Mobile fungi [μg/mm³]
- `F_i`: Insulated fungi [μg/mm³]
- `F_n`: Non-insulated fungi [μg/mm³]
- `ε_F`: Regularization constant [μg/mm³] (default 1e-4)

# Returns
- Π [-] (dimensionless ratio)

# Notes
- ε_F prevents division by zero when F_i + F_n → 0
- ε_F = 1e-4 (comparable to B_min) naturally bounds Π at physically reasonable values
- When F_i + F_n → 0: Π → F_m/ε_F remains finite (e.g., 2e-5/1e-4 = 0.2)
- Π appears in all transition rate expressions
- Nonlinear mobilization (∝ Π^δ with δ > 1) provides stability
"""
function Pi_protected(F_m::Real, F_i::Real, F_n::Real, ε_F::Real)
    F_m / (F_i + F_n + ε_F)
end

"""
    R_F(C_aq::Real, O_aq::Real, F_i::Real, F_n::Real, λ::Real, ψ::Real,
        r_F_max_T::Real, K_F::Real, L_F::Real, ν_F::Real)

Fungal carbon uptake rate (gross).

# Formula (manuscript Eq. 282)
    R_F = r_F_max(T) · C_aq/(K_F + C_aq) · O_aq/(L_F + O_aq) · (λ·F_i + F_n) · exp(ν_F·ψ)

# Arguments
- `C_aq`: Aqueous DOC concentration [μg/mm³ water]
- `O_aq`: Aqueous oxygen concentration [μg/mm³ water]
- `F_i`: Insulated fungi [μg/mm³]
- `F_n`: Non-insulated fungi [μg/mm³]
- `λ`: Reduced uptake efficiency of insulated hyphae (λ ≪ 1)
- `ψ`: Matric potential [kPa]
- `r_F_max_T`: Max specific uptake rate at current T [1/day]
- `K_F`: Half-saturation for carbon [μg/mm³]
- `L_F`: Half-saturation for oxygen [μg/mm³]
- `ν_F`: Water potential sensitivity [1/kPa] (ν_F < ν_B, fungi more drought-tolerant)

# Returns
- Uptake rate [μg-C/mm³/day]

# Notes
- Uses C_aq, NOT raw C
- F_n (primary uptaker) and λ·F_i (insulated, reduced efficiency) contribute to uptake
- F_m does not directly uptake (internal translocation only)
- r_F_max_T already includes Arrhenius factor
"""
function R_F(C_aq::Real, O_aq::Real, F_i::Real, F_n::Real, λ::Real, ψ::Real,
             r_F_max_T::Real, K_F::Real, L_F::Real, ν_F::Real)
    monod_C = C_aq / (K_F + C_aq)
    monod_O = O_aq / (L_F + O_aq)
    uptake_biomass = λ * F_i + F_n
    water_stress = exp(ν_F * ψ)

    r_F_max_T * monod_C * monod_O * uptake_biomass * water_stress
end

"""
    Y_F_const(Y_F::Real)

Constant fungal yield (simple model).

# Arguments
- `Y_F`: Constant yield [-]

# Returns
- Yield [-]

# Notes
- Simplest model: Y_F is fixed parameter
- Alternative: Y_F_uptake_dependent for Monod-like yield response
"""
function Y_F_const(Y_F::Real)
    Y_F
end

"""
    Y_F_uptake_dependent(R_F_val::Real, Y_F_max::Real, K_YF::Real)

Uptake-dependent fungal yield (alternative model).

# Formula (manuscript Eq. 327, alternative)
    Y_F = Y_F_max · R_F / (R_F + K_YF)

# Arguments
- `R_F_val`: Fungal uptake rate [μg-C/mm³/day]
- `Y_F_max`: Maximum yield [-]
- `K_YF`: Half-saturation for yield response [μg-C/mm³/day]

# Returns
- Effective yield [-]

# Notes
- Use this OR Y_F_const, not both
- Parameters specify which model via Y_F field interpretation
"""
function Y_F_uptake_dependent(R_F_val::Real, Y_F_max::Real, K_YF::Real)
    Y_F_max * R_F_val / (R_F_val + K_YF)
end

"""
    Y_F_func(Y_F_base::Real, F_i::Real, F_n::Real, F_m::Real, F_S::Real)

Fungal yield with space limitation.

# Formula
    Y_F = Y_F_base · F_S / (F_total + F_S)

where F_total = F_i + F_n + F_m.

# Arguments
- `Y_F_base`: Base fungal yield (without space limitation) [-]
- `F_i`: Insulated fungi [μg/mm³]
- `F_n`: Non-insulated fungi [μg/mm³]
- `F_m`: Mobile fungi [μg/mm³]
- `F_S`: Half-saturation for space limitation [μg/mm³]

# Returns
- Effective yield [-] (0 ≤ Y_F ≤ Y_F_base)

# Notes
- Space limitation: Y_F → 0 as F_total → ∞ (prevents unbounded growth)
- At F_total = F_S: yield reduced by 50% (space-limited)
- When F_total << F_S: Y_F ≈ Y_F_base (space-unlimited)
"""
function Y_F_func(Y_F_base::Real, F_i::Real, F_n::Real, F_m::Real, F_S::Real)
    F_total = F_i + F_n + F_m
    Y_F_base * F_S / (F_total + F_S)
end

"""
    Gamma_F(Y_F_val::Real, R_F_val::Real)

Total fungal carbon assimilation (enters F_m pool).

# Formula (manuscript Eq. 333)
    Γ_F = Y_F · R_F

# Arguments
- `Y_F_val`: Current yield [-]
- `R_F_val`: Uptake rate [μg-C/mm³/day]

# Returns
- Assimilation rate [μg-C/mm³/day]

# Notes
- All assimilated carbon enters F_m (mobile pool)
- F_m then distributes to F_i and F_n via transitions
"""
function Gamma_F(Y_F_val::Real, R_F_val::Real)
    Y_F_val * R_F_val
end

"""
    Resp_F(R_F_val::Real, Y_F_val::Real)

Direct fungal respiration (from uptake).

# Formula (manuscript Eq. 339)
    Resp_F = (1 - Y_F) · R_F

# Arguments
- `R_F_val`: Uptake rate [μg-C/mm³/day]
- `Y_F_val`: Current yield [-]

# Returns
- Respiration rate [μg-C/mm³/day]

# Notes
- Fraction (1 - Y_F) of uptake is respired
- Does NOT include conversion respiration (see Resp_F_conv)
- Contributes to S_O: O₂ consumption = α_O · (Resp_F + Resp_F_conv)
"""
function Resp_F(R_F_val::Real, Y_F_val::Real)
    (1.0 - Y_F_val) * R_F_val
end

"""
    h_Fi(F_i::Real, F_i_min::Real)

Sigmoid threshold function for insulated fungal viability.

# Formula (manuscript Eq. 389)
    h_Fi = exp(β_F·F_i) / [exp(β_F·F_i) + exp(β_F·F_i_min)]

where β_F = 50/F_i_min.

# Arguments
- `F_i`: Insulated fungal biomass [μg/mm³]
- `F_i_min`: Minimum viable concentration [μg/mm³]

# Returns
- Threshold factor [-] (0 ≤ h_Fi ≤ 1)

# Notes
- As F_i → 0: h_Fi → 0 (shuts off death)
- As F_i → ∞: h_Fi → 1 (full mortality)
- Applied to R_rec_F (fungal death)
- Analogous to h_B for bacteria
"""
function h_Fi(F_i::Real, F_i_min::Real)
    β_F = 50.0 / F_i_min
    # Numerically stable sigmoid: 1 / (1 + exp(-β_F * (F_i - F_i_min)))
    # Avoids overflow when F_i >> F_i_min
    1.0 / (1.0 + exp(-β_F * (F_i - F_i_min)))
end

"""
    R_rec_F(μ_F_T::Real, F_i::Real, F_i_min::Real)

Fungal death (recycling) rate - acts only on insulated pool.

# Formula
    R_rec_F = μ_F(T) · F_i · h_Fi

# Arguments
- `μ_F_T`: Mortality rate at current T [1/day]
- `F_i`: Insulated fungal biomass [μg/mm³]
- `F_i_min`: Minimum viable concentration [μg/mm³]

# Returns
- Death rate [μg-C/mm³/day]

# Notes
- μ_F_T already includes Arrhenius factor
- Only F_i dies (mature hyphae); F_n and F_m are sustained by network
- Returns carbon to DOC pool (contributes to R_rec in S_C)
"""
function R_rec_F(μ_F_T::Real, F_i::Real, F_i_min::Real)
    h_Fi_val = h_Fi(F_i, F_i_min)
    μ_F_T * F_i * h_Fi_val
end

"""
    fungal_transitions(F_i::Real, F_n::Real, F_m::Real, Π::Real,
                      α_i_T::Real, α_n_T::Real, β_i_T::Real, β_n_T::Real,
                      ζ_T::Real, δ::Real, η::Real, ε_F::Real)

Compute all fungal transition rates and conversion respiration.

# Returns
Named tuple with:
- `trans_i`: Net transfer to F_i from F_m [μg-C/mm³/day]
- `trans_n`: Net transfer to F_n from F_m [μg-C/mm³/day]
- `insulation`: Insulation F_n → F_i [μg-C/mm³/day]
- `Resp_F_conv`: Conversion respiration cost [μg-C/mm³/day]

# Formulas (manuscript Eqs. 360-363)
Insulation (non-insulated → insulated):
    ζ(T) · F_n

Net immobilization to F_i (mobile → insulated):
    η · (β_i(T)·Π - α_i(T)·Π^δ) · F_i

Net immobilization to F_n (mobile → non-insulated):
    η · (β_n(T)·Π - α_n(T)·Π^δ) · F_n

Conversion respiration (CRITICAL: uses abs()):
    (1-η) · |[(β_i·Π - α_i·Π^δ)·F_i + (β_n·Π - α_n·Π^δ)·F_n]|

# Arguments
- `F_i`, `F_n`, `F_m`: Three fungal pools [μg/mm³]
- `Π`: Mobile-to-immobile ratio (from Pi_protected) [-]
- `α_i_T`, `α_n_T`: Mobilization rates at current T [1/day]
- `β_i_T`, `β_n_T`: Immobilization rates at current T [1/day]
- `ζ_T`: Insulation rate at current T [1/day]
- `δ`: Mobilization exponent [-] (δ > 1)
- `η`: Conversion efficiency [-] (η < 1)
- `ε_F`: Regularization for Π [μg/mm³]

# Notes
- **CRITICAL**: Resp_F_conv uses abs() (MANUSCRIPT_CHANGES #2)
- When Π is large: immobilization dominates (β·Π > α·Π^δ)
- When Π is small: mobilization dominates (α·Π^δ > β·Π)
- Conversion cost (1-η) applies to both directions
- All rates scale with Arrhenius factor (already applied)
"""
function fungal_transitions(F_i::Real, F_n::Real, F_m::Real, Π::Real,
                           α_i_T::Real, α_n_T::Real, β_i_T::Real, β_n_T::Real,
                           ζ_T::Real, δ::Real, η::Real, ε_F::Real)
    # Insulation: F_n → F_i (no conversion cost)
    insulation = ζ_T * F_n

    # Net transition rates (positive = immobilization, negative = mobilization)
    Π_delta = Π^δ
    net_i = (β_i_T * Π - α_i_T * Π_delta) * F_i
    net_n = (β_n_T * Π - α_n_T * Π_delta) * F_n

    # Retained fractions after conversion cost
    trans_i = η * net_i
    trans_n = η * net_n

    # Conversion respiration: CRITICAL abs() for mobilization
    # When mobilization dominates (net < 0), respiration is still positive
    Resp_F_conv = (1.0 - η) * abs(net_i + net_n)

    (trans_i = trans_i, trans_n = trans_n, insulation = insulation,
     Resp_F_conv = Resp_F_conv)
end
