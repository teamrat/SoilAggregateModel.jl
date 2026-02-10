"""
    bacteria.jl

Bacterial uptake, allocation, and turnover.

All functions compute a single term from the manuscript. Temperature-dependent
rates (r_B_max(T), μ_B(T)) are passed in pre-computed from TemperatureCache.

# Dependencies
Requires: BiologicalProperties

# Critical C usage
- R_B uses C_aq (aqueous concentration), NOT raw C
- R_Bb uses C_B (constant parameter), NOT C_aq or C
- Compute C_aq = C / (θ + ρ_b·k_d) ONCE per node, pass to all functions
"""

"""
    h_B(B::Real, B_min::Real)

Sigmoid threshold function for bacterial viability.

# Formula
    h_B = exp(β·B) / [exp(β·B) + exp(β·B_min)]

where β = 50/B_min.

# Arguments
- `B`: Bacterial biomass [μg/mm³]
- `B_min`: Minimum viable concentration [μg/mm³]

# Returns
- Threshold factor [-] (0 ≤ h_B ≤ 1)

# Notes
- As B → 0: h_B → 0 (shuts off maintenance and death)
- As B → ∞: h_B → 1 (full activity)
- Applied to R_Bb (maintenance) and R_rec_B (death)
- Prevents biomass from going negative
"""
function h_B(B::Real, B_min::Real)
    β = 50.0 / B_min
    # Numerically stable sigmoid: 1 / (1 + exp(-β * (B - B_min)))
    # Avoids overflow when B >> B_min
    1.0 / (1.0 + exp(-β * (B - B_min)))
end

"""
    R_B(C_aq::Real, O_aq::Real, B::Real, ψ::Real, r_B_max_T::Real,
        K_B::Real, L_B::Real, ν_B::Real)

Bacterial carbon uptake rate (gross).

# Formula (manuscript Eq. 281)
    R_B = r_B_max(T) · C_aq/(K_B + C_aq) · O_aq/(L_B + O_aq) · B · exp(ν_B·ψ)

# Arguments
- `C_aq`: Aqueous DOC concentration [μg/mm³ water]
- `O_aq`: Aqueous oxygen concentration [μg/mm³ water]
- `B`: Bacterial biomass [μg/mm³]
- `ψ`: Matric potential [kPa] (ψ ≤ 0)
- `r_B_max_T`: Max specific uptake rate at current T [1/day]
- `K_B`: Half-saturation for carbon [μg/mm³]
- `L_B`: Half-saturation for oxygen [μg/mm³]
- `ν_B`: Water potential sensitivity [1/kPa]

# Returns
- Uptake rate [μg-C/mm³/day]

# Notes
- Uses C_aq, NOT raw C state variable
- O_aq = O·θ / (θ + K_H·θ_a) computed elsewhere
- r_B_max_T already includes Arrhenius factor
- exp(ν_B·ψ) exponentially reduces uptake in dry soil
"""
function R_B(C_aq::Real, O_aq::Real, B::Real, ψ::Real, r_B_max_T::Real,
             K_B::Real, L_B::Real, ν_B::Real)
    monod_C = C_aq / (K_B + C_aq)
    monod_O = O_aq / (L_B + O_aq)
    water_stress = exp(ν_B * ψ)

    r_B_max_T * monod_C * monod_O * B * water_stress
end

"""
    R_Bb(C_B::Real, O_aq::Real, B::Real, r_B_max_T::Real, L_B::Real, B_min::Real)

Bacterial basal maintenance metabolism.

# Formula (manuscript Eq. 296)
    R_Bb = r_B_max(T) · C_B/(K_B + C_B) · O_aq/(L_B + O_aq) · B · h_B

# Arguments
- `C_B`: Basal carbon requirement [μg/mm³] (constant parameter)
- `O_aq`: Aqueous oxygen concentration [μg/mm³ water]
- `B`: Bacterial biomass [μg/mm³]
- `r_B_max_T`: Max specific uptake rate at current T [1/day]
- `L_B`: Half-saturation for oxygen [μg/mm³]
- `B_min`: Minimum viable biomass [μg/mm³]

# Returns
- Maintenance uptake rate [μg-C/mm³/day]

# Notes
- Uses C_B (constant parameter), NOT C_aq or C
- h_B sigmoid ensures R_Bb → 0 as B → B_min
- R_Bb defines minimum carbon flux for survival
- If R_B < R_Bb, bacteria catabolize biomass
"""
function R_Bb(C_B::Real, O_aq::Real, B::Real, r_B_max_T::Real, K_B::Real,
              L_B::Real, B_min::Real)
    monod_C_B = C_B / (K_B + C_B)
    monod_O = O_aq / (L_B + O_aq)
    h_B_val = h_B(B, B_min)

    r_B_max_T * monod_C_B * monod_O * B * h_B_val
end

"""
    Y_B_func(R_diff::Real, Y_B_max::Real, K_Y::Real, B::Real, B_S::Real, ε_Y::Real)

Bacterial growth yield — substrate-limited and space-limited.

# Formula (modified from manuscript Eq. 308)
    Y_B = Y_B_max · softplus(R_diff, ε_Y) / [softplus(R_diff, ε_Y) + K_Y] · B_S / (B + B_S)

where R_diff = R_B - R_Bb.

# Arguments
- `R_diff`: Excess uptake (R_B - R_Bb) [μg-C/mm³/day]
- `Y_B_max`: Maximum yield [-]
- `K_Y`: Half-saturation for yield response [μg-C/mm³/day]
- `B`: Bacterial biomass [μg/mm³]
- `B_S`: Half-saturation for space limitation [μg/mm³]
- `ε_Y`: Softplus smoothing width [μg-C/mm³/day]

# Returns
- Effective yield [-] (0 ≤ Y_B ≤ Y_B_max)

# Notes
- Substrate limitation: Y_B → Y_B_max as R_diff → ∞ (Monod-like)
- Space limitation: Y_B → 0 as B → ∞ (prevents unbounded growth)
- At B = B_S: yield reduced by 50% (space-limited)
- softplus replaces max(0, R_diff) for smooth C∞ transition at R_diff = 0
"""
function Y_B_func(R_diff::Real, Y_B_max::Real, K_Y::Real, B::Real, B_S::Real, ε_Y::Real)
    R_diff_smooth = softplus(R_diff, ε_Y)
    Y_B_max * R_diff_smooth / (R_diff_smooth + K_Y) * B_S / (B + B_S)
end

"""
    Gamma_B(R_B_val::Real, R_Bb_val::Real, Y_B_val::Real, γ::Real, ε_Y::Real)

Bacterial biomass growth rate with smooth allocation.

# Formula (modified from manuscript Eq. 314)
    Γ_B = Y_B · softplus(R_diff, ε_Y) · (1 - γ) - softplus(-R_diff, ε_Y)

where R_diff = R_B - R_Bb.

# Arguments
- `R_B_val`: Gross uptake [μg-C/mm³/day]
- `R_Bb_val`: Maintenance uptake [μg-C/mm³/day]
- `Y_B_val`: Current yield [-]
- `γ`: EPS allocation fraction [-]
- `ε_Y`: Softplus smoothing width [μg-C/mm³/day]

# Returns
- Growth rate [μg-C/mm³/day]

# Notes
- When R_diff > 0: Γ_B ≈ Y_B·R_diff·(1-γ) (growth at yield, fraction to biomass)
- When R_diff < 0: Γ_B ≈ R_diff < 0 (biomass catabolism to meet deficit)
- softplus provides smooth C∞ transition at R_diff = 0
- Fraction γ of assimilated carbon goes to EPS (see Gamma_E)
"""
function Gamma_B(R_B_val::Real, R_Bb_val::Real, Y_B_val::Real, γ::Real, ε_Y::Real)
    R_diff = R_B_val - R_Bb_val
    Y_B_val * softplus(R_diff, ε_Y) * (1.0 - γ) - softplus(-R_diff, ε_Y)
end

"""
    Gamma_E(R_B_val::Real, R_Bb_val::Real, Y_B_val::Real, γ::Real, ε_Y::Real)

Bacterial EPS production rate with smooth allocation.

# Formula (modified from manuscript Eq. 315)
    Γ_E = Y_B · softplus(R_diff, ε_Y) · γ

where R_diff = R_B - R_Bb.

# Arguments
- `R_B_val`: Gross uptake [μg-C/mm³/day]
- `R_Bb_val`: Maintenance uptake [μg-C/mm³/day]
- `Y_B_val`: Current yield [-]
- `γ`: EPS allocation fraction [-]
- `ε_Y`: Softplus smoothing width [μg-C/mm³/day]

# Returns
- EPS production rate [μg-C/mm³/day]

# Notes
- Only produces EPS when R_diff > 0 (excess uptake)
- softplus provides smooth C∞ transition at R_diff = 0
- Fraction γ of assimilated carbon (typically 0.2-0.5)
"""
function Gamma_E(R_B_val::Real, R_Bb_val::Real, Y_B_val::Real, γ::Real, ε_Y::Real)
    R_diff = R_B_val - R_Bb_val
    Y_B_val * softplus(R_diff, ε_Y) * γ
end

"""
    Resp_B(R_Bb_val::Real, R_diff::Real, Y_B_val::Real, ε_Y::Real)

Bacterial respiration rate with smooth allocation.

# Formula (modified from manuscript Eq. 321)
    Resp_B = R_Bb + softplus(R_diff, ε_Y) · (1 - Y_B)

where R_diff = R_B - R_Bb.

# Arguments
- `R_Bb_val`: Maintenance uptake [μg-C/mm³/day]
- `R_diff`: Excess uptake (R_B - R_Bb) [μg-C/mm³/day]
- `Y_B_val`: Current yield [-]
- `ε_Y`: Softplus smoothing width [μg-C/mm³/day]

# Returns
- Respiration rate (CO₂ production) [μg-C/mm³/day]

# Notes
- Maintenance (R_Bb) is fully respired
- Growth uptake: fraction (1 - Y_B) respired, fraction Y_B assimilated
- softplus provides smooth C∞ transition at R_diff = 0
- Contributes to S_O: O₂ consumption = α_O · Resp_B
"""
function Resp_B(R_Bb_val::Real, R_diff::Real, Y_B_val::Real, ε_Y::Real)
    R_Bb_val + softplus(R_diff, ε_Y) * (1.0 - Y_B_val)
end

"""
    R_rec_B(μ_B_T::Real, B::Real, B_min::Real)

Bacterial death (recycling) rate.

# Formula
    R_rec_B = μ_B(T) · B · h_B

# Arguments
- `μ_B_T`: Mortality rate at current T [1/day]
- `B`: Bacterial biomass [μg/mm³]
- `B_min`: Minimum viable biomass [μg/mm³]

# Returns
- Death rate [μg-C/mm³/day]

# Notes
- μ_B_T already includes Arrhenius factor
- h_B sigmoid ensures death rate vanishes as B → B_min
- Returns carbon to DOC pool (contributes to R_rec in S_C)
"""
function R_rec_B(μ_B_T::Real, B::Real, B_min::Real)
    h_B_val = h_B(B, B_min)
    μ_B_T * B * h_B_val
end
