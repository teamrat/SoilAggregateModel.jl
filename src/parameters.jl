"""
    parameters.jl

Parameter structs for biological and soil properties.

All parameters follow the architecture specification (ARCHITECTURE_CLAUDE_CODE.md).
Units: μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol throughout.
"""

#═══════════════════════════════════════════════════════════════════════════════
# Biological Properties
#═══════════════════════════════════════════════════════════════════════════════

"""
    BiologicalProperties

Biological and biogeochemical parameters.

All rate constants are at reference temperature T_ref.
Activation energies control temperature dependence via Arrhenius.

Units: μg/mm³, mm, days, kPa, K, J/mol
"""
struct BiologicalProperties
    # --- Bacterial ---
    r_B_max::Float64        # Max specific uptake rate at T_ref [1/day]
    K_B::Float64            # Half-saturation for DOC [μg/mm³]
    L_B::Float64            # Half-saturation for O₂ [μg/mm³]
    ν_B::Float64            # Water potential sensitivity [1/kPa]
    Y_B_max::Float64        # Maximum growth yield [-]
    K_Y::Float64            # Half-saturation for yield [μg-C/mm³/day]
    γ::Float64              # EPS allocation fraction [-]
    C_B::Float64            # Basal carbon requirement [μg/mm³]
    μ_B::Float64            # Mortality rate at T_ref [1/day]
    B_min::Float64          # Minimum viable biomass [μg/mm³]
    Ea_B::Float64           # Activation energy [J/mol]

    # --- Fungal ---
    r_F_max::Float64        # Max specific uptake rate at T_ref [1/day]
    K_F::Float64            # Half-saturation for DOC [μg/mm³]
    L_F::Float64            # Half-saturation for O₂ [μg/mm³]
    ν_F::Float64            # Water potential sensitivity [1/kPa]
    Y_F::Float64            # Growth yield [-]
    μ_F::Float64            # Mortality rate at T_ref [1/day]
    F_i_min::Float64        # Minimum viable insulated biomass [μg/mm³]
    Ea_F::Float64           # Activation energy [J/mol] — shared by ALL fungal rates

    # --- Fungal transitions ---
    α_i::Float64            # Mobilization rate, insulated [1/day]
    α_n::Float64            # Mobilization rate, non-insulated [1/day]
    β_i::Float64            # Immobilization rate, insulated [1/day]
    β_n::Float64            # Immobilization rate, non-insulated [1/day]
    delta::Float64          # Mobilization exponent (δ > 1) [-]
    η_conv::Float64         # Conversion efficiency [-]
    ζ::Float64              # Insulation rate F_n → F_i [1/day]
    λ::Float64              # Fraction of F_n at uptake surfaces [-]
    D_Fn0::Float64          # Hyphal extension diffusivity at T_ref [mm²/day]
    D_Fm0::Float64          # Internal translocation rate at T_ref [mm²/day]
    ε_F::Float64            # Π denominator protection [μg/mm³]

    # --- EPS ---
    μ_E_max::Float64        # Max EPS degradation rate at T_ref [1/day]
    K_E::Float64            # Substrate inhibition concentration [μg/mm³]
    E_min::Float64          # Minimum EPS for h_E sigmoid [μg/mm³]
    Ea_EPS::Float64         # Activation energy [J/mol]

    # --- MAOC ---
    κ_s_ref::Float64        # Sorption rate at T_ref [1/day]
    κ_d_ref::Float64        # Desorption rate at T_ref [1/day]
    Ea_MAOC_sorb::Float64   # Activation energy, sorption [J/mol]
    Ea_MAOC_desorb::Float64 # Activation energy, desorption [J/mol]
    ε_maoc::Float64         # Softplus smoothing width [μg/mm³]

    # --- POM ---
    R_P_max::Float64        # Max dissolution rate at T_ref [μg-C/mm²/day]
    P_0::Float64            # Initial POM mass [μg-C]
    r_0::Float64            # POM radius [mm]
    θ_P::Float64            # Half-saturation water content for dissolution [-]
    L_P::Float64            # Half-saturation O₂ for dissolution [μg/mm³]
    K_B_P::Float64          # Half-saturation bacteria for dissolution [μg/mm³]
    K_F_P::Float64          # Half-saturation fungi for dissolution [μg/mm³]
    ρ_POM::Float64          # POM carbon density [μg-C/mm³]
    Ea_POM::Float64         # Activation energy [J/mol]

    # --- Oxygen ---
    α_O::Float64            # Respiratory quotient [μg-O₂/μg-C]

    # --- Reference ---
    T_ref::Float64          # Reference temperature [K]
end

"""
    BiologicalProperties(; kwargs...)

Construct BiologicalProperties with default values.

Default parameters are placeholders for testing. Actual values should come from
the manuscript or calibration.

# Keywords
All fields can be overridden via keyword arguments.
"""
function BiologicalProperties(;
    # Bacterial (defaults are placeholders)
    r_B_max = 5.0,
    K_B = 50.0,
    L_B = 0.1,
    ν_B = 0.01,
    Y_B_max = 0.5,
    K_Y = 1.0,
    γ = 0.3,
    C_B = 10.0,
    μ_B = 0.01,
    B_min = 0.1,
    Ea_B = 60_000.0,

    # Fungal
    r_F_max = 3.0,
    K_F = 100.0,
    L_F = 0.05,
    ν_F = 0.005,
    Y_F = 0.4,
    μ_F = 0.005,
    F_i_min = 0.1,
    Ea_F = 55_000.0,

    # Fungal transitions
    α_i = 0.1,
    α_n = 0.05,
    β_i = 0.2,
    β_n = 0.15,
    delta = 1.5,
    η_conv = 0.9,
    ζ = 0.01,
    λ = 0.05,           # Fraction at uptake surfaces (λ ≪ 1; FIXED: was 0.5)
    D_Fn0 = 0.01,
    D_Fm0 = 1.0,
    ε_F = 1e-10,

    # EPS
    μ_E_max = 0.1,
    K_E = 5.0,          # [μg/mm³] in C_aq units (FIXED: was 200.0 for total C scale)
    E_min = 0.1,
    Ea_EPS = 50_000.0,

    # MAOC
    κ_s_ref = 0.001,
    κ_d_ref = 0.0001,
    Ea_MAOC_sorb = 25_000.0,
    Ea_MAOC_desorb = 40_000.0,
    ε_maoc = 0.01,

    # POM
    R_P_max = 0.5,
    P_0 = 1000.0,
    r_0 = 0.1,
    θ_P = 0.1,
    L_P = 0.05,
    K_B_P = 1.0,
    K_F_P = 1.0,
    ρ_POM = 500.0,      # [μg-C/mm³] ≈ 0.5 g-C/cm³ (FIXED: was 1.0e6)
    Ea_POM = 60_000.0,

    # Oxygen
    α_O = 2.67,

    # Reference
    T_ref = 293.15
)
    BiologicalProperties(
        r_B_max, K_B, L_B, ν_B, Y_B_max, K_Y, γ, C_B, μ_B, B_min, Ea_B,
        r_F_max, K_F, L_F, ν_F, Y_F, μ_F, F_i_min, Ea_F,
        α_i, α_n, β_i, β_n, delta, η_conv, ζ, λ, D_Fn0, D_Fm0, ε_F,
        μ_E_max, K_E, E_min, Ea_EPS,
        κ_s_ref, κ_d_ref, Ea_MAOC_sorb, Ea_MAOC_desorb, ε_maoc,
        R_P_max, P_0, r_0, θ_P, L_P, K_B_P, K_F_P, ρ_POM, Ea_POM,
        α_O,
        T_ref
    )
end

#═══════════════════════════════════════════════════════════════════════════════
# Soil Properties
#═══════════════════════════════════════════════════════════════════════════════

"""
    SoilProperties

Soil physical and chemical properties.

Units: μg/mm³, mm, days, kPa, K, J/mol
"""
struct SoilProperties
    # Van Genuchten water retention
    θ_r::Float64            # Residual water content [-]
    θ_s::Float64            # Saturated water content [-]
    α_vg::Float64           # van Genuchten α [1/kPa]
    n_vg::Float64           # van Genuchten n [-]

    # EPS/fungi modification of water retention
    ω_E::Float64            # EPS effect on α (negative) [mm³/μg]
    ω_F::Float64            # Fungi effect on α (negative) [mm³/μg]

    # Equilibrium sorption
    k_d_eq::Float64         # Linear partition coefficient [mm³/μg]
    ρ_b::Float64            # Bulk density [μg/mm³]

    # MAOC capacity (Langmuir-Freundlich)
    M_max::Float64          # Maximum sorption capacity [μg/mm³]
    k_L::Float64            # Langmuir affinity [mm³/μg]
    n_LF::Float64           # Freundlich exponent [-]
    k_ma::Float64           # Mineral activity coefficient [μg-C/g-mineral]
    f_clay_silt::Float64    # Clay+silt mass fraction [-]

    # Reference diffusion at T_ref [mm²/day]
    D_C0_ref::Float64       # DOC in water
    D_O2_w_ref::Float64     # O₂ in water
    D_O2_a_ref::Float64     # O₂ in air
    D_B_rel::Float64        # Bacterial motility relative to D_C [-]

    # Aggregate stability
    k_F::Float64            # Specific binding strength [kPa/(μg/mm³)]
    χ::Float64              # Particle adhesion length scale [mm]
    a_p::Float64            # Particle radius [mm]
end

"""
    SoilProperties(; kwargs...)

Construct SoilProperties with default values.

Default parameters are typical for a medium-textured soil. Override via keyword arguments.

# Keywords
All fields can be overridden via keyword arguments.
"""
function SoilProperties(;
    # Van Genuchten (typical loamy soil)
    θ_r = 0.05,
    θ_s = 0.45,
    α_vg = 0.01,        # [1/kPa]
    n_vg = 1.4,

    # EPS/fungi effects
    ω_E = -0.001,
    ω_F = -0.0005,

    # Equilibrium sorption
    k_d_eq = 0.01,      # [mm³/μg] (FIXED: was 0.5 → retardation ~30 instead of ~1460)
    ρ_b = 1.3e3,        # [μg/mm³] = 1.3 g/cm³ (FIXED: was 1.3e6)

    # MAOC capacity
    M_max = 100.0,      # [μg/mm³]
    k_L = 0.01,         # [mm³/μg]
    n_LF = 0.8,
    k_ma = 0.05,        # [μg-C/g-mineral]
    f_clay_silt = 0.5,

    # Reference diffusion at 293.15 K
    # Convert from literature values (cm²/s) to mm²/day: multiply by 8.64e6
    D_C0_ref = 1.0e-5 * 8.64e6,   # DOC: ~1e-5 cm²/s → 0.864 mm²/day
    D_O2_w_ref = 2.0e-5 * 8.64e6, # O₂ in water: ~2e-5 cm²/s → 1.728 mm²/day
    D_O2_a_ref = 0.2 * 8.64e6,    # O₂ in air: ~0.2 cm²/s → 1728 mm²/day
    D_B_rel = 0.001,              # Bacterial motility << DOC diffusion

    # Aggregate stability
    k_F = 1.0,          # [kPa/(μg/mm³)]
    χ = 0.001,          # [mm]
    a_p = 0.01          # [mm]
)
    SoilProperties(
        θ_r, θ_s, α_vg, n_vg,
        ω_E, ω_F,
        k_d_eq, ρ_b,
        M_max, k_L, n_LF, k_ma, f_clay_silt,
        D_C0_ref, D_O2_w_ref, D_O2_a_ref, D_B_rel,
        k_F, χ, a_p
    )
end
