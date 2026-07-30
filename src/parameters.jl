"""
    parameters.jl

Parameter structs for biological and soil properties.

All parameters follow the architecture specification (docs/ARCHITECTURE.md — stale in parts; docs/BACKLOG.md item 12).
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
    ε_Y::Float64            # Softplus smoothing width for yield transition [μg-C/mm³/day]
    γ::Float64              # EPS allocation fraction [-]
    C_B::Float64            # Basal carbon requirement [μg/mm³]
    μ_B::Float64            # Mortality rate at T_ref [1/day]
    B_min::Float64          # Minimum viable biomass [μg/mm³]
    B_S::Float64            # Half-saturation for space limitation on bacterial yield [μg/mm³]
    Ea_B::Float64           # Activation energy [J/mol]

    # --- Fungal ---
    r_F_max::Float64        # Max specific uptake rate at T_ref [1/day]
    K_F::Float64            # Half-saturation for DOC [μg/mm³]
    L_F::Float64            # Half-saturation for O₂ [μg/mm³]
    ν_F::Float64            # Water potential sensitivity [1/kPa]
    Y_F::Float64            # Base growth yield [-]
    μ_F::Float64            # Mortality rate at T_ref [1/day]
    F_i_min::Float64        # Minimum viable insulated biomass [μg/mm³]
    F_n_min::Float64        # Minimum viable insulated biomass [μg/mm³]
    F_m_min::Float64        # Minimum viable insulated biomass [μg/mm³]
    F_S::Float64            # Half-saturation for space limitation on fungal yield (total F_i+F_n+F_m) [μg/mm³]
    Ea_F::Float64           # Activation energy [J/mol] — shared by ALL fungal rates

    # --- Fungal transitions ---
    # Falconer (2005) eq. 2.4 / (2008) Box 1 naming — γ(α·π^θ − β·π)·b :
    #   α = IMMOBILIZATION (mobile → sessile; the GAIN; carries the exponent)
    #   β = MOBILIZATION   (sessile → mobile; the LOSS; linear in Π)
    # These two were inverted before 2026-07, which is how the exponent came
    # to sit on the loss term. See dev_notes/falconer_answers.md §B1, §B3.
    α_i::Float64            # Immobilization rate, insulated [1/day] (Falconer α_i)
    α_n::Float64            # Immobilization rate, non-insulated [1/day] (Falconer α_n)
    β_i::Float64            # Mobilization rate, insulated [1/day] (Falconer β_i)
    β_n::Float64            # Mobilization rate, non-insulated [1/day] (Falconer β_n)
    delta::Float64          # Immobilization exponent — Falconer's θ (θ > 1) [-]
    η_conv::Float64         # Conversion efficiency — Falconer's γ [-]
    ζ::Float64              # Insulation SPLITTING FRACTION on the F_n tendency [-] — NOT a rate,
                            # and not temperature-scaled (see fungi.jl)
    λ::Float64              # Fraction of F_n at uptake surfaces [-]
    D_Fn0::Float64          # Hyphal extension diffusivity at T_ref [mm²/day]
    D_Fm0::Float64          # Internal translocation rate at T_ref [mm²/day]
    ε_F::Float64            # Π denominator protection [μg/mm³]
    K_Fm_net::Float64    # Half-saturation for network-dependent F_m translocation [μg/mm³]

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
    # Bacterial (from REFERENCE.md)
    r_B_max = 2.0,
    K_B = 1.0e-4,
    L_B = 0.00129,
    ν_B = 5.8e-4,
    Y_B_max = 0.7,
    K_Y = 3.33e-4,      # Derived: 10 * r_B_max * C_B/(K_B+C_B) * B_min = 3.33 * 1e-4 [μg/mm³/day]
    ε_Y = 3.33e-6,      # Derived: K_Y / 100, softplus smoothing width [μg-C/mm³/day]
    γ = 0.2,
    C_B = 2.0e-5,       # Derived: K_B/5
    μ_B = 0.012,
    B_min = 1.0e-4,
    B_S = 0.5,          # [μg/mm³] = 0.5 kg/m³, half-saturation for space limitation
    Ea_B = 60_000.0,

    # Fungal (from REFERENCE.md)
    r_F_max = 2.0,
    K_F = 1.0e-4,
    L_F = 0.00129,
    ν_F = 7.58e-5,
    Y_F = 0.6,
    μ_F = 0.012,
    F_i_min = 1.0e-6,
    F_n_min = 1.0e-4,
    F_m_min = 1.0e-6,
    F_S = 0.2,          # [μg/mm³] = 0.2 kg/m³, half-saturation for space limitation (total fungi)
    Ea_F = 55_000.0,

    # Fungal transitions. Falconer naming: α = immobilization (gain, carries the
    # exponent), β = mobilization (loss, linear in Π).
    #
    # These values are the MATLAB reference's, which is Falconer with θ = 1 and
    # β = 0 — a strict special case. They are working assumptions, not choices
    # made for this model. Falconer's own published values, for calibration to
    # start from (dev_notes/falconer_answers.md §F1):
    #
    #            [2005]                      [2008]              here
    #   α_n      0.2 → 0.8                   0.87 throughout     0.15
    #   α_i      0.5 (0.1 in fig 2)          0.0 / 0.01 / 0.9    0.1
    #   β_n      0.8 → 0.2                   0.0 / 0.1 / 0.9     0.0
    #   β_i      0.5 (0.9 in fig 2)          0.0 / 0.34 / 0.9    0.0
    #   θ        1.0 control, 3.0 nonlinear  1.0 throughout      1.0
    #   ζ        0.01                        0.01                0.2
    #
    # θ = 3.0 is the only published value above 1 and is the case that produces
    # the ring/aggregation patterns ([2005] figs 1c–1f); θ < 1 is never used.
    # With β = 0 the exponent has no linear loss to cross over against, so δ
    # only rescales the transition rate — δ and β are one decision, not two.
    α_i = 0.1,      # immobilization to insulated     (MATLAB: alpha_i = 0.1)
    α_n = 0.15,     # immobilization to non-insulated (MATLAB: alpha_n = 0.15)
    β_i = 0.0,      # mobilization from insulated     (MATLAB: beta_i = 0)
    β_n = 0.0,      # mobilization from non-insulated (MATLAB: beta_n = 0)
    delta = 1.0,    # Falconer's θ
    η_conv = 0.8,   # Falconer's γ
    ζ = 0.2,        # insulation SPLITTING FRACTION, constant — Falconer [2005]
                    # p. 1728 "a constant rate, ζ", held at 0.01 in every
                    # simulation of both papers. Not temperature-scaled.
    λ = 0.05,
    D_Fn0 = 0.01,
    D_Fm0 = 1.0,
    ε_F = 1e-4,
    K_Fm_net = 20.0 * 1e-4,   # 2e-3; MATLAB: 10*(Fi_min + Fn_min)    

    # EPS (from REFERENCE.md)
    μ_E_max = 0.002,
    K_E = 0.001,        # Derived: 50 * C_B = 50 * 2e-5
    E_min = 1.0e-4,
    Ea_EPS = 50_000.0,

    # MAOC (from REFERENCE.md)
    κ_s_ref = 0.1,
    κ_d_ref = 0.01,
    Ea_MAOC_sorb = 25_000.0,
    Ea_MAOC_desorb = 40_000.0,
    ε_maoc = 0.01,

    # POM (from REFERENCE.md)
    R_P_max = 1.0,
    P_0 = 1000.0,
    r_0 = 0.1,
    θ_P = 0.1,
    L_P = 0.00129,
    K_B_P = 1.0e-3,
    K_F_P = 1.0e-3,
    ρ_POM = 200.0,
    Ea_POM = 60_000.0,

    # Oxygen
    α_O = 2.2,

    # Reference
    T_ref = 293.15
)
    BiologicalProperties(
        r_B_max, K_B, L_B, ν_B, Y_B_max, K_Y, ε_Y, γ, C_B, μ_B, B_min, B_S, Ea_B,
        r_F_max, K_F, L_F, ν_F, Y_F, μ_F, F_i_min,F_n_min,F_m_min, F_S, Ea_F,
        α_i, α_n, β_i, β_n, delta, η_conv, ζ, λ, D_Fn0, D_Fm0, ε_F, K_Fm_net,
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

    # MAOC capacity (Langmuir-Freundlich).
    # There is NO M_max field. The capacity is `maoc_capacity(soil)` =
    # k_ma·f_clay_silt·ρ_b, computed in ONE place (biology/maoc.jl).
    k_L::Float64            # Langmuir affinity [mm³/μg]
    n_LF::Float64           # Freundlich exponent [-]
    k_ma::Float64           # MAOC capacity per unit clay+silt [g-C/g-mineral, dimensionless]
    f_clay_silt::Float64    # Clay+silt mass fraction [-]

    # Reference diffusion at T_ref [mm²/day]
    D_C0_ref::Float64       # DOC in water
    D_O2_w_ref::Float64     # O₂ in water
    D_O2_a_ref::Float64     # O₂ in air
    D_B_rel::Float64        # Bacterial motility relative to D_C [-]

    # Aggregate stability (see docs/REFERENCE.md §4.4)
    κ_b::Float64            # Specific binding strength per unit binder C [Pa·mm/(μg/mm³)]
    w_E::Float64            # EPS binding weight relative to F_i, per unit C [-]
    d_32::Float64           # Sauter mean particle diameter [mm]
    p_Gc::Float64           # Size dependence of the threshold, G_c ∝ (r/δ_s)^p_Gc [-]
end

"""
    SoilProperties(; kwargs...)

Construct SoilProperties with default values.

Default parameters are for sandy loam from REFERENCE.md. Override via keyword arguments.

# Keywords
All fields can be overridden via keyword arguments.
"""
function SoilProperties(;
    # Van Genuchten (sandy loam from REFERENCE.md)
    θ_r = 0.06,
    θ_s = 0.5,
    α_vg = 0.1133,      # [1/kPa]
    n_vg = 1.47,

    # EPS/fungi effects
    ω_E = -0.001,
    ω_F = -0.0005,

    # Equilibrium sorption
    k_d_eq = 0.05,      # [mm³/μg]
    ρ_b = 1.5e3,        # [μg/mm³] = 1.5 g/cm³

    # No `M_max` field — the capacity is `maoc_capacity(soil)`. See erratum 12.
    k_L = 10.0,         # [mm³/μg]
    n_LF = 0.7,
    # Low-activity minerals (Georgiou et al. 2022). DIMENSIONLESS g-C/g-mineral.
    # Was 0.48 with the unit "µg-C/g-mineral" — the same 48 mg/g, with the mg->g
    # conversion dropped, so it read 10x high in a formula that needed a mass
    # fraction. Use K_MA_HIGH_ACTIVITY for high-activity mineralogy.
    k_ma = K_MA_LOW_ACTIVITY,
    f_clay_silt = 0.40,

    # Reference diffusion at 293.15 K
    D_C0_ref = 86.4,         # DOC: from REFERENCE.md
    D_O2_w_ref = 173.7,      # O₂ in water: from REFERENCE.md (Han-Bartels)
    D_O2_a_ref = 1.52e6,     # O₂ in air: from REFERENCE.md (Chapman-Enskog)
    D_B_rel = 0.001,         # Bacterial motility << DOC diffusion

    # Aggregate stability — τ_c = (κ_b/d_32)·(F_i + w_E·E) ≥ τ_w. See §4.4/§5a.
    # BOTH κ_b AND w_E ARE FITTED. Neither is derived; the bond-counting argument
    # fixes the FORM (strength ∝ binder concentration / particle size) and not the
    # coefficients. d_32 is a measured soil property — compute it with
    # `sauter_from_texture` and override it for every real problem.
    # 0.16 replaces 0.0143869 on 2026-07-29 (0.012 ... 0.15, 0.18). The old
    # value had no measurement behind it: it was `2.25 × d_32(soil 3)`,
    # reverse-engineered to reproduce a legacy `G_c = 0.0194` that came from
    # `k_F = 2.25`, whose own derivation was dimensionally inconsistent (§26
    # erratum 11). Carrying it forward preserved a number, not a result. 0.16
    # is likewise fitted — see §5a. NOTE `G_c = τ_w·d_32/κ_b`: raising κ_b
    # LOWERS the threshold.
    κ_b = 0.16,         # [Pa·mm/(μg/mm³)] fitted
    w_E = 0.5,          # [-] fitted; was hard-coded in aggregate_radius.jl before 2026-07
    d_32 = 0.01,        # [mm] NOMINAL placeholder — override per soil
    # Size dependence of the threshold: G_c(r) = G_c(δ_s)·(r/δ_s)^p_Gc.
    # p_Gc = 0 is a threshold independent of aggregate size, which is what the
    # model did before 2026-07-29 and is therefore the default: every existing
    # result stands unchanged. p_Gc > 0 makes a larger aggregate need more
    # binder to survive. FITTED — see §4.4a. Not derived, and not zero for any
    # derived reason either; zero is a default, not a finding.
    p_Gc = 0.0          # [-] fitted
)
    SoilProperties(
        θ_r, θ_s, α_vg, n_vg,
        ω_E, ω_F,
        k_d_eq, ρ_b,
        k_L, n_LF, k_ma, f_clay_silt,
        D_C0_ref, D_O2_w_ref, D_O2_a_ref, D_B_rel,
        κ_b, w_E, d_32, p_Gc
    )
end
