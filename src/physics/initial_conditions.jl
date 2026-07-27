# initial_conditions.jl
# SOC partitioning into reactive pools for model initialization
#
# See soc_partitioning.md for full derivation and manuscript text.

"""
    InitialConditions

Describes the initial state of a soil at the start of an experiment.

These are properties of the **soil at t = 0**, not model physics parameters.
They change between experiments (e.g., de Gryze's 5 soils) but not during
a simulation.

# Fields
- `SOC`: Total soil organic carbon mass fraction [-] (e.g., 0.0183 for 1.83%)
- `s_M`: MAOC saturation ratio [-] (fraction of mineral capacity occupied)
- `f_bact`: Bacterial C as fraction of total SOC [-]
- `f_fungi`: Fungal C as fraction of total SOC [-]
- `f_eps`: EPS C as fraction of total SOC [-]
- `T_0`: Initial temperature [K]
- `ψ_0`: Initial matric potential [kPa]
- `O2_gas`: Atmospheric O₂ gas-phase density [µg/mm³]

# Notes
- SOC is mass-based (µg-C per g-soil), converted internally via ρ_b
- s_M uses texture information through M_max = k_ma × f_clay_silt × ρ_b
- Microbial fractions (f_bact, f_fungi, f_eps) are literature defaults;
  override with soil-specific measurements when available
- Typical ranges: f_bact ∈ [0.005, 0.03], f_fungi ∈ [0.005, 0.05],
  f_eps ∈ [0.001, 0.01], s_M ∈ [0.2, 0.8]

# Example
```julia
# De Gryze Soil 1: low SOC, sandy
ic1 = InitialConditions(SOC = 0.0125, s_M = 0.3)

# De Gryze Soil 3: high SOC, clayey
ic3 = InitialConditions(SOC = 0.0221, s_M = 0.5, f_fungi = 0.02)
```
"""
struct InitialConditions
    SOC::Float64            # Total SOC mass fraction [-] (e.g., 0.0183 for 1.83%)
    s_M::Float64            # MAOC saturation ratio [-]
    f_bact::Float64         # Bacterial C fraction of SOC [-]
    f_fungi::Float64        # Fungal C fraction of SOC [-]
    f_eps::Float64          # EPS C fraction of SOC [-]
    T_0::Float64            # Initial temperature [K]
    ψ_0::Float64            # Initial matric potential [kPa]
    O2_gas::Float64         # Atmospheric O₂ [µg/mm³]
end

"""
    InitialConditions(; kwargs...)

Construct InitialConditions with default values.

# Keywords
All fields can be overridden. Only `SOC` is typically required per soil.
"""
function InitialConditions(;
    SOC = 0.015,            # [-] 1.5% SOC, moderate agricultural soil
    s_M = 0.4,              # 40% of mineral capacity occupied
    f_bact = 0.01,          # 1% of SOC in bacteria
    f_fungi = 0.01,         # 1% of SOC in fungi
    f_eps = 0.005,          # 0.5% of SOC in EPS
    T_0 = 293.15,           # 20°C
    ψ_0 = -29.0,            # field capacity
    O2_gas = 0.2785         # ~21% atmospheric O₂
)
    InitialConditions(SOC, s_M, f_bact, f_fungi, f_eps, T_0, ψ_0, O2_gas)
end

"""
    partition_CM(R, s_M, M_max, k_L, k_d, θ, ρ_b, n_LF;
                 tol=1e-12, maxiter=50)

Partition residual SOC between DOC (C) and MAOC (M) at thermodynamic equilibrium.

Given a residual R = SOC_vol − (biotic pools), find C and M such that:
  1. C + M = R                                        (mass balance)
  2. M = s_M × M_max × f_LF(β × C)                   (isotherm equilibrium)

where β = k_L × k_d / (θ + ρ_b × k_d) is the lumped isotherm parameter that
converts total DOC concentration C into the equivalent aqueous → sorbed
concentration seen by the mineral surfaces, and f_LF(x) = x^n_LF / (1 + x^n_LF)
is the Langmuir-Freundlich isotherm function.

# Two regimes

The system has qualitatively different behavior depending on whether mineral
surfaces can absorb all available carbon:

**DOC-limited regime** (M_eq(R) < R): The isotherm evaluated at C = R (maximum
possible DOC) returns M_eq < R. There is more carbon than minerals can hold.
The solution has substantial DOC and M < s_M × M_max. Forward fixed-point
iteration converges reliably: guess M → compute C = R − M → evaluate M_eq →
update M. This is a contraction mapping when M_eq < R.

**Mineral-limited regime** (M_eq(R) ≥ R): The mineral capacity (s_M × M_max)
exceeds the available carbon R. Nearly all carbon should be mineral-associated,
with only trace DOC needed to maintain equilibrium. The correct solution has
M ≈ R and C ≪ R.

Forward iteration *fails* in this regime because:
  1. Guess M = R/2 → C = R/2 → isotherm wants M ≫ R → cap M at R
  2. C → C_floor → isotherm wants M ≈ 0 → M collapses
  3. Oscillation between extremes; no convergence

Instead, we solve for C directly using bisection on the residual equation:
  g(C) = C + s_M × M_max × f_LF(β × C) − R = 0

This is robust because g(C) is strictly monotonically increasing in C
(both C and M_eq increase with C), guaranteeing a unique root. Bisection
on C ∈ [C_floor, R] always converges.

# Arguments
- `R`: Residual carbon to partition [µg/mm³]
- `s_M`: MAOC saturation ratio [-] (fraction of capacity, 0 < s_M ≤ 1)
- `M_max`: Maximum MAOC capacity [µg/mm³] (= k_ma × f_clay_silt × ρ_b)
- `k_L`: Langmuir affinity constant [mm³/µg]
- `k_d`: Solid-liquid partition coefficient (equilibrium sorption) [mm³/µg]
- `θ`: Water content [-]
- `ρ_b`: Bulk density [µg/mm³]
- `n_LF`: Freundlich exponent [-] (n_LF < 1 for heterogeneous sites)

# Keyword Arguments
- `tol`: Convergence tolerance [µg/mm³] (default: 1e-12)
- `maxiter`: Maximum iterations (default: 50)

# Returns
Named tuple `(C=..., M=..., regime=...)` where:
- `C`: Total DOC concentration [µg/mm³]
- `M`: MAOC concentration [µg/mm³]
- `regime`: `:DOC_limited` or `:mineral_limited` (diagnostic)

# Notes
- C + M = R is enforced exactly (M is always computed as R − C)
- The returned M satisfies the isotherm to within `tol`
- Prints a warning if the achieved MAOC saturation differs significantly
  from the requested s_M (indicates mineral-limited regime where full
  saturation is impossible given available carbon)

# Example
```julia
# Typical silt loam: plenty of mineral capacity
R = 27.9       # µg-C/mm³ residual SOC
M_max = 461.8  # from texture
s_M = 0.6      # 60% saturation target
CM = partition_CM(R, 0.6, 461.8, 1000.0, 0.05, 0.35, 1300.0, 0.7)
# → mineral-limited: M ≈ 27.5, C ≈ 0.4 (nearly all SOC is MAOC)
# → actual saturation ≈ 27.5/461.8 = 6.0% (well below target — not enough C)
```
"""
function partition_CM(R::Real, s_M::Real, M_max::Real,
                      k_L::Real, k_d::Real, θ::Real, ρ_b::Real, n_LF::Real;
                      tol::Real=1e-12, maxiter::Int=50)
    C_floor = 1e-6  # minimum DOC [µg/mm³]

    # Lumped isotherm parameter: converts total C to effective concentration
    # seen by mineral surfaces
    β = k_L * k_d / (θ + ρ_b * k_d)

    # --- Isotherm evaluation helper ---
    function M_eq(C::Real)
        βC = β * C
        if βC > 0.0
            βC_n = βC^n_LF
            return s_M * M_max * βC_n / (1.0 + βC_n)
        else
            return 0.0
        end
    end

    # --- Regime detection ---
    # Evaluate isotherm at C = R (maximum possible DOC).
    # If M_eq(R) < R, there is excess DOC beyond what minerals can hold.
    # If M_eq(R) ≥ R, minerals want more carbon than exists → mineral-limited.
    M_eq_at_R = M_eq(R)

    if M_eq_at_R < R - C_floor
        # ═══ DOC-limited regime ═══
        # Forward fixed-point iteration: M_{k+1} = M_eq(R - M_k)
        # Contraction mapping guaranteed when M_eq(R) < R.
        M = min(0.5 * R, s_M * M_max)

        for iter in 1:maxiter
            C = R - M
            if C < C_floor
                C = C_floor
            end
            M_new = M_eq(C)
            M_new = min(M_new, R - C_floor)

            if abs(M_new - M) < tol
                C_final = R - M_new
                return (C = C_final, M = M_new, regime = :DOC_limited)
            end
            M = M_new
        end

        # Return best estimate after maxiter
        C_final = R - M
        @warn "partition_CM: DOC-limited iteration did not converge after $maxiter steps (|ΔM| > $tol)"
        return (C = C_final, M = M, regime = :DOC_limited)

    else
        # ═══ Mineral-limited regime ═══
        # Nearly all carbon goes to MAOC. Solve for the small C that
        # satisfies mass balance + isotherm simultaneously.
        #
        # Root-finding on g(C) = C + M_eq(C) − R = 0
        # g is strictly increasing (both C and M_eq increase with C),
        # so bisection on [C_floor, R] is guaranteed to converge.

        C_lo = C_floor
        C_hi = R

        # Verify bracket: g(C_lo) should be < 0, g(C_hi) should be ≥ 0
        g_lo = C_lo + M_eq(C_lo) - R
        g_hi = C_hi + M_eq(C_hi) - R
        # g_hi ≥ 0 is guaranteed since M_eq(R) ≥ R - C_floor → g(R) ≥ 0
        # g_lo < 0 is guaranteed since C_floor + M_eq(C_floor) ≪ R

        for iter in 1:maxiter
            C_mid = 0.5 * (C_lo + C_hi)
            g_mid = C_mid + M_eq(C_mid) - R

            if abs(g_mid) < tol || (C_hi - C_lo) < tol
                M_final = R - C_mid
                M_eq_full = M_eq(C_mid) / s_M  # full equilibrium (undo s_M scaling)
                s_actual = M_final / M_eq_full
                if abs(s_actual - s_M) > 0.01 * s_M
                        @info "partition_CM: mineral-limited regime. " *
                          "M=$(round(M_final, sigdigits=4)) µg/mm³, " *
                          "M_eq=$(round(M_eq_full, sigdigits=4)) µg/mm³, " *
                          "saturation=$(round(s_actual * 100, digits=1))% " *
                          "(target=$(round(s_M * 100, digits=1))%)"                end
                return (C = C_mid, M = M_final, regime = :mineral_limited)
            end

            if g_mid < 0.0
                C_lo = C_mid
            else
                C_hi = C_mid
            end
        end

        # Return best estimate
        C_mid = 0.5 * (C_lo + C_hi)
        M_final = R - C_mid
        @warn "partition_CM: mineral-limited bisection did not converge after $maxiter steps"
        return (C = C_mid, M = M_final, regime = :mineral_limited)
    end
end

"""
    create_initial_state(n::Int, bio::BiologicalProperties,
                         soil::SoilProperties, ic::InitialConditions;
                         P_0::Union{Nothing,Real}=nothing,
                         ω::Real=1.0)

Create initial state by partitioning measured SOC into model pools.

# Algorithm (see soc_partitioning.md)
1. Convert SOC to volumetric: SOC_vol = SOC × ρ_b (fraction × µg/mm³ = µg-C/mm³)
2. Biotic pools: B, F_i, E from prescribed fractions of SOC
3. F_n, F_m seeded at minimum viable concentrations
4. MAOC capacity: M_max from texture (k_ma × f_clay_silt × ρ_b)
5. C and M from iterative partition (thermodynamically consistent)
6. O₂ from gas-total conversion
7. Dilute concentrations by 1/ω for overlapping model domains
8. Validation checks

# Domain overlap dilution (ω)
When the model domain (f_domain) exceeds the physical packing cell (f_pack),
adjacent domains overlap by a factor ω = (f_domain/f_pack)³. The physical
SOC is partitioned at true concentrations (steps 1–6), then all pool
concentrations are divided by ω so that each model domain contains exactly
its proportionate share of background carbon. POM mass (P_0) is NOT diluted
— it is a lumped scalar localized at the center.

See domain_tessellation_supplemental.md for derivation and error analysis.

# Arguments
- `n::Int`: Number of grid points
- `bio::BiologicalProperties`: Biological parameters
- `soil::SoilProperties`: Soil properties (texture, sorption, water retention)
- `ic::InitialConditions`: SOC and environmental conditions at t = 0
- `P_0::Union{Nothing,Real}`: Initial POM mass [µg-C] (default: bio.P_0)
- `ω::Real`: Domain overlap factor (default: 1.0, no dilution)

# Returns
- `AggregateState`: Initialized state consistent with measured SOC

# Example
```julia
bio  = BiologicalProperties()
soil = SoilProperties(f_clay_silt = 0.42)
ic   = InitialConditions(SOC = 0.0183, s_M = 0.4)
# Heavy amendment: f_pack = 3.2, f_domain = 10 → ω = 30.5
state = create_initial_state(200, bio, soil, ic; P_0 = 500.0, ω = 30.5)
```
"""
function create_initial_state(n::Int, bio::BiologicalProperties,
                              soil::SoilProperties, ic::InitialConditions;
                              P_0::Union{Nothing,Real}=nothing,
                              ω::Real=1.0)
    state = AggregateState(n)

    # === Step 1: Volumetric SOC ===
    # SOC [-] × ρ_b [µg/mm³] = SOC_vol [µg-C/mm³]
    SOC_vol = ic.SOC * soil.ρ_b  # e.g., 0.0183 × 1500 = 27.45 µg-C/mm³

    # === Step 2: Biotic pools (analytical) ===
    B_0   = ic.f_bact  * SOC_vol
    F_i_0 = ic.f_fungi * SOC_vol   # steady-state: all fungal C in F_i
    E_0   = ic.f_eps   * SOC_vol

    # F_n and F_m: seeded at minimum viable (not from SOC fraction)
    F_n_0 = bio.F_n_min
    F_m_0 = bio.F_m_min

    # === Step 3: Water content at ψ_0 ===
    m_vg = 1.0 - 1.0 / soil.n_vg
    θ_0 = soil.θ_r + (soil.θ_s - soil.θ_r) *
          (1.0 + (-soil.α_vg * ic.ψ_0)^soil.n_vg)^(-m_vg)
    θ_a_0 = soil.θ_s - θ_0

    # === Step 4: Residual for C + M ===
    biotic_total = B_0 + F_i_0 + F_n_0 + F_m_0 + E_0
    R = SOC_vol - biotic_total

    # Validate: residual must be positive
    if R ≤ 0.0
        error("SOC partitioning failed: biotic fractions " *
              "(f_bact=$(ic.f_bact) + f_fungi=$(ic.f_fungi) + f_eps=$(ic.f_eps) " *
              "= $(ic.f_bact + ic.f_fungi + ic.f_eps)) exceed 1.0, or SOC too low.")
    end



    # === Step 5: Iterative C–M partition ===
    # M_max from texture (may differ from soil.M_max if f_clay_silt overridden)
    M_max = soil.k_ma * soil.f_clay_silt * soil.ρ_b

    CM = partition_CM(R, ic.s_M, M_max,
                      soil.k_L, soil.k_d_eq, θ_0, soil.ρ_b, soil.n_LF)

    C_0 = CM.C
    M_0 = CM.M

    # Diagnostic: report MAOC saturation
    β_diag = soil.k_L * soil.k_d_eq / (θ_0 + soil.ρ_b * soil.k_d_eq)
    βC_diag = β_diag * C_0
    βC_n_diag = βC_diag^soil.n_LF
    M_eq_full_diag = M_max * βC_n_diag / (1.0 + βC_n_diag)
    s_M_actual = M_0 / M_eq_full_diag
    @info "SOC partition: regime=$(CM.regime), " *
          "C=$(round(C_0, sigdigits=4)) µg/mm³, " *
          "M=$(round(M_0, sigdigits=4)) µg/mm³, " *
          "MAOC saturation=$(round(s_M_actual * 100, digits=1))% " *
          "(target=$(round(ic.s_M * 100, digits=1))%)"




    # === Step 6: Oxygen (gas → total soil) ===
    K_H = K_H_O2(ic.T_0)
    O_aq = ic.O2_gas / K_H

    # === Step 7: Validation at physical concentrations ===
    total_check = C_0 + B_0 + F_i_0 + F_n_0 + F_m_0 + E_0 + M_0
    rel_err = abs(total_check - SOC_vol) / SOC_vol
    if rel_err > 1e-10
        @warn "SOC partition conservation error: $(rel_err)"
    end
    if C_0 < 1e-6
        @warn "Very low initial DOC ($(C_0) µg/mm³) — SOC may be over-allocated to MAOC/biota"
    end
    if M_0 > M_max
        @warn "Initial MAOC ($(M_0)) exceeds capacity ($(M_max))"
    end

    # === Step 8: Domain overlap dilution ===
    # Divide all concentrations by ω so each model domain contains
    # exactly its proportionate share of background carbon.
    # POM (P_0) is NOT diluted — it is a lumped scalar at the center.
    # O₂ is NOT diluted — it is a boundary condition, not a carbon pool.
    if ω != 1.0
        C_0   /= ω
        B_0   /= ω
        F_i_0 /= ω
        F_n_0 /= ω
        F_m_0 /= ω
        E_0   /= ω
        M_0   /= ω
    end

    # === Step 9: Fill state ===
    state.C   .= C_0
    state.B   .= B_0
    state.F_i .= F_i_0
    state.F_n .= F_n_0
    state.F_m .= F_m_0
    state.E   .= E_0
    state.M   .= M_0
    state.O .= O_aq

    # POM (not diluted in overlapping domains)
    state.P = isnothing(P_0) ? bio.P_0 : Float64(P_0)
    state.P_0 = state.P       
    state.CO2_cumulative = 0.0

    return state
end
