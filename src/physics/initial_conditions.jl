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
- `f_bact`: Bacterial C as fraction of total SOC [-]
- `f_fungi`: Fungal C as fraction of total SOC [-]
- `f_insulated`: Fraction of native fungal C that starts insulated (F_i) [-];
  the remainder starts as non-insulated hyphae (F_n)
- `f_eps`: EPS C as fraction of total SOC [-]
- `T_0`: Initial temperature [K]
- `ψ_0`: Initial matric potential [kPa]
- `O2_gas`: Atmospheric O₂ gas-phase density [µg/mm³]

# Notes
- SOC is mass-based (µg-C per g-soil), converted internally via ρ_b
- MAOC follows from sorption equilibrium with `maoc_capacity(soil)`; its
  saturation `M/M_max` is an OUTPUT, not an input
- Microbial fractions (f_bact, f_fungi, f_eps) are literature defaults;
  override with soil-specific measurements when available
- Typical ranges: f_bact ∈ [0.005, 0.03], f_fungi ∈ [0.005, 0.05],
  f_eps ∈ [0.001, 0.01], f_insulated ∈ [0.2, 0.8]
- `f_insulated = 1.0` reproduces the pre-2026-07 behaviour exactly (all fungal
  C in F_i, F_n seeded at bio.F_n_min) and is useful as a control

# Example
```julia
# De Gryze Soil 1: low SOC, sandy
ic1 = InitialConditions(SOC = 0.0125)

# De Gryze Soil 3: high SOC, clayey
ic3 = InitialConditions(SOC = 0.0221, f_fungi = 0.02)
```
"""
struct InitialConditions
    SOC::Float64            # Total SOC mass fraction [-] (e.g., 0.0183 for 1.83%)
    f_bact::Float64         # Bacterial C fraction of SOC [-]
    f_fungi::Float64        # Fungal C fraction of SOC [-]
    f_insulated::Float64    # Fraction of fungal C starting as F_i [-]
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
    f_bact = 0.01,          # 1% of SOC in bacteria
    f_fungi = 0.01,         # 1% of SOC in fungi
    f_insulated = 0.5,      # half of native fungal C starts insulated
    f_eps = 0.005,          # 0.5% of SOC in EPS
    T_0 = 293.15,           # 20°C
    ψ_0 = -29.0,            # field capacity
    O2_gas = 0.2785         # ~21% atmospheric O₂
)
    if !(0.0 ≤ f_insulated ≤ 1.0)
        error("f_insulated must lie in [0, 1]; got $(f_insulated)")
    end
    InitialConditions(SOC, f_bact, f_fungi, f_insulated, f_eps,
                      T_0, ψ_0, O2_gas)
end

"""
    partition_CM(SOC_residual, M_max, k_L, k_d, θ, ρ_b, n_LF;
                 tol=1e-12, maxiter=50)

Partition residual SOC between DOC (C) and MAOC (M) at thermodynamic equilibrium.

Given `SOC_residual = SOC_vol − (biotic pools)`, find C and M such that:
  1. C + M = R                                        (mass balance)
  2. M = M_max × f_LF(β × C)                         (sorption equilibrium)

where β = k_L × k_d / (θ + ρ_b × k_d) is the lumped isotherm parameter that
converts total DOC concentration C into the equivalent aqueous → sorbed
concentration seen by the mineral surfaces, and f_LF(x) = x^n_LF / (1 + x^n_LF)
is the Langmuir-Freundlich isotherm function.

# Two regimes

The system has qualitatively different behavior depending on whether mineral
surfaces can absorb all available carbon:

**DOC-limited regime** (M_eq(R) < R): The isotherm evaluated at C = R (maximum
possible DOC) returns M_eq < R. There is more carbon than minerals can hold.
The solution has substantial DOC and M < M_max. Forward fixed-point
iteration converges reliably: guess M → compute C = R − M → evaluate M_eq →
update M. This is a contraction mapping when M_eq < R.

**Mineral-limited regime** (M_eq(R) ≥ R): The mineral capacity `M_max`
exceeds the available carbon R. Nearly all carbon should be mineral-associated,
with only trace DOC needed to maintain equilibrium. The correct solution has
M ≈ R and C ≪ R.

Forward iteration *fails* in this regime because:
  1. Guess M = R/2 → C = R/2 → isotherm wants M ≫ R → cap M at R
  2. C → C_floor → isotherm wants M ≈ 0 → M collapses
  3. Oscillation between extremes; no convergence

Instead, we solve for C directly using bisection on the residual equation:
  g(C) = C + M_max × f_LF(β × C) − SOC_residual = 0

This is robust because g(C) is strictly monotonically increasing in C
(both C and M_eq increase with C), guaranteeing a unique root. Bisection
on C ∈ [C_floor, R] always converges.

# Arguments
- `R`: Residual carbon to partition [µg/mm³]
- `M_max`: Maximum MAOC capacity [µg/mm³] — pass `maoc_capacity(soil)`
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
- **Saturation is an OUTPUT, not an input.** `M/M_max` is whatever sorption
  equilibrium with the available carbon produces. As the affinity `k_L` rises,
  `C → 0` and the saturation approaches `SOC_residual/M_max` — set entirely by
  measured SOC and measured texture. That makes it a prediction testable against
  reported saturation deficits (Georgiou et al. 2022), not a knob.

  Until 2026-07-29 an `s_M` argument scaled the isotherm, so the state was placed
  at `s_M` of local equilibrium and the solver — which has no `s_M` — sorbed from
  t = 0 to close the gap. Two isotherms for one quantity; see REFERENCE §26
  erratum 14.

- **On the name.** `SOC_residual`, not `R`. The manuscript uses a script R for
  total respiration and `R_*` is a rate everywhere in `src/`; a carbon stock
  called `R` collided with both.

# Example
```julia
# Typical silt loam: plenty of mineral capacity
SOC_residual = 27.9   # µg-C/mm³, SOC_vol minus the biotic pools
M_max        = 46.18  # maoc_capacity: 0.048 × 0.74 × 1300
CM = partition_CM(SOC_residual, M_max, 1000.0, 0.05, 0.35, 1300.0, 0.7)
# → mineral-limited: nearly all SOC becomes MAOC; CM.M/M_max is the
#   predicted saturation
```
"""
function partition_CM(SOC_residual::Real, M_max::Real,
                      k_L::Real, k_d::Real, θ::Real, ρ_b::Real, n_LF::Real;
                      tol::Real=1e-12, maxiter::Int=50)
    R = SOC_residual        # local shorthand; see the docstring on the name
    C_floor = 1e-6  # minimum DOC [µg/mm³]

    # Lumped isotherm parameter: converts total C to effective concentration
    # seen by mineral surfaces
    β = k_L * k_d / sorption_capacity(θ, ρ_b, k_d)

    # --- Isotherm evaluation helper ---
    # M_eq_langmuir_freundlich, not a second copy of the isotherm. β lumps
    # k_L·k_d/(θ + ρ_b·k_d) so that β·C is the k_L·C_eq the primitive expects,
    # and passing k_L = 1 avoids applying the affinity twice. CLAUDE.md §8.
    M_eq(C::Real) = C > 0.0 ?
        M_eq_langmuir_freundlich(β * C, M_max, 1.0, n_LF) : 0.0

    # --- Regime detection ---
    # Evaluate isotherm at C = R (maximum possible DOC).
    # If M_eq(R) < R, there is excess DOC beyond what minerals can hold.
    # If M_eq(R) ≥ R, minerals want more carbon than exists → mineral-limited.
    M_eq_at_R = M_eq(R)

    if M_eq_at_R < R - C_floor
        # ═══ DOC-limited regime ═══
        # Forward fixed-point iteration: M_{k+1} = M_eq(R - M_k)
        # Contraction mapping guaranteed when M_eq(R) < R.
        M = min(0.5 * R, M_max)

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
3. Fungal C split between F_i and F_n by ic.f_insulated; F_m seeded at
   minimum viable concentration
4. MAOC capacity: `maoc_capacity(soil)` = k_ma × f_clay_silt × ρ_b
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
ic   = InitialConditions(SOC = 0.0183)
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
    E_0   = ic.f_eps   * SOC_vol

    # Native fungal C is split between the insulated (F_i) and non-insulated
    # (F_n) pools by ic.f_insulated.
    #
    # This split is a FREE INITIAL CONDITION, not a steady state. With
    # F_m ≈ F_m_min the protection ratio Π = F_m/(F_i + F_n + ε_F) ≈ 0, so
    # every transition rate vanishes and *any* partition is stationary. The
    # previous code put 100% of fungal C in F_i and called it "steady state";
    # that is one arbitrary choice among many, and it is the worst one:
    #
    #   1. Uptake is throttled. R_F ∝ (λ·F_i + F_n) with λ = 0.05, so all the
    #      biomass sits in the pool with 5% uptake efficiency.
    #   2. Π is suppressed from below. F_i is in the denominator of Π, so
    #      loading F_i actively closes the only escape route.
    #   3. F_n cannot bootstrap. Its only source term,
    #         S_Fn = (1-ζ)·η·(α_n·Π^δ - β_n·Π)·F_n
    #      is proportional to F_n itself, making F_n ≈ 0 an absorbing state.
    #
    # Together these pin F_n at its floor for the whole simulation: with the
    # defaults (α_n = 0.15, β_n = 0, δ = 1, η = 0.8, ζ = 0.2) the growth rate
    # is (1-ζ)·η·α_n·Π ≈ 0.096·Π, and Π ≈ 0.03 gives an e-folding time of
    # ~350 days against a 21-day incubation. No parameter choice within a
    # defensible range escapes this; the partition has to change.
    #
    # Set f_insulated = 1.0 to recover the old behaviour exactly.
    F_fungi_0 = ic.f_fungi * SOC_vol
    F_i_0 = ic.f_insulated * F_fungi_0
    F_n_0 = max((1.0 - ic.f_insulated) * F_fungi_0, bio.F_n_min)

    # F_m: seeded at minimum viable (not from SOC fraction)
    F_m_0 = bio.F_m_min

    # === Step 3: Water content at ψ_0 ===
    # van_genuchten, not an inline copy: the inline form raised a negative base
    # to the power n_vg at saturation (ψ_0 = 0, which van_genuchten_inverse
    # returns exactly) and threw DomainError. CLAUDE.md §8.
    θ_0 = van_genuchten(ic.ψ_0, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
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
    # The one definition (biology/maoc.jl). Formerly computed inline here while
    # the solver read a separate `soil.M_max` field — see erratum 12.
    M_max = maoc_capacity(soil)

    CM = partition_CM(R, M_max,
                      soil.k_L, soil.k_d_eq, θ_0, soil.ρ_b, soil.n_LF)

    C_0 = CM.C
    M_0 = CM.M

    # Diagnostic: the PREDICTED MAOC saturation. maxlog=1 — a diameter sweep
    # calls this once per aggregate with an identical initial condition, so it
    # printed ten identical lines per run and 26 per test suite. Nothing asserts
    # on this log, so capping it is safe. Compare with reported
    # saturation deficits (Georgiou et al. 2022) — this is a model output.
    # Ceiling: as k_L rises, C -> 0 and saturation -> R/M_max, which is fixed by
    # measured SOC and measured texture alone.
    @info "SOC partition: regime=$(CM.regime), " *
          "C=$(round(C_0, sigdigits=4)) µg/mm³, " *
          "M=$(round(M_0, sigdigits=4)) µg/mm³, " *
          "MAOC saturation=$(round(100 * M_0 / M_max, digits=1))% " *
          "(ceiling at this SOC and texture: $(round(100 * R / M_max, digits=1))%), " *
          "pore-water DOC=$(round(1000 * C_aqueous(C_0, θ_0, soil), sigdigits=3)) mg/L" maxlog = 1

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
    # NO-OP at ω = 1, which is the only consistent setting — see
    # tessellation.jl and REFERENCE.md §26 erratum 13. Retained because
    # f_domain_min can still be raised for experiments, but any ω > 1 leaves the
    # state diluted while every concentration constant it is compared against is
    # not, and no rescaling fixes that.
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
