"""
Tests for physics/*.jl modules
"""

import SoilAggregateModel: TemperatureCache, AggregateState,
                           alpha_effective, van_genuchten, water_content, update_water_content!,
                           van_genuchten_inverse,
                           tortuosity_millington_quirk, D_eff_DOC, D_eff_bacteria,
                           D_eff_fungi_noninsulated, D_eff_fungi_mobile, D_eff_oxygen

@testset "Water retention" begin
    @testset "Alpha modification" begin
        soil = SoilProperties()

        # No EPS/fungi → α unchanged
        α_base = alpha_effective(0.0, 0.0, soil)
        @test α_base ≈ soil.α_vg

        # EPS increases retention (lowers α)
        α_eps = alpha_effective(10.0, 0.0, soil)
        @test α_eps < α_base  # ω_E < 0

        # Fungi increase retention
        α_fungi = alpha_effective(0.0, 5.0, soil)
        @test α_fungi < α_base  # ω_F < 0

        # Combined effect
        α_both = alpha_effective(10.0, 5.0, soil)
        @test α_both < α_eps
        @test α_both < α_fungi
    end

    @testset "Van Genuchten boundary conditions" begin
        soil = SoilProperties()

        # Saturated (ψ = 0)
        θ_sat = van_genuchten(0.0, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
        @test θ_sat ≈ soil.θ_s

        # Very dry (ψ → -∞)
        θ_dry = van_genuchten(-1e6, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
        @test θ_dry ≈ soil.θ_r rtol=0.2  # Asymptotic approach, allow 20%

        # Intermediate (field capacity ~-33 kPa)
        θ_fc = van_genuchten(-33.0, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
        @test soil.θ_r < θ_fc < soil.θ_s
    end

    @testset "Van Genuchten monotonicity" begin
        soil = SoilProperties()
        ψ_values = [-1500.0, -500.0, -100.0, -33.0, -10.0, -1.0, 0.0]

        θ_values = [van_genuchten(ψ, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
                    for ψ in ψ_values]

        # Wetter conditions (less negative ψ) → higher θ
        @test all(diff(θ_values) .> 0)
    end

    @testset "EPS/fungi increase retention" begin
        soil = SoilProperties()
        ψ = -33.0  # Field capacity

        θ_base = water_content(ψ, 0.0, 0.0, soil)
        θ_eps = water_content(ψ, 10.0, 0.0, soil)
        θ_fungi = water_content(ψ, 0.0, 5.0, soil)
        θ_both = water_content(ψ, 10.0, 5.0, soil)

        # Higher EPS/fungi → higher water content at same ψ
        @test θ_eps > θ_base
        @test θ_fungi > θ_base
        @test θ_both > θ_eps
        @test θ_both > θ_fungi
    end

    @testset "update_water_content!" begin
        n = 50
        soil = SoilProperties()

        # Create state with varying EPS and fungi
        state = AggregateState(n)
        state.E .= range(0.0, 20.0, length=n)
        state.F_i .= range(0.0, 10.0, length=n)

        # Allocate workspace arrays
        θ = zeros(n)
        θ_a = zeros(n)

        # Update
        ψ = -33.0
        update_water_content!(θ, θ_a, ψ, state, soil)

        # All values should be valid
        @test all(soil.θ_r .<= θ .<= soil.θ_s)
        @test all(0.0 .<= θ_a .<= soil.θ_s)

        # θ + θ_a = θ_s
        @test all(θ .+ θ_a .≈ soil.θ_s)

        # Higher EPS/fungi → higher θ
        @test θ[end] > θ[1]  # More EPS/fungi at end
    end

    @testset "Extreme EPS/fungi concentrations" begin
        soil = SoilProperties()
        ψ = -33.0  # Field capacity

        # Test with very high EPS and fungi (heavy biofilm)
        # With default ω_E = -0.01, ω_F = -0.02:
        # E = 100, F_i = 50 → exponent = -0.01*100 + -0.02*50 = -2.0
        # α_eff = α₀ × exp(-2.0) ≈ α₀ × 0.135
        E_extreme = 100.0
        F_extreme = 50.0

        α_eff = alpha_effective(E_extreme, F_extreme, soil)
        θ_extreme = water_content(ψ, E_extreme, F_extreme, soil)

        # α should be reduced but still positive
        @test 0.0 < α_eff < soil.α_vg
        @test α_eff ≈ soil.α_vg * exp(soil.ω_E * E_extreme + soil.ω_F * F_extreme)

        # θ should be increased (higher retention with biofilm)
        # Note: actual increase depends on ω_E, ω_F parameters
        @test θ_extreme > 0.25  # Above typical field capacity
        @test θ_extreme <= soil.θ_s

        # Air-filled porosity should be reduced but remain physical
        θ_a_extreme = soil.θ_s - θ_extreme
        @test θ_a_extreme < soil.θ_s  # Less than total porosity
        @test θ_a_extreme >= 0.0

        # CRITICAL: Verify diffusion doesn't break
        # Even with extreme biofilm, diffusion coefficients should be valid
        bio = BiologicalProperties()
        temp_cache = TemperatureCache()
        temp_cache.D_DOC_w = 0.864
        temp_cache.D_O2_w = 150.0
        temp_cache.D_O2_a = 1.73e6
        temp_cache.K_H_O = 29.0
        temp_cache.f_fun = 1.0

        D_C_extreme = D_eff_DOC(temp_cache.D_DOC_w, θ_extreme, soil.θ_s,
                                soil.ρ_b, soil.k_d_eq)
        D_O_extreme = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                   temp_cache.K_H_O, θ_extreme, θ_a_extreme, soil.θ_s)

        # Should be positive and finite (no NaN/Inf)
        # This is the CRITICAL check: model doesn't break under extreme conditions
        @test isfinite(D_C_extreme)
        @test isfinite(D_O_extreme)
        @test D_C_extreme > 0.0
        @test D_O_extreme > 0.0

        # O₂ diffusion behavior with biofilm
        # Note: Net effect depends on balance between reduced gas phase and changed aqueous phase
        D_O_baseline = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                    temp_cache.K_H_O, 0.3, 0.15, soil.θ_s)
        # Just verify both are valid - don't enforce specific ordering
        @test D_O_baseline > 0.0 && isfinite(D_O_baseline)
    end

    @testset "Tortuosity limits" begin
        soil = SoilProperties()

        # Dry soil
        τ_dry = tortuosity_millington_quirk(0.0, soil.θ_s)
        @test τ_dry == 0.0

        # Saturated: τ = θ²/θ_s^(2/3) = θ_s²/θ_s^(2/3) = θ_s^(4/3)
        τ_sat = tortuosity_millington_quirk(soil.θ_s, soil.θ_s)
        @test τ_sat ≈ soil.θ_s^(4/3)
        @test 0.3 < τ_sat < 0.4  # For θ_s ≈ 0.45

        # Monotonicity
        θ_range = range(0.0, soil.θ_s, length=10)
        τ_range = [tortuosity_millington_quirk(θ, soil.θ_s) for θ in θ_range]
        @test all(diff(τ_range) .>= 0)  # Increasing with θ
    end
end

@testset "Effective diffusion" begin
    soil = SoilProperties()
    bio = BiologicalProperties()

    # Mock temperature cache
    temp_cache = TemperatureCache()
    temp_cache.D_DOC_w = 0.864  # mm²/day at 20°C
    temp_cache.D_O2_w = 150.0
    temp_cache.D_O2_a = 1.73e6
    temp_cache.K_H_O = 29.0
    temp_cache.f_fun = 1.0  # At T_ref

    @testset "DOC diffusion with sorption" begin
        # No sorption (k_d = 0)
        soil_no_sorb = SoilProperties(k_d_eq=0.0)
        D_C_no_sorb = D_eff_DOC(temp_cache.D_DOC_w, 0.3, soil.θ_s, soil.ρ_b, 0.0)

        # With sorption
        D_C_sorb = D_eff_DOC(temp_cache.D_DOC_w, 0.3, soil.θ_s, soil.ρ_b, soil.k_d_eq)

        # Sorption retards diffusion
        @test D_C_sorb < D_C_no_sorb

        # Both should be positive and less than pure-phase
        @test 0.0 < D_C_sorb < temp_cache.D_DOC_w
        @test 0.0 < D_C_no_sorb < temp_cache.D_DOC_w
    end

    @testset "Bacterial motility" begin
        D_C = 0.1  # mm²/day
        D_B = D_eff_bacteria(D_C, soil.D_B_rel)

        # Bacteria much slower than DOC
        @test D_B < D_C
        @test D_B ≈ D_C * soil.D_B_rel
    end

    @testset "Fungal diffusion" begin
        # Non-insulated (with tortuosity)
        D_Fn_wet = D_eff_fungi_noninsulated(bio.D_Fn0, 1.0, 0.3, soil.θ_s)
        D_Fn_dry = D_eff_fungi_noninsulated(bio.D_Fn0, 1.0, 0.05, soil.θ_s)

        # Drier soil → lower diffusion
        @test D_Fn_dry < D_Fn_wet

        # Mobile: no tortuosity, but network-dependent. Translocation is
        # internal to hyphae, so D_Fm scales with local sessile biomass
        # (F_n + F_i) through a Michaelis-Menten factor with half-saturation
        # K_Fm_net (Falconer 2005 eq. 2.1; MATLAB single_aggregate_beta.m:394).
        F_dense = 100.0 * bio.K_Fm_net   # network far above half-saturation
        D_Fm = D_eff_fungi_mobile(bio.D_Fm0, 1.0, F_dense/2, F_dense/2,
                                  bio.K_Fm_net)

        @test D_Fm > 0.0
        # Approaches D_Fm0 for a dense network; no tortuosity reduction
        @test D_Fm ≈ bio.D_Fm0 rtol=0.02

        # No network → no translocation pathway
        @test D_eff_fungi_mobile(bio.D_Fm0, 1.0, 0.0, 0.0, bio.K_Fm_net) == 0.0

        # Half-maximum when F_n + F_i == K_Fm_net
        @test D_eff_fungi_mobile(bio.D_Fm0, 1.0, bio.K_Fm_net/2,
                                 bio.K_Fm_net/2, bio.K_Fm_net) ≈ 0.5 * bio.D_Fm0
    end

    @testset "Oxygen dual-phase diffusion" begin
        # Saturated (only aqueous phase)
        D_O_sat = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                temp_cache.K_H_O, soil.θ_s, 0.0, soil.θ_s)

        # Intermediate water content
        D_O_mid = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                temp_cache.K_H_O, 0.3, 0.15, soil.θ_s)

        # Dry (mostly gas phase)
        D_O_dry = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                temp_cache.K_H_O, 0.05, 0.40, soil.θ_s)

        # All should be positive
        @test D_O_sat > 0.0
        @test D_O_mid > 0.0
        @test D_O_dry > 0.0

        # Gas phase dominates when available (much faster diffusion)
        # At saturation: only aqueous phase (~150 mm²/day)
        # With air-filled pores: gas phase adds ~10⁴ faster diffusion
        @test D_O_mid > D_O_sat  # Some gas phase helps

        # CRITICAL: Magnitude test for gas-phase engagement
        # If partitioning fails, gas diffusion won't activate
        # D_O2_a/D_O2_w ≈ 11,500, so expect 2+ orders of magnitude speedup
        @test D_O_mid / D_O_sat > 100.0  # Gas phase should dominate
    end

    @testset "Water content limits on diffusion" begin
        # Dry soil (θ → 0)
        D_C_dry = D_eff_DOC(temp_cache.D_DOC_w, 0.0, soil.θ_s, soil.ρ_b, soil.k_d_eq)
        D_Fn_dry = D_eff_fungi_noninsulated(bio.D_Fn0, 1.0, 0.0, soil.θ_s)

        # Aqueous diffusion stops
        @test D_C_dry == 0.0
        @test D_Fn_dry == 0.0

        # But mobile fungi still translocate wherever a network exists:
        # D_Fm depends on network density, not on water content
        D_Fm_dry = D_eff_fungi_mobile(bio.D_Fm0, 1.0, 10*bio.K_Fm_net,
                                      10*bio.K_Fm_net, bio.K_Fm_net)
        @test D_Fm_dry > 0.0

        # And oxygen can diffuse through gas phase
        D_O_dry = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                temp_cache.K_H_O, 0.01, soil.θ_s - 0.01, soil.θ_s)
        @test D_O_dry > 0.0
    end

    # The `update_effective_diffusion!` testset, and the "Physics integration"
    # block that exercised the per-step Workspace update sequence, were archived
    # 2026-07-30 with the split solver. See _archive/split_solver_20260730/.
end

@testset "Type stability" begin
    soil = SoilProperties()
    bio = BiologicalProperties()
    temp_cache = TemperatureCache()
    temp_cache.D_DOC_w = 0.864
    temp_cache.D_O2_w = 150.0
    temp_cache.D_O2_a = 1.73e6
    temp_cache.K_H_O = 29.0
    temp_cache.f_fun = 1.0

    @testset "Water retention functions" begin
        @inferred alpha_effective(10.0, 5.0, soil)
        @inferred van_genuchten(-33.0, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
        @inferred water_content(-33.0, 10.0, 5.0, soil)
        @inferred tortuosity_millington_quirk(0.3, soil.θ_s)
    end

    @testset "Effective diffusion functions" begin
        @inferred D_eff_DOC(0.864, 0.3, soil.θ_s, soil.ρ_b, soil.k_d_eq)
        @inferred D_eff_bacteria(0.1, soil.D_B_rel)
        @inferred D_eff_fungi_noninsulated(bio.D_Fn0, 1.0, 0.3, soil.θ_s)
        @inferred D_eff_fungi_mobile(bio.D_Fm0, 1.0, 0.01, 0.01, bio.K_Fm_net)
        @inferred D_eff_oxygen(150.0, 1.73e6, 29.0, 0.3, 0.15, soil.θ_s)
    end
end


# ============================================================================
# van_genuchten_inverse
#
# This is the entry point for every experiment specified as water content or
# WFPS rather than as a potential (e.g. De Gryze 2006, "60 % WFPS"). A sign slip
# or an m/n confusion in
#       ψ = -((S_e^(-1/m) - 1)^(1/n)) / α,   m = 1 - 1/n
# would silently mis-set the water driver for a whole run.
#
# The round-trip needs no expected values at all: it is its own reference.
# ============================================================================
@testset "van Genuchten inverse" begin
    soil = SoilProperties()
    α, n_vg, θ_r, θ_s = soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s

    @testset "round-trip against van_genuchten" begin
        for frac in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 0.999)
            θ_target = θ_r + frac * (θ_s - θ_r)
            ψ = van_genuchten_inverse(θ_target, α, n_vg, θ_r, θ_s)
            @test ψ <= 0.0
            @test van_genuchten(ψ, α, n_vg, θ_r, θ_s) ≈ θ_target rtol=1e-10
        end
    end

    @testset "round-trip the other way" begin
        for ψ in (-1.0, -10.0, -33.0, -100.0, -1500.0)
            θ = van_genuchten(ψ, α, n_vg, θ_r, θ_s)
            @test van_genuchten_inverse(θ, α, n_vg, θ_r, θ_s) ≈ ψ rtol=1e-8
        end
    end

    @testset "monotone decreasing in θ" begin
        θs = [θ_r + f * (θ_s - θ_r) for f in 0.1:0.1:0.9]
        ψs = [van_genuchten_inverse(θ, α, n_vg, θ_r, θ_s) for θ in θs]
        @test all(diff(ψs) .> 0.0)      # wetter -> less negative
    end

    @testset "boundaries" begin
        @test van_genuchten_inverse(θ_s, α, n_vg, θ_r, θ_s) == 0.0
        @test van_genuchten_inverse(θ_s + 0.1, α, n_vg, θ_r, θ_s) == 0.0
        @test_throws ArgumentError van_genuchten_inverse(θ_r, α, n_vg, θ_r, θ_s)
        @test_throws ArgumentError van_genuchten_inverse(θ_r - 0.01, α, n_vg, θ_r, θ_s)
    end

    @testset "60 % WFPS — the De Gryze use case" begin
        # Not a magic number: 0.6·θ_s is the definition of 60 % water-filled
        # pore space. The assertion is the round-trip, not the value of ψ.
        θ_60 = 0.60 * θ_s
        ψ_60 = van_genuchten_inverse(θ_60, α, n_vg, θ_r, θ_s)
        @test ψ_60 < 0.0
        @test van_genuchten(ψ_60, α, n_vg, θ_r, θ_s) ≈ θ_60 rtol=1e-10
    end
end
