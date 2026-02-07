"""
Tests for physics/*.jl modules
"""

import SoilAggregateModel: TemperatureCache, AggregateState, Workspace,
                           alpha_effective, van_genuchten, water_content, update_water_content!,
                           tortuosity_millington_quirk, D_eff_DOC, D_eff_bacteria,
                           D_eff_fungi_noninsulated, D_eff_fungi_mobile, D_eff_oxygen,
                           update_effective_diffusion!

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

        # θ should approach θ_s (very high retention)
        @test θ_extreme > 0.9 * soil.θ_s
        @test θ_extreme <= soil.θ_s

        # Air-filled porosity nearly zero
        θ_a_extreme = soil.θ_s - θ_extreme
        @test θ_a_extreme < 0.1 * soil.θ_s
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

        # Gas-phase O₂ diffusion should collapse (θ_a → 0)
        # Compare to baseline with moderate water content
        D_O_baseline = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                    temp_cache.K_H_O, 0.3, 0.15, soil.θ_s)
        @test D_O_extreme < D_O_baseline  # Less gas phase
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

        # Mobile (no tortuosity)
        D_Fm = D_eff_fungi_mobile(bio.D_Fm0, 1.0)  # Mobile fungi use D_Fm0

        # D_Fm independent of water content
        @test D_Fm > 0.0
        @test D_Fm ≈ bio.D_Fm0  # No reduction from tortuosity
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

        # But mobile fungi still work
        D_Fm_dry = D_eff_fungi_mobile(bio.D_Fm0, 1.0)
        @test D_Fm_dry > 0.0

        # And oxygen can diffuse through gas phase
        D_O_dry = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                temp_cache.K_H_O, 0.01, soil.θ_s - 0.01, soil.θ_s)
        @test D_O_dry > 0.0
    end

    @testset "update_effective_diffusion!" begin
        n = 50
        ws = Workspace(n)
        bio = BiologicalProperties()

        # Set water content (varying)
        ws.θ .= range(0.1, 0.4, length=n)
        ws.θ_a .= soil.θ_s .- ws.θ

        # Update diffusion coefficients
        update_effective_diffusion!(ws, soil, bio, temp_cache)

        # All should be positive
        @test all(ws.D_C .> 0.0)
        @test all(ws.D_B .> 0.0)
        @test all(ws.D_Fn .> 0.0)
        @test all(ws.D_O .> 0.0)

        # D_Fm should be scalar and positive
        @test temp_cache.D_Fm > 0.0

        # D_B should be less than D_C (slower motility)
        @test all(ws.D_B .< ws.D_C)

        # Wetter nodes → higher aqueous diffusion
        @test ws.D_C[end] > ws.D_C[1]
        @test ws.D_Fn[end] > ws.D_Fn[1]
    end
end

@testset "Physics integration" begin
    @testset "Full workspace update sequence" begin
        n = 100
        soil = SoilProperties()
        bio = BiologicalProperties()

        # Create state with some EPS and fungi
        state = AggregateState(n)
        state.E .= 5.0
        state.F_i .= 2.0

        # Create workspace
        ws = Workspace(n)

        # Simulate timestep update sequence
        T = 293.15  # 20°C
        ψ = -33.0   # Field capacity

        # Step 1: Update temperature cache
        # (Simplified - in real code this uses update_temperature_cache!)
        ws.f_T.D_DOC_w = 0.864
        ws.f_T.D_O2_w = 150.0
        ws.f_T.D_O2_a = 1.73e6
        ws.f_T.K_H_O = 29.0
        ws.f_T.f_fun = 1.0

        # Step 2: Update water content
        update_water_content!(ws.θ, ws.θ_a, ψ, state, soil)

        # Step 3: Update effective diffusion
        update_effective_diffusion!(ws, soil, bio, ws.f_T)

        # Verify everything is reasonable
        @test all(soil.θ_r .<= ws.θ .<= soil.θ_s)
        @test all(ws.D_C .> 0.0)
        @test all(ws.D_B .> 0.0)
        @test all(ws.D_Fn .> 0.0)
        @test all(ws.D_O .> 0.0)

        # Relationships
        @test all(ws.D_B .< ws.D_C)  # Bacteria slower than DOC
        @test ws.f_T.D_Fm > 0.0       # Mobile fungi active
    end
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
        @inferred D_eff_fungi_mobile(bio.D_Fm0, 1.0)
        @inferred D_eff_oxygen(150.0, 1.73e6, 29.0, 0.3, 0.15, soil.θ_s)
    end
end
