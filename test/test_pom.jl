# test_pom.jl
# Tests for POM dissolution (carbon/pom_dissolution.jl)

using Test
import SoilAggregateModel: J_P, R_P, pom_delay_factor, sigmoid_threshold,
                           POM_DELAY_STEEPNESS

@testset "POM Dissolution" begin

    # === POM activation delay ===
    # Previously written out in both timestepper.jl and mol.jl and tested in
    # neither, which is the one place the two integrators could disagree on when
    # POM turns on.
    @testset "pom_delay_factor" begin
        # No delay is exactly 1.0, so multiplying by it cannot perturb a run
        @test pom_delay_factor(0.0, 0.0) === 1.0
        @test pom_delay_factor(7.0, 0.0) === 1.0
        @test pom_delay_factor(7.0, -1.0) === 1.0

        # Half activation exactly at t_delay
        @test pom_delay_factor(10.0, 10.0) == 0.5

        # Suppressed before, released after; monotone through the switch
        @test pom_delay_factor(0.0, 10.0) < 1.0e-4
        f = [pom_delay_factor(t, 10.0) for t in 0.0:0.5:20.0]
        @test all(diff(f) .> 0.0)

        # Documented width: 90 % activation by t ≈ 12.2 for t_delay = 10
        @test pom_delay_factor(12.2, 10.0) > 0.9
        @test pom_delay_factor(11.0, 10.0) < 0.9

        # It is the shared primitive, not a second sigmoid
        @test pom_delay_factor(6.0, 10.0) ==
              sigmoid_threshold(6.0, 10.0, POM_DELAY_STEEPNESS)
        @test POM_DELAY_STEEPNESS == 10.0
    end

    # === Test 1: Known inputs — verify formula ===
    @testset "Known inputs" begin
        # Reference values
        P = 800.0            # μg-C (80% of initial)
        P_0 = 1000.0         # μg-C
        B_0 = 5.0            # μg-C/mm³
        F_n_0 = 3.0          # μg-C/mm³
        θ_0 = 0.25           # -
        O_aq_0 = 0.1         # μg/mm³
        R_P_max_T = 0.5      # μg-C/mm²/day
        K_B_P = 1.0          # μg-C/mm³
        K_F_P = 1.0          # μg-C/mm³
        θ_P = 0.1            # -
        L_P = 0.05           # μg/mm³
        r_0 = 0.1            # mm

        # Hand-calculated expected values
        pom_factor = P / P_0                     # 0.8
        bacterial = B_0 / (K_B_P + B_0)          # 5/(1+5) = 5/6
        fungal = F_n_0 / (K_F_P + F_n_0)         # 3/(1+3) = 3/4 = 0.75
        microbial = 0.5 * (bacterial + fungal)   # 0.5 * (5/6 + 3/4) = 0.5 * 19/12
        moisture = θ_0 / (θ_P + θ_0)             # 0.25/(0.1+0.25) = 25/35 = 5/7
        oxygen = O_aq_0 / (L_P + O_aq_0)         # 0.1/(0.05+0.1) = 2/3

        # Additive formula: 0.5 * (bacterial + fungal) instead of bacterial * fungal
        expected_J_P_calc = 0.5 * 0.8 * 0.5 * (5.0/6.0 + 0.75) * (5.0/7.0) * (2.0/3.0)

        J_P_val = J_P(P, P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T, K_B_P, K_F_P, θ_P, L_P)
        @test J_P_val ≈ expected_J_P_calc rtol=1e-10

        # Test R_P relationship
        expected_R_P = 4.0 * π * r_0^2 * J_P_val
        R_P_val = R_P(J_P_val, r_0)
        @test R_P_val ≈ expected_R_P rtol=1e-12
        @test R_P_val ≈ 4.0 * π * 0.01 * J_P_val rtol=1e-12
    end

    # === Test 2: Zero POM → no flux ===
    @testset "No POM" begin
        J_P_val = J_P(0.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        @test J_P_val == 0.0

        R_P_val = R_P(J_P_val, 0.1)
        @test R_P_val == 0.0
    end

    # === Test 3: No bacteria → fungi still contribute (additive coupling) ===
    @testset "No bacteria" begin
        J_P_val = J_P(800.0, 1000.0, 0.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        # Additive: microbial_factor = 0.5 * (0 + 3/(1+3)) = 0.5 * 0.75 = 0.375
        # J_P = 0.5 * 0.8 * 0.375 * (25/35) * (2/3)
        expected = 0.5 * 0.8 * 0.375 * (25.0/35.0) * (2.0/3.0)
        @test J_P_val ≈ expected rtol=1e-10
        @test J_P_val > 0.0  # Non-zero when only fungi present
    end

    # === Test 4: No fungi → bacteria still contribute (additive coupling) ===
    @testset "No fungi" begin
        J_P_val = J_P(800.0, 1000.0, 5.0, 0.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        # Additive: microbial_factor = 0.5 * (5/(1+5) + 0) = 0.5 * 5/6 = 5/12
        # J_P = 0.5 * 0.8 * (5/12) * (25/35) * (2/3)
        expected = 0.5 * 0.8 * (5.0/12.0) * (25.0/35.0) * (2.0/3.0)
        @test J_P_val ≈ expected rtol=1e-10
        @test J_P_val > 0.0  # Non-zero when only bacteria present
    end

    # === Test 5: No water → minimal flux ===
    @testset "Dry conditions" begin
        J_P_val = J_P(800.0, 1000.0, 5.0, 3.0, 0.0, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        @test J_P_val ≈ 0.0 atol=1e-12  # 0/(θ_P + 0) = 0
    end

    # === Test 6: No oxygen → minimal flux ===
    @testset "Anoxic conditions" begin
        J_P_val = J_P(800.0, 1000.0, 5.0, 3.0, 0.25, 0.0, 0.5, 1.0, 1.0, 0.1, 0.05)
        @test J_P_val ≈ 0.0 atol=1e-12  # 0/(L_P + 0) = 0
    end

    # === Test 7: Saturating conditions → approach R_P_max_T ===
    @testset "Saturating conditions" begin
        # Very high biomass, water, oxygen
        B_0 = 1000.0         # >> K_B_P
        F_n_0 = 1000.0       # >> K_F_P
        θ_0 = 10.0           # >> θ_P
        O_aq_0 = 10.0        # >> L_P

        J_P_val = J_P(1000.0, 1000.0, B_0, F_n_0, θ_0, O_aq_0, 0.5, 1.0, 1.0, 0.1, 0.05)

        # Each Monod term → 1, so J_P → R_P_max_T * (P/P_0) = 0.5 * 1.0 = 0.5
        # (with finite saturation, 1000/(1+1000) ≈ 0.999, giving ~1.6% total deviation)
        @test J_P_val ≈ 0.5 rtol=2e-2
    end

    # === Test 8: Half-saturation behavior ===
    @testset "Half-saturation (bacteria)" begin
        # When B_0 = K_B_P, bacterial contribution = 0.5
        # With additive coupling and saturating fungi, microbial = 0.5*(0.5 + ~1.0) = 0.75
        B_0 = 1.0  # = K_B_P
        F_n_0 = 1000.0  # Saturating
        θ_0 = 10.0      # Saturating
        O_aq_0 = 10.0   # Saturating

        J_P_val = J_P(1000.0, 1000.0, B_0, F_n_0, θ_0, O_aq_0, 1.0, 1.0, 1.0, 0.1, 0.05)

        # Expected: 1.0 * 1.0 * 0.5 * (0.5 + ~1.0) * ~1.0 * ~1.0 ≈ 0.75
        # (other factors at 1000 give ~0.999 each, so ~1% total deviation)
        @test J_P_val ≈ 0.75 rtol=2e-2
    end

    @testset "Half-saturation (fungi)" begin
        # When F_n_0 = K_F_P, fungal contribution = 0.5
        # With additive coupling and saturating bacteria, microbial = 0.5*(~1.0 + 0.5) = 0.75
        B_0 = 1000.0    # Saturating
        F_n_0 = 1.0     # = K_F_P
        θ_0 = 10.0      # Saturating
        O_aq_0 = 10.0   # Saturating

        J_P_val = J_P(1000.0, 1000.0, B_0, F_n_0, θ_0, O_aq_0, 1.0, 1.0, 1.0, 0.1, 0.05)

        # Expected: 1.0 * 1.0 * 0.5 * (~1.0 + 0.5) * ~1.0 * ~1.0 ≈ 0.75
        # (other factors at 1000 give ~0.999 each, so ~1% total deviation)
        @test J_P_val ≈ 0.75 rtol=2e-2
    end

    @testset "Half-saturation (moisture)" begin
        # When θ_0 = θ_P, moisture factor = 0.5
        B_0 = 1000.0    # Saturating
        F_n_0 = 1000.0  # Saturating
        θ_0 = 0.1       # = θ_P
        O_aq_0 = 10.0   # Saturating

        J_P_val = J_P(1000.0, 1000.0, B_0, F_n_0, θ_0, O_aq_0, 1.0, 1.0, 1.0, 0.1, 0.05)

        # Expected: 1.0 * 1.0 * ~1.0 * ~1.0 * 0.5 * ~1.0 ≈ 0.5
        @test J_P_val ≈ 0.5 rtol=1e-2
    end

    @testset "Half-saturation (oxygen)" begin
        # When O_aq_0 = L_P, oxygen factor = 0.5
        B_0 = 1000.0    # Saturating
        F_n_0 = 1000.0  # Saturating
        θ_0 = 10.0      # Saturating
        O_aq_0 = 0.05   # = L_P

        J_P_val = J_P(1000.0, 1000.0, B_0, F_n_0, θ_0, O_aq_0, 1.0, 1.0, 1.0, 0.1, 0.05)

        # Expected: 1.0 * 1.0 * ~1.0 * ~1.0 * ~1.0 * 0.5 ≈ 0.5
        # (other factors at 1000 give ~0.999 each, so ~1.6% total deviation)
        @test J_P_val ≈ 0.5 rtol=2e-2
    end

    # === Test 9: R_P scaling with radius ===
    @testset "R_P scales with r²" begin
        J_P_val = 0.3  # μg-C/mm²/day

        r1 = 0.1   # mm
        r2 = 0.2   # mm

        R_P_1 = R_P(J_P_val, r1)
        R_P_2 = R_P(J_P_val, r2)

        # Should scale as r²
        @test R_P_2 / R_P_1 ≈ (r2/r1)^2 rtol=1e-12
        @test R_P_2 / R_P_1 ≈ 4.0 rtol=1e-12
    end

    # === Test 10: POM depletion → reduced flux ===
    @testset "POM depletion reduces flux" begin
        # Reference flux at full POM
        J_P_full = J_P(1000.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)

        # Half POM remaining
        J_P_half = J_P(500.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)

        # Should be exactly half (linear in P)
        @test J_P_half / J_P_full ≈ 0.5 rtol=1e-12

        # 10% POM remaining
        J_P_depleted = J_P(100.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        @test J_P_depleted / J_P_full ≈ 0.1 rtol=1e-12
    end

    # === Test 11: Type stability ===
    @testset "Type stability" begin
        @test @inferred J_P(800.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05) isa Float64
        @test @inferred R_P(0.3, 0.1) isa Float64
    end

    # === Test 12: Units consistency check ===
    @testset "Units consistency" begin
        # Realistic values with proper units
        P = 500.0            # μg-C
        P_0 = 1000.0         # μg-C
        B_0 = 2.0            # μg-C/mm³
        F_n_0 = 1.5          # μg-C/mm³
        θ_0 = 0.3            # dimensionless
        O_aq_0 = 0.08        # μg/mm³
        R_P_max_T = 0.5      # μg-C/mm²/day
        K_B_P = 1.0          # μg-C/mm³
        K_F_P = 1.0          # μg-C/mm³
        θ_P = 0.1            # dimensionless
        L_P = 0.05           # μg/mm³
        r_0 = 0.1            # mm

        J_P_val = J_P(P, P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T, K_B_P, K_F_P, θ_P, L_P)
        @test 0.0 < J_P_val < R_P_max_T  # Flux must be between 0 and max

        R_P_val = R_P(J_P_val, r_0)
        # For r_0 = 0.1 mm, surface area = 4π(0.1)² ≈ 0.1257 mm²
        # So R_P ≈ 0.1257 * J_P_val
        @test R_P_val ≈ 4.0 * π * 0.01 * J_P_val rtol=1e-12
    end

    # === Test 13: Monotonicity ===
    @testset "Monotonic increases" begin
        base_args = (1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)

        # Increasing P → increasing J_P
        J_P_1 = J_P(400.0, base_args...)
        J_P_2 = J_P(600.0, base_args...)
        J_P_3 = J_P(800.0, base_args...)
        @test J_P_1 < J_P_2 < J_P_3

        # Increasing B_0 → increasing J_P
        J_P_B1 = J_P(800.0, 1000.0, 1.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        J_P_B2 = J_P(800.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        J_P_B3 = J_P(800.0, 1000.0, 10.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        @test J_P_B1 < J_P_B2 < J_P_B3

        # Increasing θ_0 → increasing J_P
        J_P_θ1 = J_P(800.0, 1000.0, 5.0, 3.0, 0.15, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        J_P_θ2 = J_P(800.0, 1000.0, 5.0, 3.0, 0.25, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        J_P_θ3 = J_P(800.0, 1000.0, 5.0, 3.0, 0.35, 0.1, 0.5, 1.0, 1.0, 0.1, 0.05)
        @test J_P_θ1 < J_P_θ2 < J_P_θ3
    end
end
