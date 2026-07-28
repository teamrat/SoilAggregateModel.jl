# test_tessellation.jl
# Tests for src/physics/tessellation.jl
#
# `ω` scales every background carbon pool in the model, so the arithmetic here
# has to be exact.
#
# The central test is `total_POM_C / (V_soil·ρ_b) == I_input`. That identity is
# EXACT and independent of the size distribution:
#
#   Σ N_i·P_0,i = Σ [f_i·V/((4/3)π(d_i f_pack/2)³)] · (4/3)π(d_i/2)³·ρ_POM
#               = (V·ρ_POM/f_pack³) · Σ f_i
#               = V·ρ_POM·φ_POM             (since f_pack³ = 1/φ_POM)
#               = V·I_input·ρ_b
#
# so it holds for ANY diameters and ANY normalised fractions. It simultaneously
# pins the f_pack exponent, both radius/diameter halvings, and the (4/3)π in
# both places — a wrong exponent or a diameter-for-radius slip breaks it.

using Test

@testset "Domain tessellation" begin

    @testset "domain_tessellation — worked example (docstring)" begin
        # De Gryze 2006: 1.5 g wheat stems (44.3 % C) per 150 g soil.
        t = domain_tessellation(ρ_POM=200.0, I_input=4.43e-3, ρ_b=1300.0)

        @test t.φ_POM  ≈ 4.43e-3 * 1300.0 / 200.0      # 0.028795, by definition
        @test t.f_pack ≈ (1.0 / t.φ_POM)^(1/3)
        @test t.f_pack ≈ 3.26258 rtol=1e-5
        @test t.f_domain == 10.0                        # f_domain_min binds here
        @test t.ω ≈ 28.7950 rtol=1e-5
        @test t.ω ≈ (t.f_domain / t.f_pack)^3
    end

    @testset "ω is exactly 1 when packing already exceeds the minimum" begin
        # Sparse amendment: f_pack = 32.6 > f_domain_min, so domains do not
        # overlap and no correction may be applied.
        t = domain_tessellation(ρ_POM=200.0, I_input=4.43e-6, ρ_b=1300.0)
        @test t.f_pack > 10.0
        @test t.f_domain == t.f_pack
        @test t.ω == 1.0            # exact, not approximate
    end

    @testset "f_domain_min is honoured" begin
        t = domain_tessellation(ρ_POM=200.0, I_input=4.43e-3, ρ_b=1300.0,
                                f_domain_min=25.0)
        @test t.f_domain == 25.0
        @test t.ω ≈ (25.0 / t.f_pack)^3
    end

    @testset "pom_population — the amendment identity" begin
        ρ_POM   = 200.0
        I_input = 4.43e-3
        ρ_b     = 1300.0
        V       = 1.0e6                     # 1 litre
        t       = domain_tessellation(ρ_POM=ρ_POM, I_input=I_input, ρ_b=ρ_b)

        # Three unrelated size distributions. The identity is distribution-free,
        # so all three must reproduce I_input to machine precision.
        for (diams, fracs) in (
                ([0.65, 0.95, 1.25, 1.55, 1.85], [0.1, 0.2, 0.4, 0.2, 0.1]),
                ([0.3, 2.0],                     [0.7, 0.3]),
                ([1.0],                          [1.0]))

            pop = pom_population(diams, fracs, t; ρ_POM=ρ_POM, soil_volume_mm3=V)
            @test pop.total_POM_C / (V * ρ_b) ≈ I_input rtol=1e-12

            # Per-particle carbon is a sphere of POM, not of the packing cell
            for (k, d) in enumerate(diams)
                @test pop.P_0_per_particle[k] ≈ (4/3) * π * (d/2)^3 * ρ_POM
            end
            @test all(pop.N_POM .> 0.0)

            # V_pack must tile the soil exactly: Σ Nᵢ·V_pack,ᵢ = V_soil.
            # This is why the bulk subtraction equals the per-particle one.
            @test sum(pop.N_POM .* pop.V_pack) ≈ V rtol=1e-12
            for (k, d) in enumerate(diams)
                @test pop.V_pack[k] ≈ (4/3) * π * (d * t.f_pack / 2)^3
            end
        end
    end

    @testset "pom_population — counts scale as expected" begin
        t = domain_tessellation(ρ_POM=200.0, I_input=4.43e-3, ρ_b=1300.0)
        # Halving the diameter at fixed mass fraction must give 8x the particles
        a = pom_population([1.0], [1.0], t; ρ_POM=200.0)
        b = pom_population([0.5], [1.0], t; ρ_POM=200.0)
        @test b.N_POM[1] ≈ 8.0 * a.N_POM[1] rtol=1e-12
        @test b.P_0_per_particle[1] ≈ a.P_0_per_particle[1] / 8.0 rtol=1e-12
        # ...and the same total carbon
        @test b.total_POM_C ≈ a.total_POM_C rtol=1e-12
    end

    @testset "an unnormalised distribution loses mass in proportion" begin
        # A truncated distribution summing to less than 1 recovers an I_input
        # low by exactly that factor — the residual behind
        # "Total POM input 4425.1 vs expected 4430" in run_degryze.jl.
        ρ_POM, I_input, ρ_b, V = 200.0, 4.43e-3, 1300.0, 1.0e6
        t = domain_tessellation(ρ_POM=ρ_POM, I_input=I_input, ρ_b=ρ_b)

        fracs = [0.1, 0.2, 0.4, 0.2, 0.1] .* 0.99889     # deliberately short
        diams = [0.65, 0.95, 1.25, 1.55, 1.85]
        pop = pom_population(diams, fracs, t; ρ_POM=ρ_POM, soil_volume_mm3=V)

        @test pop.total_POM_C / (V * ρ_b) ≈ I_input * 0.99889 rtol=1e-12
        # ...and renormalising restores it exactly
        pop2 = pom_population(diams, fracs ./ sum(fracs), t;
                              ρ_POM=ρ_POM, soil_volume_mm3=V)
        @test pop2.total_POM_C / (V * ρ_b) ≈ I_input rtol=1e-12
    end

    @testset "log_interpolate_fraction" begin
        # Sand (53-2000 µm) split at the 250 µm sieve — De Gryze spec §0a A1
        @test log_interpolate_fraction(53, 250, 2000) ≈ log(250/53)/log(2000/53)
        @test log_interpolate_fraction(53, 250, 2000) ≈ 0.42725 rtol=1e-4
        # Midpoint in log space gives exactly one half
        @test log_interpolate_fraction(1, 10, 100) ≈ 0.5
        # Bounded by construction
        @test 0.0 < log_interpolate_fraction(53, 100, 2000) < 1.0
    end

    @testset "input validation" begin
        @test_throws ArgumentError domain_tessellation(ρ_POM=0.0, I_input=4.43e-3, ρ_b=1300.0)
        @test_throws ArgumentError domain_tessellation(ρ_POM=200.0, I_input=0.0, ρ_b=1300.0)
        @test_throws ArgumentError domain_tessellation(ρ_POM=200.0, I_input=4.43e-3, ρ_b=0.0)

        # φ_POM >= 1: the amendment cannot fill more than the soil
        @test_throws ArgumentError domain_tessellation(ρ_POM=1.0, I_input=0.5, ρ_b=1300.0)

        t = domain_tessellation(ρ_POM=200.0, I_input=4.43e-3, ρ_b=1300.0)
        @test_throws ArgumentError pom_population([1.0, 2.0], [1.0], t; ρ_POM=200.0)

        @test_throws ArgumentError log_interpolate_fraction(250, 53, 2000)   # lo > mid
        @test_throws ArgumentError log_interpolate_fraction(53, 2000, 250)   # mid > hi
        @test_throws ArgumentError log_interpolate_fraction(0, 250, 2000)    # lo <= 0
    end
end
