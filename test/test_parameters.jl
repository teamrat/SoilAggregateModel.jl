"""
Tests for parameters.jl
"""

@testset "BiologicalProperties default constructor" begin
    bio = BiologicalProperties()

    # Check a few key fields exist and have positive values
    @test bio.r_B_max > 0.0
    @test bio.K_B > 0.0
    @test bio.Y_B_max > 0.0
    @test bio.Y_B_max <= 1.0  # Yield must be ≤ 1

    @test bio.r_F_max > 0.0
    @test bio.K_F > 0.0
    @test bio.Y_F > 0.0
    @test bio.Y_F <= 1.0

    # Activation energies should be positive
    @test bio.Ea_B > 0.0
    @test bio.Ea_F > 0.0
    @test bio.Ea_EPS > 0.0
    @test bio.Ea_MAOC_sorb > 0.0
    @test bio.Ea_MAOC_desorb > 0.0
    @test bio.Ea_POM > 0.0

    # Reference temperature should be ~293 K (20°C)
    @test bio.T_ref ≈ 293.15

    # Delta should be > 1
    @test bio.delta > 1.0

    # Efficiency should be between 0 and 1
    @test 0.0 < bio.η_conv <= 1.0
    @test 0.0 < bio.γ <= 1.0

    # Protection parameters should be small
    @test bio.ε_F < 1e-6
    @test bio.ε_maoc > 0.0

    # POM parameters
    @test bio.P_0 > 0.0
    @test bio.r_0 > 0.0
    @test bio.ρ_POM > 0.0

    # Respiratory quotient
    @test bio.α_O > 0.0
end

@testset "BiologicalProperties custom values" begin
    bio = BiologicalProperties(r_B_max=10.0, K_B=25.0, T_ref=298.15)

    @test bio.r_B_max == 10.0
    @test bio.K_B == 25.0
    @test bio.T_ref == 298.15

    # Other fields should still have defaults
    @test bio.Y_B_max > 0.0
end

@testset "SoilProperties default constructor" begin
    soil = SoilProperties()

    # Van Genuchten parameters
    @test 0.0 <= soil.θ_r < soil.θ_s <= 1.0
    @test soil.α_vg > 0.0
    @test soil.n_vg > 1.0

    # EPS/fungi effects should be negative
    @test soil.ω_E < 0.0
    @test soil.ω_F < 0.0

    # Sorption parameters
    @test soil.k_d_eq > 0.0
    @test soil.ρ_b > 0.0

    # MAOC capacity
    @test soil.M_max > 0.0
    @test soil.k_L > 0.0
    @test 0.0 < soil.n_LF <= 1.0
    @test soil.k_ma > 0.0
    @test 0.0 <= soil.f_clay_silt <= 1.0

    # Diffusion coefficients should be positive
    @test soil.D_C0_ref > 0.0
    @test soil.D_O2_w_ref > 0.0
    @test soil.D_O2_a_ref > 0.0

    # O₂ in air should be much faster than in water
    @test soil.D_O2_a_ref > 100 * soil.D_O2_w_ref

    # Bacterial motility should be small
    @test 0.0 < soil.D_B_rel < 0.1

    # Aggregate stability parameters
    @test soil.k_F > 0.0
    @test soil.χ > 0.0
    @test soil.a_p > 0.0
end

@testset "SoilProperties custom values" begin
    soil = SoilProperties(θ_s=0.5, θ_r=0.08, ρ_b=1.5e6)

    @test soil.θ_s == 0.5
    @test soil.θ_r == 0.08
    @test soil.ρ_b == 1.5e6

    # Other fields should still have defaults
    @test soil.n_vg > 1.0
end

@testset "Parameter struct immutability" begin
    bio = BiologicalProperties()
    soil = SoilProperties()

    # Structs should be immutable
    @test !ismutable(bio)
    @test !ismutable(soil)
end
