"""
Tests for environment.jl
"""

import SoilAggregateModel: EnvironmentalDrivers, evaluate

@testset "EnvironmentalDrivers with constants" begin
    env = EnvironmentalDrivers(293.15, -33.0, 0.27)

    # Should be callable at any time
    @test env.T(0.0) == 293.15
    @test env.T(100.0) == 293.15
    @test env.ψ(50.0) == -33.0
    @test env.O2(1000.0) == 0.27

    # evaluate() helper
    vals = evaluate(env, 10.0)
    @test vals.T == 293.15
    @test vals.ψ == -33.0
    @test vals.O2 == 0.27
end

@testset "EnvironmentalDrivers with functions" begin
    # Temperature varies sinusoidally
    T_func = t -> 293.15 + 5.0 * sin(2π * t / 365)

    # Matric potential is constant
    # O₂ varies linearly (unrealistic, but tests function handling)
    O2_func = t -> 0.27 + 0.01 * t

    env = EnvironmentalDrivers(T_func, -33.0, O2_func)

    # Check at t=0
    @test env.T(0.0) ≈ 293.15
    @test env.ψ(0.0) == -33.0
    @test env.O2(0.0) ≈ 0.27

    # Check at t=365/4 (quarter year, max temperature)
    t_max = 365.0 / 4.0
    @test env.T(t_max) ≈ 293.15 + 5.0 * sin(π/2) ≈ 298.15

    # Check O₂ varies
    @test env.O2(10.0) ≈ 0.37

    # evaluate() helper
    vals = evaluate(env, 100.0)
    @test vals.T ≈ 293.15 + 5.0 * sin(2π * 100.0 / 365)
    @test vals.ψ == -33.0
    @test vals.O2 ≈ 0.27 + 0.01 * 100.0
end

@testset "EnvironmentalDrivers mixed types" begin
    # Mix of constant and function
    env = EnvironmentalDrivers(
        t -> 290.0 + 0.01 * t,  # function
        -50.0,                   # constant
        t -> 0.3 - 0.0001 * t   # function
    )

    vals = evaluate(env, 100.0)
    @test vals.T ≈ 291.0
    @test vals.ψ == -50.0
    @test vals.O2 ≈ 0.29
end

@testset "EnvironmentalDrivers type stability" begin
    # Constant version
    env_const = EnvironmentalDrivers(293.15, -33.0, 0.27)
    @inferred env_const.T(0.0)
    @inferred env_const.ψ(0.0)
    @inferred env_const.O2(0.0)

    # Function version
    env_func = EnvironmentalDrivers(t -> 293.15, t -> -33.0, t -> 0.27)
    @inferred env_func.T(0.0)
    @inferred env_func.ψ(0.0)
    @inferred env_func.O2(0.0)

    # evaluate() should be type-stable
    @inferred evaluate(env_const, 0.0)
    @inferred evaluate(env_func, 0.0)
end
