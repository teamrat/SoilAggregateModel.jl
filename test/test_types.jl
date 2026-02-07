"""
Tests for types.jl
"""

import SoilAggregateModel: TemperatureCache, Workspace

@testset "TemperatureCache" begin
    cache = TemperatureCache()

    # Check all fields exist
    @test isdefined(cache, :f_bac)
    @test isdefined(cache, :f_fun)
    @test isdefined(cache, :f_eps)
    @test isdefined(cache, :f_maoc_s)
    @test isdefined(cache, :f_maoc_d)
    @test isdefined(cache, :f_pom)
    @test isdefined(cache, :D_O2_w)
    @test isdefined(cache, :D_DOC_w)
    @test isdefined(cache, :D_O2_a)
    @test isdefined(cache, :D_Fm)
    @test isdefined(cache, :K_H_O)

    # Default constructor should initialize to NaN
    @test isnan(cache.f_bac)
    @test isnan(cache.K_H_O)

    # Should be mutable
    cache.f_bac = 1.5
    @test cache.f_bac == 1.5
end

@testset "AggregateState" begin
    n = 100
    state = AggregateState(n)

    # Check all vectors have correct length
    @test length(state.C) == n
    @test length(state.B) == n
    @test length(state.F_n) == n
    @test length(state.F_m) == n
    @test length(state.O) == n
    @test length(state.F_i) == n
    @test length(state.E) == n
    @test length(state.M) == n

    # Scalar fields
    @test state.P == 0.0
    @test state.CO2_cumulative == 0.0

    # Should be mutable
    state.P = 1000.0
    @test state.P == 1000.0

    state.C[1] = 5.0
    @test state.C[1] == 5.0
end

@testset "AggregateState copy" begin
    n = 50
    state = AggregateState(n)
    state.C .= 10.0
    state.B .= 1.0
    state.P = 500.0

    # Deep copy
    state_copy = copy(state)

    # Values should be equal
    @test state_copy.C == state.C
    @test state_copy.B == state.B
    @test state_copy.P == state.P

    # But arrays should be independent
    state_copy.C[1] = 20.0
    @test state_copy.C[1] == 20.0
    @test state.C[1] == 10.0  # Original unchanged

    state_copy.P = 1000.0
    @test state_copy.P == 1000.0
    @test state.P == 500.0  # Original unchanged
end

@testset "Workspace" begin
    n = 150
    ws = Workspace(n)

    # Tridiagonal arrays
    @test length(ws.lower) == n - 1
    @test length(ws.diag) == n
    @test length(ws.upper) == n - 1
    @test length(ws.rhs) == n

    # Spatially varying quantities
    @test length(ws.θ) == n
    @test length(ws.θ_a) == n
    @test length(ws.D_C) == n
    @test length(ws.D_B) == n
    @test length(ws.D_Fn) == n
    @test length(ws.D_O) == n

    # Temperature cache
    @test ws.f_T isa TemperatureCache

    # Workspace is immutable but contents are mutable
    ws.θ[1] = 0.3
    @test ws.θ[1] == 0.3
end

@testset "OutputRecord" begin
    n = 100
    state = AggregateState(n)
    state.C .= 5.0
    state.P = 800.0

    rec = OutputRecord(10.0, copy(state), 1e-14)

    @test rec.t == 10.0
    @test rec.state.C[1] == 5.0
    @test rec.state.P == 800.0
    @test rec.mass_balance_error == 1e-14

    # State should be independent copy
    state.C[1] = 100.0
    @test rec.state.C[1] == 5.0  # Record unchanged
end
