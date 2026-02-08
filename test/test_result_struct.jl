# test_result_struct.jl
# Tests for SimulationResult, GridInfo, and ParameterSet types

using Test

@testset "SimulationResult structure" begin
    @testset "GridInfo construction and conservation weights" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Constant environment
        T(t) = 293.15  # 20°C
        ψ(t) = -10.0   # kPa
        O2(t) = 0.3    # μg/mm³

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 5.0); n_grid=20)

        @test result isa SimulationResult

        # Grid
        @test result.grid isa GridInfo
        @test result.grid.n == 20
        @test result.grid.r_0 ≈ 0.1
        @test result.grid.r_max ≈ 2.0
        @test length(result.grid.r_grid) == 20
        @test length(result.grid.W) == 20
        @test result.grid.W[1] ≈ 4π * result.grid.r_grid[1]^2 * result.grid.h

        # Conservation weights are consistent
        for i in 1:result.grid.n
            @test result.grid.W[i] ≈ 4π * result.grid.r_grid[i]^2 * result.grid.h
        end

        # Grid spacing is uniform
        @test result.grid.h ≈ (result.grid.r_max - result.grid.r_0) / (result.grid.n - 1)

        # First and last points
        @test result.grid.r_grid[1] ≈ result.grid.r_0
        @test result.grid.r_grid[end] ≈ result.grid.r_max
    end

    @testset "ParameterSet preservation" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0); n_grid=15)

        # Parameters preserved
        @test result.params isa ParameterSet
        @test result.params.bio.P_0 == bio.P_0
        @test result.params.soil.θ_s == soil.θ_s
        @test result.params.bio.μ_B == bio.μ_B
        @test result.params.soil.α_vg == soil.α_vg
    end

    @testset "SimulationResult field access" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 5.0); n_grid=20)

        # Outputs
        @test result.outputs isa Vector{OutputRecord}
        @test length(result.outputs) > 0
        @test result.outputs[1].t ≈ 0.0
        @test result.outputs[end].t ≈ 5.0 rtol=1e-10

        # Diagnostics
        @test result.diagnostics isa Dict
        @test haskey(result.diagnostics, "n_steps")
        @test result.diagnostics["n_steps"] > 0
        @test haskey(result.diagnostics, "final_time")
        @test result.diagnostics["final_time"] ≈ 5.0 rtol=1e-10

        # Environment
        @test result.env isa EnvironmentalDrivers
        @test result.env.T(0.0) == 293.15
        @test result.env.ψ(0.0) == -10.0
        @test result.env.O2(0.0) == 0.3
    end

    @testset "GridInfo with different grid sizes" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        # Coarse grid
        result_10 = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0); n_grid=10)
        @test result_10.grid.n == 10
        @test length(result_10.grid.W) == 10

        # Fine grid
        result_50 = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0); n_grid=50)
        @test result_50.grid.n == 50
        @test length(result_50.grid.W) == 50

        # Different grid spacing
        @test result_10.grid.h > result_50.grid.h
    end

    @testset "Custom grid boundaries" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        # Custom r_0 and r_max
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                              n_grid=25, r_0=0.05, r_max=3.0)

        @test result.grid.r_0 ≈ 0.05
        @test result.grid.r_max ≈ 3.0
        @test result.grid.r_grid[1] ≈ 0.05
        @test result.grid.r_grid[end] ≈ 3.0
        @test result.grid.h ≈ (3.0 - 0.05) / (25 - 1)
    end
end
