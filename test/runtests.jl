"""
Test suite for SoilAggregateModel

Run with: julia --project test/runtests.jl
"""

using Test
using SoilAggregateModel

@testset "SoilAggregateModel.jl" begin
    @testset "Types" begin
        include("test_types.jl")
    end

    @testset "Parameters" begin
        include("test_parameters.jl")
    end

    @testset "Environment" begin
        include("test_environment.jl")
    end

    @testset "Temperature" begin
        include("test_temperature.jl")
    end

    @testset "Tridiagonal Solver" begin
        include("test_tridiagonal.jl")
    end

    @testset "Physics" begin
        include("test_physics.jl")
    end

    @testset "Biology" begin
        include("test_biology.jl")
    end

    @testset "POM Dissolution" begin
        include("test_pom.jl")
    end

    @testset "Crank-Nicolson Diffusion" begin
        include("test_crank_nicolson.jl")
    end

    @testset "Reaction Step" begin
        include("test_reactions.jl")
    end

    @testset "Time-stepping Integration" begin
        include("test_timestepper.jl")
    end

    @testset "User API" begin
        include("test_api.jl")
    end

    @testset "Result Structures" begin
        include("test_result_struct.jl")
    end

    @testset "Output Times" begin
        include("test_output_times.jl")
    end

    @testset "Post-processing" begin
        include("test_postprocessing.jl")
    end
end
