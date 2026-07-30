"""
Test suite for SoilAggregateModel

Run with: julia --project test/runtests.jl
"""

using Test
using SoilAggregateModel
using Logging

# Info-level logging is suppressed for the suite. `domain_tessellation` and
# `create_initial_state` each emit one diagnostic per call: useful in a
# production run, five identical repetitions here. Warn and above still reach
# the console — a solver warning during the suite is signal, not noise — which
# is why this is `ConsoleLogger(stderr, Warn)` and not `disable_logging`, since
# the latter short-circuits before any logger sees the message and would break
# the `@test_logs` assertion in test_tessellation.jl. `@test_logs` blocks
# install their own logger for their scope and are unaffected by this.
with_logger(ConsoleLogger(stderr, Logging.Warn)) do

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

        @testset "Physics" begin
            include("test_physics.jl")
        end

        @testset "Particle Size Distribution" begin
            include("test_particle_size.jl")
        end

        @testset "Domain Tessellation" begin
            include("test_tessellation.jl")
        end

        @testset "Biology" begin
            include("test_biology.jl")
        end

        @testset "POM Dissolution" begin
            include("test_pom.jl")
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

        @testset "Population Upscaling" begin
            include("test_population.jl")
        end

        @testset "Method of Lines" begin
            include("test_mol.jl")
        end

        @testset "Carbon Closure" begin
            include("test_closure.jl")
        end

        @testset "Stiff Solver" begin
            include("test_stiff.jl")
        end
    end
end
