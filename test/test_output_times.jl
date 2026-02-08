# test_output_times.jl
# Tests for user-specified output times

using Test

@testset "User-specified output times" begin
    bio = BiologicalProperties()
    soil = SoilProperties()

    # Constant environment functions
    T(t) = 293.15  # 20°C
    ψ(t) = -10.0   # kPa
    O2(t) = 0.3    # μg/mm³

    @testset "Exact landing" begin
        target = [0.0, 1.0, 5.0, 10.0, 15.0, 30.0]
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20, output_times=target)
        actual = [rec.t for rec in result.outputs]
        @test length(actual) == length(target)
        for (a, e) in zip(actual, target)
            @test abs(a - e) < 1e-10
        end
    end

    @testset "Out-of-range times filtered" begin
        result = run_aggregate(bio, soil, T, ψ, O2, (5.0, 20.0);
                              n_grid=20, output_times=[0.0, 3.0, 10.0, 15.0, 25.0])
        actual = [rec.t for rec in result.outputs]
        @test actual[1] ≈ 5.0
        @test actual[end] ≈ 20.0
        @test any(t -> abs(t - 10.0) < 1e-10, actual)
        @test any(t -> abs(t - 15.0) < 1e-10, actual)
        @test !any(t -> abs(t - 3.0) < 1e-10, actual)
        @test !any(t -> abs(t - 25.0) < 1e-10, actual)
    end

    @testset "Unordered input sorted" begin
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 10.0);
                              n_grid=20, output_times=[5.0, 2.0, 8.0])
        actual = [rec.t for rec in result.outputs]
        @test issorted(actual)
    end

    @testset "Duplicates removed" begin
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 10.0);
                              n_grid=20, output_times=[0.0, 5.0, 5.0, 10.0])
        actual = [rec.t for rec in result.outputs]
        @test length(actual) == 3  # 0, 5, 10
    end

    @testset "Dense-sparse schedule" begin
        dense_sparse = vcat(collect(0.0:0.1:1.0), collect(5.0:5.0:30.0))
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20, output_times=dense_sparse)
        @test length(result.outputs) == length(unique(dense_sparse))
    end

    @testset "Conservation at all output times" begin
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20, output_times=Float64.(0:1:30))
        for rec in result.outputs
            @test abs(rec.mass_balance_error) < 1e-12
        end
    end
end
