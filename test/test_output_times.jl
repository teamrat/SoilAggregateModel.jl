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
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20, output_times=target)
        actual = [rec.t for rec in result.outputs]
        @test length(actual) == length(target)
        for (a, e) in zip(actual, target)
            @test abs(a - e) < 1e-10
        end
    end

    @testset "Out-of-range times filtered" begin
        # The requested list IS the contract: out-of-range entries are dropped and
        # nothing is added. The archived split solver also injected t_start and
        # t_end, so the same call returned four records there and two here. That
        # divergence was invisible while both solvers existed; with one, the
        # caller's list is what comes back. Every call in paper/ already includes
        # its endpoints, so nothing downstream changes.
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (5.0, 20.0);
                              n_grid=20, output_times=[0.0, 3.0, 10.0, 15.0, 25.0])
        actual = [rec.t for rec in result.outputs]
        @test actual ≈ [10.0, 15.0]
        @test any(t -> abs(t - 10.0) < 1e-10, actual)
        @test any(t -> abs(t - 15.0) < 1e-10, actual)
        @test !any(t -> abs(t - 3.0) < 1e-10, actual)
        @test !any(t -> abs(t - 25.0) < 1e-10, actual)
    end

    @testset "Unordered input sorted" begin
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 10.0);
                              n_grid=20, output_times=[5.0, 2.0, 8.0])
        actual = [rec.t for rec in result.outputs]
        @test issorted(actual)
    end

    @testset "Duplicates removed" begin
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 10.0);
                              n_grid=20, output_times=[0.0, 5.0, 5.0, 10.0])
        actual = [rec.t for rec in result.outputs]
        @test length(actual) == 3  # 0, 5, 10
    end

    @testset "Dense-sparse schedule" begin
        dense_sparse = vcat(collect(0.0:0.1:1.0), collect(5.0:5.0:30.0))
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20, output_times=dense_sparse)
        @test length(result.outputs) == length(unique(dense_sparse))
    end

    @testset "Conservation at all output times" begin
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20, output_times=Float64.(0:1:30))
        for rec in result.outputs
            @test isnan(rec.mass_balance_error)
        end
    end
end
