# test_api.jl
# Tests for user-facing API

using Test
import SoilAggregateModel: create_initial_state, compute_total_carbon, compute_carbon_balance_error

@testset "User API" begin
    # === Test 1: Basic run with default parameters ===
    @testset "Basic run (default params)" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Constant environment
        T(t) = 293.15  # 20°C
        ψ(t) = -10.0   # kPa
        O2(t) = 0.3    # μg/mm³

        # Short 1-day run
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                              n_grid=20, dt_initial=0.01)

        # Check structure
        @test result isa SimulationResult
        @test result.grid isa GridInfo
        @test result.params isa ParameterSet
        @test result.outputs isa Vector{OutputRecord}
        @test result.diagnostics isa Dict

        # Check outputs
        @test length(result.outputs) > 0
        @test result.outputs[1].t == 0.0
        @test result.outputs[end].t ≈ 1.0 rtol=1e-10

        # Check OutputRecord structure
        rec = result.outputs[1]
        @test rec isa OutputRecord
        @test rec.t isa Float64
        @test rec.state isa AggregateState
        @test rec.mass_balance_error isa Float64

        # Check diagnostics
        @test haskey(result.diagnostics, "n_steps")
        @test haskey(result.diagnostics, "final_time")
        @test result.diagnostics["n_steps"] > 0
        @test result.diagnostics["final_time"] ≈ 1.0 rtol=1e-10
    end

    # === Test 2: Initial state creation ===
    @testset "Initial state creation" begin
        bio = BiologicalProperties()
        state = create_initial_state(10, bio)

        # Check all fields exist and are positive
        @test all(state.C .> 0.0)
        @test all(state.B .> 0.0)
        @test all(state.F_n .> 0.0)
        @test all(state.F_m .> 0.0)
        @test all(state.F_i .>= 0.0)  # F_i starts at 0 (develops over time)
        @test all(state.E .>= 0.0)   # E starts at 0 (produced by bacteria)
        @test all(state.M .>= 0.0)   # M starts at 0 (accumulates from sorption)
        @test all(state.O .> 0.0)

        # POM should match initial value
        @test state.P == bio.P_0

        # CO2 should be zero initially
        @test state.CO2_cumulative == 0.0
    end

    # === Test 3: Carbon balance diagnostic ===
    @testset "Carbon balance diagnostic" begin
        bio = BiologicalProperties()
        n = 20
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        # Create state with known total carbon
        state = create_initial_state(n, bio)

        # Initial error should be zero (all carbon in POM)
        error = compute_carbon_balance_error(state, r_grid, h, bio.P_0)

        # Compute expected: pools + POM + CO2 vs P_0
        integral = 0.0
        for i in 1:n
            C_pools = state.C[i] + state.B[i] + state.F_n[i] + state.F_m[i] +
                     state.F_i[i] + state.E[i] + state.M[i]
            W_i = 4.0 * π * r_grid[i]^2 * h
            integral += C_pools * W_i
        end
        C_total = state.P + integral + state.CO2_cumulative
        expected_error = (C_total - bio.P_0) / bio.P_0

        @test error ≈ expected_error rtol=1e-12
    end

    # === Test 4: Custom initial state ===
    @testset "Custom initial state" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Create custom initial state
        n = 10
        custom_state = AggregateState(n)
        custom_state.C .= 5.0
        custom_state.B .= 3.0
        custom_state.F_n .= 2.0
        custom_state.F_m .= 1.0
        custom_state.F_i .= 1.5
        custom_state.E .= 2.5
        custom_state.M .= 4.0
        custom_state.O .= 0.3
        custom_state.P = 500.0  # Different from bio.P_0
        custom_state.CO2_cumulative = 0.0

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 0.1);
                              n_grid=10, initial_state=custom_state)

        # Initial state should match custom state (first output)
        @test result.outputs[1].state.P == 500.0
        @test result.outputs[1].state.C[1] == 5.0
        @test result.outputs[1].state.B[1] == 3.0
    end

    # === Test 5: Time-varying environment ===
    @testset "Time-varying environment" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Time-varying temperature (warming)
        T(t) = 293.15 + 5.0 * t  # 20°C → 25°C over 1 day
        ψ(t) = -10.0 - 2.0 * t   # Drying
        O2(t) = 0.3

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                              n_grid=15)

        # Should complete without errors
        @test result.diagnostics["final_time"] ≈ 1.0 rtol=1e-10
        @test length(result.outputs) > 0

        # All states should have finite values
        for rec in result.outputs
            @test all(isfinite.(rec.state.C))
            @test all(isfinite.(rec.state.B))
            @test isfinite(rec.state.P)
            @test isfinite(rec.mass_balance_error)
        end
    end

    # === Test 6: Carbon conservation over simulation ===
    @testset "Carbon conservation (API level)" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        # Run 10-day simulation
        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 10.0);
                              n_grid=20)

        # All mass balance errors should be at machine precision
        for rec in result.outputs
            @test abs(rec.mass_balance_error) < 1e-12
        end

        # Final mass balance error
        @test abs(result.outputs[end].mass_balance_error) < 1e-12
    end

    # === Test 7: Grid resolution parameter ===
    @testset "Grid resolution" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        # Coarse grid
        result_coarse = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                                     n_grid=10)

        # Fine grid
        result_fine = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                                   n_grid=50)

        # Both should conserve carbon
        @test abs(result_coarse.outputs[end].mass_balance_error) < 1e-12
        @test abs(result_fine.outputs[end].mass_balance_error) < 1e-12

        # Results should be qualitatively similar
        # (POM should decrease in both)
        @test result_coarse.outputs[end].state.P < bio.P_0
        @test result_fine.outputs[end].state.P < bio.P_0
    end

    # === Test 8: Timestep parameters ===
    @testset "Timestep control" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        # Small timesteps (should take more steps)
        result_small = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                                    n_grid=15, dt_max=0.01)

        # Larger timesteps allowed
        result_large = run_aggregate(bio, soil, T, ψ, O2, (0.0, 1.0);
                                    n_grid=15, dt_max=0.1)

        # Small dt should take at least as many steps (or more if dt_max is limiting)
        @test result_small.diagnostics["n_steps"] >= result_large.diagnostics["n_steps"]

        # Both should reach same final time
        @test result_small.diagnostics["final_time"] ≈ 1.0 rtol=1e-10
        @test result_large.diagnostics["final_time"] ≈ 1.0 rtol=1e-10

        # Both should conserve carbon
        @test abs(result_small.outputs[end].mass_balance_error) < 1e-12
        @test abs(result_large.outputs[end].mass_balance_error) < 1e-12
    end

    # === Test 9: Non-negativity throughout simulation ===
    @testset "Non-negativity" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 10.0);
                              n_grid=20)

        # Check all outputs for non-negativity
        for rec in result.outputs
            @test all(rec.state.C .>= 0.0)
            @test all(rec.state.B .>= 0.0)
            @test all(rec.state.F_n .>= 0.0)
            @test all(rec.state.F_m .>= 0.0)
            @test all(rec.state.F_i .>= 0.0)
            @test all(rec.state.E .>= 0.0)
            @test all(rec.state.M .>= 0.0)
            @test all(rec.state.O .>= 0.0)
            @test rec.state.P >= 0.0
            @test rec.state.CO2_cumulative >= 0.0
        end
    end

    # === Test 10: CO2 accumulation ===
    @testset "CO2 accumulation" begin
        bio = BiologicalProperties()
        soil = SoilProperties()

        T(t) = 293.15
        ψ(t) = -10.0
        O2(t) = 0.3

        result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                              n_grid=20)

        # CO2 should be zero initially
        @test result.outputs[1].state.CO2_cumulative == 0.0

        # CO2 should accumulate over time
        @test result.outputs[end].state.CO2_cumulative > 0.0

        # CO2 should be monotonically increasing
        for i in 2:length(result.outputs)
            @test result.outputs[i].state.CO2_cumulative >= result.outputs[i-1].state.CO2_cumulative
        end
    end
end
