# test_timestepper.jl
# Integration tests for main time-stepping loop

using Test
import SoilAggregateModel: AggregateState, Workspace, run_simulation, adapt_timestep

@testset "Time-stepping Integration" begin
    # === Test 1: Short integration (1 day) ===
    @testset "1-day integration" begin
        # Set up grid
        n = 20
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        # Create initial state
        state = AggregateState(n)
        state.C .= 10.0
        state.B .= 2.0
        state.F_n .= 1.0
        state.F_m .= 0.5
        state.F_i .= 0.8
        state.E .= 3.0
        state.M .= 5.0
        state.O .= 0.3
        state.P = 1000.0
        state.CO2_cumulative = 0.0

        # Create temperature cache
        temp_cache = TemperatureCache()

        # Create workspace
        workspace = Workspace(
            Vector{Float64}(undef, n-1),  # lower
            Vector{Float64}(undef, n),    # diag
            Vector{Float64}(undef, n-1),  # upper
            Vector{Float64}(undef, n),    # rhs
            fill(0.3, n),                  # θ
            fill(0.2, n),                  # θ_a
            Vector{Float64}(undef, n),    # D_C
            Vector{Float64}(undef, n),    # D_B
            Vector{Float64}(undef, n),    # D_Fn
            Vector{Float64}(undef, n),    # D_Fm
            Vector{Float64}(undef, n),    # D_O
            temp_cache
        )

        # Parameters
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Environmental drivers (constant)
        T_func(t) = 298.15  # K (25°C)
        ψ_func(t) = -10.0   # kPa
        O2_func(t) = 0.3    # μg/mm³

        # Save initial state for comparison
        C_initial = copy(state.C)
        P_initial = state.P
        CO2_initial = state.CO2_cumulative

        # Run simulation
        result = run_simulation(
            state, workspace, r_grid, h, bio, soil,
            T_func, ψ_func, O2_func,
            0.0, 1.0, 0.01,  # t_start, t_end, dt_initial
            output_interval=0.25
        )

        # Check that simulation completed
        @test result.diagnostics["final_time"] ≈ 1.0 rtol=1e-10
        @test result.diagnostics["n_steps"] > 0

        # Check that outputs were saved
        @test length(result.times) >= 5  # At least start + 4 quarter-day outputs + end
        @test result.times[1] == 0.0
        @test result.times[end] ≈ 1.0 rtol=1e-10

        # Check that state evolved
        @test state.C != C_initial  # Should have changed
        @test state.P < P_initial   # POM should decrease
        @test state.CO2_cumulative > CO2_initial  # CO2 should accumulate

        # Check non-negativity
        @test all(state.C .>= 0.0)
        @test all(state.B .>= 0.0)
        @test all(state.F_n .>= 0.0)
        @test all(state.F_m .>= 0.0)
        @test all(state.F_i .>= 0.0)
        @test all(state.E .>= 0.0)
        @test all(state.M .>= 0.0)
        @test all(state.O .>= 0.0)
        @test state.P >= 0.0

        # Check finite values
        @test all(isfinite.(state.C))
        @test all(isfinite.(state.B))
        @test isfinite(state.P)
        @test isfinite(state.CO2_cumulative)
    end

    # === Test 2: Carbon conservation over time ===
    @testset "Carbon conservation (30-day)" begin
        # Set up grid
        n = 20
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        # Create initial state
        state = AggregateState(n)
        state.C .= 10.0
        state.B .= 2.0
        state.F_n .= 1.0
        state.F_m .= 0.5
        state.F_i .= 0.8
        state.E .= 3.0
        state.M .= 5.0
        state.O .= 0.3
        state.P = 1000.0
        state.CO2_cumulative = 0.0

        # Compute total carbon using conservation weights (4πr²h)
        # P + ∑(C+B+F_n+F_m+F_i+E+M)×4πr_i²h + CO2
        function total_carbon(state, r_grid, h)
            integral = 0.0
            for i in 1:length(state.C)
                C_pools = state.C[i] + state.B[i] + state.F_n[i] + state.F_m[i] +
                         state.F_i[i] + state.E[i] + state.M[i]
                integral += C_pools * 4.0 * π * r_grid[i]^2 * h
            end
            return state.P + integral + state.CO2_cumulative
        end

        C_total_initial = total_carbon(state, r_grid, h)

        # Create workspace
        temp_cache = TemperatureCache()
        workspace = Workspace(
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            fill(0.3, n), fill(0.2, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n),
            temp_cache
        )

        # Parameters
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Environmental drivers
        T_func(t) = 298.15
        ψ_func(t) = -10.0
        O2_func(t) = 0.3

        # Run 30-day simulation
        result = run_simulation(
            state, workspace, r_grid, h, bio, soil,
            T_func, ψ_func, O2_func,
            0.0, 30.0, 0.01,
            output_interval=5.0
        )

        # Compute final total carbon
        C_total_final = total_carbon(state, r_grid, h)

        # Carbon should be conserved to machine precision
        # Conservation weights W_i = 4πr_i²h ensure exact telescoping
        @test C_total_final ≈ C_total_initial rtol=1e-12

        # CO2 should have accumulated
        @test state.CO2_cumulative > 0.0

        # POM should have decreased
        @test state.P < 1000.0
    end

    # === Test 3: Adaptive timestep behavior ===
    @testset "Adaptive timestep" begin
        n = 10
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        state = AggregateState(n)
        state.C .= 10.0
        state.B .= 2.0
        state.F_n .= 1.0
        state.F_m .= 0.5
        state.F_i .= 0.8
        state.E .= 3.0
        state.M .= 5.0
        state.O .= 0.3
        state.P = 1000.0
        state.CO2_cumulative = 0.0

        temp_cache = TemperatureCache()
        workspace = Workspace(
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            fill(0.3, n), fill(0.2, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n),
            temp_cache
        )

        bio = BiologicalProperties()
        soil = SoilProperties()

        T_func(t) = 298.15
        ψ_func(t) = -10.0
        O2_func(t) = 0.3

        # Test adapt_timestep function directly
        # Test halving (max_rel_change > 0.10)
        dt = 0.01
        dt_new = adapt_timestep(0.15, dt, 1e-4, 0.1)
        @test dt_new == dt / 2.0

        # Test doubling (max_rel_change < 0.01)
        dt_new = adapt_timestep(0.005, dt, 1e-4, 0.1)
        @test dt_new == dt * 2.0

        # Test no change (0.01 ≤ max_rel_change ≤ 0.10)
        dt_new = adapt_timestep(0.05, dt, 1e-4, 0.1)
        @test dt_new == dt

        # Test bounds enforcement (lower)
        dt_new = adapt_timestep(0.20, 1e-4, 1e-4, 0.1)
        @test dt_new == 1e-4

        # Test bounds enforcement (upper)
        dt_new = adapt_timestep(0.001, 0.1, 1e-4, 0.1)
        @test dt_new == 0.1
    end

    # === Test 4: Monotonic POM decrease ===
    @testset "POM monotonic decrease" begin
        n = 10
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        state = AggregateState(n)
        state.C .= 10.0
        state.B .= 5.0  # High bacteria → fast POM dissolution
        state.F_n .= 3.0
        state.F_m .= 0.5
        state.F_i .= 0.8
        state.E .= 3.0
        state.M .= 5.0
        state.O .= 0.3
        state.P = 1000.0
        state.CO2_cumulative = 0.0

        temp_cache = TemperatureCache()
        workspace = Workspace(
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            fill(0.3, n), fill(0.2, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            Vector{Float64}(undef, n),
            temp_cache
        )

        bio = BiologicalProperties()
        soil = SoilProperties()

        T_func(t) = 298.15
        ψ_func(t) = -10.0
        O2_func(t) = 0.3

        # Run simulation with frequent outputs
        result = run_simulation(
            state, workspace, r_grid, h, bio, soil,
            T_func, ψ_func, O2_func,
            0.0, 10.0, 0.01,
            output_interval=1.0
        )

        # Check that POM decreases monotonically
        for i in 2:length(result.states)
            P_prev = result.states[i-1].P
            P_curr = result.states[i].P
            @test P_curr <= P_prev  # Should not increase
        end

        # Final POM should be less than initial (even if only slightly)
        # POM dissolution is slow with default parameters
        @test result.states[end].P < result.states[1].P
    end
end
