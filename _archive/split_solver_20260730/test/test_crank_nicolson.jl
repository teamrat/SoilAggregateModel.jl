# test_crank_nicolson.jl
# Tests for Crank-Nicolson diffusion solver with spherical geometry

using Test
import SoilAggregateModel: crank_nicolson_step!, BoundaryCondition,
                           neumann_zero, neumann_flux, dirichlet

@testset "Crank-Nicolson Spherical Diffusion" begin
    # === Test 1: Zero-flux transient (smooth initial condition) ===
    @testset "Zero-flux transient smoothing" begin
        # Start with a peaked distribution, verify it smooths out while conserving mass

        n = 50
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        # Initial condition (peak in the middle)
        r_mid = (r_0 + r_max) / 2.0
        u_initial = [10.0 * exp(-((r_grid[i] - r_mid) / 0.2)^2) for i in 1:n]
        u = copy(u_initial)

        # Constant diffusion
        D = fill(2.0, n)

        # Workspace
        lower = Vector{Float64}(undef, n-1)
        diag = Vector{Float64}(undef, n)
        upper = Vector{Float64}(undef, n-1)
        rhs = Vector{Float64}(undef, n)

        # Compute initial mass
        mass_initial = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])
        std_initial = sqrt(sum([(u[i] - 5.0)^2 for i in 1:n]) / n)

        # Run diffusion
        dt = 0.01
        for step in 1:200
            crank_nicolson_step!(u, D, dt/2, r_grid, h, lower, diag, upper, rhs,
                                neumann_zero, neumann_zero, 0.0, 0.0)
        end

        # Compute final mass and spread
        mass_final = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])
        std_final = sqrt(sum([(u[i] - 5.0)^2 for i in 1:n]) / n)

        # Mass conserved
        @test mass_final ≈ mass_initial rtol=1e-10

        # Distribution has spread out (std decreased)
        @test std_final < std_initial
    end

    # === Test 2: Neumann flux BC (functional test) ===
    @testset "Neumann flux BC functional" begin
        # Test that Neumann flux BC doesn't crash and produces sensible results

        n = 50
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        # Parameters
        D_val = 10.0
        J = 0.5

        # Initial condition
        u = fill(5.0, n)
        D = fill(D_val, n)

        # Workspace
        lower = Vector{Float64}(undef, n-1)
        diag = Vector{Float64}(undef, n)
        upper = Vector{Float64}(undef, n-1)
        rhs = Vector{Float64}(undef, n)

        # Initial mass
        mass_initial = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])

        # Run for a while
        dt = 0.01
        for step in 1:50
            crank_nicolson_step!(u, D, dt/2, r_grid, h, lower, diag, upper, rhs,
                                neumann_flux, neumann_zero, J, 0.0)
        end

        # Final mass
        mass_final = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])

        # Mass should have increased (flux is coming in)
        @test mass_final > mass_initial

        # No NaNs or Infs
        @test all(isfinite.(u))

        # Outer boundary should have zero gradient
        grad_outer = abs((u[n] - u[n-1]) / h)
        @test grad_outer < 0.01  # Small gradient
    end

    # === Test 3: Mass conservation with zero-flux BCs ===
    @testset "Mass conservation (zero-flux)" begin
        # Time-dependent diffusion with zero-flux at both boundaries
        # Total mass must be conserved exactly (to machine precision)

        n = 50
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)

        r_grid = [r_0 + i*h for i in 0:n-1]

        # Initial condition (Gaussian-like)
        r_mid = (r_0 + r_max) / 2.0
        u = [10.0 * exp(-((r_grid[i] - r_mid) / 0.3)^2) for i in 1:n]

        # Constant diffusion
        D = fill(1.0, n)

        # Compute initial mass (integrate over spherical shell)
        # M = ∫ u(r) × 4πr² dr ≈ Σ u_i × 4πr_i² × h
        mass_initial = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])

        # Workspace
        lower = Vector{Float64}(undef, n-1)
        diag = Vector{Float64}(undef, n)
        upper = Vector{Float64}(undef, n-1)
        rhs = Vector{Float64}(undef, n)

        # Run diffusion for many steps
        dt = 0.01
        for step in 1:100
            crank_nicolson_step!(u, D, dt/2, r_grid, h, lower, diag, upper, rhs,
                                neumann_zero, neumann_zero, 0.0, 0.0)
        end

        # Compute final mass
        mass_final = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])

        # Mass should be conserved to machine precision
        @test mass_final ≈ mass_initial rtol=1e-10

        # Also verify that concentration has smoothed out (diffusion occurred)
        # Standard deviation should decrease
        u_mean = mass_final / sum([4.0 * π * r_grid[i]^2 * h for i in 1:n])
        std_final = sqrt(sum([(u[i] - u_mean)^2 * 4.0 * π * r_grid[i]^2 * h for i in 1:n]) / mass_final)

        # Initial std (for Gaussian with σ=0.3)
        u_initial = [10.0 * exp(-((r_grid[i] - r_mid) / 0.3)^2) for i in 1:n]
        u_mean_initial = mass_initial / sum([4.0 * π * r_grid[i]^2 * h for i in 1:n])
        std_initial = sqrt(sum([(u_initial[i] - u_mean_initial)^2 * 4.0 * π * r_grid[i]^2 * h for i in 1:n]) / mass_initial)

        # Diffusion should reduce standard deviation
        @test std_final < std_initial
    end

    # === Test 4: Dirichlet BC (oxygen-style) ===
    @testset "Dirichlet BC at outer boundary" begin
        # Simulate oxygen diffusion: zero flux at inner, Dirichlet at outer

        n = 50
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)

        r_grid = [r_0 + i*h for i in 0:n-1]

        # Initial condition (depleted oxygen)
        u = fill(0.1, n)

        # Constant diffusion
        D = fill(10.0, n)  # Fast diffusion (like oxygen in air)

        # Boundary value
        O_amb = 0.3  # μg/mm³

        # Workspace
        lower = Vector{Float64}(undef, n-1)
        diag = Vector{Float64}(undef, n)
        upper = Vector{Float64}(undef, n-1)
        rhs = Vector{Float64}(undef, n)

        # Run to steady state
        dt = 0.01
        for step in 1:1000
            crank_nicolson_step!(u, D, dt/2, r_grid, h, lower, diag, upper, rhs,
                                neumann_zero, dirichlet, 0.0, O_amb)
        end

        # At steady state with zero inner flux and Dirichlet outer:
        # u should be uniform at O_amb (no sources/sinks)
        @test u ≈ fill(O_amb, n) rtol=1e-7

        # Specifically check boundary conditions
        @test u[n] ≈ O_amb rtol=1e-12  # Dirichlet exactly enforced
        @test abs((u[2] - u[1]) / h) < 1e-6  # Zero flux at inner
    end

    # === Test 5: Spatially varying diffusion coefficient ===
    @testset "Spatially varying D" begin
        # Test with D(r) that varies linearly with radius
        # This tests face-averaging: D_{i+1/2} = (D_i + D_{i+1})/2

        n = 50
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)

        r_grid = [r_0 + i*h for i in 0:n-1]

        # D(r) = D_min + (D_max - D_min) × (r - r_0)/(r_max - r_0)
        D_min = 0.5
        D_max = 5.0
        D = [D_min + (D_max - D_min) * (r_grid[i] - r_0) / (r_max - r_0) for i in 1:n]

        # Initial condition
        u = fill(1.0, n)

        # Workspace
        lower = Vector{Float64}(undef, n-1)
        diag = Vector{Float64}(undef, n)
        upper = Vector{Float64}(undef, n-1)
        rhs = Vector{Float64}(undef, n)

        # Run a few steps (just check it doesn't crash and conserves mass)
        mass_initial = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])

        dt = 0.01
        for step in 1:10
            crank_nicolson_step!(u, D, dt/2, r_grid, h, lower, diag, upper, rhs,
                                neumann_zero, neumann_zero, 0.0, 0.0)
        end

        mass_final = sum([u[i] * 4.0 * π * r_grid[i]^2 * h for i in 1:n])

        # Mass conservation
        @test mass_final ≈ mass_initial rtol=1e-10

        # No crashes = success!
    end

    # === Test 6: EXACT discrete flux balance at the Neumann boundary ===
    #
    # Test 2 above asserts only `mass_final > mass_initial`, which passes for
    # ANY inflow magnitude — a factor-of-2 error in the ghost
    # spacing (the source docstring says 2h, the code uses h) was invisible.
    #
    # The scheme injects exactly 4π·r_0²·J per unit time through the inner face,
    # so with a zero-flux outer boundary the mass gain is known in closed form.
    # Verified against an independent reimplementation: agreement to 7e-13.
    @testset "Neumann flux BC — exact discrete mass injection" begin
        n = 50
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        u = fill(0.1, n)
        D = fill(1.0, n)
        J = 0.5
        dt_half = 0.005
        n_steps = 50

        lower = zeros(n-1); diag = zeros(n); upper = zeros(n-1); rhs = zeros(n)
        W = [4.0 * π * r_grid[i]^2 * h for i in 1:n]
        mass_initial = sum(u .* W)

        for _ in 1:n_steps
            crank_nicolson_step!(u, D, dt_half, r_grid, h, lower, diag, upper, rhs,
                                 neumann_flux, neumann_zero, J, 0.0)
        end

        mass_final = sum(u .* W)
        expected_gain = 4.0 * π * r_0^2 * J * dt_half * n_steps

        @test mass_final - mass_initial ≈ expected_gain rtol=1e-9
    end

    # === Test 7: analytic steady state under flux-in / Dirichlet-out ===
    #
    # This is requirement (2) from the architecture review quoted verbatim in
    # the header of test_tridiagonal.jl. It was specified and never implemented.
    #
    #     u(r) = u_out + (J·r_0²/D)·(1/r − 1/r_max)
    #
    # Uses theta = 1.0 (Backward Euler) purely to reach steady state quickly —
    # the steady state solves L[u] = 0 and is independent of theta. This is also
    # the only test in the suite that exercises the theta = 1.0 path, which
    # diffusion_step.jl uses for oxygen.
    #
    # Verified against an independent reimplementation: max relative error
    # 1.256e-2 (n=50), 3.190e-3 (n=100), 7.970e-4 (n=200) — second order.
    @testset "Analytic spherical steady state (Neumann flux + Dirichlet)" begin
        r_0 = 0.1
        r_max = 2.0
        J = 0.5
        D_const = 1.0
        u_out = 0.0

        errors = Float64[]

        for n in (50, 100, 200)
            h = (r_max - r_0) / (n - 1)
            r_grid = [r_0 + i*h for i in 0:n-1]
            u = fill(u_out, n)
            D = fill(D_const, n)
            lower = zeros(n-1); diag = zeros(n); upper = zeros(n-1); rhs = zeros(n)

            # Backward Euler to steady state (converges in ~12 steps; 60 is ample)
            for _ in 1:60
                crank_nicolson_step!(u, D, 5.0, r_grid, h, lower, diag, upper, rhs,
                                     neumann_flux, dirichlet, J, u_out; theta=1.0)
            end

            analytic = [u_out + (J * r_0^2 / D_const) * (1.0/r_grid[i] - 1.0/r_max)
                        for i in 1:n]
            # Outer node is the Dirichlet value itself (analytic = 0 there)
            rel_err = maximum(abs(u[i] - analytic[i]) / analytic[i] for i in 1:n-1)
            push!(errors, rel_err)
        end

        # Absolute accuracy at the coarsest grid
        @test errors[1] < 0.02

        # Second-order spatial convergence: halving h must cut the error ~4x.
        # Anything materially below 3.5 means the operator or a boundary row is
        # not second-order accurate.
        @test errors[1] / errors[2] > 3.5
        @test errors[2] / errors[3] > 3.5
    end
end
