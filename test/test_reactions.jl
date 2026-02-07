# test_reactions.jl
# Tests for reaction step (reactions.jl + reaction_step.jl)

using Test
import SoilAggregateModel: AggregateState, Workspace, TemperatureCache, SourceTerms,
                           compute_source_terms, reaction_step!

@testset "Reaction Step" begin
    # === Test 1: compute_source_terms returns correct struct ===
    @testset "Source terms structure" begin
        # Set up state
        C = 10.0
        B = 2.0
        F_n = 1.0
        F_m = 0.5
        F_i = 0.8
        E = 3.0
        M = 5.0
        O = 0.3
        θ = 0.3
        θ_a = 0.2
        ψ = -10.0

        # Parameters
        bio = BiologicalProperties()
        soil = SoilProperties()

        # Temperature cache (set to reference temperature for simplicity)
        T = 298.15  # K
        temp_cache = TemperatureCache()
        temp_cache.f_bac = 1.0
        temp_cache.f_fun = 1.0
        temp_cache.f_eps = 1.0
        temp_cache.f_maoc_s = 1.0
        temp_cache.f_maoc_d = 1.0
        temp_cache.f_pom = 1.0
        temp_cache.D_O2_w = 1.0
        temp_cache.D_DOC_w = 0.5
        temp_cache.D_O2_a = 100.0
        temp_cache.D_Fm = 0.1
        temp_cache.K_H_O = 30.0

        # Compute source terms
        sources = compute_source_terms(C, B, F_n, F_m, F_i, E, M, O, θ, θ_a, ψ,
                                       bio, soil, temp_cache)

        # Check that all fields exist and are finite
        @test isfinite(sources.S_C)
        @test isfinite(sources.S_B)
        @test isfinite(sources.S_Fn)
        @test isfinite(sources.S_Fm)
        @test isfinite(sources.S_Fi)
        @test isfinite(sources.S_E)
        @test isfinite(sources.S_M)
        @test isfinite(sources.S_O)
        @test isfinite(sources.Resp_total)

        # Respiration should be non-negative
        @test sources.Resp_total >= 0.0
    end

    # === Test 2: Carbon conservation (via respiration) ===
    @testset "Carbon conservation" begin
        # This test verifies the key identity from test_biology.jl:
        # Carbon_in_out + Carbon_respired = 0 (to machine precision)

        C = 15.0
        B = 3.0
        F_n = 1.5
        F_m = 0.7
        F_i = 1.2
        E = 4.0
        M = 8.0
        O = 0.25
        θ = 0.35
        θ_a = 0.15
        ψ = -15.0

        bio = BiologicalProperties()
        soil = SoilProperties()

        temp_cache = TemperatureCache()
        temp_cache.f_bac = 1.0
        temp_cache.f_fun = 1.0
        temp_cache.f_eps = 1.0
        temp_cache.f_maoc_s = 1.0
        temp_cache.f_maoc_d = 1.0
        temp_cache.f_pom = 1.0
        temp_cache.D_O2_w = 1.0
        temp_cache.D_DOC_w = 0.5
        temp_cache.D_O2_a = 100.0
        temp_cache.D_Fm = 0.1
        temp_cache.K_H_O = 30.0

        sources = compute_source_terms(C, B, F_n, F_m, F_i, E, M, O, θ, θ_a, ψ,
                                       bio, soil, temp_cache)

        # Carbon in/out of pools (excluding O and respiration)
        Carbon_in_out = sources.S_C + sources.S_B + sources.S_Fn + sources.S_Fm +
                       sources.S_Fi + sources.S_E + sources.S_M

        # Carbon respired (negative of Resp_total)
        Carbon_respired = -sources.Resp_total

        # Identity: Carbon_in_out = -Resp_total
        # All carbon flow into/out of pools must equal respiration loss
        @test Carbon_in_out ≈ Carbon_respired atol=1e-12
    end

    # === Test 3: reaction_step! advances state ===
    @testset "State advancement" begin
        n = 10
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        # Create state
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
        temp_cache = TemperatureCache(
            1.0,   # f_bac
            1.0,   # f_fun
            1.0,   # f_eps
            1.0,   # f_maoc_s
            1.0,   # f_maoc_d
            1.0,   # f_pom
            1.0,   # D_O2_w
            0.5,   # D_DOC_w
            100.0, # D_O2_a
            0.1,   # D_Fm
            30.0   # K_H_O
        )

        # Create workspace with filled arrays
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
        ψ = -10.0

        # Save initial state
        C_initial = copy(state.C)
        P_initial = state.P
        CO2_initial = state.CO2_cumulative

        # Compute POM dissolution rate
        B_0 = state.B[1]
        F_n_0 = state.F_n[1]
        θ_0 = workspace.θ[1]
        θ_a_0 = workspace.θ_a[1]
        O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
        R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
        J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                     bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
        R_P_val = R_P(J_P_val, bio.r_0)

        # Perform reaction step
        dt = 0.01  # day
        reaction_step!(state, workspace, dt, r_grid, h, bio, soil, ψ, R_P_val)

        # State should have changed
        @test state.C != C_initial  # At least one element different
        @test state.P != P_initial
        @test state.CO2_cumulative > CO2_initial

        # All state variables should be non-negative
        @test all(state.C .>= 0.0)
        @test all(state.B .>= 0.0)
        @test all(state.F_n .>= 0.0)
        @test all(state.F_m .>= 0.0)
        @test all(state.F_i .>= 0.0)
        @test all(state.E .>= 0.0)
        @test all(state.M .>= 0.0)
        @test all(state.O .>= 0.0)
        @test state.P >= 0.0
    end

    # === Test 4: POM dissolution ===
    @testset "POM dissolution" begin
        n = 10
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        state = AggregateState(n)
        state.C .= 10.0
        state.B .= 5.0      # High bacteria at surface
        state.F_n .= 3.0    # High fungi at surface
        state.F_m .= 0.5
        state.F_i .= 0.8
        state.E .= 3.0
        state.M .= 5.0
        state.O .= 0.3      # Sufficient oxygen
        state.P = 1000.0
        state.CO2_cumulative = 0.0

        # Temperature cache
        temp_cache = TemperatureCache(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 100.0, 0.1, 30.0)
        workspace = Workspace(
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n), Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            fill(0.3, n), fill(0.2, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            temp_cache
        )

        bio = BiologicalProperties()
        soil = SoilProperties()
        ψ = -10.0

        P_initial = state.P

        # Compute POM dissolution rate
        B_0 = state.B[1]
        F_n_0 = state.F_n[1]
        θ_0 = workspace.θ[1]
        θ_a_0 = workspace.θ_a[1]
        O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
        R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
        J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                     bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
        R_P_val = R_P(J_P_val, bio.r_0)

        dt = 0.1
        reaction_step!(state, workspace, dt, r_grid, h, bio, soil, ψ, R_P_val)

        # POM should have decreased (dissolved)
        @test state.P < P_initial

        # POM decrease should be reasonable (not too fast)
        @test state.P > P_initial * 0.9  # Less than 10% loss in 0.1 day
    end

    # === Test 5: CO₂ accumulation ===
    @testset "CO2 accumulation" begin
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

        # Temperature cache
        temp_cache = TemperatureCache(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 100.0, 0.1, 30.0)
        workspace = Workspace(
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n), Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            fill(0.3, n), fill(0.2, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            temp_cache
        )

        bio = BiologicalProperties()
        soil = SoilProperties()
        ψ = -10.0

        # Run multiple steps
        dt = 0.01
        for step in 1:10
            # Compute POM dissolution rate
            B_0 = state.B[1]
            F_n_0 = state.F_n[1]
            θ_0 = workspace.θ[1]
            θ_a_0 = workspace.θ_a[1]
            O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
            R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
            J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                         bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
            R_P_val = R_P(J_P_val, bio.r_0)

            reaction_step!(state, workspace, dt, r_grid, h, bio, soil, ψ, R_P_val)
        end

        # CO₂ should have accumulated
        @test state.CO2_cumulative > 0.0

        # CO₂ should be finite and reasonable
        @test isfinite(state.CO2_cumulative)
    end

    # === Test 6: Non-negativity enforcement ===
    @testset "Non-negativity enforcement" begin
        n = 5
        r_0 = 0.1
        r_max = 2.0
        h = (r_max - r_0) / (n - 1)
        r_grid = [r_0 + i*h for i in 0:n-1]

        state = AggregateState(n)
        # Start with very low values
        state.C .= 0.01
        state.B .= 0.01
        state.F_n .= 0.01
        state.F_m .= 0.01
        state.F_i .= 0.01
        state.E .= 0.01
        state.M .= 0.01
        state.O .= 0.01
        state.P = 0.1
        state.CO2_cumulative = 0.0

        # Temperature cache
        temp_cache = TemperatureCache(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 100.0, 0.1, 30.0)
        workspace = Workspace(
            Vector{Float64}(undef, n-1), Vector{Float64}(undef, n), Vector{Float64}(undef, n-1), Vector{Float64}(undef, n),
            fill(0.3, n), fill(0.2, n),
            Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n),
            temp_cache
        )

        bio = BiologicalProperties()
        soil = SoilProperties()
        ψ = -10.0

        # Compute POM dissolution rate
        B_0 = state.B[1]
        F_n_0 = state.F_n[1]
        θ_0 = workspace.θ[1]
        θ_a_0 = workspace.θ_a[1]
        O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
        R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
        J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                     bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
        R_P_val = R_P(J_P_val, bio.r_0)

        # Run with large timestep (might try to make things negative)
        dt = 1.0
        reaction_step!(state, workspace, dt, r_grid, h, bio, soil, ψ, R_P_val)

        # Nothing should be negative
        @test all(state.C .>= 0.0)
        @test all(state.B .>= 0.0)
        @test all(state.F_n .>= 0.0)
        @test all(state.F_m .>= 0.0)
        @test all(state.F_i .>= 0.0)
        @test all(state.E .>= 0.0)
        @test all(state.M .>= 0.0)
        @test all(state.O .>= 0.0)
        @test state.P >= 0.0
    end
end
