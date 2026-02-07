"""
Tests for biology/*.jl modules

Strategy: Test each reaction term with known inputs, compare to hand-computed
expected values. Then test edge cases and carbon conservation.
"""

import SoilAggregateModel: TemperatureCache,
                           h_B, R_B, R_Bb, Resp_B, R_rec_B, Y_B_func, Gamma_B, Gamma_E,
                           R_F, Resp_F, h_Fi, R_rec_F, fungal_transitions,
                           Pi_protected, Y_F_const, Y_F_uptake_dependent, Gamma_F,
                           h_E, R_rec_E,
                           softplus, M_eq_langmuir_freundlich, J_M

@testset "Bacteria" begin
    @testset "h_B sigmoid" begin
        # At B = B_min: h_B = 0.5
        @test h_B(1.0, 1.0) ≈ 0.5

        # At B = 0: h_B → 0
        @test h_B(0.0, 1.0) < 0.01

        # At B >> B_min: h_B → 1
        @test h_B(10.0, 1.0) > 0.99

        # Monotonic increasing
        B_range = [0.1, 0.5, 1.0, 2.0, 5.0]
        h_vals = [h_B(B, 1.0) for B in B_range]
        @test all(diff(h_vals) .>= 0)  # >= to handle numerical precision
    end

    @testset "R_B uptake" begin
        # Known inputs
        C_aq = 10.0    # μg/mm³
        O_aq = 5.0
        B = 2.0
        ψ = -33.0      # kPa
        r_B_max_T = 0.5
        K_B = 5.0
        L_B = 2.0
        ν_B = 0.01     # POSITIVE: exp(ν·ψ) < 1 when ψ < 0 (dry)

        # Hand-compute expected
        monod_C = 10.0 / (5.0 + 10.0)  # = 2/3
        monod_O = 5.0 / (2.0 + 5.0)    # = 5/7
        water_stress = exp(0.01 * -33.0)  # = exp(-0.33) ≈ 0.7189
        expected = 0.5 * (2/3) * (5/7) * 2.0 * 0.7189  # ≈ 0.3424

        R_B_val = R_B(C_aq, O_aq, B, ψ, r_B_max_T, K_B, L_B, ν_B)
        @test R_B_val ≈ expected rtol=1e-4

        # Edge cases
        @test R_B(0.0, O_aq, B, ψ, r_B_max_T, K_B, L_B, ν_B) == 0.0  # No carbon
        @test R_B(C_aq, 0.0, B, ψ, r_B_max_T, K_B, L_B, ν_B) == 0.0  # No oxygen
        @test R_B(C_aq, O_aq, 0.0, ψ, r_B_max_T, K_B, L_B, ν_B) == 0.0  # No biomass

        # Dry soil (ψ → -∞): uptake → 0
        R_B_dry = R_B(C_aq, O_aq, B, -1500.0, r_B_max_T, K_B, L_B, ν_B)
        @test R_B_dry < 0.01 * R_B_val  # Much lower than baseline
    end

    @testset "R_Bb maintenance" begin
        C_B = 3.0
        O_aq = 5.0
        B = 2.0
        r_B_max_T = 0.5
        K_B = 5.0
        L_B = 2.0
        B_min = 1.0

        R_Bb_val = R_Bb(C_B, O_aq, B, r_B_max_T, K_B, L_B, B_min)

        # Should be positive
        @test R_Bb_val > 0.0

        # At B = B_min: h_B = 0.5 → R_Bb reduced
        R_Bb_min = R_Bb(C_B, O_aq, B_min, r_B_max_T, K_B, L_B, B_min)
        @test R_Bb_min < R_Bb_val

        # At B = 0: R_Bb → 0
        R_Bb_zero = R_Bb(C_B, O_aq, 0.0, r_B_max_T, K_B, L_B, B_min)
        @test R_Bb_zero < 0.01
    end

    @testset "Bacterial allocation" begin
        R_B_val = 1.0
        R_Bb_val = 0.4
        R_diff = R_B_val - R_Bb_val  # = 0.6
        Y_B_max = 0.6
        K_Y = 0.5
        γ = 0.3

        # Yield
        Y_B_val = Y_B_func(R_diff, Y_B_max, K_Y)
        expected_Y = 0.6 * 0.6 / (0.6 + 0.5)  # = 0.36 / 1.1 ≈ 0.3273
        @test Y_B_val ≈ expected_Y rtol=1e-4

        # Growth
        Γ_B_val = Gamma_B(R_B_val, R_Bb_val, Y_B_val, γ)
        expected_Gamma_B = Y_B_val * 0.6 * (1.0 - 0.3)  # = Y_B_val * 0.42
        @test Γ_B_val ≈ expected_Gamma_B rtol=1e-4

        # EPS production
        Γ_E_val = Gamma_E(R_B_val, R_Bb_val, Y_B_val, γ)
        expected_Gamma_E = Y_B_val * 0.6 * 0.3  # = Y_B_val * 0.18
        @test Γ_E_val ≈ expected_Gamma_E rtol=1e-4

        # Respiration
        Resp_B_val = Resp_B(R_Bb_val, R_diff, Y_B_val)
        expected_Resp = 0.4 + 0.6 * (1.0 - Y_B_val)
        @test Resp_B_val ≈ expected_Resp rtol=1e-4

        # Starvation (R_diff < 0)
        R_B_starve = 0.3
        R_diff_neg = R_B_starve - R_Bb_val  # = -0.1
        Y_B_starve = Y_B_func(R_diff_neg, Y_B_max, K_Y)
        @test Y_B_starve == 0.0  # No growth under starvation

        Γ_B_starve = Gamma_B(R_B_starve, R_Bb_val, Y_B_starve, γ)
        @test Γ_B_starve ≈ R_diff_neg  # = -0.1 (catabolism)

        Γ_E_starve = Gamma_E(R_B_starve, R_Bb_val, Y_B_starve, γ)
        @test Γ_E_starve == 0.0  # No EPS under starvation
    end

    @testset "R_rec_B death" begin
        μ_B_T = 0.02
        B = 2.0
        B_min = 1.0

        R_rec_B_val = R_rec_B(μ_B_T, B, B_min)

        # Should be positive
        @test R_rec_B_val > 0.0

        # At B = B_min: reduced by h_B = 0.5
        R_rec_at_min = R_rec_B(μ_B_T, B_min, B_min)
        @test R_rec_at_min ≈ μ_B_T * B_min * 0.5 rtol=1e-3

        # At B = 0: → 0
        @test R_rec_B(μ_B_T, 0.0, B_min) < 0.01
    end
end

@testset "Fungi" begin
    @testset "Π protection" begin
        # Normal case
        Π = Pi_protected(5.0, 3.0, 2.0, 1e-10)
        @test Π ≈ 5.0 / 5.0  # = 1.0

        # Near-zero immobile pools
        Π_protect = Pi_protected(1.0, 0.0, 0.0, 1e-10)
        @test Π_protect ≈ 1.0 / 1e-10  # Large but finite
        @test isfinite(Π_protect)
    end

    @testset "h_Fi sigmoid" begin
        # At F_i = F_i_min: h_Fi = 0.5
        @test h_Fi(1.0, 1.0) ≈ 0.5

        # At F_i = 0: h_Fi → 0
        @test h_Fi(0.0, 1.0) < 0.01

        # At F_i >> F_i_min: h_Fi → 1
        @test h_Fi(10.0, 1.0) > 0.99
    end

    @testset "R_F uptake" begin
        C_aq = 10.0
        O_aq = 5.0
        F_i = 3.0
        F_n = 2.0
        λ = 0.1
        ψ = -33.0
        r_F_max_T = 0.4
        K_F = 8.0
        L_F = 3.0
        ν_F = 0.005  # POSITIVE, smaller than ν_B (more drought-tolerant)

        # Hand-compute
        monod_C = 10.0 / (8.0 + 10.0)  # = 10/18 ≈ 0.5556
        monod_O = 5.0 / (3.0 + 5.0)    # = 5/8 = 0.625
        uptake_biomass = 3.0 + 0.1 * 2.0  # = 3.2
        water_stress = exp(0.005 * -33.0)  # = exp(-0.165)
        expected = 0.4 * monod_C * monod_O * uptake_biomass * water_stress

        R_F_val = R_F(C_aq, O_aq, F_i, F_n, λ, ψ, r_F_max_T, K_F, L_F, ν_F)
        @test R_F_val ≈ expected rtol=1e-3  # Relaxed tolerance for hand computation

        # Edge cases
        @test R_F(0.0, O_aq, F_i, F_n, λ, ψ, r_F_max_T, K_F, L_F, ν_F) == 0.0  # No carbon
        @test R_F(C_aq, 0.0, F_i, F_n, λ, ψ, r_F_max_T, K_F, L_F, ν_F) == 0.0  # No oxygen
        @test R_F(C_aq, O_aq, 0.0, 0.0, λ, ψ, r_F_max_T, K_F, L_F, ν_F) == 0.0  # No biomass
    end

    @testset "Fungal yield and allocation" begin
        # Constant yield
        Y_F = 0.4
        @test Y_F_const(Y_F) == 0.4

        # Uptake-dependent yield
        R_F_val = 0.8
        Y_F_max = 0.5
        K_YF = 0.6
        Y_F_dep = Y_F_uptake_dependent(R_F_val, Y_F_max, K_YF)
        expected_Y_F = 0.5 * 0.8 / (0.8 + 0.6)  # = 0.4 / 1.4 ≈ 0.2857
        @test Y_F_dep ≈ expected_Y_F rtol=1e-4

        # Allocation
        Γ_F_val = Gamma_F(0.4, 0.8)
        @test Γ_F_val ≈ 0.4 * 0.8  # = 0.32

        # Respiration
        Resp_F_val = Resp_F(0.8, 0.4)
        @test Resp_F_val ≈ (1.0 - 0.4) * 0.8  # = 0.48
    end

    @testset "Fungal transitions" begin
        F_i = 3.0
        F_n = 2.0
        F_m = 4.0
        ε_F = 1e-10
        Π = Pi_protected(F_m, F_i, F_n, ε_F)  # = 4/5 = 0.8

        α_i_T = 0.05
        α_n_T = 0.1
        β_i_T = 0.15
        β_n_T = 0.2
        ζ_T = 0.02
        δ = 1.5
        η = 0.9

        trans = fungal_transitions(F_i, F_n, F_m, Π, α_i_T, α_n_T, β_i_T, β_n_T,
                                   ζ_T, δ, η, ε_F)

        # Insulation: ζ·F_n
        @test trans.insulation ≈ 0.02 * 2.0  # = 0.04

        # Net transitions
        Π_delta = 0.8^1.5  # ≈ 0.7155
        net_i = (0.15 * 0.8 - 0.05 * 0.7155) * 3.0  # = (0.12 - 0.03578) * 3.0 ≈ 0.2527
        net_n = (0.2 * 0.8 - 0.1 * 0.7155) * 2.0    # = (0.16 - 0.07155) * 2.0 ≈ 0.1769

        @test trans.trans_i ≈ 0.9 * net_i rtol=1e-3
        @test trans.trans_n ≈ 0.9 * net_n rtol=1e-3

        # Conversion respiration (CRITICAL: uses abs())
        @test trans.Resp_F_conv ≈ 0.1 * abs(net_i + net_n) rtol=1e-3
        @test trans.Resp_F_conv > 0.0  # Always positive

        # Test mobilization case (Π small)
        F_m_small = 0.1
        Π_small = Pi_protected(F_m_small, F_i, F_n, ε_F)  # = 0.1/5 = 0.02
        trans_mob = fungal_transitions(F_i, F_n, F_m_small, Π_small, α_i_T, α_n_T,
                                       β_i_T, β_n_T, ζ_T, δ, η, ε_F)

        # Mobilization should dominate (α·Π^δ > β·Π when Π is small and δ > 1)
        # Resp_F_conv should still be positive due to abs()
        @test trans_mob.Resp_F_conv > 0.0
    end

    @testset "R_rec_F death" begin
        μ_F_T = 0.03
        F_i = 3.0
        F_i_min = 1.0

        R_rec_F_val = R_rec_F(μ_F_T, F_i, F_i_min)
        @test R_rec_F_val > 0.0

        # At F_i = 0: → 0
        @test R_rec_F(μ_F_T, 0.0, F_i_min) < 0.01
    end
end

@testset "EPS" begin
    @testset "h_E sigmoid" begin
        @test h_E(1.0, 1.0) ≈ 0.5
        @test h_E(0.0, 1.0) < 0.01
        @test h_E(10.0, 1.0) > 0.99
    end

    @testset "R_rec_E degradation" begin
        μ_E_max_T = 0.05
        C_aq = 5.0  # CRITICAL: uses C_aq, not C
        K_E = 10.0
        E = 8.0
        E_min = 1.0

        R_rec_E_val = R_rec_E(μ_E_max_T, C_aq, K_E, E, E_min)

        # Hand-compute
        inhibition = 10.0 / (10.0 + 5.0)  # = 2/3
        h_E_val = h_E(E, E_min)  # Should be high for E = 8.0
        expected = μ_E_max_T * inhibition * E * h_E_val

        @test R_rec_E_val ≈ expected rtol=1e-4
        @test R_rec_E_val > 0.0

        # Substrate inhibition: higher C_aq → lower R_rec_E
        R_rec_E_high_C = R_rec_E(μ_E_max_T, 50.0, K_E, E, E_min)
        @test R_rec_E_high_C < R_rec_E_val

        # Low C_aq → high degradation (scavenge EPS)
        R_rec_E_low_C = R_rec_E(μ_E_max_T, 0.1, K_E, E, E_min)
        @test R_rec_E_low_C > R_rec_E_val
    end
end

@testset "MAOC" begin
    @testset "Softplus" begin
        ε = 0.01

        # For x >> ε: softplus(x, ε) ≈ x
        @test softplus(1.0, ε) ≈ 1.0 atol=0.1*ε

        # For x << -ε: softplus(x, ε) ≈ 0
        @test softplus(-1.0, ε) < 0.1*ε

        # At x = 0: softplus(0, ε) ≈ ε·ln(2)
        @test softplus(0.0, ε) ≈ ε * log(2.0) rtol=0.1

        # Smoothness: no discontinuity
        x_vals = range(-0.1, 0.1, length=100)
        soft_vals = [softplus(x, ε) for x in x_vals]
        @test all(soft_vals .>= 0.0)

        # Numerical stability
        @test isfinite(softplus(100.0, ε))
        @test isfinite(softplus(-100.0, ε))
    end

    @testset "M_eq Langmuir-Freundlich" begin
        C_eq = 5.0  # CRITICAL: uses C_eq, not C or C_aq
        M_max = 20.0
        k_L = 0.5
        n_LF = 0.8

        M_eq_val = M_eq_langmuir_freundlich(C_eq, M_max, k_L, n_LF)

        # Hand-compute
        k_L_C = 0.5 * 5.0  # = 2.5
        k_L_C_n = 2.5^0.8  # ≈ 2.0814
        expected = 20.0 * k_L_C_n / (1.0 + k_L_C_n)  # ≈ 13.509

        @test M_eq_val ≈ expected rtol=1e-3

        # Boundary cases
        @test M_eq_langmuir_freundlich(0.0, M_max, k_L, n_LF) == 0.0
        @test M_eq_langmuir_freundlich(1000.0, M_max, k_L, n_LF) ≈ M_max rtol=0.01
    end

    @testset "J_M flux" begin
        M = 10.0
        M_eq = 15.0  # M < M_eq → sorption
        κ_s_T = 0.1
        κ_d_T = 0.05
        ε_maoc = 0.01

        J_M_val = J_M(M, M_eq, κ_s_T, κ_d_T, ε_maoc)

        # Should be positive (sorption)
        @test J_M_val > 0.0

        # Reverse case: M > M_eq → desorption
        J_M_des = J_M(20.0, 15.0, κ_s_T, κ_d_T, ε_maoc)
        @test J_M_des < 0.0

        # At equilibrium: J_M ≈ 0
        J_M_eq = J_M(15.0, 15.0, κ_s_T, κ_d_T, ε_maoc)
        @test abs(J_M_eq) < 0.1  # Small residual due to softplus

        # Smoothness across equilibrium
        M_range = range(14.0, 16.0, length=100)
        J_range = [J_M(M_val, 15.0, κ_s_T, κ_d_T, ε_maoc) for M_val in M_range]
        # Should transition smoothly from positive to negative
        @test J_range[1] > 0.0
        @test J_range[end] < 0.0
    end
end

@testset "Type stability" begin
    # Bacteria
    @inferred h_B(2.0, 1.0)
    @inferred R_B(10.0, 5.0, 2.0, -33.0, 0.5, 5.0, 2.0, 0.01)  # POSITIVE ν_B
    @inferred R_Bb(3.0, 5.0, 2.0, 0.5, 5.0, 2.0, 1.0)
    @inferred Y_B_func(0.6, 0.6, 0.5)
    @inferred Gamma_B(1.0, 0.4, 0.3, 0.3)
    @inferred Gamma_E(1.0, 0.4, 0.3, 0.3)
    @inferred Resp_B(0.4, 0.6, 0.3)
    @inferred R_rec_B(0.02, 2.0, 1.0)

    # Fungi
    @inferred Pi_protected(5.0, 3.0, 2.0, 1e-10)
    @inferred R_F(10.0, 5.0, 3.0, 2.0, 0.1, -33.0, 0.4, 8.0, 3.0, 0.005)  # POSITIVE ν_F
    @inferred Y_F_const(0.4)
    @inferred Gamma_F(0.4, 0.8)
    @inferred Resp_F(0.8, 0.4)
    @inferred h_Fi(3.0, 1.0)
    @inferred R_rec_F(0.03, 3.0, 1.0)

    # EPS
    @inferred h_E(8.0, 1.0)
    @inferred R_rec_E(0.05, 5.0, 10.0, 8.0, 1.0)

    # MAOC
    @inferred softplus(1.0, 0.01)
    @inferred M_eq_langmuir_freundlich(5.0, 20.0, 0.5, 0.8)
    @inferred J_M(10.0, 15.0, 0.1, 0.05, 0.01)
end

@testset "Carbon Conservation (THE CRITICAL TEST)" begin
    """
    This is THE most important test in the biology module.

    For any single node, carbon balance must hold:
        S_C + S_B + S_Fi + S_Fn + S_Fm + S_E + S_M = -(Resp_B + Resp_F + Resp_F_conv)

    All carbon entering or leaving the biological pools must be accounted for.
    Respiration is the only carbon sink (to CO₂).
    """

    # Set up known state
    C = 50.0
    B = 2.0
    F_i = 3.0
    F_n = 2.5
    F_m = 4.0
    E = 8.0
    M = 10.0
    O = 6.0

    # Parameters
    θ = 0.35
    θ_a = 0.10
    ψ = -33.0
    ρ_b = 1200.0  # μg/mm³
    k_d = 0.1     # mm³/μg
    K_H = 29.0

    # Biological parameters (simplified)
    r_B_max_T = 0.5
    K_B = 5.0
    L_B = 2.0
    ν_B = 0.01     # POSITIVE
    C_B = 3.0
    B_min = 1.0
    Y_B_max = 0.6
    K_Y = 0.5
    γ = 0.3
    μ_B_T = 0.02

    r_F_max_T = 0.4
    K_F = 8.0
    L_F = 3.0
    ν_F = 0.005    # POSITIVE, smaller than ν_B
    λ = 0.1
    Y_F = 0.4
    μ_F_T = 0.03
    F_i_min = 1.0

    α_i_T = 0.05
    α_n_T = 0.1
    β_i_T = 0.15
    β_n_T = 0.2
    ζ_T = 0.02
    δ = 1.5
    η = 0.9
    ε_F = 1e-10

    μ_E_max_T = 0.05
    K_E = 10.0
    E_min = 1.0

    κ_s_T = 0.1
    κ_d_T = 0.05
    ε_maoc = 0.01
    M_max = 20.0
    k_L = 0.5
    n_LF = 0.8

    α_O = 1.0  # Respiratory quotient

    # STEP 1: Compute C_aq, C_eq, O_aq ONCE (Opus ground rule #2, #3)
    C_aq = C / (θ + ρ_b * k_d)
    C_eq = k_d * C_aq
    O_aq = O * θ / (θ + K_H * θ_a)

    # STEP 2: Compute all bacterial terms
    R_B_val = R_B(C_aq, O_aq, B, ψ, r_B_max_T, K_B, L_B, ν_B)
    R_Bb_val = R_Bb(C_B, O_aq, B, r_B_max_T, K_B, L_B, B_min)
    R_diff = R_B_val - R_Bb_val
    Y_B_val = Y_B_func(R_diff, Y_B_max, K_Y)
    Γ_B_val = Gamma_B(R_B_val, R_Bb_val, Y_B_val, γ)
    Γ_E_val = Gamma_E(R_B_val, R_Bb_val, Y_B_val, γ)
    Resp_B_val = Resp_B(R_Bb_val, R_diff, Y_B_val)
    R_rec_B_val = R_rec_B(μ_B_T, B, B_min)

    # STEP 3: Compute all fungal terms
    Π = Pi_protected(F_m, F_i, F_n, ε_F)
    R_F_val = R_F(C_aq, O_aq, F_i, F_n, λ, ψ, r_F_max_T, K_F, L_F, ν_F)
    Y_F_val = Y_F_const(Y_F)
    Γ_F_val = Gamma_F(Y_F_val, R_F_val)
    Resp_F_val = Resp_F(R_F_val, Y_F_val)
    R_rec_F_val = R_rec_F(μ_F_T, F_i, F_i_min)
    trans = fungal_transitions(F_i, F_n, F_m, Π, α_i_T, α_n_T, β_i_T, β_n_T,
                               ζ_T, δ, η, ε_F)

    # STEP 4: Compute EPS term
    R_rec_E_val = R_rec_E(μ_E_max_T, C_aq, K_E, E, E_min)

    # STEP 5: Compute MAOC terms
    M_eq_val = M_eq_langmuir_freundlich(C_eq, M_max, k_L, n_LF)
    J_M_val = J_M(M, M_eq_val, κ_s_T, κ_d_T, ε_maoc)

    # STEP 6: Assemble total recycling
    R_rec = R_rec_B_val + R_rec_F_val + R_rec_E_val

    # STEP 7: Compute source terms (manuscript Section 2.4.6)
    # CRITICAL FIX: S_C MAOC coupling is -J_M (NO FACTOR)
    # The manuscript has a bug: both C and M are per bulk volume, so coupling is 1:1
    # The factor (θ + ρ_b·k_d)/k_d would over-drain C by ~1200x
    S_C = -R_B_val - R_F_val + R_rec - J_M_val
    S_B = Γ_B_val - R_rec_B_val
    S_Fn = trans.trans_n - trans.insulation  # trans.trans_n already has η in it!
    S_Fm = Γ_F_val - trans.trans_i - trans.trans_n - trans.Resp_F_conv
    S_Fi = trans.insulation + trans.trans_i - R_rec_F_val
    S_E = Γ_E_val - R_rec_E_val
    S_M = J_M_val

    # STEP 8: Total respiration
    Total_Resp = Resp_B_val + Resp_F_val + trans.Resp_F_conv

    # STEP 9: CARBON CONSERVATION CHECK
    # All carbon sources/sinks must balance with respiration
    Carbon_in_out = S_C + S_B + S_Fi + S_Fn + S_Fm + S_E + S_M
    Carbon_respired = -Total_Resp

    # With correct S_C formula (no erroneous factor), conservation should hold precisely
    @test Carbon_in_out ≈ Carbon_respired rtol=1e-10

    # Additional checks
    @test isfinite(Carbon_in_out)
    @test isfinite(Total_Resp)
    @test Total_Resp > 0.0  # Respiration should be positive

    # Verify individual balances
    # Bacterial balance: uptake = growth + EPS + respiration (NOT death - that's biomass loss)
    Bacterial_balance = R_B_val - (Γ_B_val + Γ_E_val + Resp_B_val)
    @test abs(Bacterial_balance) < 1e-10

    # Fungal balance: uptake = growth + respiration (NOT conversion cost - that's from transitions)
    # Resp_F_conv is cost of F_i ↔ F_m and F_n ↔ F_m transitions, not from new uptake
    Fungal_balance = R_F_val - (Γ_F_val + Resp_F_val)
    @test abs(Fungal_balance) < 1e-10
end
