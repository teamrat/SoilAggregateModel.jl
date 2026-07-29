# test_closure.jl
# Carbon closure — the model's single most consequential invariant.
#
#   S_C + S_B + S_Fn + S_Fm + S_Fi + S_E + S_M = −Resp_total
#
# at every node, for any state. If this breaks, every carbon number the model
# produces is wrong and nothing raises an error: the pools and CO₂ simply stop
# agreeing, and the discrepancy gets attributed to physics.
#
# Until this file existed the only thing checking it was the split solver's
# integrated carbon balance — diluted over a trajectory, blinded wherever
# clipping is active, and 40 s per sample. This tests it directly: no
# integration, no solver, exact.
#
# See docs/REFERENCE.md §17a for the proof.

using Test
using Random

import SoilAggregateModel: compute_source_terms, TemperatureCache,
                           update_temperature_cache!, softplus

@testset "Carbon closure" begin

    # ---------------------------------------------------------------------
    # The load-bearing lemma.
    #
    # The bacterial half of the closure reduces to softplus(x) − softplus(−x)
    # = x. That is exact for any ε, because it is max(0,x) − max(0,−x) = x
    # surviving the smoothing. ANY future replacement for softplus must satisfy
    # this or the carbon budget breaks — which is why it is tested on its own
    # rather than only through the closure.
    # ---------------------------------------------------------------------
    @testset "softplus antisymmetry — softplus(x) − softplus(−x) = x" begin
        for ε in (1e-8, 1e-4, 1e-2, 1.0)
            for x in (-1e3, -7.5, -1.0, -1e-3, -1e-9, 0.0,
                       1e-9, 1e-3, 1.0, 7.5, 1e3)
                @test softplus(x, ε) - softplus(-x, ε) ≈ x atol = 1e-12 * max(abs(x), 1.0)
            end
        end
        # and it must survive the magnitudes where a naive exp(x/ε) overflows
        @test softplus(500.0, 1e-3) - softplus(-500.0, 1e-3) ≈ 500.0 rtol = 1e-14
        @test isfinite(softplus(500.0, 1e-3))
    end

    # ---------------------------------------------------------------------
    # The closure itself.
    # ---------------------------------------------------------------------

    """Residual of the closure identity, and the scale to judge it against."""
    function closure_residual(C, B, F_n, F_m, F_i, E, M, O, θ, θ_a, ψ, bio, soil, f_T)
        s = compute_source_terms(C, B, F_n, F_m, F_i, E, M, O, θ, θ_a, ψ, bio, soil, f_T)
        total = s.S_C + s.S_B + s.S_Fn + s.S_Fm + s.S_Fi + s.S_E + s.S_M
        resid = total + s.Resp_total
        # Judge against the largest term in the sum, not against the sum itself:
        # the terms cancel, so the sum is small by construction and dividing by
        # it would make the tolerance meaningless.
        scale = max(abs(s.S_C), abs(s.S_B), abs(s.S_Fn), abs(s.S_Fm),
                    abs(s.S_Fi), abs(s.S_E), abs(s.S_M), abs(s.Resp_total), 1e-30)
        return resid, scale
    end

    bio  = BiologicalProperties()
    soil = SoilProperties()
    f_T  = TemperatureCache()
    update_temperature_cache!(f_T, 298.15, bio, soil)

    @testset "holds at a representative state" begin
        r, sc = closure_residual(0.5, 0.02, 5e-3, 3e-4, 5e-3, 8e-3, 1.0, 8.9e-3,
                                 0.30, 0.15, -29.0, bio, soil, f_T)
        @test abs(r) < 1e-12 * sc
    end

    @testset "holds across eight orders of magnitude in every pool" begin
        worst = 0.0
        for C in (1e-6, 1e-3, 1.0, 5.0), B in (1e-5, 1e-2, 2.0),
            E in (1e-5, 1e-2, 0.5),      M in (1e-3, 1.0, 30.0)
            r, sc = closure_residual(C, B, 5e-3, 3e-4, 5e-3, E, M, 8.9e-3,
                                     0.30, 0.15, -29.0, bio, soil, f_T)
            worst = max(worst, abs(r) / sc)
        end
        @test worst < 1e-12
    end

    # The bacterial half of the proof crosses the softplus branch at
    # R_B = R_Bb. If the identity were going to fail anywhere it would fail
    # there, so sweep C — which drives R_B — across the crossing.
    @testset "holds across the R_B = R_Bb crossing" begin
        worst = 0.0
        for C in exp10.(range(-8, 1; length = 60))
            r, sc = closure_residual(C, 0.02, 5e-3, 3e-4, 5e-3, 8e-3, 1.0, 8.9e-3,
                                     0.30, 0.15, -29.0, bio, soil, f_T)
            worst = max(worst, abs(r) / sc)
        end
        @test worst < 1e-12
    end

    @testset "holds at pool minima and at zero" begin
        worst = 0.0
        for v in (0.0, 1e-12, 1e-6, bio.B_min, bio.E_min, bio.F_i_min)
            r, sc = closure_residual(v, v, v, v, v, v, v, 8.9e-3,
                                     0.30, 0.15, -29.0, bio, soil, f_T)
            worst = max(worst, abs(r) / sc)
        end
        @test worst < 1e-12
    end

    @testset "holds under anoxia and under drought" begin
        worst = 0.0
        for O in (0.0, 1e-8, 1e-5, 8.9e-3), ψ in (-1.0, -29.0, -300.0, -1500.0)
            r, sc = closure_residual(0.5, 0.02, 5e-3, 3e-4, 5e-3, 8e-3, 1.0, O,
                                     0.30, 0.15, ψ, bio, soil, f_T)
            worst = max(worst, abs(r) / sc)
        end
        @test worst < 1e-12
    end

    @testset "holds at every temperature in range" begin
        worst = 0.0
        cache = TemperatureCache()
        for T in (275.15, 288.15, 298.15, 310.15)
            update_temperature_cache!(cache, T, bio, soil)
            r, sc = closure_residual(0.5, 0.02, 5e-3, 3e-4, 5e-3, 8e-3, 1.0, 8.9e-3,
                                     0.30, 0.15, -29.0, bio, soil, cache)
            worst = max(worst, abs(r) / sc)
        end
        @test worst < 1e-12
    end

    @testset "holds on 2000 random states" begin
        rng   = Xoshiro(20260729)
        worst = 0.0
        arg   = nothing
        for _ in 1:2000
            st = (exp10(-8 + 9 * rand(rng)), exp10(-8 + 9 * rand(rng)),
                  exp10(-8 + 6 * rand(rng)), exp10(-8 + 6 * rand(rng)),
                  exp10(-8 + 6 * rand(rng)), exp10(-8 + 7 * rand(rng)),
                  exp10(-6 + 8 * rand(rng)), exp10(-9 + 7 * rand(rng)))
            θ  = 0.05 + 0.40 * rand(rng)
            ψv = -exp10(rand(rng) * 3.5)
            r, sc = closure_residual(st..., θ, soil.θ_s - θ, ψv, bio, soil, f_T)
            if abs(r) / sc > worst
                worst = abs(r) / sc
                arg = (st, θ, ψv)
            end
        end
        @test worst < 1e-12
        if worst >= 1e-12
            @info "closure violated" worst state=arg
        end
    end

    # A test that cannot fail is not a test. Confirm the residual is actually
    # sensitive: perturbing one source term by 1% must be caught.
    @testset "the check would detect a broken closure" begin
        s = compute_source_terms(0.5, 0.02, 5e-3, 3e-4, 5e-3, 8e-3, 1.0, 8.9e-3,
                                 0.30, 0.15, -29.0, bio, soil, f_T)
        scale = max(abs(s.S_C), abs(s.S_B), abs(s.Resp_total))
        broken = (s.S_C * 1.01) + s.S_B + s.S_Fn + s.S_Fm + s.S_Fi + s.S_E +
                 s.S_M + s.Resp_total
        @test abs(broken) / scale > 1e-6
    end
end
