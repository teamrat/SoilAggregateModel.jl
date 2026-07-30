# test_stiff.jl
# The production solver, and the duplication it is allowed to have.
#
# `run_aggregate_stiff` is where every number in the paper comes from, and until
# 2026-07-30 nothing in this suite called it.

using Test

import SoilAggregateModel: compute_total_carbon

@testset "Stiff solver (production path)" begin

    bio  = BiologicalProperties()
    soil = SoilProperties()
    T(t) = 293.15
    ψ(t) = -10.0
    O2(t) = 0.3

    # The mol_laplacian/crank_nicolson pin, the cross-solver agreement check and
    # the solver dispatch went with the split solver on 2026-07-30. There is no
    # duplication left to pin: mol.jl's Laplacian is the only one, and
    # test_mol.jl checks it against the two properties that matter — it sums to
    # zero against the conservation weights, and the flux boundary delivers
    # exactly 4πr₀²J. Those are analytic, which the pin was not.

    # ── End to end ───────────────────────────────────────────────────────────
    @testset "run_aggregate_stiff runs and returns a well-formed result" begin
        times  = Float64.(0:0.5:3)
        result = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 3.0);
                                     n_grid = 30, output_times = times)

        @test result isa SimulationResult
        @test length(result.outputs) == length(times)
        @test [rec.t for rec in result.outputs] ≈ times
        @test result.diagnostics["retcode"] == "Success"

        first_, last_ = result.outputs[1].state, result.outputs[end].state
        @test last_.P < first_.P                       # POM is consumed
        @test last_.CO2_cumulative > 0.0               # and respired
        for rec in result.outputs                      # recovered CO₂ never falls
            @test rec.state.CO2_cumulative >= -1e-12
        end
        @test issorted([rec.state.CO2_cumulative for rec in result.outputs])

        for f in (:C, :B, :F_n, :F_m, :F_i, :E, :M, :O)
            @test all(isfinite, getfield(last_, f))
            @test all(>=(-1e-10), getfield(last_, f))
        end

        # NaN, not 0.0, and this is load-bearing: the stiff path recovers CO₂ by
        # difference, so its balance is zero by construction. Reporting 0.0 would
        # read as a passing check on something never checked (REFERENCE §17a).
        @test all(isnan(rec.mass_balance_error) for rec in result.outputs)
    end

    # ── Defaults ─────────────────────────────────────────────────────────────
    @testset "defaults" begin
        @test N_GRID_DEFAULT == 200

        ts = default_output_times(0.0, 3.0)
        @test ts[1] == 0.0 && ts[end] == 3.0
        @test issorted(ts)
        @test all(diff(ts) .<= 1.0 + 1e-12)
        @test default_output_times(0.0, 100.0)[2] == 1.0   # capped at 1 day

        a = run_aggregate_stiff(bio, soil, T, ψ, O2, (0.0, 2.0); n_grid = 20)
        @test length(a.outputs) == length(default_output_times(0.0, 2.0))
    end
end
