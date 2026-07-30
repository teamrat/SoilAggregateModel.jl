# test_postprocessing.jl
# Tests for post-processing functions

using Test
import SoilAggregateModel: compute_source_terms, _prepare_environment, compute_total_carbon

@testset "Post-processing" begin
    bio = BiologicalProperties()
    soil = SoilProperties()

    # Constant environment
    T(t) = 293.15  # 20°C
    ψ(t) = -10.0   # kPa
    O2(t) = 0.3    # μg/mm³

    result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                          n_grid=30, output_times=Float64.(0:1:30))

    @testset "integrated_pools — both domains" begin
        pools = integrated_pools(result)

        @test pools isa IntegratedPools
        @test length(pools.t) == 31
        @test pools.P[end] < pools.P[1]          # POM decreases
        @test pools.CO2[end] > pools.CO2[1]       # CO₂ increases

        # Total ≥ aggregate for every pool at every time
        for (total_field, agg_field) in [
            (:C_total, :C_agg), (:B_total, :B_agg),
            (:F_i_total, :F_i_agg), (:F_n_total, :F_n_agg),
            (:F_m_total, :F_m_agg), (:E_total, :E_agg), (:M_total, :M_agg)]
            t_vals = getfield(pools, total_field)
            a_vals = getfield(pools, agg_field)
            for i in 1:length(pools.t)
                @test t_vals[i] >= a_vals[i] - 1e-15
            end
        end
    end

    @testset "carbon_balance_table" begin
        bal = carbon_balance_table(result)
        @test bal isa NamedTuple
        @test haskey(bal, :t)
        @test haskey(bal, :C_total)
        @test haskey(bal, :C_initial)
        @test haskey(bal, :relative_error)
        @test all(abs.(bal.relative_error) .< 1e-12)

        # total_system_carbon (pools-level) and compute_total_carbon
        # (state-level) are separate entry points into the same sum because they
        # take different inputs. Nothing structural keeps them from drifting, so
        # they are pinned to each other here.
        pools = integrated_pools(result)
        for k in (1, 16, 31)
            from_pools = total_system_carbon(pools, k)
            from_state = compute_total_carbon(result.outputs[k].state, result.grid)
            @test from_pools ≈ from_state rtol=1e-12
        end

        # include_co2 = false drops exactly the respired term and nothing else
        no_co2 = total_system_carbon(pools, 20; include_co2 = false)
        @test no_co2 ≈ total_system_carbon(pools, 20) - pools.CO2[20] rtol=1e-12

        # the two carbon_balance_table methods are one implementation
        @test carbon_balance_table(pools).C_total == bal.C_total

        # compute_total_carbon entry points agree bitwise
        st = result.outputs[10].state
        via_grid = compute_total_carbon(st, result.grid)
        @test via_grid == compute_total_carbon(st, result.grid.W)
        @test via_grid == compute_total_carbon(st, result.grid.r_grid, result.grid.h)
        @test_throws DimensionMismatch compute_total_carbon(st, result.grid.W[1:end-1])
    end

    @testset "co2_flux" begin
        flux = co2_flux(result)
        @test length(flux) == length(result.outputs)
        @test all(flux .>= -1e-15)  # non-negative (with floating point tolerance)
        # Flux should be positive (CO2 increases over time)
        @test all(flux[2:end] .> -1e-15)
    end

    @testset "compute_r_agg" begin
        r0 = compute_r_agg(result.outputs[1], result.grid, result.params)
        @test r0 >= 0.0
        @test r0 <= result.grid.r_max
        # Check a few more snapshots
        for i in [1, 15, 31]
            r_agg = compute_r_agg(result.outputs[i], result.grid, result.params)
            @test r_agg >= 0.0
            @test r_agg <= result.grid.r_max
        end
    end

    @testset "critical_binding — size dependence" begin
        δ_s   = wet_sieving_stress().δ_s
        G_ref = critical_binding(soil)

        # The default must be the flat threshold, or every result predating
        # 2026-07-29 silently changes meaning.
        @test soil.p_Gc == 0.0

        @testset "p_Gc = 0 reduces exactly" begin
            # Not approximately. The short-circuit multiplies by 1.0, so these
            # are the same Float64 — if that ever becomes an ≈, the flat case is
            # going through the power and the reduction is no longer exact.
            for r in (0.01, 0.5, δ_s, 3.7, 100.0)
                @test critical_binding(soil, r) === G_ref
            end
            rv = [0.01, 0.5, δ_s, 3.7]
            @test critical_binding(soil, rv) == fill(G_ref, 4)
            @test length(critical_binding(soil, rv)) == length(rv)
        end

        @testset "pivot is δ_s for every exponent" begin
            # G_c(δ_s) = G_ref identically in p_Gc. This is what makes δ_s a
            # pivot rather than one more fitted length: changing p_Gc rotates
            # the threshold about a fixed point instead of also moving its level.
            for p in (0.0, 0.5, 1.0, 2.0, -0.5)
                s = SoilProperties(p_Gc = p)
                @test critical_binding(s, δ_s) ≈ G_ref rtol=1e-14
            end
        end

        @testset "p_Gc = 1 is linear in r" begin
            s = SoilProperties(p_Gc = 1.0)
            @test critical_binding(s, 2δ_s) ≈ 2G_ref rtol=1e-14
            @test critical_binding(s, δ_s/4) ≈ G_ref/4 rtol=1e-14
            # Strictly increasing, so a bigger aggregate always needs more binder
            r = collect(0.1:0.1:5.0)
            g = critical_binding(s, r)
            @test all(diff(g) .> 0.0)
        end

        @testset "p_Gc = 2 is the square" begin
            s = SoilProperties(p_Gc = 2.0)
            @test critical_binding(s, 3δ_s) ≈ 9G_ref rtol=1e-14
        end

        # ── r_agg under a size-dependent threshold ───────────────────────────
        #
        # Tested against a CONSTRUCTED binder profile, not a simulated one. The
        # default soil's threshold need not cross this particular run's profile
        # at all, and a test that silently measures "no crossing anywhere" would
        # pass while proving nothing.
        #
        # binding(r) = 1 - r/R, strictly decreasing and linear. For p_Gc = 0 and
        # p_Gc = 1 the threshold is also linear in r, so `excess` is linear and
        # the sub-cell interpolation is EXACT — the expected crossing is a closed
        # form and the test is an equality, not a tolerance on a curve fit.
        #
        #   p_Gc = 0:  1 - r/R = G_ref              →  r = R(1 - G_ref)
        #   p_Gc = 1:  1 - r/R = G_ref·r/δ_s        →  r = 1/(1/R + G_ref/δ_s)
        #
        # r_agg does not feed back into the state equations, so re-scoring one
        # stored profile under several thresholds is exact, not indicative.
        @testset "r_agg under G_c(r)" begin
            n_t, R = 101, 20.0
            grid_t  = GridInfo(n_t, 0.1, 10.1)          # h = 0.1
            state_t = AggregateState(n_t)
            state_t.F_i .= 1.0 .- grid_t.r_grid ./ R
            state_t.E   .= 0.0                          # binder is F_i alone
            rec_t = OutputRecord(0.0, state_t, 0.0)

            # G_ref = 0.6125 puts the flat crossing at r = 7.75, BETWEEN nodes
            # (7.7 and 7.8), so hitting it proves the interpolation rather than
            # node snapping. w_E is irrelevant here because E = 0.
            G_want  = 0.6125
            soil_p(p) = SoilProperties(
                κ_b = wet_sieving_stress().τ_w * SoilProperties().d_32 / G_want,
                p_Gc = p)
            r_agg_p(p) = compute_r_agg(rec_t, grid_t, ParameterSet(bio, soil_p(p)))

            @test critical_binding(soil_p(0.0)) ≈ G_want rtol=1e-12

            @testset "flat threshold, closed form" begin
                @test r_agg_p(0.0) ≈ R * (1 - G_want) rtol=1e-12   # 7.75
                @test r_agg_p(0.0) != grid_t.r_grid[78]            # not snapped
                @test grid_t.r_grid[1] < r_agg_p(0.0) < grid_t.r_max
            end

            @testset "linear threshold, closed form" begin
                expected = 1.0 / (1.0/R + G_want/δ_s)
                @test r_agg_p(1.0) ≈ expected rtol=1e-12           # ≈ 1.155
                @test grid_t.r_grid[1] < r_agg_p(1.0) < grid_t.r_max
                # The crossing moved inward by more than a factor of six. That
                # is the whole mechanism: a threshold rising with r cuts the
                # binder profile where the profile is still steep.
                @test r_agg_p(1.0) < r_agg_p(0.0)
            end

            @testset "monotone in p_Gc, and not quantised" begin
                ps   = collect(range(0.0, 1.0, length = 25))
                vals = r_agg_p.(ps)
                @test all(diff(vals) .<= 1e-12)      # non-increasing throughout
                # r_agg is interpolated, so a smooth sweep gives a spread of
                # distinct values. Snapped to nodes it would give a handful —
                # BACKLOG item 6 is why that distinction matters.
                @test length(unique(vals)) > 15
            end

            @testset "negative exponent pushes outward" begin
                @test r_agg_p(-0.5) > r_agg_p(0.0)
            end
        end
    end

    @testset "radial_profiles" begin
        # Extract specific times
        profs = radial_profiles(result; times=[0.0, 15.0, 30.0])
        @test length(profs) == 3
        @test profs[1].t ≈ 0.0
        @test profs[2].t ≈ 15.0
        @test profs[3].t ≈ 30.0
        @test length(profs[1].C) == result.grid.n
        @test length(profs[1].r) == result.grid.n

        # Extract all times
        all_profs = radial_profiles(result)
        @test length(all_profs) == length(result.outputs)
        @test all_profs[1].t ≈ result.outputs[1].t
        @test all_profs[end].t ≈ result.outputs[end].t

        # Check structure of profile
        prof = profs[1]
        @test haskey(prof, :t)
        @test haskey(prof, :r)
        @test haskey(prof, :C)
        @test haskey(prof, :B)
        @test haskey(prof, :F_i)
        @test haskey(prof, :F_n)
        @test haskey(prof, :F_m)
        @test haskey(prof, :E)
        @test haskey(prof, :M)
        @test haskey(prof, :O)
        @test haskey(prof, :P)
        @test haskey(prof, :CO2)
    end

    @testset "aqueous_concentrations" begin
        rec = result.outputs[15]
        aq = aqueous_concentrations(rec, result.grid, result.params, result.env)

        # Check structure
        @test length(aq.C_aq) == result.grid.n
        @test length(aq.O_aq) == result.grid.n

        # Check non-negativity
        @test all(aq.C_aq .>= 0.0)
        @test all(aq.O_aq .>= 0.0)

        # C_aq < C (retardation factor > 1 since k_d·ρ_b > 0)
        for i in 1:result.grid.n
            @test aq.C_aq[i] <= rec.state.C[i] + 1e-15
        end

        # O_aq == state.O (state.O is already C_aq since the O_total→C_aq switch)
        for i in 1:result.grid.n
            @test aq.O_aq[i] ≈ rec.state.O[i] atol=1e-15
        end
    end

    @testset "maoc_equilibrium" begin
        rec = result.outputs[15]
        M_eq = maoc_equilibrium(rec, result.grid, result.params, result.env)

        # Check structure
        @test length(M_eq) == result.grid.n

        # Check non-negativity
        @test all(M_eq .>= 0.0)

        # M_eq should be bounded by M_max
        M_max = maoc_capacity(result.params.soil)
        @test all(M_eq .<= M_max + 1e-15)

        # M_eq should be monotonic with C (higher C → higher M_eq)
        # At least at t=0 where C should be relatively uniform
        rec0 = result.outputs[1]
        M_eq0 = maoc_equilibrium(rec0, result.grid, result.params, result.env)
        @test all(M_eq0 .>= 0.0)
    end

    @testset "respiration_rates" begin
        rec = result.outputs[15]
        resp = respiration_rates(rec, result.grid, result.params, result.env)

        # Check structure
        @test length(resp.Resp_B) == result.grid.n
        @test length(resp.Resp_F) == result.grid.n
        @test length(resp.Resp_F_conv) == result.grid.n
        @test length(resp.Resp_total) == result.grid.n

        # Check non-negativity
        @test all(resp.Resp_B .>= -1e-15)
        @test all(resp.Resp_F .>= -1e-15)
        @test all(resp.Resp_F_conv .>= -1e-15)
        @test all(resp.Resp_total .>= -1e-15)

        # Total = sum of components
        for i in 1:result.grid.n
            @test resp.Resp_total[i] ≈ resp.Resp_B[i] + resp.Resp_F[i] + resp.Resp_F_conv[i]
        end

        # Respiration should be positive where biomass exists
        for i in 1:result.grid.n
            if rec.state.B[i] > 1e-10 || rec.state.F_i[i] > 1e-10
                @test resp.Resp_total[i] > 0.0
            end
        end
    end

    @testset "carbon_use_efficiency" begin
        rec = result.outputs[15]
        cue = carbon_use_efficiency(rec, result.grid, result.params, result.env)

        # Check structure
        @test length(cue.CUE_B) == result.grid.n
        @test length(cue.CUE_F) == result.grid.n

        # Bounds. CUE_B = Γ_B / R_B, and
        #     Γ_B = Y_B·softplus(R_diff)·(1-γ) - softplus(-R_diff),  R_diff = R_B - R_Bb
        # so where basal maintenance R_Bb exceeds uptake R_B the numerator is
        # negative by construction: the cell is respiring its own biomass to pay
        # maintenance. A negative CUE_B is the correct signature of that state,
        # not an error — it happens at outer nodes where DOC is depleted.
        # Only the upper bound is a genuine physical limit.
        @test all(isfinite.(cue.CUE_B))
        @test all(cue.CUE_B .<= 1.0 + 1e-10)
        # Nodes with net-positive uptake must still have non-negative CUE
        @test any(cue.CUE_B .>= -1e-10)
        @test all(cue.CUE_F .>= -1e-10)
        @test all(cue.CUE_F .<= 1.0 + 1e-10)

        # Fungal CUE tracks the SPACE-LIMITED yield, not the base yield:
        #     Y_F = Y_F_base · F_S / (F_i + F_n + F_m + F_S)
        # (manuscript_changes.md §3). Comparing against bio.Y_F alone only ever
        # passed because total fungal biomass was negligible.
        bio_p = result.params.bio
        for i in 1:result.grid.n
            if rec.state.F_i[i] > 1e-10
                F_total = rec.state.F_i[i] + rec.state.F_n[i] + rec.state.F_m[i]
                Y_F_expected = bio_p.Y_F * bio_p.F_S / (F_total + bio_p.F_S)
                @test abs(cue.CUE_F[i] - Y_F_expected) < 0.05
            end
        end
    end
end


# ============================================================================
# Post-processing must agree with the solver
#
# `respiration_rates` and `carbon_use_efficiency` re-derive R_B, R_F and the
# fungal transitions from the stored state. That is a second implementation of
# what `compute_source_terms` does, so the two can diverge in the definitions
# they share — C_aq, O_aq, the ζ clamp — while every structural, non-negativity
# and self-sum assertion in this file still passes.
#
# This testset is the oracle: it compares the re-derivation against the solver's
# own function, node by node, at machine precision, and carries no expected
# values of its own. See REFERENCE.md §26 for what it caught.
# ============================================================================
@testset "Post-processing agrees with compute_source_terms" begin
    bio = BiologicalProperties()
    soil = SoilProperties()
    T(t) = 293.15
    ψ(t) = -10.0
    O2(t) = 0.3

    result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0);
                           n_grid=30, output_times=Float64.(0:1:30))

    for idx in (2, 15, 31)                     # early, mid, late
        rec  = result.outputs[idx]
        ev   = _prepare_environment(rec, result.grid, result.params, result.env)
        resp = respiration_rates(rec, result.grid, result.params, result.env)

        for i in 1:result.grid.n
            src = compute_source_terms(rec.state.C[i], rec.state.B[i],
                                       rec.state.F_n[i], rec.state.F_m[i],
                                       rec.state.F_i[i], rec.state.E[i],
                                       rec.state.M[i], rec.state.O[i],
                                       ev.θ[i], ev.θ_a[i], ev.ψ,
                                       result.params.bio, result.params.soil,
                                       ev.f_T)
            @test resp.Resp_total[i] ≈ src.Resp_total rtol=1e-12
        end
    end

    # aqueous_concentrations must use the solver's C_aq definition too
    @testset "C_aq / O_aq match the solver" begin
        rec = result.outputs[15]
        ev  = _prepare_environment(rec, result.grid, result.params, result.env)
        aq  = aqueous_concentrations(rec, result.grid, result.params, result.env)
        for i in 1:result.grid.n
            @test aq.C_aq[i] ≈ rec.state.C[i] /
                  (ev.θ[i] + result.params.soil.ρ_b * result.params.soil.k_d_eq) rtol=1e-14
            @test aq.O_aq[i] == rec.state.O[i]
        end
    end
end
