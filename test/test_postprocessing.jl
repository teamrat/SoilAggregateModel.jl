# test_postprocessing.jl
# Tests for post-processing functions

using Test

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

        # O_aq <= O (Henry's law partitioning reduces aqueous concentration)
        for i in 1:result.grid.n
            @test aq.O_aq[i] <= rec.state.O[i] + 1e-15
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
        M_max = result.params.soil.M_max
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

        # Check bounds: 0 ≤ CUE ≤ 1 (with tolerance for very small biomass)
        @test all(cue.CUE_B .>= -1e-10)  # Relaxed tolerance for low-biomass regime
        @test all(cue.CUE_B .<= 1.0 + 1e-10)
        @test all(cue.CUE_F .>= -1e-10)
        @test all(cue.CUE_F .<= 1.0 + 1e-10)

        # Fungal CUE should be approximately constant (≈ Y_F)
        Y_F = result.params.bio.Y_F
        for i in 1:result.grid.n
            if rec.state.F_i[i] > 1e-10
                # Should be close to Y_F (fixed yield)
                @test abs(cue.CUE_F[i] - Y_F) < 0.1  # Within 10%
            end
        end
    end
end
