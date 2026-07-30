# test_mol.jl
# Tests for the method-of-lines formulation (src/solver/mol.jl).
#
# The load-bearing test here is the Jacobian superset check. A hand-built
# sparsity pattern that is missing an entry does not error — the Newton solve
# just converges badly or to the wrong thing. So the pattern is checked against
# a Jacobian computed independently by finite differences, entry by entry.

using Test
using SparseArrays

import SoilAggregateModel: mol_sid, mol_iP, mol_neq, MOL_NSP,
                           MOL_C, MOL_B, MOL_FN, MOL_FM, MOL_O, MOL_FI, MOL_E, MOL_M,
                           mol_rhs!, mol_laplacian, mol_jacobian_prototype,
                           state_to_vector, vector_to_state!, mol_outer_oxygen!,
                           AggregateState, TemperatureCache,
                           update_temperature_cache!, create_initial_state,
                           GridInfo

@testset "Method of lines" begin

    # ---------------------------------------------------------------- indexing
    @testset "state vector round-trips" begin
        n = 12
        st = AggregateState(n)
        for (v, base) in ((st.C, 1.0), (st.B, 2.0), (st.F_n, 3.0), (st.F_m, 4.0),
                          (st.O, 5.0), (st.F_i, 6.0), (st.E, 7.0), (st.M, 8.0))
            v .= base .+ (1:n) ./ 100
        end
        st.P = 42.0; st.P_0 = 99.0; st.CO2_cumulative = 7.5

        u = state_to_vector(st)
        @test length(u) == mol_neq(n) == 8n + 1

        back = AggregateState(n)
        back.P_0 = 99.0
        back.CO2_cumulative = -1.0        # sentinel: must survive untouched
        vector_to_state!(back, u)

        @test back.C == st.C
        @test back.B == st.B
        @test back.F_n == st.F_n
        @test back.F_m == st.F_m
        @test back.O == st.O
        @test back.F_i == st.F_i
        @test back.E == st.E
        @test back.M == st.M
        @test back.P == st.P
        # Respired carbon is recovered from the carbon balance by the caller, so
        # the vector must not carry it and must not overwrite it.
        @test back.CO2_cumulative == -1.0
        @test back.P_0 == 99.0
    end

    @testset "indices are distinct and cover the vector" begin
        n = 7
        idx = [mol_sid(i, k) for i in 1:n for k in 1:MOL_NSP]
        push!(idx, mol_iP(n))
        @test sort(idx) == collect(1:mol_neq(n))
    end

    # ------------------------------------------------------- the spatial operator
    #
    # The conservation weight is W_i = 4π r_i² h, chosen so that W_i/(r_i²h²) is
    # constant and the sum over nodes telescopes. These two tests are the reason
    # the POM→DOC transfer is exact rather than merely consistent.

    @testset "zero-flux operator conserves mass exactly" begin
        n = 40
        grid = GridInfo(n, 0.5, 5.0)
        r, h = grid.r_grid, grid.h
        u = zeros(mol_neq(n))
        for i in 1:n                                # arbitrary, non-smooth profile
            u[mol_sid(i, MOL_C)] = 1.0 + sin(3.1 * i) + 0.5 * cos(0.7 * i)
        end
        D = [0.2 + 0.05 * sin(0.4 * i) for i in 1:n]   # spatially varying

        total = 0.0
        for i in 1:n
            W = 4.0 * π * r[i]^2 * h
            total += W * mol_laplacian(u, D, i, n, r, h, MOL_C, 0.0)
        end
        scale = sum(4.0 * π * r[i]^2 * h * abs(u[mol_sid(i, MOL_C)]) for i in 1:n)
        @test abs(total) < 1e-10 * scale
    end

    @testset "flux boundary delivers exactly 4πr₀²J" begin
        n = 40
        grid = GridInfo(n, 0.5, 5.0)
        r, h = grid.r_grid, grid.h
        u = zeros(mol_neq(n))
        for i in 1:n
            u[mol_sid(i, MOL_C)] = 2.0 + cos(1.3 * i)
        end
        D = [0.2 + 0.05 * sin(0.4 * i) for i in 1:n]
        J = 3.7

        total = 0.0
        for i in 1:n
            W = 4.0 * π * r[i]^2 * h
            total += W * mol_laplacian(u, D, i, n, r, h, MOL_C, J)
        end
        expected = 4.0 * π * r[1]^2 * J        # = R_P(J, r_0)

        # The zero-flux case conserves (previous testset), so the whole of the
        # weighted sum with J is the boundary delivery.
        @test total ≈ expected rtol=1e-12

        # It must also carry no dependence on D. D_1 cancels analytically in the
        # ghost-node substitution u_0 = u_1 + h·J/D_1; if it did not, the
        # transfer out of P and the transfer into C would differ by a
        # state-dependent factor and the coupling would be conservative only
        # approximately. Perturbing D_1 by 17x must move the delivered amount
        # not at all.
        D2 = copy(D); D2[1] *= 17.0
        delivered(Dv) = sum(4.0 * π * r[i]^2 * h *
                            (mol_laplacian(u, Dv, i, n, r, h, MOL_C, J) -
                             mol_laplacian(u, Dv, i, n, r, h, MOL_C, 0.0))
                            for i in 1:n)
        @test delivered(D)  ≈ expected rtol=1e-12
        @test delivered(D2) ≈ expected rtol=1e-12
    end

    # -------------------------------------------------------- Jacobian sparsity
    @testset "prototype is a superset of the true Jacobian" begin
        n    = 20
        r_0  = 0.5
        grid = GridInfo(n, r_0, 5.0)
        bio  = BiologicalProperties()
        soil = SoilProperties()
        ic   = InitialConditions(SOC = 0.02)

        st = create_initial_state(n, bio, soil, ic; P_0 = 10.0, ω = 1.0)
        u0 = state_to_vector(st)

        f_T = TemperatureCache()
        update_temperature_cache!(f_T, 298.15, bio, soil)
        O_amb = 0.2785 / f_T.K_H_O
        mol_outer_oxygen!(u0, n, O_amb)

        p = (n = n, r_grid = grid.r_grid, h = grid.h, bio = bio, soil = soil,
             T_func = t -> 298.15, ψ_func = t -> -29.0, O2_func = t -> 0.2785,
             f_T = f_T, P_0 = st.P_0, t_delay = 0.0)

        # Every pool strictly positive, so no max(0, u) guard is sitting on its
        # kink and no coupling is hidden at the evaluation point.
        u = [max(x, 1e-4) for x in u0]

        N  = mol_neq(n)
        f0 = zeros(N); mol_rhs!(f0, u, p, 0.0)
        @test all(isfinite, f0)

        proto = mol_jacobian_prototype(n)
        @test size(proto) == (N, N)

        f1 = zeros(N)
        missing_entries = Tuple{Int,Int}[]
        for j in 1:N
            δ  = 1e-7 * max(abs(u[j]), 1e-3)
            up = copy(u); up[j] += δ
            mol_rhs!(f1, up, p, 0.0)
            for i in 1:N
                dfi = (f1[i] - f0[i]) / δ
                # threshold well above finite-difference noise, well below any
                # real coupling
                if abs(dfi) > 1e-6 * max(1.0, abs(f0[i])) && proto[i, j] == 0
                    push!(missing_entries, (i, j))
                end
            end
        end
        @test isempty(missing_entries)
        if !isempty(missing_entries)
            @info "prototype missing entries" first_ten=first(missing_entries, 10) count=length(missing_entries)
        end
    end

    @testset "prototype structure" begin
        n = 50
        N = mol_neq(n)
        S = mol_jacobian_prototype(n)

        @test all(S[d, d] == 1.0 for d in 1:N)          # full diagonal, for W = I − γJ
        @test all(nonzeros(S) .== 1.0)                  # structural, not summed duplicates

        # No dense rows. This is the property the colour count depends on: a
        # single row touching every column would put every pair of columns in
        # conflict and force one RHS evaluation per column on every Jacobian.
        # One assertion, not N of them: the max over rows is the single bit of
        # information here (CLAUDE.md §6).
        @test maximum(count(!iszero, S[i, :]) for i in 1:N) <= 3 * MOL_NSP + 2

        # Two columns may share a colour only if they share no row, so this
        # checks the property directly: columns three nodes apart must not
        # conflict.
        @test all(iszero(S[:, mol_sid(5, k)] .* S[:, mol_sid(9, k)]) for k in 1:MOL_NSP)

        # Block-tridiagonal means nonzeros per row are bounded independently of
        # n — that is the property worth asserting, not a fraction of N², which
        # is trivially satisfied once n is large.
        @test nnz(S) / N < 30
        @test nnz(mol_jacobian_prototype(200)) / mol_neq(200) < 30
        @test nnz(S) > 8 * 8 * 3 * (n - 2)               # the blocks really are there

        # No coupling beyond the immediate neighbours.
        @test all(S[mol_sid(i, MOL_C), mol_sid(j, MOL_B)] == 0.0
                  for i in 1:n, j in 1:n if abs(i - j) > 1)
    end
end
