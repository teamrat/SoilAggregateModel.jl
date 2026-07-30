# mol.jl
# Method-of-lines formulation of the same model integrated by timestepper.jl.
#
# WHY THIS EXISTS
#
# The operator-split solver advances reactions with Forward Euler and sizes its
# step from |S·Δt|/u. That criterion is a stability guard for the explicit
# reaction, so the step is bounded by the FASTEST pool anywhere in the domain.
# Measured on De Gryze soil 3, DOC at the nodes adjacent to the POM surface
# reaches |S|/u ≈ 2000 /day and the required step falls below the 1e-4 floor,
# so a 45-day run costs ~3.9e5 steps and the cost per simulated day GROWS as
# the run proceeds. That does not scale to multi-year integrations at any
# constant factor.
#
# An implicit method sizes its step from accuracy rather than stability, so the
# step grows once the fast pools reach quasi-steady state and the cost tracks
# activity rather than elapsed time.
#
# WHAT IS DELIBERATELY IDENTICAL
#
# The spatial operator was written to reproduce crank_nicolson.jl term for term —
# node-centred spherical stencil, same ARITHMETIC face average of D, same ghost
# node treatment at the flux boundary. A better discretisation is available
# (finite-volume flux form with harmonic face averages) and is NOT used, because
# then a disagreement between the two integrators could not be attributed to the
# time integration. Change the space discretisation separately, and test it
# separately.
#
# Reactions come from compute_source_terms — the same function the split solver
# calls. There is one implementation of the biology.
#
# WHAT DIFFERS, AND IT IS NOT COSMETIC
#
# 1. θ and the effective diffusivities are evaluated at the CURRENT state here.
#    timestepper.jl computes them once per step from the state at the START of
#    the step and holds them fixed across both diffusion halves and the
#    reaction. The split solver therefore lags its own coefficients by one step.
#    The two agree as Δt → 0.
#
# 2. No negativity clipping. reaction_step.jl clips each pool at zero and
#    credits the clipped carbon to CO₂ so the budget still closes. That exists
#    to catch Forward Euler overshoot; an implicit method does not overshoot the
#    same way. Any CO₂ the split solver reports that came from clipping is a
#    numerical artefact, and the difference between the two runs measures it.
#
# 3. Cumulative respired carbon is not integrated here. It is recovered from
#    the carbon balance at output times: initial total carbon minus current
#    total carbon. It is a whole-domain quantity needed at a few dozen times per
#    run, so integrating it node-by-node at every step would be machinery built
#    for a resolution nothing uses.
#
# See docs/REFERENCE.md §20a.

# State layout: node-major, species fastest. This ordering makes the Jacobian
# block-tridiagonal with 8×8 blocks; species-major would give five tridiagonal
# blocks glued by dense reaction coupling, which fills in badly under LU.
const MOL_NSP = 8
const MOL_C, MOL_B, MOL_FN, MOL_FM, MOL_O, MOL_FI, MOL_E, MOL_M = 1, 2, 3, 4, 5, 6, 7, 8

@inline mol_sid(i::Int, k::Int) = (i - 1) * MOL_NSP + k
@inline mol_iP(n::Int)   = MOL_NSP * n + 1
@inline mol_neq(n::Int)  = MOL_NSP * n + 1

"""
    state_to_vector(state::AggregateState) -> Vector{Float64}

Flatten an `AggregateState` into the MOL state vector. Inverse of
[`vector_to_state!`](@ref).
"""
function state_to_vector(state::AggregateState)
    n = length(state.C)
    u = Vector{Float64}(undef, mol_neq(n))
    @inbounds for i in 1:n
        u[mol_sid(i, MOL_C)]  = state.C[i]
        u[mol_sid(i, MOL_B)]  = state.B[i]
        u[mol_sid(i, MOL_FN)] = state.F_n[i]
        u[mol_sid(i, MOL_FM)] = state.F_m[i]
        u[mol_sid(i, MOL_O)]  = state.O[i]
        u[mol_sid(i, MOL_FI)] = state.F_i[i]
        u[mol_sid(i, MOL_E)]  = state.E[i]
        u[mol_sid(i, MOL_M)]  = state.M[i]
    end
    u[mol_iP(n)] = state.P
    return u
end

"""
    vector_to_state!(state::AggregateState, u::AbstractVector) -> AggregateState

Unpack a MOL state vector back into `state`, in place. `P_0` and
`CO2_cumulative` are not carried in the vector and are left untouched; the
caller sets `CO2_cumulative` from the carbon balance.
"""
function vector_to_state!(state::AggregateState, u::AbstractVector)
    n = length(state.C)
    @inbounds for i in 1:n
        state.C[i]   = u[mol_sid(i, MOL_C)]
        state.B[i]   = u[mol_sid(i, MOL_B)]
        state.F_n[i] = u[mol_sid(i, MOL_FN)]
        state.F_m[i] = u[mol_sid(i, MOL_FM)]
        state.O[i]   = u[mol_sid(i, MOL_O)]
        state.F_i[i] = u[mol_sid(i, MOL_FI)]
        state.E[i]   = u[mol_sid(i, MOL_E)]
        state.M[i]   = u[mol_sid(i, MOL_M)]
    end
    state.P = u[mol_iP(n)]
    return state
end

"""
    mol_laplacian(u, D, i, n, r_grid, h, k, flux_inner) -> value

Spherical diffusion operator at node `i` for the species in column `k`.

The spatial operator, matched to the conservation weights of §18:

    L_i = (1/r_i²h²)·[r_{i+½}² D_{i+½} (u_{i+1} − u_i)
                      − r_{i−½}² D_{i−½} (u_i − u_{i−1})]

with `D_{i±½}` the arithmetic mean of the adjacent nodal values.

Inner boundary uses the ghost node `u_0 = u_1 + h·J/D_1` with the face placed
AT `r_grid[1]` — the domain starts at the POM surface and does not extend below
it. Passing `flux_inner = 0` recovers the zero-flux condition. Outer boundary is
zero flux; Dirichlet species are handled by the caller, which pins the node.
"""
@inline function mol_laplacian(u, D, i::Int, n::Int, r_grid, h::Float64,
                               k::Int, flux_inner)
    @inbounds begin
    r_i = r_grid[i]
    inv_r_sq_h_sq = 1.0 / (r_i * r_i * h * h)

    if i == 1
        r_hp = (r_grid[1] + r_grid[2]) / 2.0
        D_hp = (D[1] + D[2]) / 2.0
        interior = r_hp * r_hp * D_hp * (u[mol_sid(2, k)] - u[mol_sid(1, k)])
        # ghost_diff = u_1 − u_0 = −h·J/D_1 ; the flux term is
        # −r_0²·D_1·ghost_diff = +r_0²·h·J, so D_1 cancels exactly.
        boundary = r_grid[1] * r_grid[1] * h * flux_inner
        return inv_r_sq_h_sq * (interior + boundary)
    elseif i == n
        r_hm = (r_grid[n-1] + r_grid[n]) / 2.0
        D_hm = (D[n-1] + D[n]) / 2.0
        return -inv_r_sq_h_sq * r_hm * r_hm * D_hm *
               (u[mol_sid(n, k)] - u[mol_sid(n-1, k)])
    else
        r_hp = (r_grid[i] + r_grid[i+1]) / 2.0
        r_hm = (r_grid[i-1] + r_grid[i]) / 2.0
        D_hp = (D[i] + D[i+1]) / 2.0
        D_hm = (D[i-1] + D[i]) / 2.0
        return inv_r_sq_h_sq *
               (r_hp * r_hp * D_hp * (u[mol_sid(i+1, k)] - u[mol_sid(i, k)]) -
                r_hm * r_hm * D_hm * (u[mol_sid(i, k)] - u[mol_sid(i-1, k)]))
    end
    end  # @inbounds
end

"""
    mol_rhs!(du, u, p, t)

Right-hand side of the full 8n+1 system.

`p` is a NamedTuple carrying `n`, `r_grid`, `h`, `bio`, `soil`, the three
environmental driver functions, `f_T` (a `TemperatureCache`), `P_0` and
`t_delay`.

The temperature cache is a function of `t` only, never of `u`, so it is reused
across calls as `Float64`. θ and the diffusivities DO depend on `u` and are
allocated per call at `eltype(u)` so forward-mode AD sees them; that is 7
vectors of length `n` per evaluation, and is the one knowingly unoptimised part
of this file.
"""
function mol_rhs!(du, u, p, t)
    n      = p.n
    bio    = p.bio
    soil   = p.soil
    r_grid = p.r_grid
    h      = p.h
    Tc     = p.f_T

    T      = p.T_func(t)
    ψ      = p.ψ_func(t)
    O2_gas = p.O2_func(t)

    update_temperature_cache!(Tc, T, bio, soil)
    O_amb = O2_gas / Tc.K_H_O

    θ   = similar(u, n); θ_a = similar(u, n)
    D_C = similar(u, n); D_B = similar(u, n); D_Fn = similar(u, n)
    D_Fm = similar(u, n); D_O = similar(u, n)

    # --- water content and effective diffusivities, from the current state ---
    @inbounds for i in 1:n
        E_i  = max(u[mol_sid(i, MOL_E)],  0.0)
        Fi_i = max(u[mol_sid(i, MOL_FI)], 0.0)
        Fn_i = max(u[mol_sid(i, MOL_FN)], 0.0)

        α_eff  = alpha_effective(E_i, Fi_i, soil)
        θ[i]   = van_genuchten(ψ, α_eff, soil.n_vg, soil.θ_r, soil.θ_s)
        θ_a[i] = soil.θ_s - θ[i]

        D_C[i]  = D_eff_DOC(Tc.D_DOC_w, θ[i], soil.θ_s, soil.ρ_b, soil.k_d_eq)
        D_B[i]  = D_eff_bacteria(D_C[i], soil.D_B_rel)
        D_Fn[i] = D_eff_fungi_noninsulated(bio.D_Fn0, Tc.f_fun, θ[i], soil.θ_s)
        D_Fm[i] = D_eff_fungi_mobile(bio.D_Fm0, Tc.f_fun, Fn_i, Fi_i, bio.K_Fm_net)
        D_O[i]  = D_eff_oxygen(Tc.D_O2_w, Tc.D_O2_a, Tc.K_H_O, θ[i], θ_a[i], soil.θ_s)
    end

    # --- POM dissolution: one flux, used with opposite signs so the transfer
    #     from P into C at node 1 is exact rather than merely consistent ---
    P_val    = u[mol_iP(n)]
    R_P_max_T = bio.R_P_max * Tc.f_pom * pom_delay_factor(t, p.t_delay)
    J_P_val = J_P(P_val, p.P_0,
                  max(u[mol_sid(1, MOL_B)], 0.0),
                  max(u[mol_sid(1, MOL_FN)], 0.0),
                  θ[1], max(u[mol_sid(1, MOL_O)], 0.0),
                  R_P_max_T, bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
    R_P_val = R_P(J_P_val, r_grid[1])

    # --- reactions + diffusion, node by node ---
    @inbounds for i in 1:n
        C_i  = max(u[mol_sid(i, MOL_C)],  0.0)
        B_i  = max(u[mol_sid(i, MOL_B)],  0.0)
        Fn_i = max(u[mol_sid(i, MOL_FN)], 0.0)
        Fm_i = max(u[mol_sid(i, MOL_FM)], 0.0)
        O_i  = max(u[mol_sid(i, MOL_O)],  0.0)
        Fi_i = max(u[mol_sid(i, MOL_FI)], 0.0)
        E_i  = max(u[mol_sid(i, MOL_E)],  0.0)
        M_i  = max(u[mol_sid(i, MOL_M)],  0.0)

        s = compute_source_terms(C_i, B_i, Fn_i, Fm_i, Fi_i, E_i, M_i, O_i,
                                 θ[i], θ_a[i], ψ, bio, soil, Tc)

        du[mol_sid(i, MOL_C)]  = s.S_C  + mol_laplacian(u, D_C,  i, n, r_grid, h, MOL_C,  J_P_val)
        du[mol_sid(i, MOL_B)]  = s.S_B  + mol_laplacian(u, D_B,  i, n, r_grid, h, MOL_B,  0.0)
        du[mol_sid(i, MOL_FN)] = s.S_Fn + mol_laplacian(u, D_Fn, i, n, r_grid, h, MOL_FN, 0.0)
        du[mol_sid(i, MOL_FM)] = s.S_Fm + mol_laplacian(u, D_Fm, i, n, r_grid, h, MOL_FM, 0.0)
        du[mol_sid(i, MOL_FI)] = s.S_Fi
        du[mol_sid(i, MOL_E)]  = s.S_E
        du[mol_sid(i, MOL_M)]  = s.S_M

        # Oxygen: Dirichlet at the outer node. The split solver overwrites it
        # with O_amb on every diffusion half-step, so pinning it here is the
        # same net condition, not an approximation of it. Its respiration still
        # respiration at that node still contributes to the carbon budget.
        if i == n
            du[mol_sid(n, MOL_O)] = 0.0
        else
            du[mol_sid(i, MOL_O)] = s.S_O + mol_laplacian(u, D_O, i, n, r_grid, h, MOL_O, 0.0)
        end
    end

    du[mol_iP(n)] = -R_P_val
    return nothing
end

"""
    mol_jacobian_prototype(n::Int) -> SparseMatrixCSC{Float64,Int}

Sparsity pattern of the MOL Jacobian, built from the structure of `mol_rhs!`
rather than detected by tracing.

# Why not automatic detection

`SparseConnectivityTracer` traces control flow, and `softplus`
(`math_utils.jl:30`) branches on `x > 0`. `GradientTracer` carries no primal
value, so that branch cannot be resolved and detection fails outright. The
local-tracer variant does carry a primal and would succeed, but it returns the
pattern *at that state* — and every rate law here is guarded with `max(0, u)`,
so a pool sitting at zero short-circuits its guard and the coupling never
appears. A pattern missing entries does not fail loudly: it degrades the Newton
solve silently.

The structure, by contrast, is fixed and known. This function returns a
**superset**: extra structural zeros cost a little arithmetic, missing entries
would cost correctness. `test/test_mol.jl` checks the superset property against
a Jacobian computed independently by finite differences.

# The structure

- **Full 8×8 blocks at `i−1, i, i+1`.** Reactions at node `i` read all eight
  species there (θ depends on `E` and `F_i`). The diffusion stencil reads the
  species itself at `i±1`, and the face-averaged `D` there depends on `E`,
  `F_i` and `F_n` at `i±1`. Blocking the whole 8×8 is the superset.
- **`P` into node 1**, because `J_P` enters the node-1 flux boundary term.
- **`P` row from node 1**, because `J_P` reads `B`, `F_n`, `O` and θ there.
- **Full diagonal**, so `W = I − γJ` is representable in this pattern.

The pattern needs roughly 26 colours under the column colouring that sparse
forward-mode AD uses to compress the Jacobian: columns conflict only between
nodes within distance 2, giving three node-groups × eight species. That bound is
what makes Jacobian construction cheap, and it holds only because no row of this
system is dense — a single dense row would put every pair of columns in conflict
and force one RHS evaluation per column.
"""
function mol_jacobian_prototype(n::Int)
    N  = mol_neq(n)
    iP = mol_iP(n)

    rows = Int[]
    cols = Int[]
    sizehint!(rows, 9 * MOL_NSP * MOL_NSP * n)

    @inline function block!(inode::Int, jnode::Int)
        for a in 1:MOL_NSP, b in 1:MOL_NSP
            push!(rows, mol_sid(inode, a))
            push!(cols, mol_sid(jnode, b))
        end
    end

    for i in 1:n
        for j in max(1, i - 1):min(n, i + 1)
            block!(i, j)
        end
    end

    for a in 1:MOL_NSP                 # P -> node 1  (only C, taken as the block)
        push!(rows, mol_sid(1, a)); push!(cols, iP)
        push!(rows, iP);            push!(cols, mol_sid(1, a))
    end
    push!(rows, iP); push!(cols, iP)

    for d in 1:N
        push!(rows, d); push!(cols, d)
    end

    S = sparse(rows, cols, ones(Float64, length(rows)), N, N)
    fill!(nonzeros(S), 1.0)            # duplicates summed above; a prototype is structural
    return S
end

"""
    mol_outer_oxygen!(u, n, O_amb)

Pin the outer oxygen node to the ambient aqueous concentration. `mol_rhs!` holds
`du = 0` there, so the value must be correct in the initial condition.
"""
function mol_outer_oxygen!(u::AbstractVector, n::Int, O_amb::Real)
    u[mol_sid(n, MOL_O)] = O_amb
    return u
end
