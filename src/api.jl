# api.jl
# User-facing API for running aggregate simulations.
#
# `run_aggregate` here is the Strang-split REFERENCE implementation.
# `run_aggregate_stiff` (solver/mol_solve.jl) is the default workhorse.
# See the header of solver/timestepper.jl for why both are kept.

# `run_aggregate` — the Strang-splitting reference implementation — was archived
# on 2026-07-30 together with the solver it drove, and the
# `run_aggregate(solver::Symbol, ...)` dispatcher with it: with one integrator
# there is nothing to dispatch on. `_archive/split_solver_20260730/README.md`
# gives the reasoning. `run_aggregate_stiff` is the only entry point.

# `create_initial_state_legacy` was deleted 2026-07-30. It had no caller in
# `src/`, `test/` or `paper/`, and it seeded a fixed background DOC of 1e-4
# rather than partitioning measured SOC. `create_initial_state` is the supported
# entry point and is now exported.


"""
    compute_total_carbon(state::AggregateState, W::AbstractVector{<:Real})
    compute_total_carbon(state::AggregateState, grid::GridInfo)
    compute_total_carbon(state::AggregateState, r_grid::AbstractVector{<:Real}, h::Real)

Total carbon in the system [µg-C]:

    C_total = P + ∑_i (C+B+F_n+F_m+F_i+E+M)_i · W_i + CO₂

`W` are the conservation weights (`conservation_weight`). Pass a `GridInfo` when
one is available: `grid.W` is already built, so nothing is allocated.

This is the state-level total. `total_system_carbon` computes the same quantity
from an already-integrated `IntegratedPools`. The two are separate entry points
into the same sum because they take different inputs, and they must agree.
"""
function compute_total_carbon(state::AggregateState, W::AbstractVector{<:Real})
    length(W) == length(state.C) ||
        throw(DimensionMismatch("conservation weights have length $(length(W)), state has $(length(state.C)) nodes"))
    integral = 0.0
    @inbounds for i in eachindex(state.C)
        C_pools = state.C[i] + state.B[i] + state.F_n[i] + state.F_m[i] +
                 state.F_i[i] + state.E[i] + state.M[i]
        integral += C_pools * W[i]
    end
    return state.P + integral + state.CO2_cumulative
end

compute_total_carbon(state::AggregateState, grid::GridInfo) =
    compute_total_carbon(state, grid.W)

compute_total_carbon(state::AggregateState, r_grid::AbstractVector{<:Real}, h::Real) =
    compute_total_carbon(state, conservation_weights(r_grid, h))

"""
    compute_carbon_balance_error(state, grid::GridInfo, C_initial)
    compute_carbon_balance_error(state, r_grid, h, C_initial)

Carbon mass balance error, relative: `(C_total - C_initial) / C_initial`.
"""
compute_carbon_balance_error(state::AggregateState, grid::GridInfo, C_initial::Real) =
    (compute_total_carbon(state, grid) - C_initial) / C_initial

compute_carbon_balance_error(state::AggregateState, r_grid::AbstractVector{<:Real},
                             h::Real, C_initial::Real) =
    (compute_total_carbon(state, r_grid, h) - C_initial) / C_initial

"""
    setup_run(n_grid, r_0, r_max, T_func, ψ_func, O2_func,
              bio, soil, ic, initial_state, P_0, ω) -> (grid, env, state)

The three objects both solvers need before they can start, built once here
rather than twice. `run_aggregate` and `run_aggregate_stiff` had identical
copies of this block; identical is how a divergence starts.

`initial_state` takes priority over `ic`. It is deep-copied, so a caller can
reuse one state across runs without the first run mutating it.
"""
function setup_run(n_grid::Int, r_0::Real, r_max::Real,
                   T_func, ψ_func, O2_func,
                   bio::BiologicalProperties, soil::SoilProperties,
                   ic::InitialConditions, initial_state, P_0, ω::Real)
    grid  = GridInfo(n_grid, r_0, r_max)
    env   = EnvironmentalDrivers(T_func, ψ_func, O2_func)
    state = initial_state !== nothing ? deepcopy(initial_state) :
            create_initial_state(n_grid, bio, soil, ic; P_0 = P_0, ω = ω)
    return grid, env, state
end

"""
    default_output_times(t_start, t_end) -> Vector{Float64}

Output schedule when the caller gives none: every `min(1, span/10)` days, with
the end point always included.

Both solvers use this rule. They disagreed until 2026-07-30 — the split solver
recorded on this interval while the stiff solver recorded exactly two points,
the start and the end — so the same call returned a trajectory or a pair
depending on the integrator.
"""
function default_output_times(t_start::Real, t_end::Real)
    interval = min(1.0, (Float64(t_end) - Float64(t_start)) / 10)
    interval > 0.0 || return Float64[t_start, t_end]
    ts = collect(Float64(t_start):interval:Float64(t_end))
    (isempty(ts) || ts[end] < t_end) && push!(ts, Float64(t_end))
    return ts
end
