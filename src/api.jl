# api.jl
# User-facing API for running aggregate simulations.
#
# `run_aggregate` here is the Strang-split REFERENCE implementation.
# `run_aggregate_stiff` (solver/mol_solve.jl) is the default workhorse.
# See the header of solver/timestepper.jl for why both are kept.

"""
    run_aggregate(bio::BiologicalProperties, soil::SoilProperties,
                  T_func, ψ_func, O2_func, t_span;
                  n_grid::Int=50, r_0::Real=0.1, r_max::Real=2.0,
                  initial_state=nothing, ic::Union{Nothing,InitialConditions}=nothing,
                  P_0::Union{Nothing,Real}=nothing,
                  dt_initial::Real=0.01, dt_min::Real=1e-4, dt_max::Real=0.1,
                  output_times::Vector{<:Real}=Float64[])

Run a soil aggregate biogeochemical simulation with the Strang-split solver.

!!! note "This is the reference implementation, not the default"
    [`run_aggregate_stiff`](@ref) is the default workhorse — 24x faster over
    45 days and 188x fewer steps, because an implicit method sizes its step
    from accuracy rather than from an explicit-reaction stability guard. Use
    this one to cross-check, or when you want the independent carbon-closure
    probe described in `docs/REFERENCE.md` 17a. Both call the same physics.

# Arguments
- `bio::BiologicalProperties`: Biological parameters
- `soil::SoilProperties`: Soil properties
- `T_func`: Function T(t) returning temperature [K]
- `ψ_func`: Function ψ(t) returning water potential [kPa]
- `O2_func`: Function O2(t) returning ambient O₂ concentration [μg/mm³]
- `t_span`: Tuple (t_start, t_end) [days]

# Keyword Arguments
- `n_grid::Int`: Number of radial grid points (default: 50)
- `r_0::Real`: POM radius [mm] (default: 0.1)
- `r_max::Real`: Outer aggregate radius [mm] (default: 2.0)
- `initial_state`: Pre-built AggregateState (highest priority)
- `ic::InitialConditions`: SOC-based initialization (used if initial_state is nothing)
- `P_0::Union{Nothing,Real}`: Initial POM carbon mass [μg-C] (overrides bio.P_0)
- `t_delay::Real`: POM activation delay [days] (default: 0.0). When > 0, POM
  dissolution is suppressed for ~t_delay days, allowing microbial pre-equilibration.
- `dt_initial::Real`: Initial timestep [days] (default: 0.01)
- `dt_min::Real`: Minimum timestep [days] (default: 1e-4)
- `dt_max::Real`: Maximum timestep [days] (default: 0.1)
- `output_times::Vector{<:Real}`: Specific times to save output

# State initialization priority
1. `initial_state` if provided (deepcopy)
2. `ic::InitialConditions` if provided (SOC partitioning)
3. Fallback: legacy minimal initialization (seed values)

# Returns
`SimulationResult` with fields:
- `grid::GridInfo`: Grid geometry and conservation weights
- `params::ParameterSet`: Parameters used (bio and soil)
- `env::EnvironmentalDrivers`: Environmental driver functions
- `outputs::Vector{OutputRecord}`: State snapshots with diagnostics
- `diagnostics::Dict`: Simulation diagnostics (n_steps, final_time, etc.)

# Examples

```julia
# SOC-based initialization (recommended)
bio  = BiologicalProperties()
soil = SoilProperties(f_clay_silt = 0.42)
ic   = InitialConditions(SOC = 18.3)

result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 60.0);
                       n_grid=200, ic=ic, r_0=0.5, r_max=12.5)

# Legacy initialization (backward compatible)
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 60.0))
```

# Manuscript reference
Architecture §14: Main loop, output collection
"""
function run_aggregate(bio::BiologicalProperties, soil::SoilProperties,
                      T_func, ψ_func, O2_func, t_span;
                      n_grid::Int=50, r_0::Real=0.1, r_max::Real=2.0,
                      initial_state=nothing,
                      ic::InitialConditions=InitialConditions(),
                      P_0::Union{Nothing,Real}=nothing,
                      ω::Real=1.0,
                      t_delay::Real=0.0,
                      dt_initial::Real=0.01, dt_min::Real=1e-4, dt_max::Real=0.1,
                      output_times::Vector{<:Real}=Float64[])
    # Unpack time span
    t_start, t_end = t_span

    # === Setup grid ===
    grid = GridInfo(n_grid, r_0, r_max)

    # === Environment ===
    env = EnvironmentalDrivers(T_func, ψ_func, O2_func)

    # === Initialize state (priority: initial_state > ic > legacy) ===
    if initial_state !== nothing
        state = deepcopy(initial_state)
    else
        state = create_initial_state(n_grid, bio, soil, ic; P_0=P_0, ω=ω)
    end

    # === Create workspace ===
    workspace = Workspace(n_grid)

    # === Output schedule ===
    if isempty(output_times)
        interval = min(1.0, (t_end - t_start) / 10)
        scheduled = Float64[]
    else
        interval = NaN
        scheduled = Float64.(output_times)
    end

    # === Compute initial total carbon (true reference for balance) ===
    C_total_initial = compute_total_carbon(state, grid)

    # === Run simulation ===
    result = run_simulation(
        state, workspace, grid.r_grid, grid.h, bio, soil,
        T_func, ψ_func, O2_func,
        t_start, t_end, dt_initial;
        t_delay=t_delay,
        dt_min=dt_min, dt_max=dt_max,
        output_interval=interval,
        output_times=scheduled
    )

    # === Package outputs ===
    outputs = OutputRecord[]
    for i in 1:length(result.times)
        error = compute_carbon_balance_error(result.states[i], grid, C_total_initial)
        push!(outputs, OutputRecord(result.times[i], result.states[i], error))
    end

    params = ParameterSet(bio, soil)
    diagnostics_any = Dict{String, Any}(result.diagnostics)
    return SimulationResult(grid, params, env, outputs, diagnostics_any)
end

"""
    create_initial_state_legacy(n::Int, bio::BiologicalProperties,
                                soil::SoilProperties;
                                P_0::Union{Nothing,Real}=nothing,
                                T_0::Real=293.15, ψ_0::Real=-29.0,
                                O2_gas::Real=0.2785)

Legacy initialization with seed values (backward compatible).

Uses minimum viable microbial concentrations and computes O₂ and MAOC
from first principles, but does NOT partition from measured SOC.
Use `create_initial_state(n, bio, soil, ic)` for SOC-based initialization.

# Notes
- Kept for backward compatibility with existing scripts
- New scripts should use InitialConditions-based initialization
"""
function create_initial_state_legacy(n::Int, bio::BiologicalProperties,
                                     soil::SoilProperties;
                                     P_0::Union{Nothing,Real}=nothing,
                                     T_0::Real=293.15, ψ_0::Real=-29.0,
                                     O2_gas::Real=0.2785)
    state = AggregateState(n)

    # Background DOC
    C_0 = 1e-4  # [µg/mm³]
    state.C .= C_0

    # Minimum viable microbial seed
    state.B   .= bio.B_min
    state.F_n .= bio.F_n_min
    state.F_m .= bio.F_m_min
    state.F_i .= bio.F_i_min

    # Pools that develop over time
    state.E .= 0.0

    # Water content at ψ_0
    # van_genuchten, not an inline copy: the inline form raised a negative base
    # to the power n_vg at saturation (ψ_0 = 0, which van_genuchten_inverse
    # returns exactly) and threw DomainError. CLAUDE.md §8.
    θ_0 = van_genuchten(ψ_0, soil.α_vg, soil.n_vg, soil.θ_r, soil.θ_s)
    θ_a_0 = soil.θ_s - θ_0

    # MAOC at equilibrium with initial DOC
    C_aq_0 = C_aqueous(C_0, θ_0, soil)
    C_eq_0 = soil.k_d_eq * C_aq_0
    state.M .= M_eq_langmuir_freundlich(C_eq_0, maoc_capacity(soil), soil.k_L, soil.n_LF)

    # Oxygen: total soil (dissolved + gas)
    K_H = K_H_O2(T_0)
    state.O .= O2_gas * (θ_0 + K_H * θ_a_0) / K_H

    # POM
    state.P = isnothing(P_0) ? bio.P_0 : Float64(P_0)
    state.P_0 = state.P  
    state.CO2_cumulative = 0.0

    return state
end

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
