"""
    types.jl

Core type definitions for the soil aggregate biogeochemical model.

All types follow the architecture specification (docs/ARCHITECTURE.md — stale in parts; docs/BACKLOG.md item 12).
Units: μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol throughout.

# Fill convention

Every pre-allocated array is filled with `NaN` — never zero, never `undef`.
On any correct path each array is written before it is read, so a `NaN`
surviving into a result means a path skipped its initialisation. `NaN`
propagates and fails; zero does not. A zero pool reads as empty-but-valid and a
zero diffusivity silently stops transport.
"""

#═══════════════════════════════════════════════════════════════════════════════
# Temperature Cache
#═══════════════════════════════════════════════════════════════════════════════

"""
    TemperatureCache

Stores all temperature-dependent quantities computed once per timestep.

Fields:
- `f_bac`: Arrhenius factor for bacteria (dimensionless)
- `f_fun`: Arrhenius factor for fungi (dimensionless)
- `f_eps`: Arrhenius factor for EPS degradation (dimensionless)
- `f_maoc_s`: Arrhenius factor for MAOC sorption (dimensionless)
- `f_maoc_d`: Arrhenius factor for MAOC desorption (dimensionless)
- `f_pom`: Arrhenius factor for POM dissolution (dimensionless)
- `D_O2_w`: O₂ diffusion in water [mm²/day]
- `D_DOC_w`: DOC diffusion in water [mm²/day]
- `D_O2_a`: O₂ diffusion in air [mm²/day]
- `K_H_O`: Henry's law constant for O₂ (dimensionless)
"""
mutable struct TemperatureCache
    # Arrhenius factors (dimensionless multipliers on reference rates)
    f_bac::Float64       # bacteria: μ_max_B, m_B, r_B_max
    f_fun::Float64       # fungi: μ_max_F, m_F, transitions, D_Fn0, D_Fm0
    f_eps::Float64       # EPS: μ_E_max
    f_maoc_s::Float64    # MAOC sorption: κ_s
    f_maoc_d::Float64    # MAOC desorption: κ_d
    f_pom::Float64       # POM: R_P_max

    # Pure-phase diffusion coefficients [mm²/day]
    D_O2_w::Float64
    D_DOC_w::Float64
    D_O2_a::Float64

    # Henry's law
    K_H_O::Float64
end

"""
    TemperatureCache()

Create an uninitialized TemperatureCache. All fields set to NaN.
Use `update_temperature_cache!` to populate with temperature-dependent values.
"""
function TemperatureCache()
    TemperatureCache(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
end

#═══════════════════════════════════════════════════════════════════════════════
# State Variables
#═══════════════════════════════════════════════════════════════════════════════

"""
    AggregateState

Struct-of-arrays layout for all state variables. Optimal for tridiagonal solves.

9 total state variables:
- 5 diffusing: C, B, F_n, F_m, O (each has n grid points)
- 3 immobile: F_i, E, M (each has n grid points, advanced by local ODEs)
- 1 scalar: P (POM mass)

Fields (all concentrations in μg/mm³, POM and CO2 in μg-C):
- `C`: Dissolved organic carbon [n]
- `B`: Bacteria [n]
- `F_n`: Non-insulated fungi [n]
- `F_m`: Mobile fungi [n]
- `O`: Oxygen [n]
- `F_i`: Insulated fungi [n]
- `E`: EPS (extracellular polymeric substances) [n]
- `M`: MAOC (mineral-associated organic carbon) [n]
- `P`: POM (particulate organic matter) mass [scalar, μg-C]
- `CO2_cumulative`: Total CO₂ respired [scalar, μg-C]
"""
mutable struct AggregateState
    C::Vector{Float64}      # n — dissolved organic carbon
    B::Vector{Float64}      # n — bacteria
    F_n::Vector{Float64}    # n — non-insulated fungi
    F_m::Vector{Float64}    # n — mobile fungi
    O::Vector{Float64}      # n — oxygen
    F_i::Vector{Float64}    # n — insulated fungi
    E::Vector{Float64}      # n — EPS
    M::Vector{Float64}      # n — MAOC
    P::Float64              # scalar — POM
    P_0::Float64            # Initial POM mass [µg-C] (fixed reference)
    CO2_cumulative::Float64 # diagnostic — total CO₂ respired
end

"""
    AggregateState(n::Int)

Create an AggregateState for `n` grid points, NaN-filled. Populate with
`create_initial_state`, which assigns all eight pools plus `P`, `P_0` and
`CO2_cumulative`.

`P_0` must be positive before integrating; `run_aggregate_stiff` throws otherwise.
"""
function AggregateState(n::Int)
    AggregateState(
        fill(NaN, n),  # C
        fill(NaN, n),  # B
        fill(NaN, n),  # F_n
        fill(NaN, n),  # F_m
        fill(NaN, n),  # O
        fill(NaN, n),  # F_i
        fill(NaN, n),  # E
        fill(NaN, n),  # M
        0.0,                        # P
        0.0,                        # P_0
        0.0                         # CO2_cumulative
    )
end

#═══════════════════════════════════════════════════════════════════════════════
# Workspace (Pre-Allocated Arrays)
#═══════════════════════════════════════════════════════════════════════════════

# `Workspace` was archived with the split solver on 2026-07-30. It held the
# pre-allocated per-step arrays that `timestepper.jl` reused; the stiff solver
# allocates its coefficient vectors per RHS call at `eltype(u)` so forward-mode
# AD can see them. See _archive/split_solver_20260730/.

#═══════════════════════════════════════════════════════════════════════════════
# Output
#═══════════════════════════════════════════════════════════════════════════════

"""
    OutputRecord

Snapshot of aggregate state at a specific time, with diagnostic information.

Fields:
- `t`: Time [days]
- `state`: AggregateState (deep copy of all state variables)
- `mass_balance_error`: Carbon conservation error (diagnostic)

Post-processing (aggregate radius, pool partitioning, etc.) is done on demand
from these snapshots — NOT computed during simulation.
"""
struct OutputRecord
    t::Float64
    state::AggregateState
    mass_balance_error::Float64
end

"""
    Base.copy(state::AggregateState)

Create a deep copy of an AggregateState for output recording.
"""
function Base.copy(state::AggregateState)
    AggregateState(
        copy(state.C),
        copy(state.B),
        copy(state.F_n),
        copy(state.F_m),
        copy(state.O),
        copy(state.F_i),
        copy(state.E),
        copy(state.M),
        state.P,
        state.P_0,
        state.CO2_cumulative
    )
end

#═══════════════════════════════════════════════════════════════════════════════
# Grid Information
#═══════════════════════════════════════════════════════════════════════════════

"""
    conservation_weight(r, h)
    conservation_weights(r_grid, h)

Volumetric weight of one radial node, `W = 4π r² h` [mm³].

This is the only definition of the weight. It is stencil-matched to the
spherical Laplacian discretisation: `W_i / (r_i² h²) = 4π/h` is independent of
`i`, and that cancellation is what makes the diffusive fluxes telescope exactly
so integrated pools conserve carbon. Other quadratures do not have the property
— in particular the finite-volume shell volume `(4π/3)(r₊³ − r₋³)` in
`finite_volumes.jl` is a different quantity for a different scheme and is not
interchangeable with this one.

`GridInfo` stores the vector form as `W`; use `grid.W` wherever a grid is at
hand and call `conservation_weights` only where one is not.
"""
conservation_weight(r::Real, h::Real) = 4.0 * π * r^2 * h

conservation_weights(r_grid::AbstractVector{<:Real}, h::Real) =
    [conservation_weight(r, h) for r in r_grid]

"""
    GridInfo

Immutable record of the radial grid and precomputed conservation weights.

Fields:
- `r_grid::Vector{Float64}`: Radial coordinates [mm], length n
- `h::Float64`: Grid spacing [mm]
- `r_0::Float64`: POM radius (inner boundary) [mm]
- `r_max::Float64`: Outer boundary [mm]
- `n::Int`: Number of grid points
- `W::Vector{Float64}`: Conservation weights W[i] = 4π r[i]² h [mm³], length n

The conservation weights are stencil-matched to the spherical Laplacian discretization.
All volumetric integration must use these weights.
"""
struct GridInfo
    r_grid::Vector{Float64}
    h::Float64
    r_0::Float64
    r_max::Float64
    n::Int
    W::Vector{Float64}
end

"""
    GridInfo(n::Int, r_0::Real, r_max::Real)

Construct a GridInfo with uniform spacing and precomputed conservation weights.

# Arguments
- `n::Int`: Number of grid points
- `r_0::Real`: Inner boundary (POM radius) [mm]
- `r_max::Real`: Outer boundary [mm]

# Returns
- `GridInfo`: Grid geometry with conservation weights W[i] = 4π r[i]² h
"""
function GridInfo(n::Int, r_0::Real, r_max::Real)
    h = (r_max - r_0) / (n - 1)
    r_grid = [r_0 + i * h for i in 0:n-1]
    W = conservation_weights(r_grid, h)
    GridInfo(r_grid, h, Float64(r_0), Float64(r_max), n, W)
end

