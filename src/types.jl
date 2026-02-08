"""
    types.jl

Core type definitions for the soil aggregate biogeochemical model.

All types follow the architecture specification (ARCHITECTURE_CLAUDE_CODE.md).
Units: μg/mm³ (= kg/m³), mm, days, kPa, K, J/mol throughout.
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
- `D_Fm`: Mobile fungi translocation rate [mm²/day] (spatially uniform, no tortuosity)
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
    D_Fm::Float64        # = D_Fm0_ref × f_fun (spatially uniform)

    # Henry's law
    K_H_O::Float64
end

"""
    TemperatureCache()

Create an uninitialized TemperatureCache. All fields set to NaN.
Use `update_temperature_cache!` to populate with temperature-dependent values.
"""
function TemperatureCache()
    TemperatureCache(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
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
    CO2_cumulative::Float64 # diagnostic — total CO₂ respired
end

"""
    AggregateState(n::Int)

Create an AggregateState with uninitialized vectors of length `n`.
Use `initialize_state!` to populate with initial conditions.
"""
function AggregateState(n::Int)
    AggregateState(
        Vector{Float64}(undef, n),  # C
        Vector{Float64}(undef, n),  # B
        Vector{Float64}(undef, n),  # F_n
        Vector{Float64}(undef, n),  # F_m
        Vector{Float64}(undef, n),  # O
        Vector{Float64}(undef, n),  # F_i
        Vector{Float64}(undef, n),  # E
        Vector{Float64}(undef, n),  # M
        0.0,                        # P
        0.0                         # CO2_cumulative
    )
end

#═══════════════════════════════════════════════════════════════════════════════
# Workspace (Pre-Allocated Arrays)
#═══════════════════════════════════════════════════════════════════════════════

"""
    Workspace

Pre-allocated workspace arrays. Zero allocations in the hot loop.

Fields:
- Tridiagonal system (reused for each of 5 diffusing species):
  - `lower`: Lower diagonal [n-1]
  - `diag`: Main diagonal [n]
  - `upper`: Upper diagonal [n-1]
  - `rhs`: Right-hand side / solution vector [n]

- Spatially varying quantities (updated once per timestep):
  - `θ`: Water content [-] [n]
  - `θ_a`: Air-filled porosity [-] [n]
  - `D_C`: Effective DOC diffusion [mm²/day] [n]
  - `D_B`: Effective bacterial diffusion [mm²/day] [n]
  - `D_Fn`: Effective non-insulated fungal diffusion [mm²/day] [n]
  - `D_O`: Effective oxygen diffusion [mm²/day] [n]

- Temperature cache:
  - `f_T`: TemperatureCache struct

Note: D_Fm (mobile fungi) is spatially uniform (constant) but stored as vector for solver compatibility.
"""
struct Workspace
    # Tridiagonal system (reused for each of 5 diffusing species)
    lower::Vector{Float64}   # n-1
    diag::Vector{Float64}    # n
    upper::Vector{Float64}   # n-1
    rhs::Vector{Float64}     # n

    # Spatially varying quantities (updated once per timestep)
    θ::Vector{Float64}       # n — water content
    θ_a::Vector{Float64}     # n — air-filled porosity
    D_C::Vector{Float64}     # n — effective C diffusion
    D_B::Vector{Float64}     # n — effective B diffusion
    D_Fn::Vector{Float64}    # n — effective F_n diffusion
    D_Fm::Vector{Float64}    # n — effective F_m diffusion (spatially uniform)
    D_O::Vector{Float64}     # n — effective O diffusion

    # Temperature cache
    f_T::TemperatureCache
end

"""
    Workspace(n::Int)

Create a Workspace with pre-allocated arrays for `n` grid points.
"""
function Workspace(n::Int)
    Workspace(
        Vector{Float64}(undef, n-1),  # lower
        Vector{Float64}(undef, n),    # diag
        Vector{Float64}(undef, n-1),  # upper
        Vector{Float64}(undef, n),    # rhs
        Vector{Float64}(undef, n),    # θ
        Vector{Float64}(undef, n),    # θ_a
        Vector{Float64}(undef, n),    # D_C
        Vector{Float64}(undef, n),    # D_B
        Vector{Float64}(undef, n),    # D_Fn
        Vector{Float64}(undef, n),    # D_Fm
        Vector{Float64}(undef, n),    # D_O
        TemperatureCache()             # f_T
    )
end

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
        state.CO2_cumulative
    )
end

#═══════════════════════════════════════════════════════════════════════════════
# Grid Information
#═══════════════════════════════════════════════════════════════════════════════

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
    W = [4.0 * π * r_grid[i]^2 * h for i in 1:n]
    GridInfo(r_grid, h, Float64(r_0), Float64(r_max), n, W)
end

