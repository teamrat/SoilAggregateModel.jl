# api.jl
# User-facing API for running aggregate simulations

"""
    run_aggregate(bio::BiologicalProperties, soil::SoilProperties,
                  T_func, ψ_func, O2_func, t_span;
                  n_grid::Int=50, r_0::Real=0.1, r_max::Real=2.0,
                  initial_state=nothing,
                  dt_initial::Real=0.01, dt_min::Real=1e-4, dt_max::Real=0.1,
                  output_times::Vector{<:Real}=Float64[])

Run a soil aggregate biogeochemical simulation.

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
- `initial_state`: Initial AggregateState (default: uniform pools from bio parameters)
- `dt_initial::Real`: Initial timestep [days] (default: 0.01)
- `dt_min::Real`: Minimum timestep [days] (default: 1e-4)
- `dt_max::Real`: Maximum timestep [days] (default: 0.1)
- `output_times::Vector{<:Real}`: Specific times to save output (default: regular intervals)

# Returns
`SimulationResult` with fields:
- `grid::GridInfo`: Grid geometry and conservation weights
- `params::ParameterSet`: Parameters used (bio and soil)
- `env::EnvironmentalDrivers`: Environmental driver functions
- `outputs::Vector{OutputRecord}`: State snapshots with diagnostics
- `diagnostics::Dict`: Simulation diagnostics (n_steps, final_time, etc.)

# Examples

```julia
# Constant environment
bio = BiologicalProperties()
soil = SoilProperties()
T(t) = 293.15  # 20°C
ψ(t) = -10.0   # kPa
O2(t) = 0.3    # μg/mm³

result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))

# Access outputs
for rec in result.outputs
    println("Time ", rec.t, ": POM = ", rec.state.P, " μg-C")
end

# Access grid information
println("Grid points: ", result.grid.n)
println("Conservation weights: ", result.grid.W)
```

# Manuscript reference
Architecture §14: Main loop, output collection
"""
function run_aggregate(bio::BiologicalProperties, soil::SoilProperties,
                      T_func, ψ_func, O2_func, t_span;
                      n_grid::Int=50, r_0::Real=0.1, r_max::Real=2.0,
                      initial_state=nothing,
                      dt_initial::Real=0.01, dt_min::Real=1e-4, dt_max::Real=0.1,
                      output_times::Vector{<:Real}=Float64[])
    # Unpack time span
    t_start, t_end = t_span

    # === Setup grid ===
    grid = GridInfo(n_grid, r_0, r_max)

    # === Environment ===
    env = EnvironmentalDrivers(T_func, ψ_func, O2_func)

    # === Initialize state ===
    if initial_state === nothing
        state = create_initial_state(n_grid, bio)
    else
        state = deepcopy(initial_state)
    end

    # === Create workspace ===
    workspace = Workspace(n_grid)

    # === Output schedule ===
    if isempty(output_times)
        # Interval-based output (default)
        interval = min(1.0, (t_end - t_start) / 10)
        scheduled = Float64[]
    else
        # User-specified output times
        interval = NaN  # Not used when output_times is specified
        scheduled = Float64.(output_times)
    end

    # === Compute initial total carbon (true reference for balance) ===
    C_total_initial = compute_total_carbon(state, grid.r_grid, grid.h)

    # === Run simulation ===
    result = run_simulation(
        state, workspace, grid.r_grid, grid.h, bio, soil,
        T_func, ψ_func, O2_func,
        t_start, t_end, dt_initial;
        dt_min=dt_min, dt_max=dt_max,
        output_interval=interval,
        output_times=scheduled
    )

    # === Package outputs ===
    outputs = OutputRecord[]
    for i in 1:length(result.times)
        error = compute_carbon_balance_error(result.states[i], grid.r_grid, grid.h, C_total_initial)
        push!(outputs, OutputRecord(result.times[i], result.states[i], error))
    end

    params = ParameterSet(bio, soil)
    # Convert diagnostics to Dict{String, Any} to match SimulationResult field type
    diagnostics_any = Dict{String, Any}(result.diagnostics)
    return SimulationResult(grid, params, env, outputs, diagnostics_any)
end

"""
    create_initial_state(n::Int, bio::BiologicalProperties)

Create default initial state with uniform pools.

# Arguments
- `n::Int`: Number of grid points
- `bio::BiologicalProperties`: Biological parameters for initial POM mass

# Returns
- `AggregateState`: Initial state with small uniform pools and full POM

# Notes
Physically realistic initialization for POM-driven aggregate formation:
- C: Background DOC (0.01 μg/mm³)
- B: Minimum viable bacteria (bio.B_min, ~0.1 μg/mm³)
- F_n, F_m: Small fungal seed (B_min/10, ~0.01 μg/mm³)
- F_i: No insulated fungi initially (0.0 μg/mm³, develops over time)
- E: No EPS initially (0.0 μg/mm³, produced by bacteria)
- M: No MAOC initially (0.0 μg/mm³, accumulates from DOC sorption)
- O: Ambient O₂ (0.27 μg/mm³, ~21% atmospheric)
- P: Initial POM mass from bio.P_0
- CO2_cumulative: 0.0

This initialization ensures POM dissolution is the primary carbon source,
allowing microbial populations to grow from seed values and aggregate
structure to emerge organically.
"""
function create_initial_state(n::Int, bio::BiologicalProperties)
    state = AggregateState(n)

    # Background DOC (trace concentration)
    state.C .= 0.01

    # Minimum viable microbial seed
    state.B .= bio.B_min           # Bacteria at minimum viable
    state.F_n .= bio.B_min / 10    # Small fungal seed
    state.F_m .= bio.B_min / 100   # Even smaller mobile pool
    state.F_i .= 0.0               # Insulated fungi develop over time

    # Pools that develop over time
    state.E .= 0.0    # EPS produced by bacteria
    state.M .= 0.0    # MAOC accumulates from DOC sorption

    # Ambient oxygen
    state.O .= 0.27   # ~21% atmospheric O₂

    # Initial POM from parameters
    state.P = bio.P_0

    # No CO2 yet
    state.CO2_cumulative = 0.0

    return state
end

"""
    compute_total_carbon(state::AggregateState, r_grid::Vector{Float64}, h::Real)

Compute total carbon in the system.

# Arguments
- `state::AggregateState`: Current state
- `r_grid::Vector{Float64}`: Radial grid [mm]
- `h::Real`: Grid spacing [mm]

# Returns
- Total carbon [μg-C]

# Notes
- Uses conservation weights W_i = 4πr_i²h
- C_total = P + ∑(C+B+F_n+F_m+F_i+E+M)×W_i + CO2
"""
function compute_total_carbon(state::AggregateState, r_grid::Vector{Float64}, h::Real)
    # Sum dissolved + biomass pools
    integral = 0.0
    for i in 1:length(state.C)
        C_pools = state.C[i] + state.B[i] + state.F_n[i] + state.F_m[i] +
                 state.F_i[i] + state.E[i] + state.M[i]
        # Conservation weight (matches spherical Laplacian stencil)
        W_i = 4.0 * π * r_grid[i]^2 * h
        integral += C_pools * W_i
    end

    # Total carbon = POM + dissolved/biomass + respired
    return state.P + integral + state.CO2_cumulative
end

"""
    compute_carbon_balance_error(state::AggregateState, r_grid::Vector{Float64},
                                 h::Real, C_initial::Real)

Compute carbon mass balance error for diagnostic purposes.

# Arguments
- `state::AggregateState`: Current state
- `r_grid::Vector{Float64}`: Radial grid [mm]
- `h::Real`: Grid spacing [mm]
- `C_initial::Real`: Initial total carbon [μg-C]

# Returns
- Relative error: (C_total - C_initial) / C_initial

# Notes
- Uses conservation weights W_i = 4πr_i²h
- C_total = P + ∑(C+B+F_n+F_m+F_i+E+M)×W_i + CO2
- Should be O(machine precision) for correct implementation

# Manuscript reference
Architecture §16: Validation, carbon conservation
"""
function compute_carbon_balance_error(state::AggregateState, r_grid::Vector{Float64},
                                     h::Real, C_initial::Real)
    C_total = compute_total_carbon(state, r_grid, h)
    return (C_total - C_initial) / C_initial
end
