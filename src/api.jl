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
Named tuple with:
- `times::Vector{Float64}`: Output times [days]
- `outputs::Vector{OutputRecord}`: State snapshots with diagnostics
- `diagnostics::Dict`: Simulation diagnostics (steps, rejections, etc.)

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
    h = (r_max - r_0) / (n_grid - 1)
    r_grid = [r_0 + i*h for i in 0:n_grid-1]

    # === Initialize state ===
    if initial_state === nothing
        state = create_initial_state(n_grid, bio)
    else
        state = deepcopy(initial_state)
    end

    # === Create workspace ===
    workspace = Workspace(n_grid)

    # === Determine output times ===
    if isempty(output_times)
        # Default: output every day
        output_interval = min(1.0, (t_end - t_start) / 10)
    else
        output_interval = NaN  # Will use explicit output_times
    end

    # === Compute initial total carbon (true reference for balance) ===
    C_total_initial = compute_total_carbon(state, r_grid, h)

    # === Run simulation ===
    result = run_simulation(
        state, workspace, r_grid, h, bio, soil,
        T_func, ψ_func, O2_func,
        t_start, t_end, dt_initial;
        dt_min=dt_min, dt_max=dt_max,
        output_interval=output_interval
    )

    # === Convert to OutputRecord format with diagnostics ===
    outputs = OutputRecord[]
    for i in 1:length(result.times)
        t = result.times[i]
        state_snapshot = result.states[i]

        # Compute carbon balance diagnostic
        mass_balance_error = compute_carbon_balance_error(
            state_snapshot, r_grid, h, C_total_initial
        )

        push!(outputs, OutputRecord(t, state_snapshot, mass_balance_error))
    end

    return (times=result.times, outputs=outputs, diagnostics=result.diagnostics)
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
- C, B, F_n, F_m, F_i, E, M: Small uniform values (1.0 μg/mm³)
- O: Ambient level (0.3 μg/mm³)
- P: Initial POM mass from bio.P_0
- CO2_cumulative: 0.0
"""
function create_initial_state(n::Int, bio::BiologicalProperties)
    state = AggregateState(n)

    # Small uniform pools
    state.C .= 1.0      # DOC
    state.B .= 1.0      # Bacteria
    state.F_n .= 1.0    # Non-insulated fungi
    state.F_m .= 0.1    # Mobile fungi
    state.F_i .= 0.5    # Insulated fungi
    state.E .= 1.0      # EPS
    state.M .= 1.0      # MAOC
    state.O .= 0.3      # Oxygen (ambient)

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
