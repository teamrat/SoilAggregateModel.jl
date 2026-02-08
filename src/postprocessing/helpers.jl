# helpers.jl
# Internal helper functions for post-processing

"""
    _prepare_environment(rec::OutputRecord, grid::GridInfo,
                        params::ParameterSet, env::EnvironmentalDrivers)

Internal helper that prepares environmental quantities for derived functions.

This avoids code duplication across post-processing functions that need
water content and temperature-dependent quantities.

# Returns
Named tuple with:
- `T::Float64`: Temperature at rec.t [K]
- `ψ::Float64`: Water potential at rec.t [kPa]
- `θ::Vector{Float64}`: Water content at each grid point [dimensionless]
- `θ_a::Vector{Float64}`: Air-filled porosity at each grid point [dimensionless]
- `f_T::TemperatureCache`: All temperature-dependent factors

# Usage (internal only, not exported)
```julia
local_env = _prepare_environment(rec, grid, params, env)
# Now use local_env.θ, local_env.f_T, etc.
```

# Notes
- Not exported - internal helper only
- Allocates θ and θ_a arrays (one-time cost per function call)
- Used by: aqueous_concentrations, respiration_rates, carbon_use_efficiency
"""
function _prepare_environment(rec::OutputRecord, grid::GridInfo,
                              params::ParameterSet, env::EnvironmentalDrivers)
    # Evaluate environmental forcings at snapshot time
    T = env.T(rec.t)
    ψ = env.ψ(rec.t)

    # Allocate and populate water content arrays
    θ = Vector{Float64}(undef, grid.n)
    θ_a = Vector{Float64}(undef, grid.n)
    update_water_content!(θ, θ_a, ψ, rec.state, params.soil)

    # Create and populate temperature cache
    f_T = TemperatureCache()
    update_temperature_cache!(f_T, T, params.bio, params.soil)

    return (T=T, ψ=ψ, θ=θ, θ_a=θ_a, f_T=f_T)
end
