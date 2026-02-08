# radial_profiles.jl
# Extract radial profiles from simulation results

"""
    radial_profiles(result::SimulationResult; times=nothing)

Extract radial profiles at specified times (or all output times).

# Arguments
- `result::SimulationResult`: Complete simulation output
- `times`: Optional vector of times [days] to extract. If `nothing`, extracts all output times.

# Returns
Vector of NamedTuples, each containing:
- `t::Float64`: Time [days]
- `r::Vector{Float64}`: Radial grid [mm]
- `C, B, F_i, F_n, F_m, E, M, O::Vector{Float64}`: Profiles at each node [μg/mm³]
- `P::Float64`: POM mass [μg-C] (scalar)
- `CO2::Float64`: Cumulative CO₂ [μg-C] (scalar)

# Time matching
- Requested times are matched to the nearest available output snapshot
- Duplicates are removed (if multiple requested times map to same snapshot)

# Notes
- Arrays are copied (not views) - safe to modify
- Useful for plotting or detailed analysis at specific times
- No recomputation - just extraction from stored states

# Examples
```julia
# All output times
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
profiles = radial_profiles(result)
for prof in profiles
    plot(prof.r, prof.C, label="t=\$(prof.t)")
end

# Specific times
profiles = radial_profiles(result; times=[0.0, 15.0, 30.0])
@assert length(profiles) == 3
```
"""
function radial_profiles(result::SimulationResult; times=nothing)
    if times === nothing
        # Extract all output times
        indices = 1:length(result.outputs)
    else
        # Match each requested time to nearest output
        out_times = [rec.t for rec in result.outputs]
        indices = [argmin(abs.(out_times .- t)) for t in times]
        indices = unique(indices)  # Remove duplicates
    end

    profiles = map(indices) do idx
        rec = result.outputs[idx]
        s = rec.state
        (t=rec.t, r=copy(result.grid.r_grid),
         C=copy(s.C), B=copy(s.B), F_i=copy(s.F_i), F_n=copy(s.F_n),
         F_m=copy(s.F_m), E=copy(s.E), M=copy(s.M), O=copy(s.O),
         P=s.P, CO2=s.CO2_cumulative)
    end

    return profiles
end
