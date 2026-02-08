"""
    result.jl

Result types for simulation output. These depend on types from parameters.jl
and environment.jl, so must be included after those files.
"""

#═══════════════════════════════════════════════════════════════════════════════
# Parameter Set
#═══════════════════════════════════════════════════════════════════════════════

"""
    ParameterSet

Bundles BiologicalProperties and SoilProperties. Immutable record of what was used.

Fields:
- `bio::BiologicalProperties`: Biological parameters
- `soil::SoilProperties`: Soil physical/chemical properties
"""
struct ParameterSet
    bio::BiologicalProperties
    soil::SoilProperties
end

#═══════════════════════════════════════════════════════════════════════════════
# Simulation Result
#═══════════════════════════════════════════════════════════════════════════════

"""
    SimulationResult{E<:EnvironmentalDrivers}

Complete output from run_aggregate.

Fields:
- `grid::GridInfo`: Grid geometry and conservation weights
- `params::ParameterSet`: Parameters used
- `env::E`: Environmental driver functions
- `outputs::Vector{OutputRecord}`: State snapshots at output times
- `diagnostics::Dict{String,Any}`: n_steps, n_rejected, wall_time, etc.

# Usage

Access output times via `[rec.t for rec in result.outputs]`.
Grid information is stored in `result.grid` for post-processing.
"""
struct SimulationResult{E<:EnvironmentalDrivers}
    grid::GridInfo
    params::ParameterSet
    env::E
    outputs::Vector{OutputRecord}
    diagnostics::Dict{String,Any}
end
