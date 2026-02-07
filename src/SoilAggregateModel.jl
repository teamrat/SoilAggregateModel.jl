"""
    SoilAggregateModel

Biogeochemical model for soil aggregate formation and carbon cycling.

Spherical reaction-diffusion system with Strang splitting:
- 9 state variables (5 diffusing, 3 immobile ODEs, 1 scalar ODE)
- Custom Crank-Nicolson solver (Thomas algorithm)
- Zero-allocation hot loop

# Main API
- [`run_aggregate`](@ref): Run aggregate simulation with environmental forcings
- [`BiologicalProperties`](@ref): Biological parameters
- [`SoilProperties`](@ref): Soil physical/chemical properties

# Units
μg/mm³ (≡ kg/m³), mm, days, kPa, K, J/mol throughout.
"""
module SoilAggregateModel

# Include all source files in dependency order
include("constants.jl")
include("types.jl")
include("parameters.jl")
include("environment.jl")

# Temperature dependencies
include("temperature/arrhenius.jl")
include("temperature/diffusion_pure.jl")
include("temperature/henry.jl")

# Physics
include("physics/water_retention.jl")
include("physics/effective_diffusion.jl")

# Biology
include("biology/bacteria.jl")
include("biology/fungi.jl")
include("biology/eps.jl")
include("biology/maoc.jl")

# Carbon
include("carbon/pom_dissolution.jl")

# Solver
include("solver/tridiagonal.jl")
include("solver/crank_nicolson.jl")
include("solver/finite_volumes.jl")
include("solver/reactions.jl")
include("solver/diffusion_step.jl")
include("solver/reaction_step.jl")
include("solver/workspace_updates.jl")
include("solver/timestepper.jl")

# API
include("api.jl")

# Exports
export BiologicalProperties, SoilProperties
export AggregateState, OutputRecord
export run_aggregate

end  # module SoilAggregateModel
