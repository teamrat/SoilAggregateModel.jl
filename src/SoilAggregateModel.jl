"""
    SoilAggregateModel

Biogeochemical model for soil aggregate formation and carbon cycling.

Spherical reaction-diffusion system, 9 state variables (5 diffusing, 3
immobile ODEs, 1 scalar ODE), integrated by either of two solvers over one
shared set of physics functions:

- [`run_aggregate_stiff`](@ref) — **default.** Method of lines, implicit stiff
  solver, sparse Jacobian. 45 days in 1.6 s / 2088 steps.
- [`run_aggregate`](@ref) — reference implementation. Strang splitting,
  Crank-Nicolson diffusion, Forward Euler reactions. 45 days in 38 s / 391773
  steps, because the explicit reaction bounds the step by stability rather than
  accuracy. Kept as a cross-check and for its independent carbon-closure probe.

Both call the same physics functions; only the time integration differs.
See docs/REFERENCE.md §17a and §20a.

# Main API
- [`BiologicalProperties`](@ref): Biological parameters
- [`SoilProperties`](@ref): Soil physical/chemical properties

# Units
μg/mm³ (≡ kg/m³), mm, days, kPa, K, J/mol throughout.
"""
module SoilAggregateModel

using SparseArrays: SparseMatrixCSC, sparse, nonzeros, nnz
using LinearSolve: KLUFactorization
using OrdinaryDiffEqBDF: FBDF
using OrdinaryDiffEqSDIRK: KenCarp47
using SciMLBase: ODEFunction, ODEProblem, solve

# Include all source files in dependency order
include("constants.jl")
include("types.jl")
include("parameters.jl")
include("environment.jl")
include("result.jl")

# Temperature dependencies
include("temperature/arrhenius.jl")
include("temperature/diffusion_pure.jl")
include("temperature/henry.jl")

# Physics
include("physics/particle_size.jl")
include("physics/tessellation.jl")
include("physics/water_retention.jl")
include("physics/effective_diffusion.jl")
include("physics/initial_conditions.jl")

# Math utilities
include("math_utils.jl")

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
include("solver/mol.jl")
include("solver/mol_solve.jl")

# API
include("api.jl")

# Post-processing
include("postprocessing/helpers.jl")
include("postprocessing/aggregate_radius.jl")
include("postprocessing/integration.jl")
include("postprocessing/derived.jl")
include("postprocessing/radial_profiles.jl")
include("postprocessing/population.jl")

# Exports
export BiologicalProperties, SoilProperties
export InitialConditions
export AggregateState, OutputRecord
export GridInfo, ParameterSet, SimulationResult, IntegratedPools
export run_aggregate, run_aggregate_stiff, mol_jacobian_prototype
export FBDF, KenCarp47, KLUFactorization
export sauter_mean_diameter, sauter_from_texture, TEXTURE_CLASS_DIAMETERS
export domain_tessellation, pom_population, log_interpolate_fraction
export van_genuchten, van_genuchten_inverse
export compute_r_agg, integrated_pools, carbon_balance_table, radial_profiles
export sieve_class, aggregate_mass, population_statistics
export aqueous_concentrations, maoc_equilibrium, respiration_rates, carbon_use_efficiency, co2_flux
end  # module SoilAggregateModel
