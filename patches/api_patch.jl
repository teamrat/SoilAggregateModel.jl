# api.jl — PATCH
# Replace the run_aggregate function signature and IC block.
#
# Changes:
#   1. Add ic::InitialConditions and ω::Real keywords
#   2. Remove initial_state and P_0 from signature (ic handles both)
#   3. IC construction: initial_state if provided, else build from ic
#
# ---- NEW SIGNATURE ----

function run_aggregate(bio::BiologicalProperties, soil::SoilProperties,
                      T_func, ψ_func, O2_func, t_span;
                      n_grid::Int=50, r_0::Real=0.1, r_max::Real=2.0,
                      initial_state=nothing,
                      ic::InitialConditions=InitialConditions(),
                      P_0::Union{Nothing,Real}=nothing,
                      ω::Real=1.0,
                      dt_initial::Real=0.01, dt_min::Real=1e-4, dt_max::Real=0.1,
                      output_times::Vector{<:Real}=Float64[])

    # ... (unchanged until IC block) ...

    # === Initialize state ===
    if initial_state !== nothing
        state = deepcopy(initial_state)
    else
        state = create_initial_state(n_grid, bio, soil, ic; P_0=P_0, ω=ω)
    end

    # ... (rest unchanged) ...
end
