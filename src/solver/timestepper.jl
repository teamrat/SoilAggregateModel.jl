# timestepper.jl
# Main time-stepping loop with Strang splitting and adaptive timestep

"""
    run_simulation(state::AggregateState, workspace::Workspace,
                   r_grid::Vector{Float64}, h::Real,
                   bio::BiologicalProperties, soil::SoilProperties,
                   T_func, ψ_func, O2_func,
                   t_start::Real, t_end::Real, dt_initial::Real;
                   dt_min::Real=1e-4, dt_max::Real=0.1,
                   output_interval::Real=1.0)

Run aggregate simulation with Strang splitting and adaptive timestep.

# Arguments
- `state::AggregateState`: Initial state (modified in-place)
- `workspace::Workspace`: Pre-allocated workspace
- `r_grid::Vector{Float64}`: Radial grid points [mm]
- `h::Real`: Grid spacing [mm]
- `bio::BiologicalProperties`: Biological parameters
- `soil::SoilProperties`: Soil parameters
- `T_func`: Function T(t) returning temperature [K]
- `ψ_func`: Function ψ(t) returning water potential [kPa]
- `O2_func`: Function O2(t) returning ambient O₂ [μg/mm³]
- `t_start::Real`: Start time [days]
- `t_end::Real`: End time [days]
- `dt_initial::Real`: Initial timestep [days]
- `dt_min::Real`: Minimum allowed timestep [days] (default: 1e-4)
- `dt_max::Real`: Maximum allowed timestep [days] (default: 0.1)
- `output_interval::Real`: Time between outputs [days] (default: 1.0)

# Returns
Named tuple with:
- `times::Vector{Float64}`: Output times
- `states::Vector{AggregateState}`: Saved states at output times
- `diagnostics::Dict`: Diagnostics (total steps, rejections, etc.)

# Algorithm
Strang splitting (2nd-order accurate):
1. Diffusion half-step (Δt/2)
2. Reaction full-step (Δt)
3. Diffusion half-step (Δt/2)

Adaptive timestep:
- If max(|S × Δt / u|) > 0.10 at any node → halve Δt
- If max(|S × Δt / u|) < 0.01 everywhere → double Δt
- Enforce dt_min ≤ Δt ≤ dt_max

# Manuscript reference
Architecture §2: Time integration, Strang splitting
"""
function run_simulation(state::AggregateState, workspace::Workspace,
                       r_grid::Vector{Float64}, h::Real,
                       bio::BiologicalProperties, soil::SoilProperties,
                       T_func, ψ_func, O2_func,
                       t_start::Real, t_end::Real, dt_initial::Real;
                       dt_min::Real=1e-4, dt_max::Real=0.1,
                       output_interval::Real=1.0)
    # Initialize
    t = t_start
    dt = dt_initial
    next_output = t_start + output_interval

    # Storage for outputs
    output_times = Float64[t_start]
    output_states = [deepcopy(state)]

    # Diagnostics
    n_steps = 0
    n_rejected = 0

    # Main loop
    while t < t_end
        # Don't overshoot final time
        if t + dt > t_end
            dt = t_end - t
        end

        # Don't overshoot output time
        if t + dt > next_output
            dt = next_output - t
        end

        # === Get environmental conditions ===
        T = T_func(t)
        ψ = ψ_func(t)
        O2_amb = O2_func(t)

        # === Update workspace (once per timestep) ===
        # Temperature cache
        update_temperature_cache!(workspace.f_T, T, bio, soil)

        # Water content
        update_water_content!(workspace.θ, workspace.θ_a, ψ, state, soil)

        # Effective diffusion coefficients
        update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

        # === Compute POM dissolution (once per timestep, Strang-consistent) ===
        # Compute flux density at beginning of step
        B_0 = state.B[1]
        F_n_0 = state.F_n[1]
        θ_0 = workspace.θ[1]
        θ_a_0 = workspace.θ_a[1]
        O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
        R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
        J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                     bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)

        # Convert flux density to total rate [μg-C/day]
        R_P_val = R_P(J_P_val, bio.r_0)

        # === Strang splitting ===
        # 1. Diffusion half-step (flux enters C via BC)
        diffusion_step!(state, workspace, dt/2, r_grid, h, J_P_val, O2_amb)

        # 2. Reaction full-step (P decreases by same amount)
        # Returns max_rel_change for adaptive timestep control
        max_rel_change = reaction_step!(state, workspace, dt, r_grid, h, bio, soil, ψ, R_P_val)

        # 3. Diffusion half-step
        diffusion_step!(state, workspace, dt/2, r_grid, h, J_P_val, O2_amb)

        # === Advance time ===
        t += dt
        n_steps += 1

        # === Adaptive timestep ===
        dt_new = adapt_timestep(max_rel_change, dt, dt_min, dt_max)
        if dt_new < dt && dt > dt_min
            # Would have reduced dt → this step might be inaccurate
            # For now, accept and reduce next step
            # (More sophisticated: reject and re-do with smaller dt)
            n_rejected += 1
        end
        dt = dt_new

        # === Save output if at output time ===
        if abs(t - next_output) < 1e-10 || t >= t_end
            push!(output_times, t)
            push!(output_states, deepcopy(state))
            next_output += output_interval
        end
    end

    # Return results
    diagnostics = Dict(
        "n_steps" => n_steps,
        "n_rejected" => n_rejected,
        "final_time" => t
    )

    return (times=output_times, states=output_states, diagnostics=diagnostics)
end

"""
    adapt_timestep(max_rel_change::Real, dt::Real, dt_min::Real, dt_max::Real)

Adapt timestep based on relative change from the reaction step.

# Adaptive criteria
- If max_rel_change > 0.10 → halve Δt
- If max_rel_change < 0.01 → double Δt
- Enforce dt_min ≤ Δt ≤ dt_max

# Arguments
- `max_rel_change::Real`: Maximum relative change max(|S × Δt / u|) from reaction step
- `dt::Real`: Current timestep [days]
- `dt_min::Real`: Minimum allowed timestep [days]
- `dt_max::Real`: Maximum allowed timestep [days]

# Returns
- New timestep [days]

# Notes
- Uses max_rel_change computed by reaction_step! (eliminates redundant source term computation)
- Simple halve/double logic with bounds enforcement
"""
function adapt_timestep(max_rel_change::Real, dt::Real, dt_min::Real, dt_max::Real)
    # Adaptive logic
    dt_new = dt

    if max_rel_change > 0.10
        # Too large → halve
        dt_new = dt / 2.0
    elseif max_rel_change < 0.01
        # Too small → double
        dt_new = dt * 2.0
    end

    # Enforce bounds
    dt_new = max(dt_min, min(dt_max, dt_new))

    return dt_new
end
