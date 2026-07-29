# timestepper.jl
# Strang-split time stepping. REFERENCE IMPLEMENTATION, not the default.
#
# `run_aggregate_stiff` (solver/mol_solve.jl) is the default workhorse: 24x
# faster over 45 days, 188x fewer steps, and its step GROWS through a run where
# this one is pinned at its floor with the step its own criterion demands still
# falling. Measured, REFERENCE.md 20a.
#
# This path is kept, and should not be deleted, for two reasons:
#
#   1. It computes cumulative respired carbon independently, by accumulating
#      Resp_total, where the stiff path recovers it from the carbon balance.
#      That independence is what makes its balance error an actual probe of the
#      closure identity rather than a tautology (REFERENCE.md 17a). Once
#      test/test_closure.jl is trusted this becomes a redundancy rather than the
#      only line of defence — but redundancy on the model's central invariant is
#      worth keeping.
#   2. Two independent integrators over shared physics are the only convergence
#      evidence this model has. They agree at day 45 on every reported quantity
#      (r_agg to 7 figures, CO2 to 2e-5) and disagree by 9.6% on DOC at the
#      POM-adjacent nodes, which is a finding neither could produce alone.
#
# Do not add features here without adding them to the stiff path. They share
# every physics function; they must not diverge in what they model.

"""
    run_simulation(state::AggregateState, workspace::Workspace,
                   r_grid::Vector{Float64}, h::Real,
                   bio::BiologicalProperties, soil::SoilProperties,
                   T_func, ψ_func, O2_func,
                   t_start::Real, t_end::Real, dt_initial::Real;
                   t_delay::Real=0.0,
                   dt_min::Real=1e-4, dt_max::Real=0.1,
                   output_interval::Real=1.0,
                   output_times::Vector{Float64}=Float64[])

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
- `O2_func`: Function O2(t) returning **gas-phase** O₂ density [μg/mm³].
  This is the atmospheric O₂ partial pressure converted to mass density
  via the ideal gas law: O₂_gas = x_O₂ × P_atm × M_O₂ / (R T).
  The solver converts this internally to aqueous O₂ at the outer
  boundary: O₂_amb = O₂_gas / K_H(T).
- `t_start::Real`: Start time [days]
- `t_end::Real`: End time [days]
- `dt_initial::Real`: Initial timestep [days]
- `t_delay::Real`: POM activation delay [days] (default: 0.0, no delay).
  When > 0, POM dissolution is suppressed via sigmoid
  1/(1+exp(-(t-t_delay)/(0.1·t_delay))), allowing microbial pools to
  equilibrate before POM input begins. Transition width is 10% of delay
  (e.g., t_delay=10 → 90% activation by t ≈ 12.2).
- `dt_min::Real`: Minimum allowed timestep [days] (default: 1e-4)
- `dt_max::Real`: Maximum allowed timestep [days] (default: 0.1)
- `output_interval::Real`: Time between outputs [days] (default: 1.0)
- `output_times::Vector{Float64}`: User-specified output times (default: empty, use interval)

# Returns
Named tuple with:
- `times::Vector{Float64}`: Output times
- `states::Vector{AggregateState}`: Saved states at output times
- `diagnostics::Dict`: Diagnostics (total steps, rejections, etc.)

# Algorithm

Each timestep proceeds as:

1. **Environment**: evaluate T(t), ψ(t), O₂_gas(t)
2. **Temperature**: update Arrhenius factors, Henry's law K_H, pure-phase diffusivities
3. **Water content**: update θ(r) from van Genuchten with EPS/F_i modification
4. **Effective diffusion**: recompute D_C, D_B, D_Fn, D_Fm, D_O from new θ
5. **O₂ boundary condition**: convert gas-phase O₂ to aqueous: O₂_amb = O₂_gas / K_H
6. **POM dissolution**: compute J_P and R_P at POM surface (r = r₀)
7. **Strang splitting** (2nd-order accurate):
   a. Diffusion half-step (Δt/2)
   b. Reaction full-step (Δt)
   c. Diffusion half-step (Δt/2)
8. **Adaptive timestep**: adjust Δt based on max relative change

# O₂ state variable convention

The O₂ state variable is **aqueous concentration** C_aq,O₂ [μg/mm³], not
total O₂. This choice ensures that spatially varying water content (from
EPS modification of retention) does not create spurious advection in the
diffusion equation. See O2_state_variable_change.md for derivation.

The effective diffusion coefficient D_O^eff (from D_eff_oxygen) already
accounts for dual-phase transport (aqueous + gas) normalized by the
capacity (θ + K_H θ_a). Gas-phase diffusion (~10⁴× faster than aqueous)
dominates in unsaturated pores and is fully preserved.

The reaction sink S_O = -α_O · Resp / (θ + K_H θ_a) is divided by
capacity because each unit of respiration depletes a larger reservoir
when air content is high (K_H ≈ 30). This concentrates O₂ depletion
in wet regions near the POM, driving anoxic zone formation.

# Manuscript reference
Architecture §2: Time integration, Strang splitting
"""
function run_simulation(state::AggregateState, workspace::Workspace,
                       r_grid::Vector{Float64}, h::Real,
                       bio::BiologicalProperties, soil::SoilProperties,
                       T_func, ψ_func, O2_func,
                       t_start::Real, t_end::Real, dt_initial::Real;
                       t_delay::Real=0.0,
                       dt_min::Real=1e-4, dt_max::Real=0.1,
                       output_interval::Real=1.0,
                       output_times::Vector{Float64}=Float64[])
    # POM dissolution divides by state.P_0 (J_P: pom_factor = P / P_0).
    # AggregateState(n) initialises P_0 = 0.0, so a state built by hand rather
    # than via create_initial_state silently yields Inf and then NaN through
    # every pool. Fail loudly instead of producing garbage.
    if !(state.P_0 > 0.0)
        throw(ArgumentError(
            "AggregateState.P_0 must be > 0 (got $(state.P_0)). " *
            "create_initial_state sets it automatically; if you built the " *
            "state directly, set state.P_0 = state.P before integrating."))
    end

    # Initialize
    t = t_start
    dt = dt_initial

    # Output schedule
    if !isempty(output_times)
        # User-specified: sort, deduplicate, ensure t_start and t_end included, filter to range
        scheduled = sort(unique([t_start; Float64.(output_times); t_end]))
        filter!(t -> t_start <= t <= t_end, scheduled)
        use_scheduled = true
        sched_idx = 2  # 1 is t_start, saved immediately below
    else
        use_scheduled = false
        next_output = t_start + output_interval
    end

    # Save initial state
    saved_times = Float64[t_start]
    saved_states = [deepcopy(state)]

    # Diagnostics
    n_steps = 0
    n_rejected = 0

    # Main loop
    while t < t_end
        # Don't overshoot final time
        if t + dt > t_end
            dt = t_end - t
        end

        # Don't overshoot next output target
        if use_scheduled
            if sched_idx <= length(scheduled) && t + dt > scheduled[sched_idx]
                dt = scheduled[sched_idx] - t
            end
        else
            if t + dt > next_output
                dt = next_output - t
            end
        end

        # === Get environmental conditions ===
        T = T_func(t)
        ψ = ψ_func(t)
        O2_gas = O2_func(t)

        # === Update workspace (once per timestep) ===

        # Temperature cache
        update_temperature_cache!(workspace.f_T, T, bio, soil)

        # Water content (θ changes due to EPS/F_i evolution since last step)
        update_water_content!(workspace.θ, workspace.θ_a, ψ, state, soil)

        # Effective diffusion coefficients (uses updated θ)
        update_effective_diffusion!(workspace, state, soil, bio, workspace.f_T)

        # === O₂ boundary condition: gas-phase → aqueous ===
        # O2_func returns atmospheric gas-phase density [μg/mm³].
        # State variable is aqueous C_aq = C_gas / K_H.
        O2_amb = O2_gas / workspace.f_T.K_H_O

        # === Compute POM dissolution (once per timestep, Strang-consistent) ===
        B_0 = state.B[1]
        F_n_0 = state.F_n[1]
        θ_0 = workspace.θ[1]
        θ_a_0 = workspace.θ_a[1]
        # O₂ state variable is already aqueous concentration
        O_aq_0 = state.O[1]
        R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
        # POM activation delay: sigmoid switch centered at t_delay
        # When t_delay = 0: no effect
        # When t_delay > 0: smoothly transitions from 0 → 1 around t = t_delay
        if t_delay > 0.0
            R_P_max_T *= 1.0 / (1.0 + exp(-(t - t_delay) / (0.1 * t_delay)))
        end
        J_P_val = J_P(state.P, state.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                     bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)

        # Convert flux density to total rate [μg-C/day]
        # CRITICAL: Use actual POM radius from grid, not bio.r_0 default
        R_P_val = R_P(J_P_val, r_grid[1])

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
            n_rejected += 1
        end
        dt = dt_new

        # === Save output at target times ===
        if use_scheduled
            if sched_idx <= length(scheduled) && abs(t - scheduled[sched_idx]) < 1e-10
                push!(saved_times, t)
                push!(saved_states, deepcopy(state))
                sched_idx += 1
            end
        else
            if abs(t - next_output) < 1e-10 || t >= t_end
                push!(saved_times, t)
                push!(saved_states, deepcopy(state))
                next_output += output_interval
            end
        end
    end

    # Return results
    diagnostics = Dict(
        "n_steps" => n_steps,
        "n_rejected" => n_rejected,
        "final_time" => t
    )

    return (times=saved_times, states=saved_states, diagnostics=diagnostics)
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
