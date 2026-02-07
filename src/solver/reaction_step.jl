# reaction_step.jl
# Reaction step for all nodes + POM + CO₂ accumulation

"""
    reaction_step!(state::AggregateState, workspace::Workspace, dt::Real,
                   r_grid::Vector{Float64}, h::Real,
                   bio::BiologicalProperties, soil::SoilProperties,
                   ψ::Real, R_P_dissolution::Real) -> Float64

Perform one reaction step for the full aggregate.

Advances all 8 spatial state variables at each node, plus POM scalar and CO₂ cumulative.
Uses Forward Euler time integration.

# Arguments
- `state::AggregateState`: State variables, updated in-place
- `workspace::Workspace`: Pre-allocated workspace (for θ, θ_a, temp_cache)
- `dt::Real`: Time step [days]
- `r_grid::Vector{Float64}`: Radial grid points [mm] [n]
- `h::Real`: Grid spacing [mm]
- `bio::BiologicalProperties`: Biological parameters
- `soil::SoilProperties`: Soil parameters
- `ψ::Real`: Water potential [kPa] (uniform in space for now)
- `R_P_dissolution::Real`: POM dissolution rate [μg-C/day] (pre-computed for Strang consistency)

# Returns
- `max_rel_change::Float64`: Maximum relative change max(|S × Δt / u|) over all nodes/species

# Notes
- Advances all nodes independently (no spatial coupling in reactions)
- POM dissolution depends on state at first grid node (r = r_POM)
- CO₂ accumulation sums respiration over all nodes × volume elements
- Enforces non-negativity after update
- Zero allocation (all computations in-place)
- Returns max_rel_change for adaptive timestep control (eliminates redundant computation)

# Integration method
Forward Euler: u^{n+1} = u^n + Δt × S(u^n)

# Manuscript reference
Architecture §4: Reaction step, POM dissolution, CO₂ accumulation
"""
function reaction_step!(state::AggregateState, workspace::Workspace, dt::Real,
                       r_grid::Vector{Float64}, h::Real,
                       bio::BiologicalProperties, soil::SoilProperties,
                       ψ::Real, R_P_dissolution::Real)
    n = length(state.C)
    max_rel_change = 0.0
    threshold = 1e-6  # Avoid division by tiny values

    # === Loop over all grid nodes ===
    @inbounds for i in 1:n
        # Get state at this node
        C_i = state.C[i]
        B_i = state.B[i]
        F_n_i = state.F_n[i]
        F_m_i = state.F_m[i]
        F_i_i = state.F_i[i]
        E_i = state.E[i]
        M_i = state.M[i]
        O_i = state.O[i]
        θ_i = workspace.θ[i]
        θ_a_i = workspace.θ_a[i]

        # Compute source terms
        sources = compute_source_terms(C_i, B_i, F_n_i, F_m_i, F_i_i, E_i, M_i, O_i,
                                      θ_i, θ_a_i, ψ, bio, soil, workspace.f_T)

        # Compute relative changes for adaptive timestep (before applying updates)
        rel_C = abs(sources.S_C * dt / max(C_i, threshold))
        rel_B = abs(sources.S_B * dt / max(B_i, threshold))
        rel_Fn = abs(sources.S_Fn * dt / max(F_n_i, threshold))
        rel_Fm = abs(sources.S_Fm * dt / max(F_m_i, threshold))
        rel_Fi = abs(sources.S_Fi * dt / max(F_i_i, threshold))
        rel_E = abs(sources.S_E * dt / max(E_i, threshold))
        rel_M = abs(sources.S_M * dt / max(M_i, threshold))
        rel_O = abs(sources.S_O * dt / max(O_i, threshold))

        # Track maximum over all species at this node
        node_max = max(rel_C, rel_B, rel_Fn, rel_Fm, rel_Fi, rel_E, rel_M, rel_O)
        max_rel_change = max(max_rel_change, node_max)

        # Forward Euler update
        state.C[i] += dt * sources.S_C
        state.B[i] += dt * sources.S_B
        state.F_n[i] += dt * sources.S_Fn
        state.F_m[i] += dt * sources.S_Fm
        state.F_i[i] += dt * sources.S_Fi
        state.E[i] += dt * sources.S_E
        state.M[i] += dt * sources.S_M
        state.O[i] += dt * sources.S_O

        # Conservation weight for this node (matches spherical Laplacian stencil)
        # W_i = 4πr_i²h ensures W_i/(r_i²h²) = 4π/h is constant → exact telescoping
        volume_i = 4.0 * π * r_grid[i]^2 * h

        # Enforce non-negativity — track clipped carbon for exact conservation
        # Any carbon clipped to zero is redirected to CO₂ (lost biomass/substrate)
        clip_carbon = 0.0

        # Carbon species (NOT oxygen)
        if state.C[i] < 0.0
            clip_carbon += abs(state.C[i]) * volume_i
            state.C[i] = 0.0
        end
        if state.B[i] < 0.0
            clip_carbon += abs(state.B[i]) * volume_i
            state.B[i] = 0.0
        end
        if state.F_n[i] < 0.0
            clip_carbon += abs(state.F_n[i]) * volume_i
            state.F_n[i] = 0.0
        end
        if state.F_m[i] < 0.0
            clip_carbon += abs(state.F_m[i]) * volume_i
            state.F_m[i] = 0.0
        end
        if state.F_i[i] < 0.0
            clip_carbon += abs(state.F_i[i]) * volume_i
            state.F_i[i] = 0.0
        end
        if state.E[i] < 0.0
            clip_carbon += abs(state.E[i]) * volume_i
            state.E[i] = 0.0
        end
        if state.M[i] < 0.0
            clip_carbon += abs(state.M[i]) * volume_i
            state.M[i] = 0.0
        end

        # Oxygen (NOT carbon — just clip)
        state.O[i] = max(0.0, state.O[i])

        # Accumulate CO₂ from respiration + clipped carbon
        state.CO2_cumulative += dt * sources.Resp_total * volume_i + clip_carbon
    end

    # === POM dissolution (using pre-computed rate from caller) ===
    # R_P_dissolution is computed in timestepper ONCE per step (Strang-consistent)
    # This ensures the flux entering C (via diffusion BC) matches the amount leaving P
    state.P -= dt * R_P_dissolution

    # Enforce non-negativity — track clipped carbon
    if state.P < 0.0
        state.CO2_cumulative += abs(state.P)  # POM is scalar, not per-volume
        state.P = 0.0
    end

    return max_rel_change
end
