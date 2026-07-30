# diffusion_step.jl
# Diffusion half-step for all 5 diffusing species

"""
    diffusion_step!(state, workspace, dt_half, r_grid, h, J_P_flux, O_amb)

Perform one diffusion half-step for all 5 diffusing species.

Calls Crank-Nicolson solver for each of the 5 species:
- C (dissolved organic carbon): Neumann flux at inner (POM dissolution), zero flux at outer
- B (bacteria): Zero flux at both boundaries
- F_n (non-insulated fungi): Zero flux at both boundaries
- F_m (mobile fungi): Zero flux at both boundaries
- O (oxygen): Zero flux at inner, Dirichlet at outer (ambient O₂)

# Arguments
- `state::AggregateState`: State variables, updated in-place
- `workspace::Workspace`: Pre-allocated workspace arrays
- `dt_half::Real`: Half time step [days] (typically Δt/2 for Strang splitting)
- `r_grid::Vector{Float64}`: Radial grid points [mm] [n]
- `h::Real`: Grid spacing [mm]
- `J_P_flux::Real`: POM dissolution flux density [μg-C/mm²/day] (for C inner BC)
- `O_amb::Real`: Ambient oxygen concentration [μg/mm³] (for O outer BC)

# Notes
- Uses pre-allocated workspace arrays (zero allocation)
- Updates state variables in-place
- 5 tridiagonal solves, each O(n)
- Total complexity: O(5n) per half-step

# Manuscript reference
Architecture §2, §3.4: Strang splitting, boundary conditions
"""
function diffusion_step!(state::AggregateState, workspace::Workspace,
                        dt_half::Real, r_grid::Vector{Float64}, h::Real,
                        J_P_flux::Real, O_amb::Real)
    # === 1. C (dissolved organic carbon) ===
    # Inner: Neumann flux from POM dissolution (J_P_flux)
    # Outer: Zero flux
    crank_nicolson_step!(state.C, workspace.D_C, dt_half, r_grid, h,
                        workspace.lower, workspace.diag, workspace.upper, workspace.rhs,
                        neumann_flux, neumann_zero, J_P_flux, 0.0)

    # === 2. B (bacteria) ===
    # Inner: Zero flux
    # Outer: Zero flux
    crank_nicolson_step!(state.B, workspace.D_B, dt_half, r_grid, h,
                        workspace.lower, workspace.diag, workspace.upper, workspace.rhs,
                        neumann_zero, neumann_zero, 0.0, 0.0)

    # === 3. F_n (non-insulated fungi) ===
    # Inner: Zero flux
    # Outer: Zero flux
    crank_nicolson_step!(state.F_n, workspace.D_Fn, dt_half, r_grid, h,
                        workspace.lower, workspace.diag, workspace.upper, workspace.rhs,
                        neumann_zero, neumann_zero, 0.0, 0.0)

    # === 4. F_m (mobile fungi) ===
    # Inner: Zero flux
    # Outer: Zero flux
    # Note: D_Fm is now spatially varying (network-dependent, see effective_diffusion.jl)
    crank_nicolson_step!(state.F_m, workspace.D_Fm, dt_half, r_grid, h,
                        workspace.lower, workspace.diag, workspace.upper, workspace.rhs,
                        neumann_zero, neumann_zero, 0.0, 0.0)

    # === 5. O (oxygen) ===
    # Inner: Zero flux
    # Outer: Dirichlet (ambient concentration)
    # NOTE: Backward Euler (theta=1.0) is required here.
    # O₂ diffuses ~100-1000× faster than other species (gas-phase transport),
    # giving Fo = D_eff × dt_half / h² ≈ 42,600 >> 1 for typical parameters.
    # CN (theta=0.5) with Fo >> 1 acts as a "reflection about steady-state":
    #   CN(u)[1] ≈ 2×O2_amb - u[1]
    # The Strang split D(dt/2)→R(dt)→D(dt/2) then INVERTS the reaction sign:
    #   net effect = u_old + δ (grows) instead of u_old - δ (depletes).
    # Backward Euler (theta=1.0) with Fo >> 1 projects to steady state, which
    # is physically correct since O₂ equilibrates in ~minutes, far faster than dt.
    crank_nicolson_step!(state.O, workspace.D_O, dt_half, r_grid, h,
                        workspace.lower, workspace.diag, workspace.upper, workspace.rhs,
                        neumann_zero, dirichlet, 0.0, O_amb; theta=1.0)

    nothing
end
