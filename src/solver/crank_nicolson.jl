# crank_nicolson.jl
# Crank-Nicolson diffusion half-step for one species on spherical domain

"""
    BoundaryCondition

Type of boundary condition at domain edges.

- `:neumann_zero`: Zero flux (∂u/∂r = 0), implemented with ghost node
- `:neumann_flux`: Specified flux J [μg/mm²/day], implemented with ghost node
- `:dirichlet`: Specified value
"""
@enum BoundaryCondition begin
    neumann_zero
    neumann_flux
    dirichlet
end

"""
    crank_nicolson_step!(u, D, dt_half, r_grid, h, lower, diag, upper, rhs,
                         bc_inner, bc_outer, flux_inner=0.0, value_outer=0.0)

Perform one Crank-Nicolson diffusion half-step for a single species.

Solves the spherical diffusion equation with spatially-varying diffusion coefficient D(r):
    ∂u/∂t = (1/r²) × ∂/∂r[r² × D × ∂u/∂r]

Uses Crank-Nicolson time integration (second-order accurate):
    u^{n+1} - u^n = (Δt/4) × [L[u^n] + L[u^{n+1}]]

where L is the spherical Laplacian operator discretized as:
    L_i = (1/r_i² h²) × [r_{i+1/2}² D_{i+1/2} (u_{i+1} - u_i)
                          - r_{i-1/2}² D_{i-1/2} (u_i - u_{i-1})]

# Arguments
- `u::Vector{Float64}`: State vector [n], updated in-place
- `D::Vector{Float64}`: Effective diffusion coefficient [mm²/day] at each node [n]
- `dt_half::Real`: Half time step [days] (usually Δt/2 for Strang splitting)
- `r_grid::Vector{Float64}`: Radial grid points [mm] [n]
- `h::Real`: Grid spacing [mm]
- `lower::Vector{Float64}`: Pre-allocated lower diagonal [n-1], modified in-place
- `diag::Vector{Float64}`: Pre-allocated main diagonal [n], modified in-place
- `upper::Vector{Float64}`: Pre-allocated upper diagonal [n-1], modified in-place
- `rhs::Vector{Float64}`: Pre-allocated RHS vector [n], modified in-place
- `bc_inner::BoundaryCondition`: Inner boundary condition type
- `bc_outer::BoundaryCondition`: Outer boundary condition type
- `flux_inner::Real=0.0`: Flux at inner boundary [μg/mm²/day] (if bc_inner == neumann_flux)
- `value_outer::Real=0.0`: Value at outer boundary (if bc_outer == dirichlet)

# Boundary conditions
Inner boundary (i=0):
- `:neumann_zero`: ∂u/∂r = 0 → ghost node u[-1] = u[0]
- `:neumann_flux`: −D₀ ∂u/∂r = flux_inner → ghost node u[-1] = u[0] + 2h×flux_inner/D₀

Outer boundary (i=n-1):
- `:neumann_zero`: ∂u/∂r = 0 → ghost node u[n] = u[n-1]
- `:dirichlet`: u[n-1] = value_outer

# Notes
- Updates `u` in-place (zero allocation)
- Uses pre-allocated workspace arrays (lower, diag, upper, rhs)
- Face-averaged diffusion coefficients: D_{i+1/2} = (D_i + D_{i+1})/2
- Solves via Thomas algorithm in O(n) time

# Manuscript reference
Architecture §3.2-3.4: Spherical Laplacian, Crank-Nicolson, boundary conditions
"""
function crank_nicolson_step!(u::Vector{Float64}, D::Vector{Float64},
                              dt_half::Real, r_grid::Vector{Float64}, h::Real,
                              lower::Vector{Float64}, diag::Vector{Float64},
                              upper::Vector{Float64}, rhs::Vector{Float64},
                              bc_inner::BoundaryCondition, bc_outer::BoundaryCondition,
                              flux_inner::Real=0.0, value_outer::Real=0.0)
    n = length(u)

    # === Step 1: Compute explicit Laplacian L[u^n] and build RHS ===
    @inbounds for i in 1:n
        r_i = r_grid[i]
        r_i_sq = r_i * r_i
        inv_r_sq_h_sq = 1.0 / (r_i_sq * h * h)

        if i == 1
            # Inner boundary
            r_half_plus = (r_grid[1] + r_grid[2]) / 2.0
            D_half_plus = (D[1] + D[2]) / 2.0

            if bc_inner == neumann_zero
                # Ghost node: u[0] = u[1] (zero flux)
                # L[1] = (1/r₁² h²) × [r_{3/2}² D_{3/2} (u[2] - u[1]) - r_{1/2}² D_{1/2} (u[1] - u[0])]
                #      = (1/r₁² h²) × [r_{3/2}² D_{3/2} (u[2] - u[1]) - 0]  (u[1] = u[0])
                L_explicit = inv_r_sq_h_sq * r_half_plus * r_half_plus * D_half_plus * (u[2] - u[1])
            elseif bc_inner == neumann_flux
                # Ghost node: u[0] = u[1] + h×J/D[1] where J = flux_inner
                # L[1] = (1/r₁² h²) × [r_{3/2}² D_{3/2} (u[2] - u[1]) - r₁² D₁ (u[1] - u[0])]
                # CRITICAL: Face is AT r_0 = r_grid[1], ghost at r_0-h, distance = h
                r_half_minus = r_grid[1]  # = r_0 (face location)
                D_half_minus = D[1]       # At boundary, use D[1]
                ghost_diff = u[1] - (u[1] + h*flux_inner/D[1])  # u[1] - u[0]
                ghost_diff = -h*flux_inner/D[1]

                L_explicit = inv_r_sq_h_sq * (r_half_plus * r_half_plus * D_half_plus * (u[2] - u[1]) -
                                              r_half_minus * r_half_minus * D_half_minus * ghost_diff)
            else
                error("Dirichlet BC not supported at inner boundary")
            end
        elseif i == n
            # Outer boundary
            if bc_outer == neumann_zero
                # Ghost node: u[n] = u[n-1] (zero flux)
                # L[n] = (1/r_n² h²) × [0 - r_{n-1/2}² D_{n-1/2} (u[n-1] - u[n-2])]
                r_half_minus = (r_grid[n-1] + r_grid[n]) / 2.0
                D_half_minus = (D[n-1] + D[n]) / 2.0
                L_explicit = -inv_r_sq_h_sq * r_half_minus * r_half_minus * D_half_minus * (u[n] - u[n-1])
            elseif bc_outer == dirichlet
                # For Dirichlet, the last equation is simply u[n] = value_outer
                # So L_explicit is not used, we'll handle this specially below
                L_explicit = 0.0
            else
                error("Neumann flux BC not supported at outer boundary")
            end
        else
            # Interior nodes
            r_half_plus = (r_grid[i] + r_grid[i+1]) / 2.0
            r_half_minus = (r_grid[i-1] + r_grid[i]) / 2.0
            D_half_plus = (D[i] + D[i+1]) / 2.0
            D_half_minus = (D[i-1] + D[i]) / 2.0

            L_explicit = inv_r_sq_h_sq * (r_half_plus * r_half_plus * D_half_plus * (u[i+1] - u[i]) -
                                          r_half_minus * r_half_minus * D_half_minus * (u[i] - u[i-1]))
        end

        # RHS: u^n + (dt/4) × L[u^n]
        # Note: dt_half = Δt/2, so dt/4 = dt_half/2
        rhs[i] = u[i] + 0.5 * dt_half * L_explicit
    end

    # === Step 2: Build implicit tridiagonal matrix ===
    # LHS: u^{n+1} - (dt/4) × L[u^{n+1}] = rhs
    # Rearranges to: [I - (dt/4)×L] × u^{n+1} = rhs
    # This gives tridiagonal system: a_i × u_{i-1} + b_i × u_i + c_i × u_{i+1} = rhs_i

    @inbounds for i in 1:n
        r_i = r_grid[i]
        r_i_sq = r_i * r_i
        inv_r_sq_h_sq = 1.0 / (r_i_sq * h * h)
        # Note: dt_half = Δt/2, so dt/4 = dt_half/2
        coeff = -0.5 * dt_half * inv_r_sq_h_sq

        if i == 1
            # Inner boundary
            r_half_plus = (r_grid[1] + r_grid[2]) / 2.0
            D_half_plus = (D[1] + D[2]) / 2.0

            if bc_inner == neumann_zero
                # Implicit stencil with ghost node u[0] = u[1]:
                # a[1]×u[0] + b[1]×u[1] + c[1]×u[2]
                # = a[1]×u[1] + b[1]×u[1] + c[1]×u[2]  (substitute u[0] = u[1])
                # = (a[1] + b[1])×u[1] + c[1]×u[2]
                # But a[1] would be for the r_{1/2}² D_{1/2} (u[1] - u[0]) term, which is zero

                c_1 = coeff * r_half_plus * r_half_plus * D_half_plus  # Coefficient of u[2]
                a_1 = 0.0  # No lower term (ghost = current)
                b_1 = 1.0 - c_1 - a_1

                diag[1] = b_1
                upper[1] = c_1
                # lower[1] not used (first equation has no lower)
                # rhs[1] already set above

            elseif bc_inner == neumann_flux
                # Ghost node: u[0] = u[1] + h×J/D[1]
                # Standard stencil: a[1]×u[0] + b[1]×u[1] + c[1]×u[2] = rhs
                # Substitute: a[1]×(u[1] + h×J/D[1]) + b[1]×u[1] + c[1]×u[2] = rhs
                # Rearrange: (a[1] + b[1])×u[1] + c[1]×u[2] = rhs - a[1]×h×J/D[1]
                # Note: a[1] < 0, so subtracting it adds a positive correction

                # CRITICAL: Face is AT r_0 = r_grid[1], ghost at r_0-h, distance = h
                # Domain starts at r_0, does not extend below POM
                r_half_minus = r_grid[1]  # = r_0 (face location)
                D_half_minus = D[1]

                a_1 = coeff * r_half_minus * r_half_minus * D_half_minus  # Coefficient of u[0] (negative)
                c_1 = coeff * r_half_plus * r_half_plus * D_half_plus      # Coefficient of u[2]
                b_1 = 1.0 - a_1 - c_1                                      # Coefficient of u[1]

                diag[1] = a_1 + b_1  # Combine diagonal with ghost contribution
                upper[1] = c_1
                rhs[1] -= a_1 * h * flux_inner / D[1]  # Subtract (a_1 < 0, so this adds)
            end

        elseif i == n
            # Outer boundary
            if bc_outer == neumann_zero
                # Ghost node: u[n] = u[n-1]
                # Standard stencil: a[n]×u[n-1] + b[n]×u[n] + c[n]×u[n+1]
                # Substitute: a[n]×u[n-1] + b[n]×u[n] + c[n]×u[n]  (u[n+1] = u[n])
                #           = a[n]×u[n-1] + (b[n] + c[n])×u[n]

                r_half_minus = (r_grid[n-1] + r_grid[n]) / 2.0
                D_half_minus = (D[n-1] + D[n]) / 2.0

                a_n = coeff * r_half_minus * r_half_minus * D_half_minus  # Coefficient of u[n-1]
                c_n = 0.0  # No upper term (ghost = current)
                b_n = 1.0 - a_n - c_n

                lower[n-1] = a_n
                diag[n] = b_n + c_n  # Combine diagonal with ghost contribution
                # upper[n-1] not used (last equation has no upper)
                # rhs[n] already set above

            elseif bc_outer == dirichlet
                # Dirichlet BC: simply set u[n] = value_outer
                # Row n: 0×u[n-1] + 1×u[n] + 0×u[n+1] = value_outer
                lower[n-1] = 0.0
                diag[n] = 1.0
                # upper[n-1] not used
                rhs[n] = value_outer
            end

        else
            # Interior nodes
            r_half_plus = (r_grid[i] + r_grid[i+1]) / 2.0
            r_half_minus = (r_grid[i-1] + r_grid[i]) / 2.0
            D_half_plus = (D[i] + D[i+1]) / 2.0
            D_half_minus = (D[i-1] + D[i]) / 2.0

            a_i = coeff * r_half_minus * r_half_minus * D_half_minus  # Coefficient of u[i-1]
            c_i = coeff * r_half_plus * r_half_plus * D_half_plus      # Coefficient of u[i+1]
            b_i = 1.0 - a_i - c_i                                       # Coefficient of u[i]

            lower[i-1] = a_i
            diag[i] = b_i
            upper[i] = c_i
            # rhs[i] already set above
        end
    end

    # === Step 3: Solve tridiagonal system (u is overwritten with solution) ===
    thomas!(lower, diag, upper, rhs)

    # Copy solution back to u
    @inbounds for i in 1:n
        u[i] = rhs[i]
    end

    nothing
end
