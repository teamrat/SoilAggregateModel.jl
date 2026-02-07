# Test Requirements for crank_nicolson.jl

**Source**: Architecture review by Claude Opus 4.6

The Thomas algorithm in `tridiagonal.jl` is **correct**. However, the current tests only validate **Cartesian geometry with Dirichlet BCs**. The actual aggregate model uses **spherical geometry with Neumann BCs**, which requires additional tests.

---

## Tests to Add in test_crank_nicolson.jl

### 1. Spherical Steady-State Test

**Problem**: Solve ∇²u = 0 in spherical shell [r₀, r_max] with boundary conditions:
- u(r₀) = 1 (Dirichlet at inner boundary)
- u(r_max) = 0 (Dirichlet at outer boundary)

**Analytical solution**:
```
u(r) = r₀(r_max - r) / [r(r_max - r₀)]
```

**Purpose**: Validates the spherical Laplacian discretization with r-dependent coefficients:
```
L_i = (1/r_i² h²) × [r_{i+1/2}² D_{i+1/2} (u_{i+1} - u_i)
                     - r_{i-1/2}² D_{i-1/2} (u_i - u_{i-1})]
```

**Acceptance**: Solution matches analytical to rtol < 1e-8

---

### 2. Neumann Boundary Condition Test

**Problem**: Solve diffusion in spherical shell with:
- Inner boundary: Constant flux J [μg/mm²/day] (Neumann, ghost node)
- Outer boundary: Zero flux (Neumann, ghost node)
- Steady state

**Analytical solution**:
```
u(r) = u₀ + (J × r₀²/D) × (1/r - 1/r_max)
```

**Ghost node implementation**:
- Inner (flux BC): u[0] = u[1] - (2h × J/D)
- Outer (zero-flux): u[n+1] = u[n]

These modify the first and last rows of the tridiagonal system.

**Purpose**: Validates proper Neumann BC implementation for the POM flux BC (inner) and zero-flux BCs (outer) used throughout the model.

**Acceptance**:
- Solution matches analytical to rtol < 1e-8
- Mass conservation to rtol < 1e-10

---

### 3. Mass Conservation with Zero-Flux BCs

**Problem**: Time-dependent diffusion with zero-flux Neumann BCs at both boundaries.

**Expected**: Total mass conserved to machine precision (rtol < 1e-10).

**Current status**: The Crank-Nicolson test in `test_tridiagonal.jl` shows ~0.034% error because it uses "implicit" Neumann BCs (no boundary modification) rather than proper ghost nodes.

**Purpose**: Confirms ghost-node implementation preserves conservative form.

---

## Why This Matters

1. **Spherical geometry**: The model uses 1/r² × d/dr(r² × ...), NOT d²/dx². Non-constant coefficients change the tridiagonal structure.

2. **POM flux BC**: The inner boundary has a source flux from POM dissolution, implemented as a Neumann BC. If this is wrong, carbon enters the system incorrectly.

3. **Mass conservation**: The architecture targets machine-precision carbon balance. Current ~0.03% error would accumulate to 0.3% over 10 years, violating the conservation requirement.

4. **Anoxia prediction**: Zero-flux at outer boundaries is critical for predicting O₂ depletion. If flux leaks out, the model will underpredict anoxic zone extent.

---

## Implementation Notes

**Ghost node for Neumann BC at i=1** (inner boundary with flux J):
```julia
# Standard stencil: a[1]×u[0] + b[1]×u[1] + c[1]×u[2] = d[1]
# Ghost node: u[0] = u[1] - 2h×J/D
# Substitute into first row:
#   a[1]×(u[1] - 2h×J/D) + b[1]×u[1] + c[1]×u[2] = d[1]
#   (a[1] + b[1])×u[1] + c[1]×u[2] = d[1] + a[1]×2h×J/D
# Modified coefficients:
lower[1] = (not used, first row has no lower)
diag[1] = b[1] + a[1]  # Combine diagonal and ghost contribution
upper[1] = c[1]         # Unchanged
rhs[1] = d[1] + a[1]×2h×J/D
```

**Ghost node for zero-flux at i=n** (outer boundary):
```julia
# Standard stencil: a[n]×u[n-1] + b[n]×u[n] + c[n]×u[n+1] = d[n]
# Ghost node: u[n+1] = u[n] (zero flux)
# Substitute:
#   a[n]×u[n-1] + b[n]×u[n] + c[n]×u[n] = d[n]
#   a[n]×u[n-1] + (b[n] + c[n])×u[n] = d[n]
# Modified coefficients:
lower[n-1] = a[n]       # Unchanged
diag[n] = b[n] + c[n]   # Combine diagonal and ghost contribution
upper[n-1] = (not used, last row has no upper)
rhs[n] = d[n]           # Unchanged for zero flux
```

---

## Current Test Gap

`test_tridiagonal.jl` line 170-172 uses loose tolerance (rtol=0.1) because the BC implementation is incomplete. This is **documented but not fixed**. The proper fix goes in `crank_nicolson.jl` when that module is built.
