"""
    tridiagonal.jl

Thomas algorithm for solving tridiagonal systems.

Solves Ax = d where A is tridiagonal:
    b[1]   c[1]    0     ...     0
    a[2]   b[2]   c[2]   ...     0
     0     a[3]   b[3]   ...     0
    ...    ...    ...    ...    ...
     0      0      0    a[n]   b[n]

The algorithm is O(n) and in-place, overwriting the RHS vector with the solution.
"""

"""
    thomas!(lower, diag, upper, rhs)

Solve tridiagonal system Ax = rhs using Thomas algorithm (in-place).

# Arguments
- `lower`: Lower diagonal [n-1] (a[2], a[3], ..., a[n])
- `diag`: Main diagonal [n] (b[1], b[2], ..., b[n])
- `upper`: Upper diagonal [n-1] (c[1], c[2], ..., c[n-1])
- `rhs`: Right-hand side [n], **overwritten with solution** on exit

# Returns
- Nothing (modifies `rhs` in-place)

# Algorithm
Two-pass Thomas algorithm:
1. Forward elimination: eliminate lower diagonal
2. Backward substitution: solve for x from bottom up

# Notes
- **Mutates `diag`, `upper`, and `rhs`** (working storage)
- Zero allocations
- Assumes system is well-conditioned (diagonally dominant)
- For repeated solves with same matrix, use `thomas_factorize!` + `thomas_solve!`

# Examples

```julia
# Solve simple tridiagonal system
lower = [1.0, 1.0, 1.0]
diag  = [2.0, 2.0, 2.0, 2.0]
upper = [1.0, 1.0, 1.0]
rhs   = [1.0, 0.0, 0.0, 1.0]

thomas!(lower, diag, upper, rhs)
# rhs now contains solution
```

# Performance
- Time: O(n)
- Space: O(1) - no allocations
- For n=1500: ~3 μs
"""
function thomas!(lower::AbstractVector, diag::AbstractVector,
                 upper::AbstractVector, rhs::AbstractVector)
    n = length(diag)

    # Forward elimination
    @inbounds for i in 2:n
        factor = lower[i-1] / diag[i-1]
        diag[i] -= factor * upper[i-1]
        rhs[i] -= factor * rhs[i-1]
    end

    # Backward substitution
    @inbounds rhs[n] /= diag[n]
    @inbounds for i in (n-1):-1:1
        rhs[i] = (rhs[i] - upper[i] * rhs[i+1]) / diag[i]
    end

    nothing
end

"""
    thomas(lower, diag, upper, rhs)

Solve tridiagonal system Ax = rhs (non-mutating version).

Returns the solution without modifying input arrays.
Allocates temporary storage - use `thomas!` for performance-critical code.

# Arguments
- `lower`: Lower diagonal [n-1]
- `diag`: Main diagonal [n]
- `upper`: Upper diagonal [n-1]
- `rhs`: Right-hand side [n]

# Returns
- Solution vector x [n]

# Examples

```julia
x = thomas(lower, diag, upper, rhs)
# Original arrays unchanged
```
"""
function thomas(lower::AbstractVector, diag::AbstractVector,
                upper::AbstractVector, rhs::AbstractVector)
    # Copy arrays to avoid mutation
    lower_copy = copy(lower)
    diag_copy = copy(diag)
    upper_copy = copy(upper)
    rhs_copy = copy(rhs)

    thomas!(lower_copy, diag_copy, upper_copy, rhs_copy)

    return rhs_copy
end

"""
    thomas_factorize!(lower, diag, upper)

Factor tridiagonal matrix in-place for repeated solves with same matrix.

Performs forward elimination step, modifying arrays to store the factorization.
Use with `thomas_solve!` for repeated RHS.

# Arguments
- `lower`: Lower diagonal [n-1], **overwritten with multipliers**
- `diag`: Main diagonal [n], **overwritten with factored diagonal**
- `upper`: Upper diagonal [n-1], **unchanged (original values preserved)**

# Returns
- Nothing (modifies `lower` and `diag` in-place)

# Usage

```julia
# Factor once
thomas_factorize!(lower, diag, upper)

# Solve multiple times with different RHS
for rhs in rhs_collection
    thomas_solve!(lower, diag, upper, rhs)  # uses factored matrix
end
```

# Notes
- After factorization:
  - `lower[i-1]` contains multiplier m[i] = original_lower[i-1] / diag[i-1]
  - `diag[i]` contains modified diagonal
  - `upper[i]` unchanged (original values)
- Useful when solving multiple systems with same matrix (e.g., different timesteps)
"""
function thomas_factorize!(lower::AbstractVector, diag::AbstractVector,
                           upper::AbstractVector)
    n = length(diag)

    # Forward elimination - store multipliers in lower
    @inbounds for i in 2:n
        lower[i-1] /= diag[i-1]  # Store multiplier: m[i] = a[i] / b[i-1]
        diag[i] -= lower[i-1] * upper[i-1]  # Update diagonal: b[i] = b[i] - m[i]*c[i-1]
    end

    nothing
end

"""
    thomas_solve!(lower_factored, diag_factored, upper_original, rhs)

Solve factorized tridiagonal system (in-place).

Use after `thomas_factorize!` to solve with factored matrix.

# Arguments
- `lower_factored`: Multipliers from `thomas_factorize!` [n-1]
- `diag_factored`: Factored diagonal from `thomas_factorize!` [n]
- `upper_original`: Original upper diagonal (unchanged by factorization) [n-1]
- `rhs`: Right-hand side [n], **overwritten with solution**

# Returns
- Nothing (modifies `rhs` in-place)

# Examples

```julia
# Factor once
thomas_factorize!(lower, diag, upper)

# Solve multiple RHS
rhs1 = [1.0, 0.0, 0.0, 1.0]
thomas_solve!(lower, diag, upper, rhs1)  # rhs1 now contains solution

rhs2 = [0.0, 1.0, 1.0, 0.0]
thomas_solve!(lower, diag, upper, rhs2)  # rhs2 now contains solution
```
"""
function thomas_solve!(lower_factored::AbstractVector, diag_factored::AbstractVector,
                       upper_original::AbstractVector, rhs::AbstractVector)
    n = length(diag_factored)

    # Forward substitution (apply multipliers to RHS)
    @inbounds for i in 2:n
        rhs[i] -= lower_factored[i-1] * rhs[i-1]
    end

    # Backward substitution
    @inbounds rhs[n] /= diag_factored[n]
    @inbounds for i in (n-1):-1:1
        rhs[i] = (rhs[i] - upper_original[i] * rhs[i+1]) / diag_factored[i]
    end

    nothing
end

"""
    is_diagonally_dominant(lower, diag, upper)

Check if tridiagonal matrix is diagonally dominant.

A tridiagonal matrix is diagonally dominant if:
    |b[i]| ≥ |a[i]| + |c[i]| for all i

# Arguments
- `lower`: Lower diagonal [n-1]
- `diag`: Main diagonal [n]
- `upper`: Upper diagonal [n-1]

# Returns
- `true` if diagonally dominant, `false` otherwise

# Notes
- Diagonal dominance guarantees Thomas algorithm stability
- Crank-Nicolson diffusion matrices are diagonally dominant
"""
function is_diagonally_dominant(lower::AbstractVector, diag::AbstractVector,
                                upper::AbstractVector)
    n = length(diag)

    # First row: |b[1]| ≥ |c[1]|
    abs(diag[1]) >= abs(upper[1]) || return false

    # Interior rows: |b[i]| ≥ |a[i]| + |c[i]|
    @inbounds for i in 2:(n-1)
        abs(diag[i]) >= abs(lower[i-1]) + abs(upper[i]) || return false
    end

    # Last row: |b[n]| ≥ |a[n]|
    abs(diag[n]) >= abs(lower[n-1]) || return false

    return true
end
