"""
    math_utils.jl

Shared mathematical utility functions.
"""

"""
    softplus(x::Real, ε::Real)

Smooth approximation to max(0, x) via softplus: ε·ln(1 + exp(x/ε)).

Numerically stable implementation:
    x > 0:  x + ε·ln(1 + exp(-x/ε))
    x ≤ 0:  ε·ln(1 + exp(x/ε))

# Arguments
- `x`: Input value (can be any real)
- `ε`: Smoothing width [same units as x]

# Returns
- Smoothed value (always ≥ 0)

# Notes
- As ε → 0: softplus(x, ε) → max(0, x)
- For |x| >> ε: softplus(x, ε) ≈ max(0, x) to within ε
- Numerically stable form avoids overflow for large |x|
- Used for smooth transitions in yield, allocation, and MAOC sorption
"""
function softplus(x::Real, ε::Real)
    if x > 0.0
        # Stable for large positive x
        x + ε * log(1.0 + exp(-x / ε))
    else
        # Stable for large negative x
        ε * log(1.0 + exp(x / ε))
    end
end
