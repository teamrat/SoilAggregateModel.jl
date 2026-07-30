"""
    math_utils.jl

Shared mathematical utility functions.
"""

"""
    sigmoid_threshold(x::Real, x_min::Real, steepness::Real)

Smooth 0 → 1 switch centred on `x_min`:

    1 / (1 + exp(-β·(x - x_min))),     β = steepness / x_min

`steepness` is dimensionless (it is β·x_min), so the width of the transition is
a fixed fraction of `x_min` whatever the units of `x`. Returns exactly 0.5 at
`x = x_min`, → 0 for `x ≪ x_min`, → 1 for `x ≫ x_min`.

This form is the numerically stable one: the algebraically equivalent
`exp(βx) / [exp(βx) + exp(βx_min)]` written in the manuscript overflows for
`x ≫ x_min`, while this does not.

The single definition behind `h_B`, `h_E`, `h_Fi` (steepness `SIGMOID_STEEPNESS`)
and the POM activation delay (`POM_DELAY_STEEPNESS`). Those four had the same
three lines written out four times.

# Arguments
- `x`: the quantity being switched on
- `x_min`: switch centre, same units as `x`; must be > 0
- `steepness`: β·x_min [-]

# Returns
- Factor in (0, 1)
"""
function sigmoid_threshold(x::Real, x_min::Real, steepness::Real)
    β = steepness / x_min
    1.0 / (1.0 + exp(-β * (x - x_min)))
end

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
