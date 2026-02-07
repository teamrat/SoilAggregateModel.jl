"""
    environment.jl

Environmental drivers (temperature, matric potential, ambient O₂).

Parametric types ensure Julia specializes: constant functions are as fast as literal constants.
"""

#═══════════════════════════════════════════════════════════════════════════════
# Environmental Drivers
#═══════════════════════════════════════════════════════════════════════════════

"""
    EnvironmentalDrivers{FT, Fψ, FO}

Externally imposed environmental conditions as functions of time.

Type parameters ensure Julia specializes on function types:
- Constant values are as fast as literals
- Time-varying functions are efficiently dispatched
- Interpolated data tables work seamlessly

Fields:
- `T`: Temperature function T(t) → K
- `ψ`: Matric potential function ψ(t) → kPa
- `O2`: Ambient O₂ concentration function O2(t) → μg/mm³

# Examples

```julia
# All constant
env = EnvironmentalDrivers(293.15, -33.0, 0.27)

# Temperature varies sinusoidally, others constant
env = EnvironmentalDrivers(
    t -> 293.15 + 5.0 * sin(2π * t / 365),
    -33.0,
    0.27
)

# From data using Interpolations.jl
using Interpolations
T_interp = linear_interpolation(t_data, T_data)
ψ_interp = linear_interpolation(t_data, ψ_data)
env = EnvironmentalDrivers(T_interp, ψ_interp, 0.27)
```

Use the outer constructor (defined below) to automatically wrap constants as functions.
"""
struct EnvironmentalDrivers{FT, Fψ, FO}
    T::FT       # T(t) → temperature [K]
    ψ::Fψ       # ψ(t) → matric potential [kPa]
    O2::FO      # O2(t) → boundary O₂ concentration [μg/mm³]

    # Inner constructor: only accepts function types (suppresses auto-generated constructor)
    function EnvironmentalDrivers{FT, Fψ, FO}(T::FT, ψ::Fψ, O2::FO) where {FT, Fψ, FO}
        new{FT, Fψ, FO}(T, ψ, O2)
    end
end

"""
    EnvironmentalDrivers(T, ψ, O2)

Construct EnvironmentalDrivers from constants or functions.

Scalar arguments are automatically wrapped as constant functions `t -> value`.
Function arguments are used as-is.

# Arguments
- `T`: Temperature [K] (constant or function of time)
- `ψ`: Matric potential [kPa] (constant or function of time)
- `O2`: Ambient O₂ [μg/mm³] (constant or function of time)

# Examples

```julia
# All constant
env = EnvironmentalDrivers(293.15, -33.0, 0.27)

# Mixed constant and function
env = EnvironmentalDrivers(
    t -> 293.15 + 10.0 * sin(2π * t / 365),  # seasonal temperature
    -33.0,                                    # constant matric potential
    0.27                                      # constant O₂
)
```
"""
function EnvironmentalDrivers(T, ψ, O2)
    # Wrap scalars as constant functions
    T_func  = _wrap_constant(T)
    ψ_func  = _wrap_constant(ψ)
    O2_func = _wrap_constant(O2)

    # Construct with explicit type parameters using inner constructor
    EnvironmentalDrivers{typeof(T_func), typeof(ψ_func), typeof(O2_func)}(T_func, ψ_func, O2_func)
end

"""
    _wrap_constant(x)

Wrap a scalar as a constant function, or return a function unchanged.

Internal helper for EnvironmentalDrivers constructor.
"""
_wrap_constant(x::Number) = t -> x
_wrap_constant(f) = f  # Already a function or callable object

"""
    evaluate(env::EnvironmentalDrivers, t::Real)

Evaluate all environmental drivers at time `t`.

Returns a NamedTuple `(T, ψ, O2)` with current values.

# Arguments
- `env`: EnvironmentalDrivers instance
- `t`: Time [days]

# Returns
- Named tuple with fields `T` [K], `ψ` [kPa], `O2` [μg/mm³]

# Example

```julia
env = EnvironmentalDrivers(293.15, -33.0, 0.27)
vals = evaluate(env, 10.0)  # t = 10 days
# vals.T = 293.15, vals.ψ = -33.0, vals.O2 = 0.27
```
"""
function evaluate(env::EnvironmentalDrivers, t::Real)
    (T = env.T(t), ψ = env.ψ(t), O2 = env.O2(t))
end
