"""
    arrhenius.jl

Arrhenius temperature dependence for biological rates.

Symbol note: 𝓔_a (mathcal E) is activation energy, NOT E (EPS state variable).

# Dependencies
Requires: R_GAS from constants.jl (must be included before this file)
"""

"""
    arrhenius(Ea::Real, T::Real, T_ref::Real)

Compute Arrhenius temperature factor.

Returns the dimensionless multiplier f(T) such that:
    k(T) = k_ref × f(T)

where:
    f(T) = exp[𝓔_a/R · (1/T_ref - 1/T)]

# Arguments
- `Ea`: Activation energy [J/mol]
- `T`: Current temperature [K]
- `T_ref`: Reference temperature [K]

# Returns
- Dimensionless factor (= 1.0 at T = T_ref)

# Examples

```julia
# Bacterial metabolism with Ea = 60 kJ/mol
f = arrhenius(60_000.0, 298.15, 293.15)
# f > 1.0 because T > T_ref (reaction faster at higher temperature)

# At reference temperature
f = arrhenius(60_000.0, 293.15, 293.15)
# f = 1.0 exactly
```

# Notes
- Higher Ea → stronger temperature sensitivity
- f > 1 when T > T_ref (rate increases)
- f < 1 when T < T_ref (rate decreases)
- Valid for T > 0 K
"""
function arrhenius(Ea::Real, T::Real, T_ref::Real)
    exp(Ea / R_GAS * (1.0 / T_ref - 1.0 / T))
end

"""
    arrhenius_ratio(Ea::Real, T1::Real, T2::Real)

Compute ratio of Arrhenius factors at two temperatures.

Returns k(T2) / k(T1) without needing a reference temperature.

# Arguments
- `Ea`: Activation energy [J/mol]
- `T1`: First temperature [K]
- `T2`: Second temperature [K]

# Returns
- Dimensionless ratio k(T2) / k(T1)

# Examples

```julia
# Q10 approximation: ratio over 10 K
ratio = arrhenius_ratio(60_000.0, 293.15, 303.15)
# ratio ≈ 2.3 for Ea = 60 kJ/mol
```
"""
function arrhenius_ratio(Ea::Real, T1::Real, T2::Real)
    exp(Ea / R_GAS * (1.0 / T1 - 1.0 / T2))
end

"""
    q10_equivalent(Ea::Real, T::Real)

Compute the effective Q10 at a given temperature.

Q10 is the factor by which a rate increases over a 10 K temperature increase.

# Arguments
- `Ea`: Activation energy [J/mol]
- `T`: Temperature [K]

# Returns
- Q10 value (typically 2-3 for biological processes)

# Examples

```julia
# Bacterial metabolism at 20°C
q10 = q10_equivalent(60_000.0, 293.15)
# q10 ≈ 2.3
```

# Notes
- Q10 is not constant but depends on temperature
- Higher Ea → higher Q10
- Q10 decreases slightly with increasing temperature for fixed Ea
"""
function q10_equivalent(Ea::Real, T::Real)
    arrhenius_ratio(Ea, T, T + 10.0)
end
