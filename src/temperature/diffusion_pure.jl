"""
    diffusion_pure.jl

Pure-phase diffusion coefficients as functions of temperature.

Three mechanisms:
1. Stokes-Einstein + VFT (DOC in water)
2. Han & Bartels empirical (O₂ in water)
3. Chapman-Enskog kinetic theory (O₂ in air)

Units: All diffusion coefficients in mm²/day.

# Dependencies
Requires: R_GAS from constants.jl (must be included before this file)
"""

"""
    water_viscosity(T::Real)

Water viscosity via Vogel-Fulcher-Tammann equation.

Returns dynamic viscosity [mPa·s] = [cP].

# Formula
    ln[η(T) / mPa·s] = -3.7188 + 578.919 / (T - 137.546)

# Arguments
- `T`: Temperature [K]

# Returns
- Dynamic viscosity [mPa·s]

# Validity
- Valid 273 K ≤ T ≤ 373 K (0°C to 100°C)
- Accurate to within 2% over this range

# Examples

```julia
η = water_viscosity(293.15)  # 20°C
# η ≈ 1.002 mPa·s

η = water_viscosity(273.15)  # 0°C
# η ≈ 1.787 mPa·s
```

# References
- Vogel (1921), Fulcher (1925), Tammann & Hesse (1926)
- NIST Chemistry WebBook
"""
function water_viscosity(T::Real)
    exp(-3.7188 + 578.919 / (T - 137.546))
end

"""
    D_DOC_water(T::Real, D_ref::Real, T_ref::Real)

DOC diffusion in water via Stokes-Einstein relation.

Assumes diffusivity scales with T/η (temperature/viscosity).

# Formula
    D(T) = D_ref × (T/T_ref) × [η(T_ref) / η(T)]

# Arguments
- `T`: Temperature [K]
- `D_ref`: Reference diffusion coefficient [mm²/day]
- `T_ref`: Reference temperature [K]

# Returns
- Diffusion coefficient [mm²/day]

# Examples

```julia
# Typical DOC: D_ref ≈ 0.864 mm²/day at 20°C
D = D_DOC_water(298.15, 0.864, 293.15)
# D ≈ 1.03 mm²/day (warmer → faster diffusion)
```

# Notes
- Valid for dissolved organic molecules with MW ~200-1000 Da
- Stokes-Einstein is exact for spherical particles in continuum limit
- Viscosity from VFT equation (valid 273-373 K)
"""
function D_DOC_water(T::Real, D_ref::Real, T_ref::Real)
    η_ref = water_viscosity(T_ref)
    η_T   = water_viscosity(T)
    D_ref * (T / T_ref) * (η_ref / η_T)
end

"""
    D_O2_water(T::Real)

O₂ diffusion in water via Han & Bartels (1996) empirical correlation.

Returns diffusion coefficient [mm²/day].

# Formula
    log₁₀[D / cm²s⁻¹] = -4.410 + 773.8/T - (506.4/T)²

Converted to mm²/day by multiplying cm²/s by 8.64 × 10⁶.

# Arguments
- `T`: Temperature [K]

# Returns
- Diffusion coefficient [mm²/day]

# Validity
- Valid 273 K ≤ T ≤ 313 K (0°C to 40°C)
- Accurate to within 3% over this range

# Examples

```julia
D = D_O2_water(293.15)  # 20°C
# D ≈ 1.73 mm²/day

D = D_O2_water(298.15)  # 25°C
# D ≈ 2.11 mm²/day
```

# References
- Han & Bartels (1996), J. Chem. Eng. Data 41:1024-1027
- DOI: 10.1021/je960093x
"""
function D_O2_water(T::Real)
    # Compute in cm²/s
    T_inv = 1.0 / T
    log10_D_cm2s = -4.410 + 773.8 * T_inv - (506.4 * T_inv)^2
    D_cm2s = 10.0^log10_D_cm2s

    # Convert to mm²/day
    # 1 cm²/s = 10² mm²/s × 86400 s/day = 8.64e6 mm²/day
    D_cm2s * 8.64e6
end

"""
    D_O2_air(T::Real, D_ref::Real, T_ref::Real)

O₂ diffusion in air via Chapman-Enskog kinetic theory.

# Formula
    D(T) = D_ref × (T / T_ref)^1.75

# Arguments
- `T`: Temperature [K]
- `D_ref`: Reference diffusion coefficient [mm²/day]
- `T_ref`: Reference temperature [K]

# Returns
- Diffusion coefficient [mm²/day]

# Examples

```julia
# Typical O₂ in air: D_ref ≈ 1.73e6 mm²/day at 20°C
D = D_O2_air(298.15, 1.73e6, 293.15)
# D ≈ 1.87e6 mm²/day
```

# Notes
- Chapman-Enskog theory for dilute gases
- Exponent 1.75 assumes hard-sphere collisions with ω = 0.5
- Valid for ideal gas behavior (1 atm, T > 200 K)
- O₂ diffusion in air ~1000× faster than in water

# References
- Chapman & Cowling (1970), The Mathematical Theory of Non-uniform Gases
- Bird, Stewart, & Lightfoot (2002), Transport Phenomena
"""
function D_O2_air(T::Real, D_ref::Real, T_ref::Real)
    D_ref * (T / T_ref)^1.75
end

"""
    D_fungal_translocation(T::Real, D_ref::Real, Ea::Real, T_ref::Real)

Mobile fungi translocation rate (internal hyphal transport).

Unlike aqueous diffusion, fungal translocation is an active biological process
and follows Arrhenius temperature dependence (like other fungal rates).

# Formula
    D_Fm(T) = D_ref × exp[Ea/R · (1/T_ref - 1/T)]

# Arguments
- `T`: Temperature [K]
- `D_ref`: Reference translocation rate [mm²/day]
- `Ea`: Activation energy for fungal processes [J/mol]
- `T_ref`: Reference temperature [K]

# Returns
- Translocation rate [mm²/day]

# Examples

```julia
# Ea_F = 55,000 J/mol (from BiologicalProperties)
D = D_fungal_translocation(298.15, 1.0, 55_000.0, 293.15)
# D ≈ 1.35 mm²/day (biological process, faster at higher T)
```

# Notes
- This is NOT tortuosity-limited (occurs within hyphal network)
- Shares activation energy with other fungal rates (Ea_F)
- Compute using arrhenius() from arrhenius.jl:
  D_Fm(T) = D_Fm_ref × arrhenius(Ea_F, T, T_ref)
"""
function D_fungal_translocation(T::Real, D_ref::Real, Ea::Real, T_ref::Real)
    # Uses R_GAS from constants.jl
    D_ref * exp(Ea / R_GAS * (1.0 / T_ref - 1.0 / T))
end
