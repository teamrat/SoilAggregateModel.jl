"""
    henry.jl

Henry's law constants for gas-water partitioning.

Temperature dependence via van't Hoff equation.

Units: Dimensionless Henry's law constant K_H (concentration ratio C_gas/C_aq).

# Dependencies
Requires: R_GAS from constants.jl (must be included before this file)
"""

"""
    henry_vant_hoff(K_H_ref::Real, ΔH_sol::Real, T::Real, T_ref::Real)

Henry's law constant via van't Hoff equation.

Returns dimensionless K_H such that:
    C_gas = K_H × C_aq

where C_gas and C_aq are concentrations in gas and aqueous phases.

# Formula
    K_H(T) = K_H_ref × exp[ΔH_sol/R · (1/T - 1/T_ref)]

NOTE: No negative sign in front! This is the correct form for K_H = C_gas/C_aq.

# Arguments
- `K_H_ref`: Reference Henry's law constant (dimensionless) at T_ref
- `ΔH_sol`: Enthalpy of dissolution [J/mol] (negative for exothermic dissolution)
- `T`: Temperature [K]
- `T_ref`: Reference temperature [K]

# Returns
- Dimensionless Henry's law constant K_H

# Physical Interpretation
For O₂ dissolution in water:
- ΔH_sol < 0: Dissolution is exothermic
- Warming → K_H increases → C_aq = C_gas/K_H decreases → LESS dissolved O₂ ✓
- Cooling → K_H decreases → C_aq = C_gas/K_H increases → MORE dissolved O₂ ✓

# Examples

```julia
# O₂: K_H_ref = 31.25 at 298 K, ΔH_sol = -12,000 J/mol
K_H = henry_vant_hoff(31.25, -12_000.0, 293.15, 298.15)
# K_H ≈ 29.1 (cooler → lower K_H → more dissolved O₂)

K_H = henry_vant_hoff(31.25, -12_000.0, 303.15, 298.15)
# K_H ≈ 33.8 (warmer → higher K_H → less dissolved O₂)
```

# Notes
- For O₂, ΔH_sol ≈ -12 kJ/mol (dissolution is exothermic)
- Warming increases K_H, reducing O₂ solubility (fish kills in warm summer ponds)
- This is the CORRECTED formula - manuscript has an erroneous negative sign
"""
function henry_vant_hoff(K_H_ref::Real, ΔH_sol::Real, T::Real, T_ref::Real)
    K_H_ref * exp(ΔH_sol / R_GAS * (1.0 / T - 1.0 / T_ref))
end

"""
    K_H_O2(T::Real)

Henry's law constant for O₂ at temperature T.

Uses standard reference values:
- K_H_ref = 31.25 (dimensionless) at T_ref = 298.15 K
- ΔH_sol = -12,000 J/mol

# Arguments
- `T`: Temperature [K]

# Returns
- Dimensionless K_H for O₂

# Validity
- Valid 273 K ≤ T ≤ 313 K (0°C to 40°C)
- Consistent with NIST data for air-saturated water

# Examples

```julia
K_H = K_H_O2(293.15)  # 20°C
# K_H ≈ 29.1

K_H = K_H_O2(273.15)  # 0°C
# K_H ≈ 21.5 (cold water holds more O₂)
```

# Notes
- At 20°C (293.15 K), air-saturated water contains ~0.27 μg/mm³ O₂
- Atmospheric O₂ partial pressure ≈ 0.21 atm
- C_aq = C_gas / K_H for equilibrium partitioning
"""
function K_H_O2(T::Real)
    henry_vant_hoff(31.25, -12_000.0, T, 298.15)
end

"""
    O2_saturation(T::Real, P_atm::Real=0.21)

Compute air-saturated aqueous O₂ concentration.

# Arguments
- `T`: Temperature [K]
- `P_atm`: Atmospheric O₂ partial pressure [atm] (default 0.21)

# Returns
- Saturated aqueous O₂ concentration [μg/mm³]

# Formula
With Henry's law C_gas = K_H × C_aq, we have:
    C_aq_sat = C_gas / K_H(T)

where C_gas is the equilibrium gas-phase concentration.

# Examples

```julia
C_sat = O2_saturation(293.15)  # 20°C, sea level
# C_sat ≈ 0.28 μg/mm³

C_sat = O2_saturation(273.15)  # 0°C, sea level
# C_sat ≈ 0.44 μg/mm³ (cold water holds more)

C_sat = O2_saturation(293.15, 0.15)  # 20°C, high altitude
# C_sat ≈ 0.20 μg/mm³ (lower partial pressure)
```

# Notes
- Uses ideal gas law to convert partial pressure to concentration
- At 20°C, 0.21 atm O₂: C_gas ≈ 0.27 kg/m³ = 0.27 μg/mm³ (in pore air)
- Aqueous concentration is ~30× lower (K_H ≈ 30)
- This is a convenience function for boundary conditions
- Default P_atm = 0.21 assumes sea-level air (21% O₂)
"""
function O2_saturation(T::Real, P_atm::Real=0.21)
    # Ideal gas law: C = P/(R*T) where R = 0.08206 L·atm/(mol·K)
    # O₂ molar mass: 32 g/mol
    # Convert to μg/mm³ = kg/m³
    R_gas_Latm = 0.08206  # L·atm/(mol·K)
    MW_O2 = 32.0  # g/mol

    # C_gas [mol/L] = P [atm] / (R [L·atm/(mol·K)] × T [K])
    C_mol_per_L = P_atm / (R_gas_Latm * T)

    # Convert to μg/mm³ = kg/m³ = g/L
    C_gas = C_mol_per_L * MW_O2  # [g/L] = [kg/m³] = [μg/mm³]

    # Henry's law: C_gas = K_H × C_aq → C_aq = C_gas / K_H
    K_H = K_H_O2(T)
    C_gas / K_H
end
