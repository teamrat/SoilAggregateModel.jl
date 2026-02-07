"""
    water_retention.jl

Water retention via modified van Genuchten equation.

EPS and fungi modify the retention curve by changing the α parameter.

# Dependencies
Requires: SoilProperties from parameters.jl
"""

"""
    alpha_effective(E::Real, F_i::Real, soil::SoilProperties)

Compute effective van Genuchten α parameter with EPS and fungal modification.

# Formula
    α_eff = α_0 × exp(ω_E × E + ω_F × F_i)

where:
- α_0: Base van Genuchten α [1/kPa]
- ω_E: EPS effect (negative, increases retention) [mm³/μg]
- ω_F: Fungal effect (negative, increases retention) [mm³/μg]
- E: EPS concentration [μg/mm³]
- F_i: Insulated fungal biomass [μg/mm³]

# Arguments
- `E`: EPS concentration [μg/mm³]
- `F_i`: Insulated fungi concentration [μg/mm³]
- `soil`: SoilProperties struct

# Returns
- Effective α [1/kPa]

# Examples

```julia
soil = SoilProperties()
α_eff = alpha_effective(10.0, 5.0, soil)
# Higher EPS/fungi → lower α_eff → higher retention
```

# Notes
- ω_E, ω_F < 0 → exp(...) < 1 → α_eff < α_0
- Lower α → flatter retention curve → holds more water at given |ψ|
- EPS and fungi create micropores that retain water
"""
function alpha_effective(E::Real, F_i::Real, soil::SoilProperties)
    soil.α_vg * exp(soil.ω_E * E + soil.ω_F * F_i)
end

"""
    van_genuchten(ψ::Real, α::Real, n::Real, θ_r::Real, θ_s::Real)

Van Genuchten water retention model.

# Formula
    θ = θ_r + (θ_s - θ_r) / [1 + (α|ψ|)^n]^m

where m = 1 - 1/n.

# Arguments
- `ψ`: Matric potential [kPa] (ψ ≤ 0 for unsaturated)
- `α`: van Genuchten α parameter [1/kPa]
- `n`: van Genuchten n parameter [-]
- `θ_r`: Residual water content [-]
- `θ_s`: Saturated water content [-]

# Returns
- Volumetric water content θ [-]

# Examples

```julia
θ = van_genuchten(-33.0, 0.01, 1.4, 0.05, 0.45)
# At field capacity (-33 kPa), typical loam
```

# Notes
- For ψ = 0 (saturated): θ = θ_s
- For ψ → -∞ (dry): θ → θ_r
- Typical field capacity: ψ ≈ -33 kPa
- Permanent wilting point: ψ ≈ -1500 kPa
"""
function van_genuchten(ψ::Real, α::Real, n::Real, θ_r::Real, θ_s::Real)
    if ψ >= 0.0
        # Saturated
        return θ_s
    else
        m = 1.0 - 1.0 / n
        θ_r + (θ_s - θ_r) / (1.0 + (α * abs(ψ))^n)^m
    end
end

"""
    water_content(ψ::Real, E::Real, F_i::Real, soil::SoilProperties)

Compute water content with EPS/fungal modification.

Combines modified α with van Genuchten equation.

# Arguments
- `ψ`: Matric potential [kPa]
- `E`: EPS concentration [μg/mm³]
- `F_i`: Insulated fungi concentration [μg/mm³]
- `soil`: SoilProperties struct

# Returns
- Volumetric water content θ [-]

# Examples

```julia
soil = SoilProperties()
θ = water_content(-33.0, 10.0, 5.0, soil)
# Higher EPS/fungi → higher θ at same ψ
```
"""
function water_content(ψ::Real, E::Real, F_i::Real, soil::SoilProperties)
    α_eff = alpha_effective(E, F_i, soil)
    van_genuchten(ψ, α_eff, soil.n_vg, soil.θ_r, soil.θ_s)
end

"""
    update_water_content!(θ::AbstractVector, θ_a::AbstractVector,
                         ψ::Real, state::AggregateState, soil::SoilProperties)

Update water content and air-filled porosity arrays (in-place).

Computes θ(r) and θ_a(r) for all grid points based on current EPS and fungi distributions.

# Arguments
- `θ`: Water content array [n], **overwritten**
- `θ_a`: Air-filled porosity array [n], **overwritten**
- `ψ`: Current matric potential [kPa] (spatially uniform)
- `state`: AggregateState with current E and F_i profiles
- `soil`: SoilProperties struct

# Returns
- Nothing (modifies θ and θ_a in-place)

# Notes
- Air-filled porosity: θ_a = θ_s - θ
- Called once per timestep before diffusion calculation
- Zero allocations (uses pre-allocated workspace arrays)

# Examples

```julia
# In timestepper
ψ = env.ψ(t)
update_water_content!(ws.θ, ws.θ_a, ψ, solver.state, solver.soil)
```
"""
function update_water_content!(θ::AbstractVector, θ_a::AbstractVector,
                               ψ::Real, state::AggregateState, soil::SoilProperties)
    n = length(θ)

    @inbounds for i in 1:n
        # Compute effective α with local EPS and fungi
        α_eff = alpha_effective(state.E[i], state.F_i[i], soil)

        # Van Genuchten
        θ[i] = van_genuchten(ψ, α_eff, soil.n_vg, soil.θ_r, soil.θ_s)

        # Air-filled porosity
        θ_a[i] = soil.θ_s - θ[i]
    end

    nothing
end

"""
    tortuosity_millington_quirk(θ::Real, θ_s::Real)

Millington-Quirk tortuosity factor for aqueous diffusion.

# Formula
    τ(θ) = θ² / θ_s^(2/3)

# Arguments
- `θ`: Volumetric water content [-]
- `θ_s`: Saturated water content [-]

# Returns
- Tortuosity factor [-] (0 ≤ τ ≤ θ_s^(1/3))

# Examples

```julia
τ = tortuosity_millington_quirk(0.3, 0.45)
# τ ≈ 0.133 (diffusion reduced by factor of ~7.5)
```

# Notes
- At saturation (θ = θ_s): τ = θ_s^(1/3) ≈ 0.77 for typical soil
- Dry soil (θ → 0): τ → 0 (diffusion stops)
- Accounts for pore connectivity and path tortuosity
- Valid for aqueous diffusion only (not gas phase)
"""
function tortuosity_millington_quirk(θ::Real, θ_s::Real)
    θ^2 / θ_s^(2/3)
end
