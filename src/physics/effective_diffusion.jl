"""
    effective_diffusion.jl

Effective diffusion coefficients for all species.

Combines pure-phase diffusion (temperature-dependent), tortuosity (water content),
and species-specific modifications (sorption, dual-phase transport).

# Dependencies
Requires: TemperatureCache, SoilProperties, tortuosity_millington_quirk
"""

"""
    D_eff_DOC(D_DOC_w::Real, θ::Real, θ_s::Real, ρ_b::Real, k_d::Real)

Effective diffusion coefficient for dissolved organic carbon.

# Formula
    D_C = D_DOC_w × τ(θ) × [θ/(θ + ρ_b·k_d)]

where:
- τ(θ) = θ²/θ_s^(2/3) (Millington-Quirk tortuosity)
- θ/(θ + ρ_b·k_d) = retardation factor (sorption reduces mobility)

# Arguments
- `D_DOC_w`: DOC diffusion in pure water [mm²/day]
- `θ`: Water content [-]
- `θ_s`: Saturated water content [-]
- `ρ_b`: Bulk density [μg/mm³]
- `k_d`: Partition coefficient [mm³/μg]

# Returns
- Effective D_C [mm²/day]

# Notes
- Sorption retardation: C partitions between solution and solid phase
- Higher k_d → more sorption → slower diffusion
- At saturation with no sorption: D_C = D_DOC_w × θ_s^(1/3) ≈ 0.77 D_DOC_w
"""
function D_eff_DOC(D_DOC_w::Real, θ::Real, θ_s::Real, ρ_b::Real, k_d::Real)
    τ = tortuosity_millington_quirk(θ, θ_s)
    retardation = θ / (θ + ρ_b * k_d)
    D_DOC_w * τ * retardation
end

"""
    D_eff_bacteria(D_C::Real, D_B_rel::Real)

Effective diffusion coefficient for bacteria (chemotaxis).

# Formula
    D_B = D_B_rel × D_C

where D_B_rel ≪ 1 (bacterial motility is much slower than DOC diffusion).

# Arguments
- `D_C`: Effective DOC diffusion [mm²/day]
- `D_B_rel`: Relative bacterial motility [-] (typically ~0.001)

# Returns
- Effective D_B [mm²/day]

# Notes
- Bacteria move via chemotaxis toward DOC gradients
- Movement is slow compared to molecular diffusion
- D_B_rel captures swimming speed and chemotactic efficiency
- Typical value: D_B_rel ≈ 0.001 → D_B ≈ 0.001 × D_C
"""
function D_eff_bacteria(D_C::Real, D_B_rel::Real)
    D_B_rel * D_C
end

"""
    D_eff_fungi_noninsulated(D_Fn0_ref::Real, f_Arrhenius::Real, θ::Real, θ_s::Real)

Effective diffusion coefficient for non-insulated fungi (hyphal extension).

# Formula
    D_Fn = D_Fn0_ref × f_Arrhenius × τ(θ)

where:
- D_Fn0_ref: Reference hyphal extension rate [mm²/day]
- f_Arrhenius: Temperature factor from Arrhenius equation
- τ(θ): Millington-Quirk tortuosity

# Arguments
- `D_Fn0_ref`: Reference diffusion at T_ref [mm²/day]
- `f_Arrhenius`: Arrhenius temperature factor [-]
- `θ`: Water content [-]
- `θ_s`: Saturated water content [-]

# Returns
- Effective D_Fn [mm²/day]

# Notes
- Hyphal extension occurs at hyphal tips through pore space
- Temperature-dependent (biological process)
- Tortuosity-limited (hyphae follow pore connectivity)
- Dry soil (θ → 0) stops hyphal extension
"""
function D_eff_fungi_noninsulated(D_Fn0_ref::Real, f_Arrhenius::Real, θ::Real, θ_s::Real)
    τ = tortuosity_millington_quirk(θ, θ_s)
    D_Fn0_ref * f_Arrhenius * τ
end

"""
    D_eff_fungi_mobile(D_Fm0_ref::Real, f_Arrhenius::Real)

Effective diffusion coefficient for mobile fungi (internal translocation).

# Formula
    D_Fm = D_Fm0_ref × f_Arrhenius

No tortuosity factor! Internal translocation occurs within the hyphal network,
not through pore space.

# Arguments
- `D_Fm0_ref`: Reference translocation rate at T_ref [mm²/day]
- `f_Arrhenius`: Arrhenius temperature factor [-]

# Returns
- Effective D_Fm [mm²/day] (scalar, spatially uniform)

# Notes
- Internal transport through existing hyphae
- NOT limited by water content or pore connectivity
- Mobile fungi can transport through dry regions
- Temperature-dependent (active biological process)
- This is a key mechanism for nutrient scavenging in dry zones
"""
function D_eff_fungi_mobile(D_Fm0_ref::Real, f_Arrhenius::Real)
    D_Fm0_ref * f_Arrhenius
end

"""
    D_eff_oxygen(D_O2_w::Real, D_O2_a::Real, K_H::Real, θ::Real, θ_a::Real, θ_s::Real)

Effective diffusion coefficient for oxygen (dual-phase transport).

# Formula
    D_O = D_w × [θ/(θ + K_H·θ_a)] × [θ²/θ_s^(2/3)]
        + D_a × [K_H·θ_a/(θ + K_H·θ_a)] × [θ_a^(10/3)/θ_s²]

where:
- First term: aqueous phase diffusion (Millington-Quirk tortuosity)
- Second term: gas phase diffusion (stronger tortuosity dependence)
- Partitioning: θ/(θ + K_H·θ_a) = fraction in aqueous phase

# Arguments
- `D_O2_w`: O₂ diffusion in pure water [mm²/day]
- `D_O2_a`: O₂ diffusion in air [mm²/day]
- `K_H`: Henry's law constant [-] (C_gas = K_H × C_aq)
- `θ`: Water content [-]
- `θ_a`: Air-filled porosity [-]
- `θ_s`: Saturated water content [-]

# Returns
- Effective D_O [mm²/day]

# Notes
- O₂ diffuses through both water-filled and air-filled pores
- Gas-phase diffusion is ~10⁴ faster than aqueous, but path is more tortuous
- At saturation (θ_a → 0): only aqueous diffusion
- Dry soil (θ → 0, θ_a → θ_s): only gas diffusion
- Typical wet soil: gas phase dominates despite lower connectivity
"""
function D_eff_oxygen(D_O2_w::Real, D_O2_a::Real, K_H::Real, θ::Real, θ_a::Real, θ_s::Real)
    # Partition fractions
    f_aq = θ / (θ + K_H * θ_a)
    f_gas = K_H * θ_a / (θ + K_H * θ_a)

    # Tortuosity factors
    τ_aq = θ^2 / θ_s^(2/3)
    τ_gas = θ_a^(10/3) / θ_s^2

    # Dual-phase diffusion
    D_O2_w * f_aq * τ_aq + D_O2_a * f_gas * τ_gas
end

"""
    update_effective_diffusion!(workspace::Workspace, soil::SoilProperties,
                                bio::BiologicalProperties, temp_cache::TemperatureCache)

Update all effective diffusion coefficients in workspace (in-place).

Computes D_C, D_B, D_Fn, D_Fm, D_O for all grid points.
Called once per timestep after water content update.

# Arguments
- `workspace`: Workspace with θ, θ_a, D_C, D_B, D_Fn, D_Fm, D_O, f_T
- `soil`: SoilProperties
- `bio`: BiologicalProperties (for D_Fn0, D_Fm0)
- `temp_cache`: TemperatureCache with current temperature-dependent values

# Returns
- Nothing (modifies workspace.D_C, D_B, D_Fn, D_Fm, D_O, and temp_cache.D_Fm in-place)

# Notes
- Requires θ and θ_a to be updated first (call update_water_content! before this)
- Zero allocations (uses pre-allocated workspace arrays)
- D_Fm is spatially uniform (constant value at all nodes) but stored as vector for solver

# Examples

```julia
# In timestepper
T = env.T(t)
update_temperature_cache!(ws.f_T, T, bio, soil)
update_water_content!(ws.θ, ws.θ_a, ψ, state, soil)
update_effective_diffusion!(ws, soil, bio, ws.f_T)
# Now ws.D_C, D_B, D_Fn, D_O are ready for diffusion step
```
"""
function update_effective_diffusion!(workspace::Workspace, soil::SoilProperties,
                                     bio::BiologicalProperties, temp_cache::TemperatureCache)
    n = length(workspace.θ)

    # Mobile fungi (spatially uniform, no tortuosity)
    D_Fm_val = D_eff_fungi_mobile(bio.D_Fm0, temp_cache.f_fun)
    temp_cache.D_Fm = D_Fm_val  # Store in cache for convenience

    # Spatially varying diffusion coefficients
    @inbounds for i in 1:n
        θ_i = workspace.θ[i]
        θ_a_i = workspace.θ_a[i]

        # DOC with sorption retardation
        workspace.D_C[i] = D_eff_DOC(temp_cache.D_DOC_w, θ_i, soil.θ_s,
                                      soil.ρ_b, soil.k_d_eq)

        # Bacteria (chemotaxis proportional to DOC gradient)
        workspace.D_B[i] = D_eff_bacteria(workspace.D_C[i], soil.D_B_rel)

        # Non-insulated fungi (hyphal extension)
        workspace.D_Fn[i] = D_eff_fungi_noninsulated(bio.D_Fn0, temp_cache.f_fun,
                                                       θ_i, soil.θ_s)

        # Mobile fungi (spatially uniform)
        workspace.D_Fm[i] = D_Fm_val

        # Oxygen (dual-phase: aqueous + gas)
        workspace.D_O[i] = D_eff_oxygen(temp_cache.D_O2_w, temp_cache.D_O2_a,
                                         temp_cache.K_H_O, θ_i, θ_a_i, soil.θ_s)
    end

    nothing
end
