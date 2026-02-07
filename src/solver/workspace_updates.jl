# workspace_updates.jl
# Helper functions to update workspace arrays once per timestep

"""
    update_temperature_cache!(cache::TemperatureCache, T::Real,
                              bio::BiologicalProperties, soil::SoilProperties)

Update all temperature-dependent values in the cache.

# Arguments
- `cache::TemperatureCache`: Cache to update (modified in-place)
- `T::Real`: Temperature [K]
- `bio::BiologicalProperties`: Biological parameters
- `soil::SoilProperties`: Soil parameters (for constants/reference values)

# Notes
- Arrhenius factors for biological rates
- Pure-phase diffusion coefficients (Stokes-Einstein, Chapman-Enskog)
- Henry's law constant for O₂
- Should be called once per timestep (not per sub-step in Strang splitting)

# Manuscript reference
Architecture §6, §7: Temperature dependencies
"""
function update_temperature_cache!(cache::TemperatureCache, T::Real,
                                  bio::BiologicalProperties, soil::SoilProperties)
    # Use reference temperature from biological parameters
    T_ref = bio.T_ref

    # === Arrhenius factors for biological rates ===
    cache.f_bac = arrhenius(bio.Ea_B, T, T_ref)
    cache.f_fun = arrhenius(bio.Ea_F, T, T_ref)
    cache.f_eps = arrhenius(bio.Ea_EPS, T, T_ref)
    cache.f_maoc_s = arrhenius(bio.Ea_MAOC_sorb, T, T_ref)
    cache.f_maoc_d = arrhenius(bio.Ea_MAOC_desorb, T, T_ref)
    cache.f_pom = arrhenius(bio.Ea_POM, T, T_ref)

    # === Pure-phase diffusion coefficients ===
    # DOC: Stokes-Einstein with VFT viscosity
    # Use reference diffusion from soil parameters (already in mm²/day)
    cache.D_DOC_w = D_DOC_water(T, soil.D_C0_ref, T_ref)

    # O₂ in water: Stokes-Einstein (no reference value needed - empirical formula)
    cache.D_O2_w = D_O2_water(T)

    # O₂ in air: Chapman-Enskog
    # Use reference diffusion from soil parameters (already in mm²/day)
    cache.D_O2_a = D_O2_air(T, soil.D_O2_a_ref, T_ref)

    # Mobile fungi: temperature-dependent translocation (uses Ea_F)
    cache.D_Fm = D_fungal_translocation(T, bio.D_Fm0, bio.Ea_F, T_ref)

    # === Henry's law constant for O₂ ===
    cache.K_H_O = K_H_O2(T)

    nothing
end
