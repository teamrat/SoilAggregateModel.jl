# particle_size.jl
# Representative particle diameters from a particle-size distribution.
#
# The aggregate stability criterion needs ONE length that stands for the whole
# particle-size distribution. Which mean is correct is not a matter of taste: it
# is fixed by what the length is used for. The criterion counts bonds across a
# failure surface, bonds live on particle surfaces, and the number of bonds per
# unit area therefore scales with the interfacial area per unit volume. The mean
# that preserves the volume-to-surface ratio is the Sauter mean diameter
# (Sauter, 1926), conventionally d₃₂.
#
# See docs/REFERENCE.md §4.4.

"""
    sauter_mean_diameter(mass_fractions, diameters) -> Float64

Sauter mean diameter d₃₂ — the diameter of a sphere with the same volume-to-
surface-area ratio as the distribution.

    d₃₂ = 1 / Σᵢ (fᵢ / dᵢ)

for mass (equivalently volume, at constant particle density) fractions `fᵢ`
summing to 1 and representative diameters `dᵢ`.

# Why this mean

For spheres the surface area per unit volume of solid is `S_v = 6 Σᵢ fᵢ/dᵢ`, so
`d₃₂ = 6/S_v`. Any quantity proportional to interfacial area — the number of
inter-particle contacts per unit area of a failure surface, hence the binding
strength at fixed binder concentration — scales as `1/d₃₂`. A geometric or
arithmetic mean does **not** have this property: both are dominated by the
coarse fraction, which contributes the fewest contacts per unit mass.

# Arguments
- `mass_fractions`: fractions summing to 1 (tolerance 1e-6)
- `diameters`: representative diameter of each class, same units as the result,
  all strictly positive

# Returns
`d₃₂` in the units of `diameters`.

# Sensitivity
`d₃₂` is dominated by the finest class, because that class contributes the
largest `fᵢ/dᵢ`. The representative diameter chosen for the finest class is
therefore the most consequential input, and it should be stated wherever this
function is used. Ordering between soils is robust to that choice; the spread
between them is not.

# Example
```julia
# Sand / silt / clay with geometric class midpoints [mm]
sauter_mean_diameter([0.44, 0.45, 0.11], [0.316, 0.010, 0.001])   # 0.00639 mm
```

# References
- Sauter, J. (1926). *Die Größenbestimmung der in Gemischnebeln von
  Verbrennungskraftmaschinen vorhandenen Brennstoffteilchen.* VDI-Verlag,
  Berlin.
"""
function sauter_mean_diameter(mass_fractions::AbstractVector,
                              diameters::AbstractVector)
    length(mass_fractions) == length(diameters) || throw(DimensionMismatch(
        "mass_fractions has $(length(mass_fractions)) entries but diameters has $(length(diameters))"))
    isempty(mass_fractions) && throw(ArgumentError("mass_fractions must not be empty"))
    all(f -> f >= 0.0, mass_fractions) || throw(ArgumentError(
        "mass_fractions must be non-negative (got $(mass_fractions))"))
    all(d -> d > 0.0, diameters) || throw(ArgumentError(
        "diameters must be strictly positive (got $(diameters))"))
    isapprox(sum(mass_fractions), 1.0; atol = 1e-6) || throw(ArgumentError(
        "mass_fractions must sum to 1 (got $(sum(mass_fractions)))"))

    inv_d = sum(f / d for (f, d) in zip(mass_fractions, diameters))
    inv_d > 0.0 || throw(ArgumentError("degenerate distribution: Σ f/d = 0"))
    return 1.0 / inv_d
end

"""
Geometric midpoints of the USDA texture classes, in mm.

- `clay`  = 0.001  (class is < 0.002 mm; 0.001 is the midpoint of 0.0005–0.002)
- `silt`  = 0.010  (0.002–0.05 mm, geometric midpoint 0.01)
- `sand`  = 0.316  (0.05–2.0 mm, geometric midpoint 0.316)

These are **assumptions**, not measurements: a texture triple gives three mass
fractions, not a distribution. The clay value matters most (see
[`sauter_mean_diameter`](@ref)) because the class has no measured lower bound.
Any problem with a measured particle-size distribution should pass that
distribution directly rather than going through the texture triple.
"""
const TEXTURE_CLASS_DIAMETERS = (clay = 0.001, silt = 0.010, sand = 0.316)

"""
    sauter_from_texture(sand, silt, clay; d_sand, d_silt, d_clay) -> Float64

Sauter mean diameter [mm] from a USDA texture triple.

Convenience wrapper on [`sauter_mean_diameter`](@ref) using
[`TEXTURE_CLASS_DIAMETERS`](@ref) by default. The three fractions must sum to 1.

The class representative diameters are keyword arguments precisely so that the
assumption is visible at every call site and can be varied in a sensitivity
analysis.

# Example
```julia
sauter_from_texture(0.44, 0.45, 0.11)   # 0.00639 mm — De Gryze soil 3
```
"""
function sauter_from_texture(sand::Real, silt::Real, clay::Real;
                             d_sand::Real = TEXTURE_CLASS_DIAMETERS.sand,
                             d_silt::Real = TEXTURE_CLASS_DIAMETERS.silt,
                             d_clay::Real = TEXTURE_CLASS_DIAMETERS.clay)
    return sauter_mean_diameter([sand, silt, clay], [d_sand, d_silt, d_clay])
end
