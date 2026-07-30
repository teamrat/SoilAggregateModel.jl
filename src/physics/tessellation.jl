"""
    tessellation.jl

Domain tessellation and population upscaling.

The solver integrates ONE spherical domain around ONE POM particle. To compare
against a bulk-soil measurement, two things must be established:

1. **How large the model domain is relative to the POM particle**, and what
   correction that implies for background carbon (`ω`).
2. **How many such particles a unit volume of soil contains** (`N_POM`).

Both follow from measured quantities — the amendment rate, the POM carbon
density and the soil bulk density — so neither belongs in an experiment script.

# Dependencies
None. Pure geometry.
"""

"""
    domain_tessellation(; ρ_POM, I_input, ρ_b, f_domain_min=10.0)

Domain geometry and the background-carbon overlap correction.

Each POM particle of diameter `d` owns a cell of diameter `f_pack·d`, where
`f_pack` is fixed by requiring the cells to fill the soil exactly:

    φ_POM  = I_input · ρ_b / ρ_POM          POM volume fraction of bulk soil
    f_pack = (1/φ_POM)^(1/3)                cell diameter / POM diameter

The model domain extends to `f_domain·d`. **It must equal `f_pack·d`**, which
is the equivalent-sphere radius of the particle's own Voronoi cell. Zero flux at
that boundary is the correct symmetry condition between neighbouring particles,
and it makes

    ω = (f_domain / f_pack)³ = 1

# Why ω must be 1

`ω > 1` means the domain encloses more soil than the particle owns, and the
convention then divides background carbon (DOC, bacteria, all three fungal
pools, EPS, MAOC) by `ω` so that summing over particles gives the right
per-litre total.

**That is only valid for terms linear in concentration.** Every saturating term
in the model compares a state variable against a fixed concentration — `K_B`,
`K_F`, `C_B`, `B_S`, `F_S`, `B_min`, `F_i_min`, `E_min`, `K_E`, `K_B_P`,
`K_F_P`, `K_Fm_net`, `1/k_L`, `M_max`, `1/ω_E` — and none of those was divided.
Traced at ω = 30.3 on De Gryze soil 3 (REFERENCE.md §26 erratum 13): the O₂ sink
was 45× too weak, so near-POM anoxia never formed; the MAOC equilibrium left 17×
too little DOC; POM dissolution ran 1.75× too slow; bacterial maintenance and
mortality were 6.4× too low with `B` parked at 0.97·`B_min`; fungal yield was
2.35× too high and fungal respiration 44 % too low.

There is exactly one rescaling of constants that restores invariance — divide
every constant with units of carbon concentration by `ω` and multiply every
constant with units of inverse concentration by `ω`. That is a change of the
unit of carbon mass, not a modelling choice, and **it is blocked by POM**: the C
equation needs `R_P_max/ω` while the P equation needs `R_P_max` unchanged, and
reconciling them would require diluting `P_0`, which is precisely what makes the
amendment inventory correct. The two demands are irreconcilable, so no choice of
constants makes the oversized domain consistent.

**Overlap is the intended design.** The default is `f_domain_min = 10`, so
`ω ≈ 30`. `domain_tessellation` records the geometry with an `@info` rather than a warning — the configuration is deliberate. What remains
settled: the oversized domain is kept, because the outer boundary must stay clear
of the near-POM gradients — shrinking it to `f_pack` puts a zero-flux surface
1.5 mm from a DOC front that travels 3.2 mm — and the `1/ω` is applied once, to
the background pools at initialisation. It is NOT applied again at the summation
stage: dilution and volume-inflation already cancel there, so a second factor
would double-correct. Manuscript Appendix (Domain Tessellation and Overlap
Correction) has the derivation and the proof; REFERENCE.md §26 erratum 13 records
how this was briefly and wrongly reopened.

# The resolution argument was backwards

`f_domain_min = 10.0` was justified as needed "to resolve the near-POM
gradients". The grid is uniform, `h = (r_max − r_0)/(n − 1)`, so at fixed
`n_grid` a **smaller** domain is finer. For d = 1.325 mm at n = 200:
`f_domain = 10` gives h = 0.0300 mm; `f_domain = f_pack = 3.21` gives
h = 0.0074 mm — 4× finer, on the correct domain.

# Arguments
- `ρ_POM`: POM carbon density [µg-C/mm³]
- `I_input`: amendment rate [µg-C per µg-soil, dimensionless]
- `ρ_b`: soil bulk density [µg/mm³]
- `f_domain_min`: floor on the domain factor [-]. Default 10; see the ω note above.

# Returns
NamedTuple `(φ_POM, f_pack, f_domain, ω)`.

# Example
```julia
# De Gryze 2006: 1.5 g wheat stems (44.3 % C) per 150 g soil
tess = domain_tessellation(ρ_POM = 200.0, I_input = 4.43e-3, ρ_b = 1300.0)
tess.ω        # 28.8  (overlap factor; background pools start at 1/ω)
tess.f_domain # 10.0
```
"""
function domain_tessellation(; ρ_POM::Real, I_input::Real, ρ_b::Real,
                             f_domain_min::Real = 10.0)
    ρ_POM > 0 || throw(ArgumentError("ρ_POM must be > 0 (got $(ρ_POM))"))
    I_input > 0 || throw(ArgumentError("I_input must be > 0 (got $(I_input))"))
    ρ_b > 0 || throw(ArgumentError("ρ_b must be > 0 (got $(ρ_b))"))

    C_input_vol = I_input * ρ_b            # [µg-C/mm³]
    φ_POM       = C_input_vol / ρ_POM      # [-]
    φ_POM < 1.0 || throw(ArgumentError(
        "POM volume fraction φ_POM = $(φ_POM) ≥ 1 — the amendment cannot fill " *
        "more than the soil. Check ρ_POM, I_input and ρ_b."))

    f_pack   = (1.0 / φ_POM)^(1/3)
    f_domain = max(f_pack, f_domain_min)
    ω        = (f_domain / f_pack)^3

    # Recorded, not warned about: an overlapping domain is the
    # intended design (field scenarios with repeated input need model volume to
    # exceed physical volume), so a standing @warn on the chosen configuration
    # would be wallpaper. What is worth seeing is the actual geometry, because a
    # run at an unintended ω is otherwise invisible. The consistency question ω
    # raises — the state is diluted, the concentration constants are not — is
    # bounded by 1/ω; see REFERENCE.md §26 erratum 13.
    ω > 1.0 + 1e-12 && @info(
        "domain_tessellation: overlapping domains, ω = $(round(ω, digits=2)). " *
        "Background carbon is diluted by ω; see REFERENCE.md §26 erratum 13.",
        f_pack, f_domain, ω)

    (φ_POM = φ_POM, f_pack = f_pack, f_domain = f_domain, ω = ω)
end

"""
    pom_population(diam_mm, f_POM, tess; ρ_POM, soil_volume_mm3=1.0e6)

Number of POM particles and their carbon, per unit volume of soil.

Each size class `i` owns packing cells of volume `(4/3)π(dᵢ·f_pack/2)³`, so the
number of particles in that class is its share `f_POM[i]` of the soil volume
divided by the cell volume.

    N_POM[i] = f_POM[i] · V_soil / [(4/3)π(dᵢ·f_pack/2)³]
    P_0[i]   = (4/3)π(dᵢ/2)³ · ρ_POM

# Arguments
- `diam_mm`: POM diameters, one per size class [mm]
- `f_POM`: mass fraction in each class, **which must sum to 1**. The amendment
  identity below is exact only if it does; a distribution truncated at the bin
  edges (e.g. a Normal clipped at 0.5 and 2.0 mm) loses mass in direct
  proportion to the shortfall, and the recovered `I_input` is low by the same
  factor. Renormalise before calling.
- `tess`: result of [`domain_tessellation`](@ref)
- `ρ_POM`: POM carbon density [µg-C/mm³]
- `soil_volume_mm3`: reference soil volume [mm³] (default 1 litre)

# Returns
NamedTuple `(N_POM, P_0_per_particle, total_POM_C, V_pack)` — counts [-], carbon
per particle [µg-C], total amendment carbon in the reference volume [µg-C], and
the packing-cell volume owned by one particle of each class [mm³].

`V_pack` is returned because it is the soil each particle owns: an aggregate
whose volume approaches `V_pack` has consumed its entire share, and beyond that
the independent-domain picture has failed. Pass it to `population_statistics`
as `cell_volume_mm3` to have that checked every run.

# Verification
`total_POM_C / (soil_volume_mm3 · ρ_b · 1e-6)` must reproduce the experiment's
stated amendment rate in µg-C per g soil. If it does not, the tessellation
inputs disagree with the experiment.
"""
function pom_population(diam_mm::AbstractVector, f_POM::AbstractVector, tess;
                        ρ_POM::Real, soil_volume_mm3::Real = 1.0e6)
    length(diam_mm) == length(f_POM) || throw(ArgumentError(
        "diam_mm has $(length(diam_mm)) entries but f_POM has $(length(f_POM))"))

    V_pack = [(4.0/3.0) * π * (d * tess.f_pack / 2.0)^3 for d in diam_mm]
    N_POM  = f_POM .* soil_volume_mm3 ./ V_pack
    P_0    = [(4.0/3.0) * π * (d / 2.0)^3 * ρ_POM for d in diam_mm]

    (N_POM = N_POM, P_0_per_particle = P_0, total_POM_C = sum(N_POM .* P_0),
     V_pack = V_pack)
end

"""
    log_interpolate_fraction(lo, mid, hi)

Fraction of a size class `[lo, hi]` that falls below `mid`, assuming equal mass
per logarithmic interval.

    f = ln(mid/lo) / ln(hi/lo)

Used where a measured size fraction spans two reporting classes and the split is
not reported. It is an assumption about particle-size distributions, made
without reference to any observable, and must be recorded as such wherever it is
applied.

# Example
```julia
# Sand (53–2000 µm) split at the 250 µm sieve
log_interpolate_fraction(53, 250, 2000)   # 0.4272
```
"""
function log_interpolate_fraction(lo::Real, mid::Real, hi::Real)
    (0 < lo < mid < hi) || throw(ArgumentError(
        "require 0 < lo < mid < hi (got $lo, $mid, $hi)"))
    log(mid / lo) / log(hi / lo)
end
