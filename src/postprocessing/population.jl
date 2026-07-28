# population.jl
# Population upscaling of single-aggregate results.
#
# The solver integrates ONE domain around ONE POM particle. A wet-sieving
# measurement reports a whole soil sample. This file holds the arithmetic that
# links the two: what one aggregate weighs, which sieve class it lands in, and
# how a population of them adds up to a sample statistic.
#
# Everything here operates on plain arrays and returns plain arrays. Table
# construction, column naming and file output are packaging and stay with the
# experiment script.
#
# Companion: src/physics/tessellation.jl supplies `n_dist` (particle counts) and
# the ω convention. See docs/REFERENCE.md §5b and §5c.

"""
    sieve_class(D, sieve_sizes) -> Int

Index of the sieve class containing diameter `D`, ascending.

`k` sieve apertures define `k+1` classes. `sieve_sizes` must be ascending. Class
1 is everything below the first aperture; class `k+1` is everything at or above
the last. The convention is `D < aperture` passes, matching wet sieving where
material smaller than the mesh is washed through.

# Example
```julia
sieve_class(0.03, [0.053, 0.25, 2.0])   # 1  (<53 µm)
sieve_class(2.0,  [0.053, 0.25, 2.0])   # 4  (≥2 mm)
```
"""
function sieve_class(D::Real, sieve_sizes::AbstractVector)
    k = 1
    for e in sieve_sizes
        D < e && break
        k += 1
    end
    return k
end

"""
    aggregate_mass(D_agg, r_0, P_C; ρ_b, f_C_POM) -> Float64

Dry mass [µg] of one aggregate: a mineral shell around a POM core.

    m = ρ_b·(4/3)π·max((D_agg/2)³ − r_0³, 0)   mineral shell
      + P_C / f_C_POM                           remaining residue

The core term is what makes retained mass time-dependent. Aggregate *volume*
does not: `r_agg` can sit still for days while the residue inside it is being
respired away. Weighting a population by volume therefore both overstates the
mass of POM-rich aggregates and discards the only time dependence the retained
mass has.

`max(..., 0)` guards the case `D_agg/2 < r_0`; `compute_r_agg` floors the
radius at the POM surface, so it should not arise.

# Arguments
- `D_agg`: aggregate diameter [mm]
- `r_0`: POM core radius [mm]
- `P_C`: POM carbon remaining in this aggregate [µg-C]

# Keyword Arguments
- `ρ_b`: bulk density [µg/mm³] of the **amended** sample — soil plus whatever
  residue has been mixed into it. Taking it as the amended value rather than
  the bare soil is a deliberate simplification: it makes `ρ_b·V` the total
  sample mass, so no separate residue bookkeeping is needed anywhere
  downstream. The amendment is ~1 % of soil mass, so the two readings differ
  by ~1 %. See REFERENCE.md §5c and the experiment's own notes.
- `f_C_POM`: carbon mass fraction of the residue [-], converting µg-C to dry
  mass. Experiment-specific — there is no default.
"""
function aggregate_mass(D_agg::Real, r_0::Real, P_C::Real;
                        ρ_b::Real, f_C_POM::Real)
    f_C_POM > 0 || throw(ArgumentError("f_C_POM must be > 0 (got $(f_C_POM))"))
    shell = ρ_b * (4.0/3.0) * π * max((D_agg / 2.0)^3 - r_0^3, 0.0)
    core  = P_C / f_C_POM
    return shell + core
end

"""
    population_statistics(D_agg, POM, CO2_cum, CO2_flux, r_0, n_dist;
                          sieve_sizes, mineral_class_fractions,
                          class_nominal_mm, soil_volume_mm3,
                          ρ_b, f_C_POM) -> NamedTuple

Sample-level statistics from a size-resolved population of single-aggregate runs.

Each of `D_agg`, `POM`, `CO2_cum`, `CO2_flux` is a `n_sizes × n_times` matrix —
one row per POM size class, one column per output time. `r_0` and `n_dist` are
length-`n_sizes` vectors: the POM core radius [mm] and the number of particles
of that class in `soil_volume_mm3`.

# Returns

NamedTuple of length-`n_times` vectors, plus one matrix:

| Field | Meaning |
|---|---|
| `MWD_agg_only` | mass-weighted mean diameter over aggregates only [mm] |
| `f_agg` | aggregate share of total sample **mass** [-]. **Not clamped** — a value > 1 means the model has produced more aggregate than there is sample, and every class share is then meaningless |
| `f_agg_vol` | aggregate share of sample **volume** [-] — geometry check |
| `POM_mass_frac` | residue share of retained coarse mass [-] |
| `shell_thickness_mm` | mass-weighted mean shell thickness `r_agg − r_0` [mm] — the continuum-validity check, see below |
| `cell_occupancy` | `n_sizes × n_times`, aggregate volume ÷ the packing-cell volume that size class owns [-]. `NaN` unless `cell_volume_mm3` is given. **A value ≥ 1 means that class has consumed more soil than it owns** |
| `CO2_total` | population cumulative CO₂ [µg-C] in `soil_volume_mm3` |
| `CO2_flux_total` | population CO₂ rate [µg-C/day] |
| `class_pct` | `(k+1) × n_times` matrix, % of sample mass per sieve class |
| `MWD_fixed_weight` | fixed-weight MWD [mm], or `NaN` if no nominals given |

# Physics

**Mean weight diameter over aggregates.**

    MWD_agg(t) = Σᵢ nᵢ·mᵢ(t)·Dᵢ(t) / Σᵢ nᵢ·mᵢ(t)

with `mᵢ` from [`aggregate_mass`](@ref). This is the mean over POM-centred
aggregates and **is not** the whole-sample statistic a sieving paper reports:
it excludes the unaggregated matrix entirely.

**Sieve classes as a share of the whole sample.** Aggregate mass is binned by
[`sieve_class`](@ref) and divided by the sample mass `ρ_b·soil_volume_mm3`. The
model predicts aggregates, not the matrix's own particle-size distribution, so
reporting classes as a share of the whole sample requires that distribution to
be **supplied** as `mineral_class_fractions` — `k+1` fractions summing to 1,
ascending. The unaggregated remainder `(1 − f_agg)` is distributed across the
classes with those fractions, after which the columns sum to 100 %.

Supply nothing and the class columns hold aggregate mass alone; the continuous
outputs (`MWD_agg_only`, `f_agg`, `POM_mass_frac`) need no such input and are
always produced.

**The classes are the assay's, not the textbook's.** `mineral_class_fractions`
is the unaggregated matrix binned by **this experiment's sieve series**, not by
sand/silt/clay. The algorithm is universal; only the cutoffs are
problem-specific. Mapping a measured texture onto a particular sieve series is
the experiment's job, not this function's.

**Why the remainder keeps the bulk distribution.** Given no information about
which particles a growing aggregate takes up, the least-assuming closure is that
the matrix is well mixed and material is incorporated in proportion to its
relative abundance — so the remainder's composition is unchanged and only its
total falls. This is a statement about what we do **not** track, not a claim
about the physics of uptake.

**Where that closure breaks, and how to see it.** The model is a continuum: a
homogeneous shell of bulk-density matrix growing around a POM core. The
class accounting is discrete: it speaks of particles of named sizes. The two are
only compatible while the shell is much thicker than the particles it is said to
contain — a 0.6 mm aggregate cannot engulf a 0.5 mm grain. `shell_thickness_mm`
is returned so this can be checked every run: **compare it against the coarse
end of `sieve_sizes`.** Where the shell is comparable to or thinner than the
coarse classes, the composition of the remainder is not defensible in detail,
even though the total mass balance still holds exactly.

No size-restricted uptake rule is applied, deliberately. Restricting absorption
to particles finer than the aggregate would mix a discrete particle model into a
continuum one, which is the incoherence above rather than a fix for it.

**Fixed-weight MWD.** Given `class_nominal_mm`, one nominal diameter per class:

    MWD = Σₖ nominalₖ · class_pctₖ / 100

This is the form of De Gryze et al. (2006) eq. (1). It is a fixed-weight sum
over class shares, not a mean of measured diameters: it saturates at the largest
nominal and is bounded below by the smallest, whatever the aggregates do.

**Per-class occupancy.** The bulk formula subtracts total shell mass from total
sample mass, which is exactly equivalent to subtracting per particle and summing
over classes — the packing cells tile the soil, so `Σ Nᵢ·ρ_b·V_pack,ᵢ = ρ_b·V`.
The equivalence is why no per-class loop is needed. But it also means a single
size class overrunning its own share is silently absorbed by the others. Pass
`cell_volume_mm3` (from `pom_population`) and `cell_occupancy` exposes that.

**No ω anywhere in this function.** `Σ nᵢ·CO₂ᵢ` and the mass sums are built from
physical particle counts and physical diameters, so they are already per-sample
totals. `ω` corrects background-carbon *concentration* inside an oversized
domain and is applied once, at initialization. See §5b.

# Diagnostics
`f_agg_vol` above ~0.6 means the aggregates interpenetrate and the
independent-domain assumption has broken down. Treat that as a signal to
revisit `r_agg` or the particle number density, not as a valid volume fraction.
"""
function population_statistics(D_agg::AbstractMatrix, POM::AbstractMatrix,
                               CO2_cum::AbstractMatrix, CO2_flux::AbstractMatrix,
                               r_0::AbstractVector, n_dist::AbstractVector;
                               sieve_sizes::AbstractVector = Float64[],
                               mineral_class_fractions::Union{Nothing,AbstractVector} = nothing,
                               class_nominal_mm::Union{Nothing,AbstractVector} = nothing,
                               cell_volume_mm3::Union{Nothing,AbstractVector} = nothing,
                               soil_volume_mm3::Real = 1.0e6,
                               ρ_b::Real, f_C_POM::Real)

    n_sizes, n_times = size(D_agg)

    # --- Shape validation -------------------------------------------------
    size(POM)      == (n_sizes, n_times) || throw(DimensionMismatch("POM is $(size(POM)), expected $((n_sizes, n_times))"))
    size(CO2_cum)  == (n_sizes, n_times) || throw(DimensionMismatch("CO2_cum is $(size(CO2_cum)), expected $((n_sizes, n_times))"))
    size(CO2_flux) == (n_sizes, n_times) || throw(DimensionMismatch("CO2_flux is $(size(CO2_flux)), expected $((n_sizes, n_times))"))
    length(r_0)    == n_sizes || throw(DimensionMismatch("r_0 has $(length(r_0)) entries, expected $(n_sizes)"))
    length(n_dist) == n_sizes || throw(DimensionMismatch("n_dist has $(length(n_dist)) entries, expected $(n_sizes)"))
    cell_volume_mm3 === nothing || length(cell_volume_mm3) == n_sizes || throw(DimensionMismatch("cell_volume_mm3 has $(length(cell_volume_mm3)) entries, expected $(n_sizes)"))
    ρ_b > 0 || throw(ArgumentError("ρ_b must be > 0 (got $(ρ_b))"))
    soil_volume_mm3 > 0 || throw(ArgumentError("soil_volume_mm3 must be > 0 (got $(soil_volume_mm3))"))
    issorted(sieve_sizes) || throw(ArgumentError("sieve_sizes must be ascending (got $(sieve_sizes))"))

    n_class = length(sieve_sizes) + 1
    report_classes = mineral_class_fractions !== nothing
    if report_classes
        length(mineral_class_fractions) == n_class || throw(ArgumentError("mineral_class_fractions must have length(sieve_sizes)+1 = $(n_class), got $(length(mineral_class_fractions))"))
        isapprox(sum(mineral_class_fractions), 1.0; atol = 1e-6) || throw(ArgumentError("mineral_class_fractions must sum to 1 (got $(sum(mineral_class_fractions)))"))
    end
    if class_nominal_mm !== nothing
        report_classes || throw(ArgumentError("class_nominal_mm requires mineral_class_fractions: without the mineral distribution the class shares are not whole-sample shares, and a fixed-weight sum over them is not comparable with a published MWD"))
        length(class_nominal_mm) == n_class || throw(ArgumentError("class_nominal_mm must have length(sieve_sizes)+1 = $(n_class), got $(length(class_nominal_mm))"))
    end

    # --- Allocate ---------------------------------------------------------
    MWD_agg    = Vector{Float64}(undef, n_times)
    f_agg_vec  = Vector{Float64}(undef, n_times)
    f_vol_vec  = Vector{Float64}(undef, n_times)
    core_frac  = Vector{Float64}(undef, n_times)
    shell_mm   = Vector{Float64}(undef, n_times)
    CO2_total  = Vector{Float64}(undef, n_times)
    flux_total = Vector{Float64}(undef, n_times)
    class_pct  = Matrix{Float64}(undef, n_class, n_times)
    occupancy  = fill(NaN, n_sizes, n_times)
    MWD_fixed  = Vector{Float64}(undef, n_times)

    soil_mass = ρ_b * soil_volume_mm3      # dry mass of the reference sample [µg]

    m_agg  = Vector{Float64}(undef, n_sizes)
    m_core = Vector{Float64}(undef, n_sizes)
    cls    = Vector{Float64}(undef, n_class)

    for t in 1:n_times
        total_mass  = 0.0
        total_vol   = 0.0
        total_core  = 0.0
        sum_mD      = 0.0
        sum_m_shell = 0.0
        fill!(cls, 0.0)

        for i in 1:n_sizes
            D          = D_agg[i, t]
            m_core[i]  = POM[i, t] / f_C_POM
            m_agg[i]   = n_dist[i] * aggregate_mass(D, r_0[i], POM[i, t];
                                                    ρ_b = ρ_b, f_C_POM = f_C_POM)
            total_mass += m_agg[i]
            total_core += n_dist[i] * m_core[i]
            total_vol  += n_dist[i] * (4.0/3.0) * π * (D / 2.0)^3
            sum_mD     += m_agg[i] * D
            sum_m_shell += m_agg[i] * max(D / 2.0 - r_0[i], 0.0)
            if cell_volume_mm3 !== nothing
                occupancy[i, t] = (4.0/3.0) * π * (D / 2.0)^3 / cell_volume_mm3[i]
            end
            cls[sieve_class(D, sieve_sizes)] += m_agg[i]
        end

        MWD_agg[t]   = total_mass > 0.0 ? sum_mD / total_mass : 0.0
        core_frac[t] = total_mass > 0.0 ? total_core / total_mass : 0.0
        shell_mm[t]  = total_mass > 0.0 ? sum_m_shell / total_mass : 0.0
        f_vol_vec[t] = total_vol / soil_volume_mm3

        # REPORTED UNCLAMPED. f_agg > 1 means the aggregates outweigh the
        # sample — a model-invalid state, not something to hide behind a clamp.
        # Only the mineral top-up below is clamped, so the class columns stay
        # non-negative and the disagreement stays visible in `f_agg` itself.
        f_agg = total_mass / soil_mass
        f_agg_vec[t] = f_agg
        if f_agg > 1.0
            @warn "population_statistics: aggregate mass exceeds sample mass" f_agg time_index=t maxlog=1
        end

        cls ./= soil_mass
        if report_classes
            f_matrix = max(1.0 - f_agg, 0.0)
            for k in 1:n_class
                cls[k] += f_matrix * mineral_class_fractions[k]
            end
        end
        class_pct[:, t] = 100.0 .* cls

        MWD_fixed[t] = class_nominal_mm === nothing ? NaN :
            sum(class_nominal_mm[k] * class_pct[k, t] for k in 1:n_class) / 100.0

        CO2_total[t]  = sum(n_dist[i] * CO2_cum[i, t]  for i in 1:n_sizes)
        flux_total[t] = sum(n_dist[i] * CO2_flux[i, t] for i in 1:n_sizes)
    end

    return (MWD_agg_only     = MWD_agg,
            f_agg            = f_agg_vec,
            f_agg_vol        = f_vol_vec,
            POM_mass_frac    = core_frac,
            shell_thickness_mm = shell_mm,
            CO2_total        = CO2_total,
            CO2_flux_total   = flux_total,
            class_pct        = class_pct,
            cell_occupancy   = occupancy,
            MWD_fixed_weight = MWD_fixed)
end
