# postprocess_dataframe.jl
# DataFrame generators and population-sweep driver for simulation results.
#
# This file is PACKAGING, not model logic. It turns SimulationResult objects and
# sweep results into tables. The quantities themselves are computed in the
# package:
#
#   integrated_pools, co2_flux, radial_profiles  -> src/postprocessing/
#   population_statistics, aggregate_mass,       -> src/postprocessing/population.jl
#   sieve_class
#   domain_tessellation, pom_population          -> src/physics/tessellation.jl
#
# It lives outside src/ only because the package does not depend on DataFrames
# or CSV. Nothing that a second experiment could reuse belongs here — if a new
# quantity is being defined rather than tabulated, it goes in src/.
using DataFrames

"""
    result_to_dataframe(result::SimulationResult) -> DataFrame

Convert a SimulationResult into a comprehensive time-series DataFrame.

This is the primary postprocessing function. It computes everything needed
for analysis, plotting, calibration, and comparison with data — so that
scripts never need to reimplement integration or normalization logic.

# Columns returned

**Time and geometry:**
- `time_days`: Simulation time [days]
- `r_agg_mm`: Aggregate radius [mm] (from binding criterion)
- `D_agg_mm`: Aggregate diameter [mm] (= 2 × r_agg)

**Total domain pools** (integrated r_0 → r_max) [μg-C]:
- `POM`, `C_total`, `B_total`, `Fi_total`, `Fn_total`, `Fm_total`,
  `F_total`, `E_total`, `M_total`, `CO2_cumulative`

**Aggregate domain pools** (integrated r_0 → r_agg) [μg-C]:
- `C_agg`, `B_agg`, `Fi_agg`, `Fn_agg`, `Fm_agg`, `F_agg`,
  `E_agg`, `M_agg`

**Derived quantities:**
- `CO2_flux`: Instantaneous CO₂ rate [μg-C/day] (backward finite difference)
- `C_total_all`: Total system carbon [μg-C] (POM + all pools + CO₂)
- `C_balance_error`: Relative conservation error [-]

**Normalized pools** (per mille of initial total C):
- `POM_permille`, `C_permille`, `B_permille`, `F_permille`,
  `E_permille`, `M_permille`, `CO2_permille`

# Arguments
- `result::SimulationResult`: Output from `run_aggregate()`

# Returns
- `DataFrame`: One row per output time

# Notes
- Calls `integrated_pools()` internally (computes r_agg, both domain integrations)
- CO₂ flux via backward finite difference on CO2_cumulative
- Normalization denominator: total C at t=0 (POM + all spatially integrated pools)
- This is the only function scripts need for standard analysis

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 60.0); output_times=0:1:60)
df = result_to_dataframe(result)
CSV.write("output.csv", df)

# Everything is ready for analysis:
plot(df.time_days, df.CO2_permille)       # cumulative respiration
plot(df.time_days, df.D_agg_mm)           # aggregate growth
filter(row -> row.time_days == 21, df)    # extract day 21 for data comparison
```
"""
function result_to_dataframe(result::SimulationResult)
    # --- Core integration (does both total and aggregate domain) ---
    pools = integrated_pools(result)

    # --- CO₂ flux (backward finite difference) ---
    flux = co2_flux(result)

    # --- Total fungi (convenience) ---
    F_total_vec  = pools.F_i_total .+ pools.F_n_total .+ pools.F_m_total
    F_agg_vec    = pools.F_i_agg   .+ pools.F_n_agg   .+ pools.F_m_agg

    # --- Total system carbon and conservation error ---
    bal = carbon_balance_table(pools)
    C_total_all     = bal.C_total
    C_balance_error = bal.relative_error

    # --- Normalization (per mille of initial total C) ---
    # Denominator: all carbon at t=0 (excluding CO₂ which starts at 0)
    C_init_no_co2 = total_system_carbon(pools, 1; include_co2 = false)

    permille = x -> x ./ C_init_no_co2 .* 1000.0

    # --- Build DataFrame ---
    return DataFrame(
        # Time and geometry
        time_days        = pools.t,
        r_agg_mm         = pools.r_agg,
        D_agg_mm         = 2.0 .* pools.r_agg,

        # Total domain pools [μg-C]
        POM              = pools.P,
        C_total          = pools.C_total,
        B_total          = pools.B_total,
        Fi_total         = pools.F_i_total,
        Fn_total         = pools.F_n_total,
        Fm_total         = pools.F_m_total,
        F_total          = F_total_vec,
        E_total          = pools.E_total,
        M_total          = pools.M_total,
        CO2_cumulative   = pools.CO2,

        # Aggregate domain pools [μg-C]
        C_agg            = pools.C_agg,
        B_agg            = pools.B_agg,
        Fi_agg           = pools.F_i_agg,
        Fn_agg           = pools.F_n_agg,
        Fm_agg           = pools.F_m_agg,
        F_agg            = F_agg_vec,
        E_agg            = pools.E_agg,
        M_agg            = pools.M_agg,

        # Derived
        CO2_flux         = flux,
        C_system_total   = C_total_all,
        C_balance_error  = C_balance_error,

        # Normalized (per mille of initial total C)
        POM_permille     = permille(pools.P),
        C_permille       = permille(pools.C_total),
        B_permille       = permille(pools.B_total),
        F_permille       = permille(F_total_vec),
        E_permille       = permille(pools.E_total),
        M_permille       = permille(pools.M_total),
        CO2_permille     = permille(pools.CO2),
    )
end


"""
    spatial_snapshots(result::SimulationResult; times::Vector{<:Real}) -> DataFrame

Extract spatial profiles at specified times into a long-format DataFrame.

Each row is one (time, radius) pair with all state variables. This format
is convenient for faceted plotting and filtering in R or Python.

# Columns returned
- `time_days`: Snapshot time [days]
- `radius_mm`: Radial position [mm]
- `node`: Grid node index (1-based)
- `C`, `B`, `F_i`, `F_n`, `F_m`, `E`, `M`, `O`: Concentrations [μg/mm³]
- `P`: POM mass at this time [μg-C] (scalar, repeated for convenience)
- `CO2`: Cumulative CO₂ at this time [μg-C] (scalar, repeated)

# Arguments
- `result::SimulationResult`: Output from `run_aggregate()`
- `times::Vector{<:Real}`: Times [days] at which to extract profiles.
  Each is matched to the nearest available output snapshot.

# Returns
- `DataFrame`: Long format, n_times × n_grid rows

# Notes
- Requested times are snapped to nearest output; duplicates removed
- Can produce large DataFrames (e.g., 10 times × 200 nodes = 2000 rows)
- Use sparingly — select a handful of diagnostic times, not all outputs

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 60.0))
df_sp = spatial_snapshots(result; times=[0.0, 7.0, 14.0, 21.0, 60.0])
CSV.write("spatial.csv", df_sp)

# In R: ggplot(df_sp, aes(radius_mm, C, color=factor(time_days))) + geom_line()
```
"""
function spatial_snapshots(result::SimulationResult; times::Vector{<:Real})
    profiles = radial_profiles(result; times=times)
    r = result.grid.r_grid
    n_nodes = result.grid.n

    rows = DataFrame[]
    for prof in profiles
        push!(rows, DataFrame(
            time_days  = fill(prof.t, n_nodes),
            radius_mm  = prof.r,
            node       = 1:n_nodes,
            C          = prof.C,
            B          = prof.B,
            F_i        = prof.F_i,
            F_n        = prof.F_n,
            F_m        = prof.F_m,
            E          = prof.E,
            M          = prof.M,
            O          = prof.O,
            P          = fill(prof.P, n_nodes),
            CO2        = fill(prof.CO2, n_nodes),
        ))
    end

    return vcat(rows...)
end


"""
    run_diameter_sweep(diam_all, bio, soil, T_func, ψ_func, O2_func;
                       t_max, output_times, n_grid=200,
                       domain_factor=25.0, ρ_POM=50.0,
                       dt_max=0.1, dt_min=1e-4) -> DataFrame

Run a population of single-aggregate simulations across multiple POM diameters.

Each POM diameter defines one aggregate lifecycle:
  - r_0    = diam / 2           (POM radius)
  - r_max  = diam × domain_factor / 2  (domain extent)
  - P_0    = (4/3)π r₀³ × ρ_POM       (initial POM mass from sphere geometry)

Results are postprocessed via `result_to_dataframe()` and concatenated
into one long DataFrame with diameter metadata prepended as the leading columns.

# Arguments
- `diam_all::Vector`: POM diameters [mm]
- `bio::BiologicalProperties`: Biological parameters (shared across all sizes)
- `soil::SoilProperties`: Soil parameters (shared across all sizes)
- `T_func`: Temperature forcing T(t) [K]
- `ψ_func`: Water potential forcing ψ(t) [kPa]
- `O2_func`: Ambient O₂ forcing O2(t) [μg/mm³]

# Keyword Arguments
- `t_max::Real`: Simulation duration [days]
- `output_times::Vector`: Times at which to save output [days]
- `n_grid::Int`: Radial grid points per aggregate (default: 200)
- `domain_factor::Real`: r_max = diam × domain_factor / 2 (default: 25.0)
- `ρ_POM::Real`: POM carbon density [μg-C/mm³] (default: 50.0)
- `dt_max::Real`: Maximum timestep [days] (default: 0.1)
- `dt_min::Real`: Minimum timestep [days] (default: 1e-4)

# Returns
- `DataFrame`: All diameters concatenated, with leading columns:
  `diam_mm`, `r_0_mm`, `r_max_mm`, `P_0_ugC`, `n_steps`, `wall_time`,
  followed by all columns from `result_to_dataframe()`

# Example
```julia
diam_all = [0.1, 0.5, 1.0, 2.0]
bio = BiologicalProperties()
soil = SoilProperties()

df = run_diameter_sweep(diam_all, bio, soil,
                        t -> 293.15, t -> -29.0, t -> 0.278;
                        t_max=60.0, output_times=collect(0:1:60))

# Filter to one diameter:
df_1mm = filter(row -> row.diam_mm == 1.0, df)

# All diameters at day 21:
df_day21 = filter(row -> row.time_days ≈ 21.0, df)
```
"""
function run_diameter_sweep(diam_all, bio, soil, T_func, ψ_func, O2_func;
                            t_max, output_times,
                            n_grid::Int=200, domain_factor::Real=25.0,
                            ρ_POM::Real=50.0,
                            dt_max::Real=0.1, dt_min::Real=1e-4,
                            ic=nothing,
                            ω::Real=1.0,
                            snap_times::Vector{Float64}=Float64[],
                            solver::Symbol=:stiff,
                            quiet::Bool=false,
                            output_dir::String="")

    solver in (:stiff, :split) || throw(ArgumentError(
        "solver must be :stiff (default) or :split (reference implementation), got $(solver)"))

    all_timeseries = DataFrame[]
    all_snapshots  = DataFrame[]
    save_snaps = !isempty(snap_times)

    if !quiet
        println("Running $(length(diam_all)) simulations...")
        save_snaps && println("  Spatial snapshots at t = $(snap_times)")
        println("-"^72)
    end

    total_start = time()

    for (idx, diam) in enumerate(diam_all)
        elapsed_start = time()

        # Geometry: POM radius, domain extent, initial mass
        r_0   = diam / 2.0
        r_max = diam * domain_factor / 2.0
        P_0   = (4.0/3.0) * π * r_0^3 * ρ_POM

        # Run single-aggregate simulation. :stiff is the default workhorse;
        # :split is the reference implementation kept for cross-checking and
        # for its independent carbon-closure probe (REFERENCE.md §17a, §20a).
        result = if solver === :stiff
            run_aggregate_stiff(bio, soil, T_func, ψ_func, O2_func, (0.0, t_max);
                                n_grid=n_grid, r_0=r_0, r_max=r_max,
                                ic=ic, P_0=P_0, ω=ω,
                                output_times=output_times)
        else
            run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, t_max);
                          n_grid=n_grid, r_0=r_0, r_max=r_max,
                          ic=ic, P_0=P_0, ω=ω,
                          dt_max=dt_max, dt_min=dt_min,
                          output_times=output_times)
        end

        elapsed = time() - elapsed_start
        n_steps = get(result.diagnostics, "n_accept", result.diagnostics["n_steps"])

        # Postprocess into DataFrame
        df = result_to_dataframe(result)
        # NaN for :stiff by design — carbon closure there is structural, not
        # measured. REFERENCE.md §17a.
        max_error = maximum(abs.(df.C_balance_error))

        # Prepend metadata columns so they appear first in the output.
        # insertcols!(df, pos, pairs...) inserts columns at the given
        # position; a scalar value is broadcast to fill all rows.
        insertcols!(df, 1,
            :diam_mm   => diam,
            :r_0_mm    => r_0,
            :r_max_mm  => r_max,
            :P_0_ugC   => P_0,
            :n_steps   => n_steps,
            :wall_time => elapsed
        )

        push!(all_timeseries, df)

        # Spatial snapshots (if requested)
        if save_snaps
            df_sp = spatial_snapshots(result; times=snap_times)
            insertcols!(df_sp, 1, :diam_mm => diam)
            push!(all_snapshots, df_sp)

            # Save per-diameter CSV if output_dir provided
            if !isempty(output_dir)
                sp_dir = joinpath(output_dir, "spatial_profiles")
                isdir(sp_dir) || mkpath(sp_dir)
                fname = "spatial_$(lpad(string(idx), 3, '0')).csv"
                CSV.write(joinpath(sp_dir, fname), df_sp)
            end
        end

        # Progress
        Printf_diam = lpad(string(round(diam, digits=3)), 6)
        quiet || println("  [$(lpad(idx, 2))/$(length(diam_all))]  D=$(Printf_diam) mm  " *
                "P₀=$(round(P_0, digits=3)) µg-C  " *
                "steps=$(n_steps)  " *
                "Δt=$(round(elapsed, digits=1))s  " *
                "err=$(round(max_error, sigdigits=2))")
    end

    total_elapsed = time() - total_start
    if !quiet
        println("-"^72)
        println("Total wall time: $(round(total_elapsed, digits=1)) s")
        println()
    end

    # Return snapshots as second value if requested
    if save_snaps
        return vcat(all_timeseries...), vcat(all_snapshots...)
    else
        return vcat(all_timeseries...)
    end
end


# ============================================================
# Population-level outputs — DataFrame adapter
# ============================================================
# The arithmetic lives in src/postprocessing/population.jl. What is left here is
# packaging: pull the per-diameter time series out of the sweep DataFrame, hand
# plain matrices to `population_statistics`, and name the resulting columns.

"""
    population_outputs(df_sweep::DataFrame, n_dist::Vector{Float64};
                       sieve_sizes, mineral_class_fractions,
                       class_nominal_mm, class_labels,
                       soil_volume_mm3, ρ_b, f_C_POM) -> DataFrame

DataFrame wrapper around [`population_statistics`](@ref).

This function does no physics. It extracts `D_agg_mm`, `POM`, `CO2_cumulative`,
`CO2_flux` and `r_0_mm` from a `run_diameter_sweep()` result into matrices,
calls `population_statistics`, and builds the output table. Every definition —
aggregate mass, sieve binning, the fixed-weight MWD, the ω convention — is
documented at the `src/` function and in `docs/REFERENCE.md` §5c.

# Arguments
- `df_sweep::DataFrame`: output of `run_diameter_sweep()`. Must carry
  `diam_mm`, `time_days`, `D_agg_mm`, `POM`, `CO2_cumulative`, `CO2_flux`,
  `r_0_mm`.
- `n_dist::Vector{Float64}`: POM particle count per size class, in
  `soil_volume_mm3`. Ascending by diameter, matching `sort(unique(diam_mm))`.

# Keyword Arguments
Passed straight through to `population_statistics`, except `class_labels`:

- `sieve_sizes`: sieve apertures [mm], ascending. `k` sieves → `k+1` classes.
- `mineral_class_fractions`: the soil's own mineral mass split across those
  classes. Required for whole-sample class shares; see the `src/` docstring.
- `class_nominal_mm`: nominal diameter per class, enabling a fixed-weight MWD.
- `cell_volume_mm3`: packing-cell volume per size class, from `pom_population`
  (`V_pack`). Enables the `max_cell_occ` column.
- `class_labels`: names for the class columns, ascending. The sieve series is a
  property of the assay, so the assay names its own classes. Without labels the
  columns are `pct_class1..N`.
- `soil_volume_mm3`, `ρ_b`, `f_C_POM`: reference sample volume [mm³], bulk
  density [µg/mm³] and residue carbon fraction [-].

There is **no `ω` argument.** ω rescales background-carbon concentration inside
an oversized domain and is applied once, at initialization. The sums here are
built from physical counts and physical diameters and are already per-sample
totals; dividing by ω again would understate the amendment by that factor.

# Returns
`DataFrame`, one row per output time: `time_days`, `MWD_agg_only`, `f_agg`,
`f_agg_vol`, `POM_mass_frac`, `shell_mm`, `CO2_total`, `CO2_flux_total`, one
`pct_<label>` column per sieve class in ascending order, and
`MWD_fixed_weight_mm` when `class_nominal_mm` is supplied.

The `pct_*` columns sum to 100 when `mineral_class_fractions` is given and to
`f_agg·100` when it is not. They are non-overlapping ranges, not cumulative
"above sieve X" fractions — the running sum from the top gives those.
"""
function population_outputs(df_sweep::DataFrame, n_dist::Vector{Float64};
                            sieve_sizes::Vector{Float64}=Float64[],
                            mineral_class_fractions::Union{Nothing,Vector{Float64}}=nothing,
                            class_nominal_mm::Union{Nothing,Vector{Float64}}=nothing,
                            class_labels::Union{Nothing,Vector{String}}=nothing,
                            cell_volume_mm3::Union{Nothing,Vector{Float64}}=nothing,
                            soil_volume_mm3::Real=1.0e6,
                            ρ_b::Real=1300.0,
                            f_C_POM::Real=0.443)

    diams   = sort(unique(df_sweep.diam_mm))
    n_sizes = length(diams)
    @assert length(n_dist) == n_sizes "n_dist has $(length(n_dist)) entries but sweep has $(n_sizes) diameters"

    # Common time grid, taken from the first diameter
    times   = df_sweep[df_sweep.diam_mm .== diams[1], :time_days]
    n_times = length(times)

    n_class = length(sieve_sizes) + 1
    if class_labels !== nothing
        @assert length(class_labels) == n_class "class_labels must have length(sieve_sizes)+1 = $(n_class)"
    end

    # --- Column extraction: [size, time] matrices -------------------------
    D_agg    = Matrix{Float64}(undef, n_sizes, n_times)
    POM      = Matrix{Float64}(undef, n_sizes, n_times)
    CO2_cum  = Matrix{Float64}(undef, n_sizes, n_times)
    CO2_flux = Matrix{Float64}(undef, n_sizes, n_times)
    r_0      = Vector{Float64}(undef, n_sizes)

    for (i, d) in enumerate(diams)
        mask = df_sweep.diam_mm .== d
        D_agg[i, :]    = df_sweep[mask, :D_agg_mm]
        POM[i, :]      = df_sweep[mask, :POM]
        CO2_cum[i, :]  = df_sweep[mask, :CO2_cumulative]
        CO2_flux[i, :] = df_sweep[mask, :CO2_flux]
        r_0[i]         = first(df_sweep[mask, :r_0_mm])
    end

    # --- Physics ----------------------------------------------------------
    stats = population_statistics(D_agg, POM, CO2_cum, CO2_flux, r_0, n_dist;
                                  sieve_sizes             = sieve_sizes,
                                  mineral_class_fractions = mineral_class_fractions,
                                  class_nominal_mm        = class_nominal_mm,
                                  cell_volume_mm3         = cell_volume_mm3,
                                  soil_volume_mm3         = soil_volume_mm3,
                                  ρ_b                     = ρ_b,
                                  f_C_POM                 = f_C_POM)

    # --- Packaging --------------------------------------------------------
    df = DataFrame(
        time_days      = times,
        MWD_agg_only   = stats.MWD_agg_only,    # mass-weighted mean over aggregates
        f_agg          = stats.f_agg,           # MASS fraction of sample that is aggregate
        f_agg_vol      = stats.f_agg_vol,       # volume fraction (geometry check, keep < 0.6)
        POM_mass_frac  = stats.POM_mass_frac,   # residue share of retained coarse mass
        shell_mm       = stats.shell_thickness_mm,  # continuum check: compare vs sieve_sizes
        # Worst per-class packing-cell occupancy. >= 1 means some size class has
        # consumed more soil than it owns; the bulk mass balance would hide it.
        max_cell_occ   = [maximum(view(stats.cell_occupancy, :, t)) for t in 1:length(times)],
        CO2_total      = stats.CO2_total,
        CO2_flux_total = stats.CO2_flux_total,
    )

    for k in 1:n_class
        nm = class_labels === nothing ? "pct_class$(k)" : "pct_" * class_labels[k]
        df[!, Symbol(nm)] = stats.class_pct[k, :]
    end
    if class_nominal_mm !== nothing
        df[!, :MWD_fixed_weight_mm] = stats.MWD_fixed_weight
    end

    return df
end
