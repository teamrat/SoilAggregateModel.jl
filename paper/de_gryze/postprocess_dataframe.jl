# postprocess_dataframe.jl
# DataFrame generators and population-sweep driver for simulation results
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
    n = length(pools.t)

    # --- CO₂ flux (backward finite difference) ---
    flux = co2_flux(result)

    # --- Total fungi (convenience) ---
    F_total_vec  = pools.F_i_total .+ pools.F_n_total .+ pools.F_m_total
    F_agg_vec    = pools.F_i_agg   .+ pools.F_n_agg   .+ pools.F_m_agg

    # --- Total system carbon and conservation error ---
    C_total_all = Vector{Float64}(undef, n)
    for i in 1:n
        C_total_all[i] = pools.P[i] + pools.CO2[i] +
                         pools.C_total[i] + pools.B_total[i] +
                         pools.F_i_total[i] + pools.F_n_total[i] + pools.F_m_total[i] +
                         pools.E_total[i] + pools.M_total[i]
    end
    C_initial = C_total_all[1]
    C_balance_error = (C_total_all .- C_initial) ./ C_initial

    # --- Normalization (per mille of initial total C) ---
    # Denominator: all carbon at t=0 (excluding CO₂ which starts at 0)
    C_init_no_co2 = pools.P[1] + pools.C_total[1] + pools.B_total[1] +
                    pools.F_i_total[1] + pools.F_n_total[1] + pools.F_m_total[1] +
                    pools.E_total[1] + pools.M_total[1]

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
                            output_dir::String="")

    all_timeseries = DataFrame[]
    all_snapshots  = DataFrame[]
    save_snaps = !isempty(snap_times)

    println("Running $(length(diam_all)) simulations...")
    if save_snaps
        println("  Spatial snapshots at t = $(snap_times)")
    end
    println("-"^72)

    total_start = time()

    for (idx, diam) in enumerate(diam_all)
        elapsed_start = time()

        # Geometry: POM radius, domain extent, initial mass
        r_0   = diam / 2.0
        r_max = diam * domain_factor / 2.0
        P_0   = (4.0/3.0) * π * r_0^3 * ρ_POM

        # Run single-aggregate simulation
        result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, t_max);
                               n_grid=n_grid, r_0=r_0, r_max=r_max,
                               ic=ic, P_0=P_0, ω=ω,
                               dt_max=dt_max, dt_min=dt_min,
                               output_times=output_times)

        elapsed = time() - elapsed_start
        n_steps = result.diagnostics["n_steps"]

        # Postprocess into DataFrame
        df = result_to_dataframe(result)
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
        println("  [$(lpad(idx, 2))/$(length(diam_all))]  D=$(Printf_diam) mm  " *
                "P₀=$(round(P_0, digits=3)) µg-C  " *
                "steps=$(n_steps)  " *
                "Δt=$(round(elapsed, digits=1))s  " *
                "err=$(round(max_error, sigdigits=2))")
    end

    total_elapsed = time() - total_start
    println("-"^72)
    println("Total wall time: $(round(total_elapsed, digits=1)) s")
    println()

    # Return snapshots as second value if requested
    if save_snaps
        return vcat(all_timeseries...), vcat(all_snapshots...)
    else
        return vcat(all_timeseries...)
    end
end


# ============================================================
# Population-level outputs
# ============================================================
# These functions take the combined DataFrame from run_diameter_sweep()
# and a POM number distribution to produce population-weighted quantities
# that can be compared directly with experimental data.

"""
    population_outputs(df_sweep::DataFrame, n_dist::Vector{Float64};
                       ω::Real=1.0,
                       sieve_sizes::Vector{Float64}=Float64[]) -> DataFrame

Compute population-weighted outputs from a diameter sweep.

Given the per-diameter time-series from `run_diameter_sweep()` and a POM
number distribution, computes:
  - Mean Weight Diameter (MWD)
  - Total CO₂ flux and cumulative CO₂ (weighted by POM count, corrected for overlap)
  - Wet aggregate stability (fraction of soil volume in aggregates above
    specified sieve sizes)

# Arguments
- `df_sweep::DataFrame`: Output from `run_diameter_sweep()`.
  Must contain columns: `diam_mm`, `time_days`, `D_agg_mm`, `CO2_cumulative`, `CO2_flux`
- `n_dist::Vector{Float64}`: Number of POM particles per size class.
  Length must match the number of unique diameters in `df_sweep`.
  This is a *number* distribution, not mass — each element is the count
  of POM fragments of that diameter class.

# Keyword Arguments
- `ω::Real`: Domain overlap correction factor (default: 1.0 = no overlap).
  When model domains are larger than packing cells (f_domain > f_pack),
  ω = (f_domain/f_pack)³ corrects for double-counted background SOC.
  See `domain_tessellation.md` for derivation.
- `sieve_sizes::Vector{Float64}`: Sieve apertures [mm] for wet aggregate
  stability calculation (default: empty = skip). Common values: [0.25, 0.5, 1.0, 2.0].

# Returns
- `DataFrame` with one row per output time, columns:
  - `time_days`: Time [days]
  - `MWD_mm`: Mean weight diameter [mm]
  - `CO2_total`: Population total cumulative CO₂ [μg-C] (overlap-corrected)
  - `CO2_flux_total`: Population total CO₂ flux [μg-C/day] (overlap-corrected)
  - `WAS_X_XXmm`: (optional) Fraction of aggregate volume above each sieve size [-]

# Physics

**Mean Weight Diameter:**
  MWD(t) = Σᵢ nᵢ · mᵢ · Dᵢ(t) / Σᵢ nᵢ · mᵢ

  where nᵢ = number of POM particles in size class i,
        mᵢ = mass of one aggregate of size class i ∝ D_agg,i(t)³,
        Dᵢ(t) = aggregate diameter at time t.

  This is the mass-weighted mean of aggregate diameters — the standard
  experimental measure from wet sieving.

**Wet Aggregate Stability (WAS):**
  WAS(s, t) = Σᵢ nᵢ · Vᵢ(t) · 𝟙[Dᵢ(t) > s] / Σᵢ nᵢ · Vᵢ(t)

  where Vᵢ(t) = (π/6)·Dᵢ(t)³ is the aggregate volume,
        s is the sieve aperture, and 𝟙 is the indicator function.

  MWD and WAS are mass-fraction-based and do not require ω correction.

**Total CO₂:**
  CO₂_total(t) = Σᵢ nᵢ · CO₂ᵢ(t) / ω

  Count-weighted sum of per-aggregate cumulative respiration, corrected
  for domain overlap when ω > 1.

# Example
```julia
# Domain tessellation with overlap correction
f_pack = 5.0   # from experimental POM input
f_domain = 10.0  # minimum for numerical resolution
ω = (f_domain / f_pack)^3  # = 8.0

df_pop = population_outputs(df_sweep, N_POM; ω=ω, sieve_sizes=[0.25, 1.0])
```
"""
function population_outputs(df_sweep::DataFrame, n_dist::Vector{Float64};
                            ω::Real=1.0,
                            sieve_sizes::Vector{Float64}=Float64[])

    # --- Validate inputs ---
    diams = sort(unique(df_sweep.diam_mm))
    n_sizes = length(diams)
    @assert length(n_dist) == n_sizes "n_dist has $(length(n_dist)) entries but sweep has $(n_sizes) diameters"

    # Get the common time grid from the first diameter
    times = df_sweep[df_sweep.diam_mm .== diams[1], :time_days]
    n_times = length(times)

    # --- Extract per-diameter time-series into arrays ---
    # D_agg[i, t], CO2[i, t], CO2_flux[i, t]
    D_agg    = Matrix{Float64}(undef, n_sizes, n_times)
    CO2_cum  = Matrix{Float64}(undef, n_sizes, n_times)
    CO2_flux = Matrix{Float64}(undef, n_sizes, n_times)

    for (i, d) in enumerate(diams)
        mask = df_sweep.diam_mm .== d
        D_agg[i, :]    = df_sweep[mask, :D_agg_mm]
        CO2_cum[i, :]  = df_sweep[mask, :CO2_cumulative]
        CO2_flux[i, :] = df_sweep[mask, :CO2_flux]
    end

    # --- Compute population outputs at each time ---
    MWD        = Vector{Float64}(undef, n_times)
    CO2_total  = Vector{Float64}(undef, n_times)
    flux_total = Vector{Float64}(undef, n_times)

    # Preallocate WAS columns
    WAS = Matrix{Float64}(undef, length(sieve_sizes), n_times)

    for t in 1:n_times
        # Aggregate volumes: V_i = (π/6) · D_i³
        V = [(π / 6.0) * D_agg[i, t]^3 for i in 1:n_sizes]

        # Mass of each size class (proportional to n_i × V_i)
        mass = [n_dist[i] * V[i] for i in 1:n_sizes]
        total_mass = sum(mass)

        # MWD = Σ(n_i · V_i · D_i) / Σ(n_i · V_i)
        if total_mass > 0.0
            MWD[t] = sum(mass[i] * D_agg[i, t] for i in 1:n_sizes) / total_mass
        else
            MWD[t] = 0.0
        end

        # Total CO₂ = Σ(n_i · CO₂_i) / ω  (overlap correction)
        CO2_total[t]  = sum(n_dist[i] * CO2_cum[i, t]  for i in 1:n_sizes)
        flux_total[t] = sum(n_dist[i] * CO2_flux[i, t] for i in 1:n_sizes)

        # WAS for each sieve size
        for (s_idx, sieve) in enumerate(sieve_sizes)
            # Volume above sieve: sum of n_i × V_i where D_i > sieve
            # init=0.0 handles case where no aggregates exceed the sieve size
            vol_above = sum(mass[i] for i in 1:n_sizes if D_agg[i, t] > sieve; init=0.0)
            WAS[s_idx, t] = total_mass > 0.0 ? vol_above / total_mass : 0.0
        end
    end

    # --- Build output DataFrame ---
    df = DataFrame(
        time_days      = times,
        MWD_mm         = MWD,
        CO2_total      = CO2_total,
        CO2_flux_total = flux_total,
    )

    # Add WAS columns with descriptive names: WAS_0_25mm, WAS_1_00mm, etc.
    for (s_idx, sieve) in enumerate(sieve_sizes)
        col_name = Symbol("WAS_$(replace(string(round(sieve, digits=2)), '.' => '_'))mm")
        df[!, col_name] = WAS[s_idx, :]
    end

    return df
end
