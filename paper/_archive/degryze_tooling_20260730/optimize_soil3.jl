"""
De Gryze (2006) — Parameter Optimization for Soil 3
=====================================================

Fits cumulative CO₂ data for Soil 3 using the Phase 1 calibration strategy:
  4 core parameters: R_P_max, r_B_max, μ_B, γ
  
Additional secondary parameters can be toggled on.

Strategy:
  1. Latin Hypercube Sampling (LHS) for global exploration
  2. Nelder-Mead refinement from best LHS point
  3. Optional: local sensitivity analysis at optimum

Target data (Soil 3, µg-C/g-soil):
  Day  0:     0.0
  Day  5:   481.954
  Day 10:  1263.633
  Day 15:  1766.159
  Day 21:  2139.010
"""

## ============================================================
## Setup
## ============================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Revise
using SoilAggregateModel
using DataFrames
using CSV
using Distributions
using Printf
using Random

include(joinpath(@__DIR__, "postprocess_dataframe.jl"))

Random.seed!(42)

## ============================================================
## Target data — Soil 3
## ============================================================

data_times = [0.0, 5.0, 10.0, 15.0, 21.0]   # days
data_CO2   = [0.0, 481.954, 1263.633, 1766.159, 2139.010]  # µg-C/g-soil

## ============================================================
## Fixed experimental setup (identical to run_degryze.jl)
## ============================================================

T_const  = 298.15   # 25 °C — De Gryze 2006 p.237 (was 293.15; diverged from run_degryze.jl)
ψ_const  = -29.0
O2_frac  = 0.21
M_O2     = 0.032
P_atm    = 101000.0
R_gas    = 8.314
O2_const = O2_frac * P_atm * M_O2 / (R_gas * T_const)

T_func  = t -> T_const
ψ_func  = t -> ψ_const
O2_func = t -> O2_const

# POM size distribution
POM_mean  = 1.25
POM_sigma = 0.23
bin_edges = collect(0.5:0.3:2.0)
diam_all  = [(bin_edges[i] + bin_edges[i+1]) / 2.0 for i in 1:length(bin_edges)-1]
F_POM = cdf.(Normal(POM_mean, POM_sigma), bin_edges)
f_POM = diff(F_POM)

# Domain tessellation and population — src/physics/tessellation.jl
# (was a copy of the block in run_degryze.jl; the two had already diverged)
ρ_POM   = 200.0
I_input = 4.43e-3
ρ_bulk  = 1300.0
soil_V  = 100.0^3

tess = domain_tessellation(ρ_POM=ρ_POM, I_input=I_input, ρ_b=ρ_bulk)
pop  = pom_population(diam_all, f_POM, tess; ρ_POM=ρ_POM, soil_volume_mm3=soil_V)

C_input_vol   = I_input * ρ_bulk
φ_POM         = tess.φ_POM
f_pack        = tess.f_pack
f_domain      = tess.f_domain
ω             = tess.ω
domain_factor = tess.f_domain
N_POM         = pop.N_POM
soil_mass_per_L = soil_V * ρ_bulk * 1e-6  # grams soil per liter

# Soil 3: SOC = 2.21% (Table 1 in de Gryze)
ic_soil3 = InitialConditions(
    SOC   = 0.0221,
    f_bact  = 0.01,
    f_fungi = 0.01,
    f_eps   = 0.005,
    T_0   = T_const,
    ψ_0   = ψ_const,
    O2_gas = O2_const
)

# Grid and time — use coarser grid for optimization speed
n_grid_opt = 80        # coarse for optimization (vs 200 for production)
t_max_opt  = 22.0      # only need to reach day 21
dt_output  = 0.5       # coarser output for speed
output_times_opt = collect(0.0:dt_output:t_max_opt)

## ============================================================
## Fixed soil properties (measured for de Gryze silt loam)
## ============================================================

soil_fixed = SoilProperties(
    k_L = 1000.0,
    D_B_rel = 0.00001,
    ρ_b = ρ_bulk,
    f_clay_silt = 0.74
)

## ============================================================
## Parameter space definition
## ============================================================
# Each entry: (name, lo, hi, default, log_scale?)
# log_scale=true means we sample in log space (better for rates)

PARAM_DEFS = [
    # Phase 1: CO₂ shape
    (:R_P_max,  0.05,  5.0,   0.5,  true),   # POM dissolution rate
    (:r_B_max,  0.2,   5.0,   1.0,  true),   # Bacterial max uptake
    (:μ_B,      0.003, 0.1,   0.02, true),   # Bacterial mortality
    (:γ,        0.05,  0.50,  0.2,  false),  # EPS allocation fraction

    # Phase 1b: secondary CO₂ controls (toggle on as needed)
     (:r_F_max,  0.2,   5.0,   1.0,  true),
     (:Y_B_max,  0.3,   0.7,   0.5,  false),
     (:B_S,      0.05,  5.0,   0.5,  true),
]

# Extract active parameter info
param_names = [p[1] for p in PARAM_DEFS]
param_lo    = [p[2] for p in PARAM_DEFS]
param_hi    = [p[3] for p in PARAM_DEFS]
param_def   = [p[4] for p in PARAM_DEFS]
param_log   = [p[5] for p in PARAM_DEFS]
n_params    = length(PARAM_DEFS)

println("="^72)
println("Parameter Optimization — Soil 3 CO₂")
println("="^72)
println("Active parameters: $(n_params)")
for (i, p) in enumerate(PARAM_DEFS)
    scale = p[5] ? "log" : "lin"
    @printf("  [%d] %-10s  range [%.4g, %.4g]  default=%.4g  (%s)\n",
            i, p[1], p[2], p[3], p[4], scale)
end
println()

## ============================================================
## Helper: encode/decode parameter vectors in [0,1]^n
## ============================================================

function decode_params(x::Vector{Float64})
    """Map from [0,1]^n to physical parameter values."""
    θ = Vector{Float64}(undef, n_params)
    for i in 1:n_params
        if param_log[i]
            log_lo = log10(param_lo[i])
            log_hi = log10(param_hi[i])
            θ[i] = 10.0^(log_lo + x[i] * (log_hi - log_lo))
        else
            θ[i] = param_lo[i] + x[i] * (param_hi[i] - param_lo[i])
        end
    end
    return θ
end

function encode_params(θ::Vector{Float64})
    """Map from physical parameter values to [0,1]^n."""
    x = Vector{Float64}(undef, n_params)
    for i in 1:n_params
        if param_log[i]
            log_lo = log10(param_lo[i])
            log_hi = log10(param_hi[i])
            x[i] = (log10(θ[i]) - log_lo) / (log_hi - log_lo)
        else
            x[i] = (θ[i] - param_lo[i]) / (param_hi[i] - param_lo[i])
        end
    end
    return x
end

## ============================================================
## Helper: build BiologicalProperties from parameter vector
## ============================================================

function make_bio(θ::Vector{Float64})
    """Construct BiologicalProperties overriding active parameters."""
    # Start from defaults (same as run_degryze.jl)
    kwargs = Dict{Symbol,Float64}(
        # Defaults from run_degryze.jl
        :κ_s_ref => 0.01,
        :κ_d_ref => 0.005,
        :F_i_min => 1e-6,
        :F_n_min => 2e-4,
        :F_m_min => 1e-6,
        :D_Fn0   => 0.00001,
        :D_Fm0   => 0.001,
        :r_B_max => 1.0,
        :r_F_max => 1.0,
        :R_P_max => 0.5,
        :C_B     => 5.0e-5,
        :μ_B     => 0.02,
        :μ_F     => 0.02,
    )

    # Override with active calibration parameters
    for (i, name) in enumerate(param_names)
        kwargs[name] = θ[i]
    end

    BiologicalProperties(; kwargs...)
end

## ============================================================
## Forward model: parameters → CO₂ at data times
## ============================================================

function forward_model(θ::Vector{Float64}; verbose::Bool=false)
    """
    Run the full population sweep and return model CO₂ at data times.
    
    Returns:
      (co2_at_data_times, wall_time, max_balance_error)
    or nothing on failure.
    """
    bio = make_bio(θ)

    t_start = time()

    local df_summary
    try
        df_summary, _ = run_diameter_sweep(
            diam_all, bio, soil_fixed, T_func, ψ_func, O2_func;
            t_max=t_max_opt, output_times=output_times_opt,
            n_grid=n_grid_opt, domain_factor=domain_factor,
            ρ_POM=ρ_POM, ic=ic_soil3, ω=ω,
            snap_times=[0.0]  # minimal snapshots
        )
    catch e
        if verbose
            println("  ⚠ Simulation failed: $(e)")
        end
        return nothing
    end

    wall = time() - t_start

    # Check conservation
    max_err = maximum(abs.(df_summary.C_balance_error))
    if max_err > 1e-6
        if verbose
            @printf("  ⚠ Conservation error: %.2e\n", max_err)
        end
        return nothing
    end

    # Population-level CO₂
    df_pop = population_outputs(df_summary, N_POM; sieve_sizes=Float64[],
                                ρ_b=ρ_bulk, f_C_POM=0.443)
    co2_model_per_g = df_pop.CO2_total ./ soil_mass_per_L

    # Interpolate to data times
    model_times = df_pop.time_days
    co2_at_data = Float64[]
    for t_d in data_times
        idx = argmin(abs.(model_times .- t_d))
        push!(co2_at_data, co2_model_per_g[idx])
    end

    return (co2_at_data, wall, max_err)
end

## ============================================================
## Cost function
## ============================================================

function cost_function(x::Vector{Float64}; verbose::Bool=false)
    """
    Compute RMSE between model and data CO₂ (µg-C/g-soil).
    
    x is in [0,1]^n (unit hypercube).
    Returns Inf on failure.
    """
    θ = decode_params(x)
    result = forward_model(θ; verbose=verbose)

    if result === nothing
        return Inf
    end

    co2_model, wall, err = result

    # RMSE over data points (skip t=0 which is trivially 0)
    residuals = co2_model[2:end] .- data_CO2[2:end]
    rmse = sqrt(sum(residuals.^2) / length(residuals))

    if verbose
        @printf("  RMSE=%.1f  wall=%.1fs  err=%.1e  params=", rmse, wall, err)
        for (i, name) in enumerate(param_names)
            @printf(" %s=%.4g", name, θ[i])
        end
        println()
    end

    return rmse
end

## ============================================================
## Phase 1: Latin Hypercube Sampling
## ============================================================

function latin_hypercube(n_samples::Int, n_dims::Int)
    """Generate LHS design in [0,1]^n_dims."""
    lhs = Matrix{Float64}(undef, n_samples, n_dims)
    for j in 1:n_dims
        perm = randperm(n_samples)
        for i in 1:n_samples
            lhs[i, j] = (perm[i] - rand()) / n_samples
        end
    end
    return lhs
end

function run_lhs(n_samples::Int)
    println("="^72)
    println("Phase 1: Latin Hypercube Sampling ($(n_samples) samples)")
    println("="^72)

    lhs = latin_hypercube(n_samples, n_params)
    results = Vector{Tuple{Float64, Vector{Float64}, Vector{Float64}}}()

    t_total = time()

    for i in 1:n_samples
        x = lhs[i, :]
        θ = decode_params(x)

        @printf("[%3d/%d] ", i, n_samples)
        for (j, name) in enumerate(param_names)
            @printf("%s=%.3g ", name, θ[j])
        end

        rmse = cost_function(x)

        if rmse == Inf
            println("→ FAILED")
        else
            @printf("→ RMSE=%.1f\n", rmse)
            push!(results, (rmse, x, θ))
        end
    end

    elapsed = time() - t_total
    n_valid = length(results)
    println()
    println("-"^72)
    @printf("LHS complete: %d/%d valid (%.1f%%), total %.1f s\n",
            n_valid, n_samples, 100*n_valid/n_samples, elapsed)

    if n_valid == 0
        println("ERROR: No valid parameter sets found!")
        return nothing
    end

    # Sort by RMSE
    sort!(results, by=x->x[1])

    println()
    println("Top 10 parameter sets:")
    println("-"^72)
    for k in 1:min(10, n_valid)
        rmse, x, θ = results[k]
        @printf("  [%2d] RMSE=%7.1f ", k, rmse)
        for (j, name) in enumerate(param_names)
            @printf(" %s=%.4g", name, θ[j])
        end
        println()
    end

    return results
end

## ============================================================
## Phase 2: Nelder-Mead refinement
## ============================================================

function nelder_mead(f, x0::Vector{Float64};
                     max_iter::Int=200, tol_f::Float64=1.0,
                     tol_x::Float64=1e-4,
                     α::Float64=1.0, β_nm::Float64=0.5,
                     γ_nm::Float64=2.0, δ::Float64=0.5,
                     initial_scale::Float64=0.08)
    """
    Nelder-Mead simplex optimization in [0,1]^n.
    
    Clamped to unit hypercube. Designed for expensive black-box functions.
    """
    n = length(x0)
    
    # Initialize simplex
    simplex = Vector{Vector{Float64}}(undef, n+1)
    simplex[1] = clamp.(x0, 0.0, 1.0)
    for i in 1:n
        p = copy(x0)
        p[i] = clamp(p[i] + initial_scale, 0.0, 1.0)
        if p[i] ≈ x0[i]  # at boundary, go the other direction
            p[i] = clamp(p[i] - 2*initial_scale, 0.0, 1.0)
        end
        simplex[i+1] = p
    end
    
    # Evaluate all vertices
    fvals = [f(s) for s in simplex]
    n_evals = n + 1
    
    println("  Nelder-Mead: initial best RMSE = $(round(minimum(fvals), digits=1))")
    
    for iter in 1:max_iter
        # Sort
        order = sortperm(fvals)
        simplex = simplex[order]
        fvals = fvals[order]
        
        f_best = fvals[1]
        f_worst = fvals[end]
        f_second_worst = fvals[end-1]
        
        # Convergence check
        f_range = f_worst - f_best
        x_range = maximum(maximum(abs.(simplex[i] .- simplex[1])) for i in 2:n+1)
        
        if f_range < tol_f && x_range < tol_x
            println("  Converged at iter $(iter): RMSE=$(round(f_best, digits=1)), " *
                    "Δf=$(round(f_range, digits=2)), Δx=$(round(x_range, sigdigits=2))")
            break
        end
        
        # Centroid (excluding worst)
        centroid = zeros(n)
        for i in 1:n
            centroid .+= simplex[i]
        end
        centroid ./= n
        
        # Reflection
        x_r = clamp.(centroid .+ α .* (centroid .- simplex[end]), 0.0, 1.0)
        f_r = f(x_r)
        n_evals += 1
        
        if f_best ≤ f_r < f_second_worst
            # Accept reflection
            simplex[end] = x_r
            fvals[end] = f_r
        elseif f_r < f_best
            # Try expansion
            x_e = clamp.(centroid .+ γ_nm .* (x_r .- centroid), 0.0, 1.0)
            f_e = f(x_e)
            n_evals += 1
            if f_e < f_r
                simplex[end] = x_e
                fvals[end] = f_e
            else
                simplex[end] = x_r
                fvals[end] = f_r
            end
        else
            # Contraction
            if f_r < f_worst
                x_c = clamp.(centroid .+ β_nm .* (x_r .- centroid), 0.0, 1.0)
            else
                x_c = clamp.(centroid .+ β_nm .* (simplex[end] .- centroid), 0.0, 1.0)
            end
            f_c = f(x_c)
            n_evals += 1
            if f_c < min(f_r, f_worst)
                simplex[end] = x_c
                fvals[end] = f_c
            else
                # Shrink toward best
                for i in 2:n+1
                    simplex[i] = clamp.(simplex[1] .+ δ .* (simplex[i] .- simplex[1]), 0.0, 1.0)
                    fvals[i] = f(simplex[i])
                    n_evals += 1
                end
            end
        end
        
        # Progress every 10 iterations
        if iter % 10 == 0
            best_θ = decode_params(simplex[sortperm(fvals)[1]])
            @printf("  iter %3d: RMSE=%.1f  evals=%d ", iter, minimum(fvals), n_evals)
            for (j, name) in enumerate(param_names)
                @printf(" %s=%.3g", name, best_θ[j])
            end
            println()
        end
    end
    
    # Return best
    best_idx = argmin(fvals)
    return simplex[best_idx], fvals[best_idx], n_evals
end

## ============================================================
## Phase 3: Diagnostics at optimum
## ============================================================

function run_diagnostics(θ_opt::Vector{Float64})
    println("="^72)
    println("Diagnostics at optimum")
    println("="^72)
    
    # Full-resolution run at optimum
    println("\nRunning production-resolution (n_grid=200) at optimum...")
    bio_opt = make_bio(θ_opt)
    
    output_times_full = collect(0.0:0.125:22.0)
    df_summary, _ = run_diameter_sweep(
        diam_all, bio_opt, soil_fixed, T_func, ψ_func, O2_func;
        t_max=22.0, output_times=output_times_full,
        n_grid=200, domain_factor=domain_factor,
        ρ_POM=ρ_POM, ic=ic_soil3, ω=ω,
        snap_times=[0.0]
    )
    
    df_pop = population_outputs(df_summary, N_POM; sieve_sizes=[0.25, 0.5, 1.0, 2.0],
                                ρ_b=ρ_bulk, f_C_POM=0.443)
    co2_model = df_pop.CO2_total ./ soil_mass_per_L
    
    # Interpolate to data times
    model_times = df_pop.time_days
    println("\nModel vs Data (Soil 3):")
    println("-"^50)
    @printf("  %6s  %12s  %12s  %10s\n", "Day", "Data", "Model", "Residual")
    println("-"^50)
    for (k, t_d) in enumerate(data_times)
        idx = argmin(abs.(model_times .- t_d))
        model_val = co2_model[idx]
        resid = model_val - data_CO2[k]
        @printf("  %6.0f  %12.1f  %12.1f  %+10.1f\n", t_d, data_CO2[k], model_val, resid)
    end
    
    # Local sensitivity: ±10% perturbation of each parameter
    println("\nLocal sensitivity (±10% perturbation):")
    println("-"^60)
    @printf("  %-10s  %10s  %10s  %10s\n", "Parameter", "RMSE(-10%)", "RMSE(base)", "RMSE(+10%)")
    println("-"^60)
    
    x_opt = encode_params(θ_opt)
    rmse_base = cost_function(x_opt)
    
    for i in 1:n_params
        θ_lo = copy(θ_opt)
        θ_hi = copy(θ_opt)
        θ_lo[i] *= 0.9
        θ_hi[i] *= 1.1
        # Clamp to bounds
        θ_lo[i] = clamp(θ_lo[i], param_lo[i], param_hi[i])
        θ_hi[i] = clamp(θ_hi[i], param_lo[i], param_hi[i])
        
        rmse_lo = cost_function(encode_params(θ_lo))
        rmse_hi = cost_function(encode_params(θ_hi))
        @printf("  %-10s  %10.1f  %10.1f  %10.1f\n",
                param_names[i], rmse_lo, rmse_base, rmse_hi)
    end
    
    return df_pop
end

## ============================================================
## Main execution
## ============================================================

println()
println("Testing forward model at defaults...")
θ_default = Float64[p[4] for p in PARAM_DEFS]
result_default = forward_model(θ_default; verbose=true)
if result_default !== nothing
    co2_def, wall_def, _ = result_default
    println("\nDefault CO₂ at data times: $(round.(co2_def, digits=1))")
    println("Target  CO₂ at data times: $(round.(data_CO2, digits=1))")
    resid = co2_def[2:end] .- data_CO2[2:end]
    println("Default RMSE: $(round(sqrt(sum(resid.^2)/length(resid)), digits=1))")
    println("Wall time per sweep: $(round(wall_def, digits=1)) s")
else
    println("Default parameters failed — check model setup")
end

## --- LHS exploration ---
println()
n_lhs = 200   # ≈ 15 per dimension; ~15 min at ~15s/eval
lhs_results = run_lhs(n_lhs)

if lhs_results !== nothing && length(lhs_results) > 0
    # Best from LHS
    best_rmse_lhs, best_x_lhs, best_θ_lhs = lhs_results[1]
    
    ## --- Nelder-Mead refinement ---
    println()
    println("="^72)
    println("Phase 2: Nelder-Mead refinement from best LHS point")
    println("="^72)
    
    x_opt, rmse_opt, n_evals = nelder_mead(cost_function, best_x_lhs;
                                            max_iter=300, tol_f=5.0, tol_x=0.005)
    θ_opt = decode_params(x_opt)
    
    println()
    println("="^72)
    println("Optimization result")
    println("="^72)
    @printf("  RMSE: %.1f µg-C/g-soil\n", rmse_opt)
    println("  Parameters:")
    for (i, name) in enumerate(param_names)
        @printf("    %-10s = %.6g\n", name, θ_opt[i])
    end
    
    ## --- Diagnostics ---
    println()
    df_pop_opt = run_diagnostics(θ_opt)
    
    ## --- Save results ---
    output_dir = joinpath(@__DIR__, "output")
    isdir(output_dir) || mkpath(output_dir)
    
    # Save optimized parameters
    open(joinpath(output_dir, "optimized_params_soil3.txt"), "w") do io
        println(io, "# De Gryze Soil 3 — Optimized Parameters")
        println(io, "# RMSE: $(round(rmse_opt, digits=1)) µg-C/g-soil")
        println(io, "#")
        for (i, name) in enumerate(param_names)
            @printf(io, "%-10s = %.6g\n", name, θ_opt[i])
        end
    end
    
    # Save LHS results for analysis
    lhs_df = DataFrame(
        rmse    = [r[1] for r in lhs_results],
    )
    for (j, name) in enumerate(param_names)
        lhs_df[!, name] = [r[3][j] for r in lhs_results]
    end
    CSV.write(joinpath(output_dir, "lhs_results_soil3.csv"), lhs_df)
    
    println()
    println("✓ Results saved to $(output_dir)/")
    println("  optimized_params_soil3.txt")
    println("  lhs_results_soil3.csv")
end

println()
println("="^72)
println("✓ Optimization complete")
println("="^72)
