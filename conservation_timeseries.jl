"""
Conservation Error Accumulation Over Time

Track how conservation error grows over 365 days to distinguish:
- Linear growth (systematic leak in physics)
- Quadratic growth (accumulation of numerical errors)
"""

using SoilAggregateModel

# Standard parameters and conditions
bio = BiologicalProperties()
soil = SoilProperties()
T_func = t -> 293.15  # K (20°C)
ψ_func = t -> -33.0   # kPa
O2_func = t -> 0.21   # 21% O₂

# Run 365-day simulation with monthly outputs
println("Running 365-day simulation with monthly outputs...")
output_times = collect(0.0:30.0:365.0)
result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 365.0);
                       output_times=output_times)

println("Simulation complete: $(result.diagnostics["n_steps"]) timesteps")
println()

# Compute carbon balance for each output time
println("="^70)
println("Conservation Error Accumulation Over Time")
println("="^70)
println()

bal = carbon_balance_table(result)
pools = integrated_pools(result)

# Print header
println("Time    | Total C    | P+Pools    | CO₂        | Error      | Error %")
println("(days)  | (μg-C)     | (μg-C)     | (μg-C)     | (μg-C)     | ")
println("-"^70)

# Print each time point
for i in 1:length(bal.t)
    t = bal.t[i]
    C_total = bal.C_total[i]
    C_initial = bal.C_initial

    # Components
    P = pools.P[i]
    CO2 = pools.CO2[i]
    C_pools = pools.C_total[i] + pools.B_total[i] + pools.F_i_total[i] +
              pools.F_n_total[i] + pools.F_m_total[i] + pools.E_total[i] + pools.M_total[i]
    P_plus_pools = P + C_pools

    error = C_total - C_initial
    rel_error_pct = bal.relative_error[i] * 100

    println("$(rpad(round(Int, t), 7)) | ",
            "$(rpad(round(C_total, digits=2), 10)) | ",
            "$(rpad(round(P_plus_pools, digits=2), 10)) | ",
            "$(rpad(round(CO2, digits=2), 10)) | ",
            "$(rpad(round(error, digits=4), 10)) | ",
            "$(round(rel_error_pct, digits=4))%")
end

println()
println("="^70)

# Analyze growth pattern
errors = bal.relative_error
times = bal.t

println()
println("Error Growth Analysis:")
println("-"^70)
println("Initial error (t=0):     $(round(errors[1]*100, digits=6))%")
println("Final error (t=365):     $(round(errors[end]*100, digits=4))%")
println()

# Fit linear and quadratic models to log10(abs(error))
if length(errors) > 3 && all(errors .!= 0)
    # Simple linear regression for error vs time
    n = length(times)
    mean_t = sum(times) / n
    mean_e = sum(errors) / n
    slope = sum((times .- mean_t) .* (errors .- mean_e)) / sum((times .- mean_t).^2)

    println("Linear growth rate: $(round(slope*365*100, digits=4))% per year")
    println()

    if slope > 0
        println("Interpretation:")
        if slope * 365 < 0.01  # < 1% per year
            println("  Very slow accumulation - likely numerical rounding")
        elseif slope * 365 < 0.05  # < 5% per year
            println("  Moderate accumulation - check integration tolerances")
        else
            println("  ⚠️  Fast accumulation - likely systematic physics leak")
        end
    end
end

println()
println("="^70)

# Component breakdown at key times
println()
println("Carbon Component Breakdown:")
println("-"^70)

# Get integrated pools
pools = integrated_pools(result)

# Print at t=0, 180d, 365d
key_indices = [1, findfirst(x -> abs(x - 180.0) < 1, bal.t), length(bal.t)]
key_indices = filter(!isnothing, key_indices)

for idx in key_indices
    t = bal.t[idx]
    println()
    println("t = $(round(Int, t)) days:")

    # Components
    P = pools.P[idx]
    CO2 = pools.CO2[idx]
    C_total = pools.C_total[idx]
    B_total = pools.B_total[idx]
    F_i = pools.F_i_total[idx]
    F_n = pools.F_n_total[idx]
    F_m = pools.F_m_total[idx]
    F_total = F_i + F_n + F_m
    E_total = pools.E_total[idx]
    M_total = pools.M_total[idx]

    # Total carbon in system
    carbon_in = P + C_total + B_total + F_total + E_total + M_total

    println("  POM (P):        $(round(P, digits=2)) μg-C")
    println("  Pools:")
    println("    DOC (C):      $(round(C_total, digits=2)) μg-C")
    println("    Bacteria (B): $(round(B_total, digits=2)) μg-C")
    println("    Fungi (F):    $(round(F_total, digits=2)) μg-C")
    println("    EPS (E):      $(round(E_total, digits=2)) μg-C")
    println("    MAOC (M):     $(round(M_total, digits=2)) μg-C")
    println("  Total in system: $(round(carbon_in, digits=2)) μg-C")
    println("  CO₂ respired:    $(round(CO2, digits=2)) μg-C")

    # Balance
    absolute_error = bal.C_total[idx] - bal.C_initial
    println("  Balance: $(round(absolute_error, digits=4)) μg-C ($(round(bal.relative_error[idx]*100, digits=4))%)")
end

println()
println("="^70)
