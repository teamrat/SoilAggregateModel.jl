"""
Test conservation with the full adaptive timestepper (run_aggregate)
to see if the 0.014%/year leak appears in the adaptive stepping logic.
"""

using SoilAggregateModel

# Standard setup
bio = BiologicalProperties()
soil = SoilProperties()

# Environmental conditions
T_func = t -> 293.15  # K
ψ_func = t -> -33.0   # kPa
O2_func = t -> 0.21   # 21% O₂

println("Running 365-day simulation with adaptive timestepper...")
println()

# Run with adaptive timestepper
result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 365.0);
                       output_times=collect(0.0:30.0:365.0))

println("Simulation complete: $(result.diagnostics["n_steps"]) timesteps")
println()

# Compute conservation
bal = carbon_balance_table(result)

println("="^70)
println("Conservation Results (Adaptive Timestepper)")
println("="^70)
println()
println("Time (days) | Total C (μg-C) | Error (μg-C) | Rel Error (%)")
println("-"^70)

for i in 1:length(bal.t)
    t = bal.t[i]
    C_total = bal.C_total[i]
    error = C_total - bal.C_initial
    rel_pct = bal.relative_error[i] * 100

    println("$(rpad(round(t, digits=1), 11)) | ",
            "$(rpad(round(C_total, digits=4), 14)) | ",
            "$(rpad(round(error, digits=6), 12)) | ",
            "$(round(rel_pct, digits=6))%")
end

println()
println("="^70)
println("Final error: $(round(bal.relative_error[end] * 100, digits=6))%")
println("Expected for 30 days at 0.014%/year: $(round(0.014 * 30/365, digits=6))%")
println("="^70)
