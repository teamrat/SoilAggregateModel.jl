## Setup
using Pkg
Pkg.activate(@__DIR__() * "/..")

using SoilAggregateModel

bio = BiologicalProperties()
soil = SoilProperties()

T_func = t -> 293.15
ψ_func = t -> -33.0
O2_func = t -> 0.21

# --- Test 1: Baseline (current dt_max=0.1) ---
println("Test 1: Baseline dt_max=0.1")
t1 = @elapsed result1 = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 365.0);
                                       output_times=[0.0, 365.0])
println("  Steps: $(result1.diagnostics["n_steps"])")
println("  Wall time: $(round(t1, digits=1))s")
bal1 = carbon_balance_table(result1)
println("  Conservation: $(bal1.relative_error[end])")
println()

# --- Test 2: Relaxed dt_max=0.5 ---
println("Test 2: dt_max=0.5")
t2 = @elapsed result2 = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 365.0);
                                       output_times=[0.0, 365.0], dt_max=0.5)
println("  Steps: $(result2.diagnostics["n_steps"])")
println("  Wall time: $(round(t2, digits=1))s")
bal2 = carbon_balance_table(result2)
println("  Conservation: $(bal2.relative_error[end])")
println()

# --- Test 3: Relaxed dt_max=1.0 ---
println("Test 3: dt_max=1.0")
t3 = @elapsed result3 = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 365.0);
                                       output_times=[0.0, 365.0], dt_max=1.0)
println("  Steps: $(result3.diagnostics["n_steps"])")
println("  Wall time: $(round(t3, digits=1))s")
bal3 = carbon_balance_table(result3)
println("  Conservation: $(bal3.relative_error[end])")
println()

# --- Test 4: 5-year run at dt_max=0.1 (feasibility check) ---
println("Test 4: 5-year baseline")
t4 = @elapsed result4 = run_aggregate(bio, soil, T_func, ψ_func, O2_func, (0.0, 1825.0);
                                       output_times=[0.0, 1825.0])
println("  Steps: $(result4.diagnostics["n_steps"])")
println("  Wall time: $(round(t4, digits=1))s")
bal4 = carbon_balance_table(result4)
println("  Conservation: $(bal4.relative_error[end])")
println()

# --- Summary ---
println("="^60)
println("Summary:")
println("  Baseline 1yr:  $(result1.diagnostics["n_steps"]) steps, $(round(t1,digits=1))s")
println("  dt_max=0.5:    $(result2.diagnostics["n_steps"]) steps, $(round(t2,digits=1))s")
println("  dt_max=1.0:    $(result3.diagnostics["n_steps"]) steps, $(round(t3,digits=1))s")
println("  Baseline 5yr:  $(result4.diagnostics["n_steps"]) steps, $(round(t4,digits=1))s")
println()
println("Target: 400 aggregates × 10yr < 1hr on M3 Ultra")
println("Per-aggregate 10yr budget: ~9s (serial) or ~90s (32-core parallel)")
println("Current 5yr rate: $(round(t4,digits=1))s → 10yr ≈ $(round(2*t4,digits=1))s")
