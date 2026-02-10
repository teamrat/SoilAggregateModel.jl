"""
Conservation Ablation Testing

Systematically disable subsystems to isolate the source of the 1% conservation error.
Strategy: If disabling a subsystem restores conservation, that subsystem contains the leak.
"""

using SoilAggregateModel

function run_conservation_test(label::String; bio_kwargs=Dict{Symbol,Any}(), soil_kwargs=Dict{Symbol,Any}(), tspan=(0.0, 30.0))
    # Merge overrides with defaults
    bio = BiologicalProperties(; bio_kwargs...)
    soil = SoilProperties(; soil_kwargs...)

    # Standard conditions (same as Theme 1)
    # Wrap constants in functions for time-varying API
    T_func = t -> 293.15  # K (20°C)
    ψ_func = t -> -33.0   # kPa
    O2_func = t -> 0.21   # 21% O₂

    # Run simulation
    result = run_aggregate(bio, soil, T_func, ψ_func, O2_func, tspan)

    # Compute conservation error
    bal = carbon_balance_table(result)
    max_error = maximum(abs.(bal.relative_error))
    conservation_pct = max_error * 100

    # Report
    println("$label:")
    println("  Conservation error: $(round(conservation_pct, digits=4))%")
    println("  Steps: $(result.diagnostics["n_steps"])")
    println()

    return (result, max_error)
end

println("="^70)
println("Conservation Ablation Tests — 30-day simulations")
println("="^70)
println()

# Test 1: Full model (baseline — should show ~1% error)
run_conservation_test("1. Full model (baseline)")

# Test 2: No MAOC (disable sorption/desorption)
run_conservation_test("2. No MAOC",
    bio_kwargs=Dict(:κ_s_ref => 0.0, :κ_d_ref => 0.0))

# Test 3: No fungi (disable all fungal processes)
run_conservation_test("3. No fungi",
    bio_kwargs=Dict(:r_F_max => 0.0, :μ_F => 0.0, :ζ => 0.0, :α_i => 0.0, :α_n => 0.0, :β_i => 0.0, :β_n => 0.0))

# Test 4: No EPS (disable bacterial EPS production)
run_conservation_test("4. No EPS",
    bio_kwargs=Dict(:γ => 0.0, :μ_E_max => 0.0))

# Test 5: No POM dissolution (closed system — no external carbon source)
run_conservation_test("5. No POM",
    bio_kwargs=Dict(:R_P_max => 0.0))

# Test 6: No space limitation (disable space-limited yield)
run_conservation_test("6. No space limit",
    bio_kwargs=Dict(:B_S => 1e10, :F_S => 1e10))

# Test 7: Abiotic only (disable all biology)
run_conservation_test("7. Abiotic only",
    bio_kwargs=Dict(
        :r_B_max => 0.0,
        :r_F_max => 0.0,
        :μ_B => 0.0,
        :μ_F => 0.0,
        :μ_E_max => 0.0,
        :R_P_max => 0.0,
        :κ_s_ref => 0.0,
        :κ_d_ref => 0.0,
        :γ => 0.0
    ))

println("="^70)
println("Interpretation:")
println("  If a test shows ~0% error, the disabled subsystem contains the leak.")
println("  If all tests still show ~1% error, the leak is in diffusion/integration.")
println("="^70)
