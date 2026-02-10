using SoilAggregateModel
import SoilAggregateModel: diffusion_step!, reaction_step!, GridInfo
import SoilAggregateModel: AggregateState, Workspace, TemperatureCache
import SoilAggregateModel: J_P, R_P
import SoilAggregateModel: update_temperature_cache!, update_water_content!, update_effective_diffusion!
import SoilAggregateModel: create_initial_state

bio = BiologicalProperties()
soil = SoilProperties()
n = 100
grid = GridInfo(n, 0.25, 5.0)

# Use the standard initial state (has seed biomass)
state = create_initial_state(n, bio)

# Full workspace initialization
workspace = Workspace(n)
T = 293.15
ψ = -33.0
update_temperature_cache!(workspace.f_T, T, bio, soil)
update_water_content!(workspace.θ, workspace.θ_a, ψ, state, soil)
update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

# Helper
function total_carbon(state, grid)
    total = state.P + state.CO2_cumulative
    for i in 1:grid.n
        total += (state.C[i] + state.B[i] + state.F_i[i] + state.F_n[i] +
                  state.F_m[i] + state.E[i] + state.M[i]) * grid.W[i]
    end
    return total
end

dt = 0.001  # Small fixed timestep
O2_amb = 0.21

max_diff1_err = 0.0
max_rxn_err = 0.0
max_diff2_err = 0.0
max_total_err = 0.0

C_initial = total_carbon(state, grid)
println("Initial total carbon: $C_initial")

n_steps = 10000  # 10 days at dt=0.001

for step in 1:n_steps
    global max_diff1_err, max_rxn_err, max_diff2_err, max_total_err

    # Compute POM flux
    B_0 = state.B[1]
    F_n_0 = state.F_n[1]
    θ_0 = workspace.θ[1]
    θ_a_0 = workspace.θ_a[1]
    O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
    R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
    J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                  bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
    R_P_val = R_P(J_P_val, grid.r_grid[1])

    expected_flux_half = J_P_val * 4.0 * π * grid.r_grid[1]^2 * (dt/2)

    C0 = total_carbon(state, grid)

    # Diffusion half-step 1
    diffusion_step!(state, workspace, dt/2, grid.r_grid, grid.h, J_P_val, O2_amb)
    C1 = total_carbon(state, grid)
    d1_err = (C1 - C0) - expected_flux_half

    # Reaction step (should conserve exactly: P decrease + pool changes + CO₂ = 0)
    reaction_step!(state, workspace, dt, grid.r_grid, grid.h, bio, soil, ψ, R_P_val)
    C2 = total_carbon(state, grid)
    r_err = C2 - C1

    # Diffusion half-step 2
    diffusion_step!(state, workspace, dt/2, grid.r_grid, grid.h, J_P_val, O2_amb)
    C3 = total_carbon(state, grid)
    d2_err = (C3 - C2) - expected_flux_half

    # Update water content (once per step, as in real timestepper)
    update_water_content!(workspace.θ, workspace.θ_a, ψ, state, soil)
    update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

    max_diff1_err = max(max_diff1_err, abs(d1_err))
    max_rxn_err = max(max_rxn_err, abs(r_err))
    max_diff2_err = max(max_diff2_err, abs(d2_err))

    total_err = C3 - C_initial
    max_total_err = max(max_total_err, abs(total_err))

    if step in [1, 5, 10, 100, 1000, 2000, 4000, 6000, 8000, 10000]
        expected_pom_decrease = -R_P_val * dt
        diff_from_expected = r_err - expected_pom_decrease
        println("Step $step (t=$(round(step*dt, digits=3))d):")
        println("  rxn_err = $(r_err)")
        println("  expected = $(expected_pom_decrease)")
        println("  diff = $(diff_from_expected)")
        println()
    end
end

println()
println("="^60)
println("Maximum absolute errors across $n_steps steps:")
println("  Diff1:  $(max_diff1_err)")
println("  Rxn:    $(max_rxn_err)")
println("  Diff2:  $(max_diff2_err)")
println("  Total:  $(max_total_err)")
println()
println("Final total carbon: $(total_carbon(state, grid))")
println("Initial total carbon: $C_initial")
println("Net error: $(total_carbon(state, grid) - C_initial)")
println("Relative:  $((total_carbon(state, grid) - C_initial) / C_initial)")
