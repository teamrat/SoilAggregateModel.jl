# Fixed dt test for 365 days
# This will take 365,000 steps — feasible
using SoilAggregateModel
import SoilAggregateModel: diffusion_step!, reaction_step!, GridInfo
import SoilAggregateModel: create_initial_state, Workspace
import SoilAggregateModel: J_P, R_P
import SoilAggregateModel: update_temperature_cache!, update_water_content!, update_effective_diffusion!

bio = BiologicalProperties()
soil = SoilProperties()
n = 100
grid = GridInfo(n, 0.25, 5.0)
state = create_initial_state(n, bio)
workspace = Workspace(n)

T = 293.15; ψ = -33.0; O2_amb = 0.21
update_temperature_cache!(workspace.f_T, T, bio, soil)
update_water_content!(workspace.θ, workspace.θ_a, ψ, state, soil)
update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

function total_carbon(state, grid)
    total = state.P + state.CO2_cumulative
    for i in 1:grid.n
        total += (state.C[i] + state.B[i] + state.F_i[i] + state.F_n[i] +
                  state.F_m[i] + state.E[i] + state.M[i]) * grid.W[i]
    end
    return total
end

C_initial = total_carbon(state, grid)
dt = 0.001
n_steps = 365_000  # 365 days

println("Running 365 days with fixed dt=$dt ($n_steps steps)...")
t = 0.0
for step in 1:n_steps
    global t

    B_0 = state.B[1]; F_n_0 = state.F_n[1]
    θ_0 = workspace.θ[1]; θ_a_0 = workspace.θ_a[1]
    O_aq_0 = state.O[1] * θ_0 / (θ_0 + workspace.f_T.K_H_O * θ_a_0)
    R_P_max_T = bio.R_P_max * workspace.f_T.f_pom
    J_P_val = J_P(state.P, bio.P_0, B_0, F_n_0, θ_0, O_aq_0, R_P_max_T,
                  bio.K_B_P, bio.K_F_P, bio.θ_P, bio.L_P)
    R_P_val = R_P(J_P_val, grid.r_grid[1])

    diffusion_step!(state, workspace, dt/2, grid.r_grid, grid.h, J_P_val, O2_amb)
    reaction_step!(state, workspace, dt, grid.r_grid, grid.h, bio, soil, ψ, R_P_val)
    diffusion_step!(state, workspace, dt/2, grid.r_grid, grid.h, J_P_val, O2_amb)

    update_water_content!(workspace.θ, workspace.θ_a, ψ, state, soil)
    update_effective_diffusion!(workspace, soil, bio, workspace.f_T)

    t += dt

    if step % 30_000 == 0
        C_now = total_carbon(state, grid)
        err = (C_now - C_initial) / C_initial * 100
        println("  t=$(round(t, digits=0))d: error=$(err)%")
    end
end

C_final = total_carbon(state, grid)
rel_err = (C_final - C_initial) / C_initial * 100
println("\nFinal: $(rel_err)% (adaptive gave 0.014354%)")
