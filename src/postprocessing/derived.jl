# derived.jl
# Derived quantities from simulation results

"""
    aqueous_concentrations(record::OutputRecord, grid::GridInfo,
                          params::ParameterSet, env::EnvironmentalDrivers)

Aqueous-phase concentrations at each node.

Computes:
- C_aq[i] = C[i] / (1 + k_d·ρ_b)
- O_aq[i] = O[i]·θ[i] / (θ[i] + K_H(T)·θ_a[i])

# Arguments
- `record::OutputRecord`: Snapshot at a single time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Parameter bundle
- `env::EnvironmentalDrivers`: Time-dependent forcings

# Returns
NamedTuple with:
- `C_aq::Vector{Float64}`: Aqueous DOC [μg/mm³ water]
- `O_aq::Vector{Float64}`: Aqueous O₂ [μg/mm³ water]

# Notes
- C_aq accounts for sorption retardation (bulk-to-aqueous conversion)
- O_aq accounts for gas-water partitioning via Henry's law
- Uses K_H(T) from temperature-dependent Henry's law constant
- Recomputes θ, θ_a at record.t using environmental forcings

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
aq = aqueous_concentrations(result.outputs[15], result.grid, result.params, result.env)
plot(result.grid.r_grid, aq.C_aq, label="Aqueous DOC")
```
"""
function aqueous_concentrations(record::OutputRecord, grid::GridInfo,
                               params::ParameterSet, env::EnvironmentalDrivers)
    n = grid.n
    soil = params.soil

    # Recompute environmental conditions at this time
    env_vals = _prepare_environment(record, grid, params, env)
    θ = env_vals.θ
    θ_a = env_vals.θ_a
    K_H = env_vals.f_T.K_H_O

    # Allocate output vectors
    C_aq = Vector{Float64}(undef, n)
    O_aq = Vector{Float64}(undef, n)

    # Retardation factor (spatially uniform)
    retardation = 1.0 + soil.k_d_eq * soil.ρ_b

    # Compute at each node
    for i in 1:n
        C_aq[i] = record.state.C[i] / retardation
        O_aq[i] = record.state.O[i] * θ[i] / (θ[i] + K_H * θ_a[i])
    end

    return (C_aq=C_aq, O_aq=O_aq)
end

"""
    maoc_equilibrium(record::OutputRecord, grid::GridInfo,
                    params::ParameterSet, env::EnvironmentalDrivers)

Equilibrium MAOC concentration at each node from Langmuir-Freundlich isotherm.

Computes:
    M_eq[i] = M_max · (k_L·C_eq[i])^n_LF / [1 + (k_L·C_eq[i])^n_LF]

where C_eq = k_d·C_aq is the equilibrium sorbed concentration.

# Arguments
- `record::OutputRecord`: Snapshot at a single time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Parameter bundle
- `env::EnvironmentalDrivers`: Time-dependent forcings

# Returns
- `M_eq::Vector{Float64}`: Equilibrium MAOC at each node [μg-C/mm³]

# Notes
- Recomputes C_aq first (via aqueous_concentrations)
- Then computes C_eq = k_d·C_aq (equilibrium sorbed concentration)
- Uses Langmuir-Freundlich isotherm from maoc.jl
- Actual M may differ from M_eq (sorption/desorption kinetics)

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
M_eq = maoc_equilibrium(result.outputs[15], result.grid, result.params, result.env)
plot(result.grid.r_grid, [M_eq, result.outputs[15].state.M], label=["M_eq" "M_actual"])
```
"""
function maoc_equilibrium(record::OutputRecord, grid::GridInfo,
                         params::ParameterSet, env::EnvironmentalDrivers)
    n = grid.n
    soil = params.soil

    # Get aqueous concentrations first
    aq = aqueous_concentrations(record, grid, params, env)
    C_aq = aq.C_aq

    # Allocate output
    M_eq = Vector{Float64}(undef, n)

    # Compute equilibrium MAOC at each node
    for i in 1:n
        C_eq = soil.k_d_eq * C_aq[i]
        M_eq[i] = M_eq_langmuir_freundlich(C_eq, soil.M_max, soil.k_L, soil.n_LF)
    end

    return M_eq
end

"""
    respiration_rates(record::OutputRecord, grid::GridInfo,
                     params::ParameterSet, env::EnvironmentalDrivers)

Instantaneous respiration rates at each node.

Returns respiration components:
- Resp_B: Bacterial respiration = (1-Y_B)·R_B + m_B·B
- Resp_F: Fungal uptake respiration = (1-Y_F)·R_F
- Resp_F_conv: Fungal conversion respiration (from transitions)
- Resp_total: Total = Resp_B + Resp_F + Resp_F_conv

# Arguments
- `record::OutputRecord`: Snapshot at a single time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Parameter bundle
- `env::EnvironmentalDrivers`: Time-dependent forcings

# Returns
NamedTuple with:
- `Resp_B::Vector{Float64}`: Bacterial respiration [μg-C/mm³/day]
- `Resp_F::Vector{Float64}`: Fungal uptake respiration [μg-C/mm³/day]
- `Resp_F_conv::Vector{Float64}`: Fungal conversion respiration [μg-C/mm³/day]
- `Resp_total::Vector{Float64}`: Total respiration [μg-C/mm³/day]

# Notes
- Calls compute_source_terms at each node (reuses existing tested functions)
- Recomputes environmental conditions at record.t
- All rates are instantaneous (not time-integrated)

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
resp = respiration_rates(result.outputs[15], result.grid, result.params, result.env)
plot(result.grid.r_grid, resp.Resp_total, label="Total respiration")
```
"""
function respiration_rates(record::OutputRecord, grid::GridInfo,
                          params::ParameterSet, env::EnvironmentalDrivers)
    n = grid.n
    bio = params.bio
    soil = params.soil

    # Recompute environmental conditions at this time
    env_vals = _prepare_environment(record, grid, params, env)
    θ = env_vals.θ
    θ_a = env_vals.θ_a
    ψ = env_vals.ψ
    f_T = env_vals.f_T

    # Allocate output vectors
    Resp_B_vec = Vector{Float64}(undef, n)
    Resp_F_vec = Vector{Float64}(undef, n)
    Resp_F_conv_vec = Vector{Float64}(undef, n)
    Resp_total_vec = Vector{Float64}(undef, n)

    # Compute at each node using existing source term functions
    for i in 1:n
        # Extract state at this node
        C = record.state.C[i]
        B = record.state.B[i]
        F_n = record.state.F_n[i]
        F_m = record.state.F_m[i]
        F_i = record.state.F_i[i]
        E = record.state.E[i]
        M = record.state.M[i]
        O = record.state.O[i]

        # Compute all source terms (includes all respiration components)
        src = compute_source_terms(C, B, F_n, F_m, F_i, E, M, O,
                                   θ[i], θ_a[i], ψ, bio, soil, f_T)

        # Extract respiration components from compute_source_terms
        # We need to recompute individual components since SourceTerms only stores Resp_total
        # Let me check what compute_source_terms actually computes...

        # From reactions.jl, I can see:
        # Resp_B_val = Resp_B(R_Bb_val, R_diff, Y_B_val)
        # Resp_F_val = Resp_F(R_F_val, Y_F_val)
        # Resp_total_val = Resp_B_val + Resp_F_val + trans.Resp_F_conv

        # But SourceTerms struct only has Resp_total, not the components
        # I need to recompute the components manually

        # Compute C_aq, O_aq
        C_aq = C / (θ[i] + soil.ρ_b * soil.k_d_eq)
        O_aq = O * θ[i] / (θ[i] + f_T.K_H_O * θ_a[i])

        # Bacterial respiration
        r_B_max_T = bio.r_B_max * f_T.f_bac
        R_B_val = R_B(C_aq, O_aq, B, ψ, r_B_max_T, bio.K_B, bio.L_B, bio.ν_B)
        R_Bb_val = R_Bb(bio.C_B, O_aq, B, r_B_max_T, bio.K_B, bio.L_B, bio.B_min)
        R_diff = R_B_val - R_Bb_val
        Y_B_val = Y_B_func(R_diff, bio.Y_B_max, bio.K_Y)
        Resp_B_val = Resp_B(R_Bb_val, R_diff, Y_B_val)

        # Fungal respiration
        r_F_max_T = bio.r_F_max * f_T.f_fun
        R_F_val = R_F(C_aq, O_aq, F_i, F_n, bio.λ, ψ, r_F_max_T, bio.K_F, bio.L_F, bio.ν_F)
        Y_F_val = Y_F_const(bio.Y_F)
        Resp_F_val = Resp_F(R_F_val, Y_F_val)

        # Fungal conversion respiration (from transitions)
        μ_F_T = bio.μ_F * f_T.f_fun
        α_i_T = bio.α_i * f_T.f_fun
        α_n_T = bio.α_n * f_T.f_fun
        β_i_T = bio.β_i * f_T.f_fun
        β_n_T = bio.β_n * f_T.f_fun
        ζ_T = bio.ζ * f_T.f_fun
        Π_val = Pi_protected(F_m, F_i, F_n, bio.ε_F)
        trans = fungal_transitions(F_i, F_n, F_m, Π_val, α_i_T, α_n_T, β_i_T, β_n_T,
                                   ζ_T, bio.delta, bio.η_conv, bio.ε_F)

        # Store results
        Resp_B_vec[i] = Resp_B_val
        Resp_F_vec[i] = Resp_F_val
        Resp_F_conv_vec[i] = trans.Resp_F_conv
        Resp_total_vec[i] = Resp_B_val + Resp_F_val + trans.Resp_F_conv
    end

    return (Resp_B=Resp_B_vec, Resp_F=Resp_F_vec,
            Resp_F_conv=Resp_F_conv_vec, Resp_total=Resp_total_vec)
end

"""
    carbon_use_efficiency(record::OutputRecord, grid::GridInfo,
                         params::ParameterSet, env::EnvironmentalDrivers)

Carbon use efficiency (CUE) at each node.

Returns:
- CUE_B = Γ_B / R_B (bacterial CUE)
- CUE_F = Γ_F / R_F (fungal CUE)

where Γ is anabolic carbon allocation (growth) and R is total uptake.

# Arguments
- `record::OutputRecord`: Snapshot at a single time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Parameter bundle
- `env::EnvironmentalDrivers`: Time-dependent forcings

# Returns
NamedTuple with:
- `CUE_B::Vector{Float64}`: Bacterial carbon use efficiency [-]
- `CUE_F::Vector{Float64}`: Fungal carbon use efficiency [-]

# Notes
- CUE = (growth) / (uptake) = fraction of uptake allocated to anabolism
- For bacteria: CUE_B varies with substrate availability (flexible yield Y_B)
- For fungi: CUE_F is approximately constant (fixed yield Y_F)
- When R ≈ 0: CUE → 0 (avoid division by zero)
- Recomputes environmental conditions at record.t

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
cue = carbon_use_efficiency(result.outputs[15], result.grid, result.params, result.env)
plot(result.grid.r_grid, [cue.CUE_B, cue.CUE_F], label=["Bacteria" "Fungi"])
```
"""
function carbon_use_efficiency(record::OutputRecord, grid::GridInfo,
                              params::ParameterSet, env::EnvironmentalDrivers)
    n = grid.n
    bio = params.bio
    soil = params.soil

    # Recompute environmental conditions at this time
    env_vals = _prepare_environment(record, grid, params, env)
    θ = env_vals.θ
    θ_a = env_vals.θ_a
    ψ = env_vals.ψ
    f_T = env_vals.f_T

    # Allocate output vectors
    CUE_B_vec = Vector{Float64}(undef, n)
    CUE_F_vec = Vector{Float64}(undef, n)

    # Compute at each node
    for i in 1:n
        # Extract state at this node
        C = record.state.C[i]
        B = record.state.B[i]
        F_n = record.state.F_n[i]
        F_m = record.state.F_m[i]
        F_i = record.state.F_i[i]
        O = record.state.O[i]

        # Compute C_aq, O_aq
        C_aq = C / (θ[i] + soil.ρ_b * soil.k_d_eq)
        O_aq = O * θ[i] / (θ[i] + f_T.K_H_O * θ_a[i])

        # Bacterial CUE
        r_B_max_T = bio.r_B_max * f_T.f_bac
        R_B_val = R_B(C_aq, O_aq, B, ψ, r_B_max_T, bio.K_B, bio.L_B, bio.ν_B)
        R_Bb_val = R_Bb(bio.C_B, O_aq, B, r_B_max_T, bio.K_B, bio.L_B, bio.B_min)
        R_diff = R_B_val - R_Bb_val
        Y_B_val = Y_B_func(R_diff, bio.Y_B_max, bio.K_Y)
        Γ_B_val = Gamma_B(R_B_val, R_Bb_val, Y_B_val, bio.γ)

        # CUE_B = Γ_B / R_B (avoid division by zero)
        if R_B_val > 1e-15
            CUE_B_vec[i] = Γ_B_val / R_B_val
        else
            CUE_B_vec[i] = 0.0
        end

        # Fungal CUE
        r_F_max_T = bio.r_F_max * f_T.f_fun
        R_F_val = R_F(C_aq, O_aq, F_i, F_n, bio.λ, ψ, r_F_max_T, bio.K_F, bio.L_F, bio.ν_F)
        Y_F_val = Y_F_const(bio.Y_F)
        Γ_F_val = Gamma_F(Y_F_val, R_F_val)

        # CUE_F = Γ_F / R_F (avoid division by zero)
        if R_F_val > 1e-15
            CUE_F_vec[i] = Γ_F_val / R_F_val
        else
            CUE_F_vec[i] = 0.0
        end
    end

    return (CUE_B=CUE_B_vec, CUE_F=CUE_F_vec)
end

"""
    co2_flux(result::SimulationResult)

Instantaneous CO₂ flux [μg-C/day] at each output time.

Computed using backward finite difference:
    flux[i] = (CO₂[i] - CO₂[i-1]) / (t[i] - t[i-1])

# Arguments
- `result::SimulationResult`: Complete simulation output

# Returns
- `flux::Vector{Float64}`: CO₂ flux at each output time [μg-C/day]
  - flux[1] = flux[2] (extrapolation for first point)
  - flux[i] for i ≥ 2 from backward difference

# Notes
- Pure finite differences on stored CO2_cumulative
- No environmental recomputation needed
- Non-negative except for floating point errors

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
flux = co2_flux(result)
plot([result.outputs[i].t for i in 1:length(result.outputs)], flux)
```
"""
function co2_flux(result::SimulationResult)
    n = length(result.outputs)
    flux = Vector{Float64}(undef, n)

    for i in 2:n
        dt = result.outputs[i].t - result.outputs[i-1].t
        dCO2 = result.outputs[i].state.CO2_cumulative - result.outputs[i-1].state.CO2_cumulative
        flux[i] = dt > 0 ? dCO2 / dt : 0.0
    end

    # Extrapolate first point
    flux[1] = n > 1 ? flux[2] : 0.0

    return flux
end
