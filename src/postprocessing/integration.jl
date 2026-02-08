# integration.jl
# Volumetric integration of carbon pools

"""
    IntegratedPools

Volumetrically integrated carbon pools at each output time.

# Fields
- `t::Vector{Float64}`: Time [days]
- `P::Vector{Float64}`: POM [μg-C] (scalar, same in both domains)
- `CO2::Vector{Float64}`: Cumulative CO₂ [μg-C]
- `r_agg::Vector{Float64}`: Aggregate radius [mm]

**Total domain** (entire r_0 to r_max):
- `C_total, B_total, F_i_total, F_n_total, F_m_total, E_total, M_total::Vector{Float64}`: [μg-C]

**Aggregate domain** (r_0 to r_agg only):
- `C_agg, B_agg, F_i_agg, F_n_agg, F_m_agg, E_agg, M_agg::Vector{Float64}`: [μg-C]
"""
struct IntegratedPools
    t::Vector{Float64}
    P::Vector{Float64}
    CO2::Vector{Float64}
    r_agg::Vector{Float64}
    # Total domain
    C_total::Vector{Float64}
    B_total::Vector{Float64}
    F_i_total::Vector{Float64}
    F_n_total::Vector{Float64}
    F_m_total::Vector{Float64}
    E_total::Vector{Float64}
    M_total::Vector{Float64}
    # Aggregate domain
    C_agg::Vector{Float64}
    B_agg::Vector{Float64}
    F_i_agg::Vector{Float64}
    F_n_agg::Vector{Float64}
    F_m_agg::Vector{Float64}
    E_agg::Vector{Float64}
    M_agg::Vector{Float64}
end

"""
    integrated_pools(result::SimulationResult)

Volumetrically integrated carbon pools at each output time.

Returns `IntegratedPools` struct with both total-domain and aggregate-domain integrations.

# Integration method
Uses conservation weights grid.W[i] = 4πr²h (precomputed in GridInfo).
- Total domain: sum over all nodes (r_0 to r_max)
- Aggregate domain: sum only nodes where r_grid[i] ≤ r_agg(t)

# Arguments
- `result::SimulationResult`: Complete simulation output

# Returns
- `IntegratedPools`: Struct with time series of integrated pools

# Notes
- r_agg computed at each output time using compute_r_agg
- P and CO2 are scalars, same in both domains
- All other pools integrated using conservation weights

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
pools = integrated_pools(result)
plot(pools.t, pools.C_total, label="Total DOC")
plot!(pools.t, pools.C_agg, label="Aggregate DOC")
```
"""
function integrated_pools(result::SimulationResult)
    n_out = length(result.outputs)
    grid = result.grid
    params = result.params

    # Allocate output vectors
    t_vec = Vector{Float64}(undef, n_out)
    P_vec = Vector{Float64}(undef, n_out)
    CO2_vec = Vector{Float64}(undef, n_out)
    r_agg_vec = Vector{Float64}(undef, n_out)

    # Allocate pool vectors (both domains, always)
    C_total = Vector{Float64}(undef, n_out)
    B_total = Vector{Float64}(undef, n_out)
    F_i_total = Vector{Float64}(undef, n_out)
    F_n_total = Vector{Float64}(undef, n_out)
    F_m_total = Vector{Float64}(undef, n_out)
    E_total = Vector{Float64}(undef, n_out)
    M_total = Vector{Float64}(undef, n_out)

    C_agg = Vector{Float64}(undef, n_out)
    B_agg = Vector{Float64}(undef, n_out)
    F_i_agg = Vector{Float64}(undef, n_out)
    F_n_agg = Vector{Float64}(undef, n_out)
    F_m_agg = Vector{Float64}(undef, n_out)
    E_agg = Vector{Float64}(undef, n_out)
    M_agg = Vector{Float64}(undef, n_out)

    for k in 1:n_out
        rec = result.outputs[k]
        t_vec[k] = rec.t
        P_vec[k] = rec.state.P
        CO2_vec[k] = rec.state.CO2_cumulative

        # Compute r_agg at this time
        r_agg_val = compute_r_agg(rec, grid, params)
        r_agg_vec[k] = r_agg_val

        # Initialize sums
        C_tot = 0.0
        B_tot = 0.0
        Fi_tot = 0.0
        Fn_tot = 0.0
        Fm_tot = 0.0
        E_tot = 0.0
        M_tot = 0.0

        C_ag = 0.0
        B_ag = 0.0
        Fi_ag = 0.0
        Fn_ag = 0.0
        Fm_ag = 0.0
        E_ag = 0.0
        M_ag = 0.0

        # Integrate over all nodes
        for i in 1:grid.n
            W = grid.W[i]

            # Total domain (all nodes)
            C_tot += rec.state.C[i] * W
            B_tot += rec.state.B[i] * W
            Fi_tot += rec.state.F_i[i] * W
            Fn_tot += rec.state.F_n[i] * W
            Fm_tot += rec.state.F_m[i] * W
            E_tot += rec.state.E[i] * W
            M_tot += rec.state.M[i] * W

            # Aggregate domain (only nodes within r_agg)
            if grid.r_grid[i] <= r_agg_val
                C_ag += rec.state.C[i] * W
                B_ag += rec.state.B[i] * W
                Fi_ag += rec.state.F_i[i] * W
                Fn_ag += rec.state.F_n[i] * W
                Fm_ag += rec.state.F_m[i] * W
                E_ag += rec.state.E[i] * W
                M_ag += rec.state.M[i] * W
            end
        end

        # Store results
        C_total[k] = C_tot
        B_total[k] = B_tot
        F_i_total[k] = Fi_tot
        F_n_total[k] = Fn_tot
        F_m_total[k] = Fm_tot
        E_total[k] = E_tot
        M_total[k] = M_tot

        C_agg[k] = C_ag
        B_agg[k] = B_ag
        F_i_agg[k] = Fi_ag
        F_n_agg[k] = Fn_ag
        F_m_agg[k] = Fm_ag
        E_agg[k] = E_ag
        M_agg[k] = M_ag
    end

    return IntegratedPools(
        t_vec, P_vec, CO2_vec, r_agg_vec,
        C_total, B_total, F_i_total, F_n_total, F_m_total, E_total, M_total,
        C_agg, B_agg, F_i_agg, F_n_agg, F_m_agg, E_agg, M_agg
    )
end

"""
    carbon_balance_table(result::SimulationResult)

Carbon conservation check at each output time.

Computes total carbon = P + Σ(pools×W) + CO₂ and checks against initial value.

# Returns
NamedTuple with:
- `t::Vector{Float64}`: Time [days]
- `C_total::Vector{Float64}`: Total carbon at each time [μg-C]
- `C_initial::Float64`: Reference value from first output [μg-C]
- `relative_error::Vector{Float64}`: (C_total[i] - C_initial) / C_initial

# Notes
- Builds on integrated_pools (uses :total domain)
- Should show |relative_error| < 1e-12 for correct implementation
- Any systematic drift indicates numerical issues

# Example
```julia
result = run_aggregate(bio, soil, T, ψ, O2, (0.0, 30.0))
bal = carbon_balance_table(result)
@assert all(abs.(bal.relative_error) .< 1e-12)
```
"""
function carbon_balance_table(result::SimulationResult)
    pools = integrated_pools(result)
    n = length(pools.t)
    C_tot = Vector{Float64}(undef, n)

    for i in 1:n
        C_tot[i] = pools.P[i] + pools.CO2[i] +
                   pools.C_total[i] + pools.B_total[i] +
                   pools.F_i_total[i] + pools.F_n_total[i] + pools.F_m_total[i] +
                   pools.E_total[i] + pools.M_total[i]
    end

    C_initial = C_tot[1]
    relative_error = [(C_tot[i] - C_initial) / C_initial for i in 1:n]

    return (t=pools.t, C_total=C_tot, C_initial=C_initial, relative_error=relative_error)
end
