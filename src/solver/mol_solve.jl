# mol_solve.jl
# Stiff-solver entry point. Mirrors run_aggregate so every existing
# post-processing path works unchanged on the result.
#
# The numerical argument for this path is in mol.jl; docs/REFERENCE.md §20a.

"""
    run_aggregate_stiff(bio, soil, T_func, ψ_func, O2_func, t_span; kwargs...)

Integrate the same model as [`run_aggregate`](@ref) with an implicit stiff
solver over the method-of-lines system in `mol.jl`, and return the same
`SimulationResult`.

Keyword arguments match `run_aggregate` where they mean the same thing
(`n_grid`, `r_0`, `r_max`, `ic`, `initial_state`, `P_0`, `ω`, `t_delay`,
`output_times`). The step-control arguments differ, because an implicit solver
controls error rather than step size:

- `reltol::Real = 1e-6`, `abstol::Real = 1e-10`: solver tolerances. `abstol`
  matters here — the pools span roughly eight orders of magnitude, so a single
  scalar either over-resolves MAOC or under-resolves DOC. `abstol_scale`
  builds a per-state vector instead; see below.
- `abstol_scale::Union{Nothing,Real} = 1e-8`: if given, `abstol` becomes a
  vector, `max(abstol, abstol_scale · |u₀ᵢ|)` per state. Set to `nothing` to
  use the scalar.
- `alg = nothing`: solver. Defaults to `FBDF` with a sparse Jacobian and a KLU
  factorisation, which is the documented choice for systems of this size. Pass
  `KenCarp47(...)` if FBDF's Newton iterations fail — the BDF family is
  documented to tolerate less stiffness than the ESDIRK family.
- `sparse_jac::Bool = true`: hand the solver a sparse Jacobian pattern. A dense
  Jacobian on 1602 states is a ~1.4e9-flop factorisation per step; the pattern
  here is block-tridiagonal. Built structurally by
  [`mol_jacobian_prototype`](@ref), not detected by tracing — see its docstring.
- `jac_prototype = nothing`: override the pattern. Must be a superset of the
  true one; a pattern missing entries degrades the Newton solve silently.
- `maxiters::Int = 10^7`: the solver default of 1e5 is too low for a 45-day run.

# Returns
`SimulationResult`, identical in structure to `run_aggregate`, so
`integrated_pools`, `compute_r_agg`, `result_to_dataframe` and the rest apply
without change.

`CO2_cumulative` in each output state is recovered from the carbon balance —
initial total carbon minus current total carbon — not integrated during the
solve. It is a whole-domain quantity wanted at output times, so carrying it
through the solver as a per-node quadrature would build machinery for a
resolution nothing uses, and its dense Jacobian row would cost one RHS
evaluation per column on every Jacobian.

The consequence is that `mass_balance_error` is `NaN` for this path rather than
a number: the recovery makes the balance identically zero, so reporting it would
assert a check that was never performed. `diagnostics["co2_monotonic"]` carries
the check that IS available — respiration cannot be negative, so the recovered
quantity cannot fall between output times.

# On comparing the two solvers
They discretise space identically and share `compute_source_terms`, so a
difference between them is a difference in time integration, coefficient
lagging, or negativity clipping — see the header of `mol.jl`, which lists all
three. Any comparison should be made against the split solver run at a `dt_min`
small enough to be converged, not against its production settings.
"""
function run_aggregate_stiff(bio::BiologicalProperties, soil::SoilProperties,
                             T_func, ψ_func, O2_func, t_span;
                             n_grid::Int = 200, r_0::Real = 0.1, r_max::Real = 2.0,
                             initial_state = nothing,
                             ic::InitialConditions = InitialConditions(),
                             P_0::Union{Nothing,Real} = nothing,
                             ω::Real = 1.0,
                             t_delay::Real = 0.0,
                             reltol::Real = 1e-6,
                             abstol::Real = 1e-10,
                             abstol_scale::Union{Nothing,Real} = 1e-8,
                             alg = nothing,
                             sparse_jac::Bool = true,
                             jac_prototype = nothing,
                             maxiters::Int = 10^7,
                             output_times::Vector{<:Real} = Float64[],
                             verbose::Bool = false)

    t_start, t_end = t_span
    grid = GridInfo(n_grid, r_0, r_max)
    env  = EnvironmentalDrivers(T_func, ψ_func, O2_func)

    state = initial_state !== nothing ? deepcopy(initial_state) :
            create_initial_state(n_grid, bio, soil, ic; P_0 = P_0, ω = ω)

    if !(state.P_0 > 0.0)
        throw(ArgumentError(
            "AggregateState.P_0 must be > 0 (got $(state.P_0)). " *
            "create_initial_state sets it automatically."))
    end

    f_T = TemperatureCache()
    update_temperature_cache!(f_T, T_func(t_start), bio, soil)
    O_amb = O2_func(t_start) / f_T.K_H_O

    u0 = state_to_vector(state)
    mol_outer_oxygen!(u0, n_grid, O_amb)

    p = (n = n_grid, r_grid = grid.r_grid, h = grid.h,
         bio = bio, soil = soil,
         T_func = T_func, ψ_func = ψ_func, O2_func = O2_func,
         f_T = f_T, P_0 = state.P_0, t_delay = Float64(t_delay))

    # --- Jacobian sparsity ------------------------------------------------
    # Structural, not traced. See mol_jacobian_prototype for why: softplus
    # branches on its argument, so the tracer cannot resolve it, and the
    # local-tracer variant that could would return a state-dependent pattern
    # that the max(0, u) guards can silently truncate.
    fode = if sparse_jac
        proto = jac_prototype === nothing ? mol_jacobian_prototype(n_grid) : jac_prototype
        if verbose
            @info "MOL Jacobian" n_states=length(u0) nnz=nnz(proto) nonzero_fraction=nnz(proto)/length(proto)
        end
        ODEFunction{true}(mol_rhs!; jac_prototype = proto)
    else
        ODEFunction{true}(mol_rhs!)
    end

    algorithm = alg === nothing ?
        (sparse_jac ? FBDF(linsolve = KLUFactorization()) : FBDF()) : alg

    at = abstol_scale === nothing ? abstol :
         [max(abstol, abstol_scale * abs(x)) for x in u0]

    scheduled = isempty(output_times) ? Float64[] :
                sort(unique(Float64.(output_times)))
    filter!(x -> t_start <= x <= t_end, scheduled)
    isempty(scheduled) && (scheduled = Float64[t_start, t_end])

    prob = ODEProblem(fode, u0, (Float64(t_start), Float64(t_end)), p)
    sol  = solve(prob, algorithm;
                 reltol = reltol, abstol = at,
                 saveat = scheduled, maxiters = maxiters,
                 save_everystep = false)

    # Compared as a string so this does not depend on which package re-exports
    # ReturnCode. A non-Success code here is the signal that matters: it means
    # the reported trajectory is short or unconverged, not merely slow.
    if string(sol.retcode) != "Success"
        @warn "run_aggregate_stiff: solver did not report success" retcode=string(sol.retcode) t_final=sol.t[end]
    end

    # Reference total, with the respired term zeroed so it is pools + POM only.
    state.CO2_cumulative = 0.0
    C_total_initial = compute_total_carbon(state, grid.r_grid, grid.h)

    # Cumulative respired carbon is recovered here rather than integrated during
    # the solve: it is what has left the pools. `mass_balance_error` is set to
    # NaN because that recovery makes the balance identically zero — reporting
    # 0.0 would read as a passing check when nothing was checked. The genuine
    # check available at this resolution is monotonicity: respiration is
    # non-negative, so the recovered quantity can never fall.
    outputs = OutputRecord[]
    co2_prev = -Inf
    co2_monotonic = true
    for (idx, t) in enumerate(sol.t)
        st = deepcopy(state)
        vector_to_state!(st, sol.u[idx])
        st.CO2_cumulative = 0.0
        st.CO2_cumulative = C_total_initial - compute_total_carbon(st, grid.r_grid, grid.h)
        if st.CO2_cumulative < co2_prev
            co2_monotonic = false
            @warn "run_aggregate_stiff: recovered respired carbon decreased — carbon is being created somewhere" t previous=co2_prev now=st.CO2_cumulative
        end
        co2_prev = st.CO2_cumulative
        push!(outputs, OutputRecord(t, st, NaN))
    end

    diagnostics = Dict{String,Any}(
        "n_steps"     => sol.stats.naccept + sol.stats.nreject,
        "n_accept"    => sol.stats.naccept,
        "n_reject"    => sol.stats.nreject,
        "n_f"         => sol.stats.nf,
        "n_jac"       => sol.stats.njacs,
        "final_time"  => sol.t[end],
        "retcode"     => string(sol.retcode),
        "algorithm"   => string(nameof(typeof(algorithm))),
        "sparse_jac"  => sparse_jac,
        "reltol"      => reltol,
        "co2_monotonic" => co2_monotonic,
    )

    return SimulationResult(grid, ParameterSet(bio, soil), env, outputs, diagnostics)
end
