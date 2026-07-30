# setup.jl — shared configuration for everything in fitting/
#
# The MODEL is not configured here. Soil, initial condition, parameters, POM
# distribution, environment and domain tessellation all come from
# ../degryze_config.jl, which run_degryze.jl includes as well. This file adds
# only what fitting needs on top of it: the fitting window, one model
# evaluation, and the measured data (CLAUDE.md §8).
#
# Days 0-21. Day 21 is where the data ends; extrapolation happens after the
# parameters are fixed, not during fitting.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using SoilAggregateModel
using DataFrames, CSV, Printf, Distributions, Logging

include(joinpath(@__DIR__, "..", "postprocess_dataframe.jl"))
include(joinpath(@__DIR__, "..", "degryze_config.jl"))

# Aliases. The scan and sensitivity scripts perturb copies of a baseline, so the
# baseline carries a `0` suffix to keep "the reference" and "the perturbed"
# visually distinct at the call site.
const SOIL0 = SOIL
const BIO0  = BIO
const TF    = T_FUNC
const PSI   = PSI_FUNC
const O2F   = O2_FUNC

# ── Fitting window ───────────────────────────────────────────────────────────
# 0.25 d output spacing: every measured day is an integer, so `at_time` finds an
# exact match and never interpolates.
const T_END        = 21.0
const OUTPUT_TIMES = collect(0.0:0.25:T_END)

# ── Rebuilding a parameter struct with one field changed ─────────────────────

"""
    with_field(x, name, val)

Copy of `x` with field `name` set to `val`, rebuilt through the keyword
constructor so any field the constructor would otherwise derive is passed
explicitly.
"""
function with_field(x::T, name::Symbol, val) where {T}
    kw = NamedTuple(n => (n === name ? val : getfield(x, n)) for n in fieldnames(T))
    return T(; kw...)
end

# A no-change rebuild must reproduce the original exactly, or every perturbation
# built on this is measuring the rebuild instead of the parameter.
for (obj, nm) in ((BIO0, :r_B_max), (SOIL0, :κ_b))
    rebuilt = with_field(obj, nm, getfield(obj, nm))
    all(getfield(rebuilt, n) === getfield(obj, n) for n in fieldnames(typeof(obj))) ||
        error("with_field is not identity for $(typeof(obj))")
end

# ── One model evaluation ─────────────────────────────────────────────────────

"""
    run_model(bio, soil) -> DataFrame

Population-level time series for the five-diameter sweep. Columns of interest:
`time_days`, `MWD_fixed_weight_mm`, `CO2_total`, `pct_*`.
"""
function run_model(bio::BiologicalProperties, soil::SoilProperties)
    # snap_times omitted, so this returns the summary DataFrame alone.
    df_summary = run_diameter_sweep(DIAM_ALL, bio, soil, TF, PSI, O2F;
                                    t_max = T_END, output_times = OUTPUT_TIMES,
                                    n_grid = N_GRID, domain_factor = TESS.f_domain,
                                    ρ_POM = ρ_POM, ic = IC, ω = TESS.ω,
                                    solver = SOLVER_USED, quiet = true)
    return population_outputs(df_summary, POP.N_POM;
                              sieve_sizes = DEGRYZE_SIEVES,
                              mineral_class_fractions = MINERAL_CLASSES,
                              class_nominal_mm = DEGRYZE_CLASS_NOMINALS,
                              class_labels = DEGRYZE_CLASS_LABELS,
                              cell_volume_mm3 = POP.V_pack,
                              soil_volume_mm3 = SOIL_V,
                              ρ_b = ρ_BULK, f_C_POM = F_C_POM)
end

"""
    at_time(df, col, t) -> value

Value of `col` at time `t`. `OUTPUT_TIMES` has 0.25 d spacing and every measured
day is an integer, so an exact match exists; this asserts that rather than
silently returning a neighbour.
"""
function at_time(df::DataFrame, col::Symbol, t::Real)
    k = argmin(abs.(df.time_days .- t))
    abs(df.time_days[k] - t) < 1e-6 ||
        error("no output at t=$(t) (nearest $(df.time_days[k])); OUTPUT_TIMES must cover every measured day")
    return df[k, col]
end

# ── Measured data, soil 3 ────────────────────────────────────────────────────
#
# MWD: day 0 EXCLUDED, and now for a second reason. No fitted parameter reaches
# it — the model is flat over days 0-1 — and since 2026-07-29 the day-0 value is
# what the mineral sand split was INVERTED FROM (§0a A1'). Scoring it would score
# an identity. Absolute MWD is still what the later days are compared on.
#
# CO2: day 0 excluded because it is identically zero on both sides and
# ln(0/0) is undefined. Not a modelling choice.

const DATA_DIR = joinpath(@__DIR__, "..", "data")

const MWD_DAYS, MWD_OBS = let
    df = CSV.read(joinpath(DATA_DIR, "degryze2006.csv"), DataFrame;
                  header = 2, missingstring = "")
    col = Symbol("Soil$(SOIL_ID)")
    keep = [i for i in 1:nrow(df) if df[i, 1] > 0.0 && !ismissing(df[i, col])]
    (Float64.(df[keep, 1]), Float64.(df[keep, col]) .* 1000.0)   # µm
end

const CO2_DAYS, CO2_OBS = let
    df = CSV.read(joinpath(DATA_DIR, "degryze_CO2_2006.csv"), DataFrame)
    t  = eltype(df[:, 1]) <: AbstractString ? parse.(Float64, df[:, 1]) : Float64.(df[:, 1])
    col = Symbol("Soil_$(SOIL_ID)")
    keep = [i for i in 1:nrow(df) if t[i] > 0.0 && !ismissing(df[i, col])]
    (t[keep], Float64.(df[keep, col]))                            # µg-C/g-soil
end

"""
    obs_at(days, values, t) -> value

The measurement at day `t`. Errors rather than returning a neighbour, so a
script that asks for a day the data does not have fails loudly.
"""
function obs_at(days::AbstractVector, values::AbstractVector, t::Real)
    k = findfirst(d -> abs(d - t) < 1e-6, days)
    k === nothing && error("no measurement at day $(t); have $(days)")
    return values[k]
end
