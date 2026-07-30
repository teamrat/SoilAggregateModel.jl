# De Gryze tooling, archived 2026-07-30

Everything here lived under `paper/de_gryze/` and did not contribute to the
forward run (`run_degryze.jl`). Archived rather than patched, on the owner's
instruction: keep outside `src/` only what the main De Gryze run needs.

Nothing here should be quoted. All of it predates or bypasses the current
configuration, and the numbers in the headers were measured at parameters that
have since changed (`k_L` 1000 → 25000, `κ_d_ref` 0.001 → 1e-4, `R_P_max`
2.5 → 1.5, `P_ATM` 101000 → 101325).

## optimize_soil3.jl  (663 lines)

LHS + Nelder–Mead fit of soil 3 to cumulative CO₂. Superseded: it forked the
whole configuration, and it hardcoded the soil-3 CO₂ series as literals when
those exact numbers are column `Soil_3` of `data/degryze_CO2_2006.csv` — the
measured data duplicated into a script.

The plan that replaces it: get the fit visually close by hand first, then build
a reusable fitting routine central enough for every example to use, tested
against De Gryze. That is a separate future project, not a repair of this file.

## fitting/  (5 files, 667 lines)

`setup.jl` (config + one model evaluation + observed data), `loss.jl` (the
objective), `scan.jl`, `scan_transport.jl`, `sensitivity.jl` (±10 % elasticities
of day-21 MWD and CO₂). Architecturally sound — all five took their model
configuration from `degryze_config.jl` — but calibration is the future project,
so they are not part of the current run. `setup.jl` is the piece worth reading
first when that project starts: `run_model`, `with_field` (with an identity
self-check), `at_time`, and the data readers are all reusable.

`sensitivity.jl`'s 2026-07-29 result is worth remembering: the six activation
energies, ν_B and ν_F are collinear with the rates they multiply when T and ψ
are constant — redundant, not inert. Hold them fixed in any fit.

## diagnostics/  (4 files, 740 lines)

`smoke_stiff.jl` — does the stiff solver load, compile, run. Permanently
answered. `bench_stiff.jl` — does the implicit step grow; the evidence for
switching solvers, and the switch is made. `converge_doc.jl` — which solver is
right about DOC at day 45; its premise ("9.6 % disagreement on C") was measured
at `k_L = 1000` against today's 25000, the one parameter that moves DOC, so the
claim it supports is unsupported until re-derived. `diagnose_speed.jl` — the
07-30 profiler that closed the audit's cost arm; its numbers are quoted in
`docs/_archive/AUDIT_20260729.md` §B, so the result survives the archive. It is the one file here
that is current and can be lifted straight back out if a profile is wanted.

## legacy/  (2 R scripts)

Superseded plotting. The R analysis layer gets built once the Julia side is
settled.
