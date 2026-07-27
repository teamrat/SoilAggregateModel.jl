# Instructions for Claude Code: Working with SoilAggregateModel.jl

## The Core Problem

You (Claude Code) have a tendency to **guess** struct field names, function signatures, and API conventions rather than inspecting the actual source code. This has already caused bugs — for example, using `grid.r` when the actual field is `grid.r_grid`, calling fungi pools "mantle" and "protected" when the model's terminology is "insulated" and "mobile", and passing O₂ as molar fraction (0.21) when the API expects concentration in μg/mm³ (≈0.27).

These are not minor style issues. In a scientific model, wrong names produce misleading outputs that could propagate into published figures. Wrong units produce wrong physics. **Never guess. Always look.**

---

## Rules

### 1. Inspect Before You Write

Before writing any script that uses `SoilAggregateModel`, read the relevant source files first. At minimum:

- **`src/types.jl`** — all struct definitions, field names, and field types
- **`src/parameters.jl`** — `BiologicalProperties` and `SoilProperties` constructors and default values
- **`src/api.jl`** — the `run_aggregate` function signature and keyword arguments
- **`src/result.jl`** — the `SimulationResult` and `OutputRecord` structs

Do not write a single line of code that references a struct field until you have read the struct definition. Do not call `run_aggregate` until you have read its signature.

### 2. Use the Model's Terminology Exactly

The manuscript and codebase use precise terminology for every pool and variable. Use these names — do not invent synonyms, "friendlier" names, or analogies.

| Symbol | Variable | Meaning | NOT this |
|--------|----------|---------|----------|
| C | `state.C` | Dissolved organic carbon (total: aqueous + sorbed) | "substrate", "soil carbon" |
| B | `state.B` | Bacterial biomass | "microbes" (too general) |
| F_i | `state.F_i` | Insulated fungal hyphae (protected within aggregate) | "protected fungi" without specifying insulated |
| F_n | `state.F_n` | Non-insulated fungal hyphae (exposed, active at periphery) | "active fungi" |
| F_m | `state.F_m` | Mobile fungal biomass (mantle zone, exploratory) | "mantle", "protected" — these are WRONG for F_m |
| E | `state.E` | Extracellular polymeric substances (EPS) | "exudates", "glue" |
| M | `state.M` | Mineral-associated organic carbon (MAOC) | "stable C", "humus" |
| O | `state.O` | Total oxygen (dissolved + gaseous) | "air" |
| P | `state.P` | Particulate organic matter (POM) — scalar, not spatially distributed | "litter", "detritus" |
| CO2_cumulative | `state.CO2_cumulative` | Cumulative respired carbon — scalar | "respiration" (that's a rate, not cumulative) |

### 3. Use Precomputed Quantities — Don't Reinvent the Math

The package provides conservation weights, post-processing functions, and integration routines that are **consistent with the solver's internal accounting**. Reimplementing these by hand risks subtle inconsistencies.

**Integration**: Use `grid.W` (precomputed conservation weights = 4π r² h, stencil-matched):
```julia
# CORRECT
total = sum(state.C[j] * grid.W[j] for j in 1:grid.n)

# WRONG — recomputes weights, may not match solver's stencil
total = sum(state.C[j] * 4π * grid.r_grid[j]^2 * grid.h for j in 1:grid.n)
```

**Post-processing**: If functions like `integrated_pools`, `carbon_balance_table`, `radial_profiles`, `aqueous_concentrations`, or `co2_flux` exist in `src/postprocessing/`, use them. Don't rewrite the logic in a script.

### 4. Check Units — The Model Has a Strict Unit System

Everything in the model uses: **μg/mm³, mm, days, kPa, K, J/mol**. No exceptions.

Common traps:
- **O₂ ambient**: The API expects concentration in μg/mm³ (≈0.27 at 20°C, 21% atmospheric), NOT molar fraction (0.21). Read how `run_aggregate` handles this — it may compute concentration internally from molar fraction, or it may expect the concentration directly.
- **Temperature**: Always in Kelvin (293.15, not 20).
- **Matric potential**: In kPa (−33.0, not −0.033 MPa or −330 cm H₂O).

### 5. Check Whether Arguments Are Scalars or Functions

Read the `run_aggregate` signature to determine whether environmental drivers (T, ψ, O₂) are passed as:
- Scalar values (constant throughout simulation), or
- Functions of time `t -> value`

Do not guess. The answer is in `src/api.jl`.

### 6. Don't Fabricate Post-Processing That Already Exists

Before writing any analysis code (CUE calculations, flux derivations, mass balance checks), search the `src/postprocessing/` directory for existing implementations. The Phase 1 infrastructure was specifically built to avoid ad hoc reimplementation in scripts.

If the function doesn't exist yet, say so explicitly: "This package doesn't appear to export a CUE function — here's a manual calculation, but it should be verified against the manuscript equations."

### 7. When Uncertain, State It

If you cannot find a field name, function signature, or expected unit by reading the source code, **say so**. Do not silently guess and produce code that compiles but does the wrong thing. A clear "I wasn't able to determine the field name for X — please check `types.jl`" is infinitely better than inventing `grid.r` when it's actually `grid.r_grid`.

---

## Pre-Flight Checklist

Before delivering any script that uses SoilAggregateModel, verify:

- [ ] Every struct field name matches `types.jl` exactly
- [ ] `run_aggregate` call matches the signature in `api.jl` exactly
- [ ] Units are correct (μg/mm³, mm, days, kPa, K)
- [ ] Integration uses `grid.W`, not hand-rolled weights
- [ ] Pool names in CSV headers and comments match the model's terminology
- [ ] Post-processing functions from the package are used where available
- [ ] Any quantity you computed manually is flagged as "manual — verify against package functions"

---

## Key Source Files to Read

| What you need | Read this file |
|---------------|---------------|
| Struct definitions, field names | `src/types.jl` |
| Parameter defaults and constructors | `src/parameters.jl` |
| `run_aggregate` signature | `src/api.jl` |
| Result/output structs | `src/result.jl` |
| Available post-processing | `src/postprocessing/*.jl` |
| Reaction equations | `src/solver/reactions.jl` |
| Conservation weights, grid | `src/solver/finite_volumes.jl` |
| Physical constants, unit conventions | `src/constants.jl` |

Read the relevant files **before** writing code. Not after you get an error.
