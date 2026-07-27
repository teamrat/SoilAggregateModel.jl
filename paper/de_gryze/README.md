# De Gryze et al. (2006) — Aggregate Formation Incubation

**Reference:** De Gryze, S., Six, J., Brits, C., & Merckx, R. (2006). A quantification of short‐term macroaggregate dynamics: influences of wheat residue input and texture. *European Journal of Soil Science*, 57, 235–246.

## Experiment

- 4-week (28-day) lab incubation with wheat residue amendment
- 5 soil textures (sandy loam → clay)
- Measurements: mean weight diameter (MWD) and cumulative CO₂ at multiple time points
- Data extends to day 21

## Model Setup

- 23 POM diameter classes (0.025 → 2.2 mm) representing a population of residue fragments
- Constant environment: ψ = −29 kPa, T = 20°C, 21% O₂
- 60-day simulation (extends past data window for post-peak dynamics)
- Each POM size class → one `run_aggregate()` call

## Files

| File | Description |
|------|-------------|
| `run_degryze.jl` | Main simulation script (forward run) |
| `degryze2006.csv` | MWD data (mm) for 5 soils |
| `degryze_CO2_2006.csv` | Cumulative CO₂ data for 5 soils |
| `output/` | Generated outputs (git-ignored) |

## MATLAB Correspondence

| MATLAB | Julia |
|--------|-------|
| `de_gryze_test.m` | `run_degryze.jl` |
| `param.m` | `BiologicalProperties()` + `SoilProperties()` |
| `single_aggregate_beta.m` | `run_aggregate()` |

## Status

- [x] Forward simulation script
- [ ] Parameter matching to MATLAB defaults
- [ ] Unit conversion for CO₂ comparison
- [ ] MWD calculation from aggregate diameters
- [ ] Calibration against 5 soils
- [ ] Figures for manuscript
