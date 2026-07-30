# Measured data — De Gryze et al. (2006)

De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R. (2006). *Water
repellence and soil aggregate dynamics in a loamy grassland soil as affected by
texture.* **European Journal of Soil Science** 57, 235–246.
PDF: `references/De Gryze et al EJSS 2006.pdf`.

**These files are what the model is compared against. They are never a model
input.** Measured inputs (texture, bulk density, SOC) live in
`../degryze_soils.jl`.

---

## The experiment (p. 236–237)

Five soils from one fertilised sheep meadow at K.U. Leuven, Lovenjoel, Belgium
(50°51′N, 4°47′E; Gleyic and Haplic Luvisols), sampled autumn 2003 along a
natural texture gradient over a short distance — so texture varies while parent
material, climate and management do not. That is what makes the five soils a
texture comparison rather than five unrelated sites.

Soils were taken from 0–10 cm, air-dried, crushed through a 250 µm sieve.
Material > 2 mm was discarded. Material > 250 µm (native POM and sand) was
dispersed in sodium hexametaphosphate, its organic matter removed with 10 %
H₂O₂ and combustion at 500 °C, then recombined with the < 250 µm material in the
original proportion and homogenised. **So the incubation started from soil with
no aggregates > 250 µm and no native POM** — every macroaggregate measured
during the incubation was formed during it.

Incubation: 150 g soil + 1.5 g wheat stems (*Triticum aestivum* L., 44.3 % C,
0.17 % N, C:N 261, cut to 0.5–2 mm), 60 % water-filled pore space, compressed to
each soil's field bulk density in 5 cm stainless-steel cylinders, 21 days at
25 °C in glass canning jars. Jars flushed with compressed air daily, so O₂ stays
ambient. 14 samples per soil; destructive sampling for aggregates.

---

## `degryze2006.csv` — aggregate mean weight diameter

**Units: mm.** Columns `Days, Soil1 … Soil5`. Line endings are **CR only**
(classic Mac) — `wc -l` reports 0. `CSV.read(...; header=2, missingstring="")`
handles it; row 1 is a citation line, not a header.

MWD is computed by the paper's equation (1), p. 237:

    MWD = (5000[A%] + 1125[B%] + 151.5[C%] + 26.5[D%]) / 100     [µm]

over the four wet-sieved fractions — `[A%]` > 2000 µm (large macroaggregates),
`[B%]` 250–2000 (small macroaggregates), `[C%]` 53–250 (microaggregates),
`[D%]` < 53 (silt and clay). The weights are **fixed constants, not class
midpoints**, so MWD saturates at 5000 µm and is a weighted sieve statistic
rather than a physical diameter. The model reproduces this same formula rather
than computing a mean of its own aggregate sizes — otherwise the two would not
be the same quantity.

Wet sieving follows Elliott (1986): 50 g soil submerged 5 min in deionised
water, then three sieves consecutively. Macroaggregate fractions by moving the
sieve 2 cm in and out of the water, 50 cycles in 2 minutes. The microaggregate
fraction is washed vigorously until the wash water runs clear — an adjustment
the authors made to remove clay completely from these loamy soils. All fractions
oven-dried at 50 °C. Day 0 is six replicates per soil.

| days present | soils with data |
|---|---|
| 0, 1, 4, 7, 13, 21 | all five |
| 2 | soils 1, 2 |
| 10 | soils 2, 3 |
| 16 | soils 1, 2, 3 |

**The Methods state destructive sampling at days 0, 1, 4, 7, 13 and 21 — six
dates. This file has nine.** The three extra dates are partial. The Methods
also disagree with the figures and with the regression `n` reported in Table 3;
see `../degryze_EJSS_2006_spec.md` §5.1. This file is the better guide, but the
provenance of days 2, 10 and 16 has not been established against the figures and
should be before they carry weight in a fit.

Cross-check: soil 3 reads 0.199 mm at day 1 and 1.014 at day 21, against
Table 2's 198 and 1008 µm. Table 2's values are also in `../degryze_soils.jl`
as `DEGRYZE_OBSERVED`, together with Table 3's fitted formation rates.

---

## `degryze_CO2_2006.csv` — cumulative respiration

**Units: µg-C per g soil, cumulative.** Columns `Time, Soil_1 … Soil_5`, five
time points: 0, 5, 10, 15, 21 days. The time column is read as text (`05`) and
parsed to Float64 in `run_degryze.jl`.

**The paper measured CO₂ every 3 hours** by GC-8A gas chromatograph coupled to
an automated system (Swerts *et al.*, 1995). This file is a five-point
digitisation of the published figure, not that series. It is adequate for a
cumulative comparison and would not support anything asking about respiration
*rate* — the model's flux panel has no measured counterpart for that reason.

Both sides of the CO₂ comparison are **total** soil respiration. There was no
unamended control, so residue-derived CO₂ cannot be separated from turnover of
native soil carbon, and no partition factor is applied on the model side.
See `../degryze_EJSS_2006_spec.md` §0a A3. Soil 5 respires 119 % of the
carbohydrate carbon added, which is the direct evidence that the measurement
includes native SOC.

---

## Soil properties (Table 1, p. 236)

In `../degryze_soils.jl`, not here. Reproduced for orientation; "Loam" is the
paper's column heading for the **silt** fraction.

| Soil | pH | CEC | Org C % | C/N | Sand | Silt | Clay | Texture | BD |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 4.4 | 5.6 | 1.78 | 11.5 | 53 | 40 | 7 | Sandy loam | 1.37 |
| 2 | 4.6 | 6.1 | 1.72 | 11.9 | 33 | 57 | 10 | Silt loam | 1.37 |
| 3 | 4.8 | 6.8 | 2.14 | 9.5 | 44 | 45 | 11 | Loam | 1.37 |
| 4 | 4.2 | 6.6 | 2.33 | 11.1 | 44 | 43 | 13 | Loam | 1.42 |
| 5 | 5.3 | 17.4 | 3.10 | 9.7 | 21 | 52 | 27 | Clay loam | 1.28 |

CEC in cmol_c/kg, BD in g/cm³. Soil 5's CEC is 2.6× the others — a
mineralogical difference the sand/silt/clay triple does not capture, and the
one place a per-soil free parameter could silently absorb the texture signal
during fitting. Note also that no soil is a clay: clay runs 7–27 %.
