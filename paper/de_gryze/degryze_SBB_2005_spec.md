# De Gryze SBB 2005 — incubation replication spec (+ citation untangle)

**Role in this project: the FITTING set** (six residue levels, an unamended
control, an 11-point time course, and the only fungal observable available).

**Written:** 2026-07-27
**Source:** De Gryze, S., Six, J., Brits, C. & Merckx, R. (2005) *A quantification
of short-term macroaggregate dynamics: influences of wheat residue input and
texture.* Soil Biology & Biochemistry **37**, 55–66.
doi:10.1016/j.soilbio.2004.07.024

**Companion document:** `degryze_EJSS_2006_spec.md` — the other incubation (the
validation set), including the MWD definition and the diff against
`run_degryze.jl`. This document supersedes its §7 and resolves its §5.3.

---

## 1. Three different De Gryze papers, and which is which

The README citation is a **chimera** — it welds paper A's title and authors onto
paper B's journal, volume and pages, with paper C's year. Resolving it lands you
on a third paper, which is what happened.

### A — the one you just added ✅

> **De Gryze, S., Six, J., Brits, C. & Merckx, R. (2005).** *A quantification of
> short-term macroaggregate dynamics: influences of wheat residue input and
> texture.* **Soil Biology & Biochemistry 37, 55–66.**
> doi:10.1016/j.soilbio.2004.07.024

**Verified from the PDF's own title page and running header** — `Soil Biology &
Biochemistry 37 (2005) 55–66`. So the volume is **37**, not 35. The "35, 55–66"
printed in the 2006 paper's reference list is a **typo in the published source**,
which is what I flagged last turn. The file you added is the right paper.

### B — the source of your CSV data ✅

> **De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R. (2006).**
> *Water repellence and soil aggregate dynamics in a loamy grassland soil as
> affected by texture.* **European Journal of Soil Science 57, 235–246.**
> doi:10.1111/j.1365-2389.2005.00733.x

Note the **five-author list** — Jassogne and Bossuyt are on this one and not on A.
This is where `degryze2006.csv` and `degryze_CO2_2006.csv` come from; I verified
the numbers against its Table 2 and Figure 3 last turn and they match exactly.

### C — the paper the bad citation resolves to ❌

**[UNVERIFIED — I could not open it; Wiley returned 403.]** My best candidate for
what you landed on is:

> De Gryze, S., Six, J. & Merckx, R. (2006). *Quantifying water-stable soil
> aggregate turnover and its implication for soil organic matter dynamics in a
> model study.* **European Journal of Soil Science**, doi:10.1111/j.1365-2389.2005.00760.x

It fits the trap perfectly: **same three-author subset** (De Gryze, Six, Merckx),
**same year (2006)**, **same journal**, and a title about aggregate *turnover* —
so a search on the README's author list + year + journal + "aggregate dynamics"
resolves here rather than to B. I could not confirm from the paywalled page that
it uses rare-earth-oxide tracers, so I'm taking your description of its contents
at face value rather than asserting it myself.

### Fix for `paper/de_gryze/README.md`

The current header block should read:

```markdown
**Data source:** De Gryze, S., Jassogne, L., Bossuyt, H., Six, J. & Merckx, R.
(2006). Water repellence and soil aggregate dynamics in a loamy grassland soil
as affected by texture. *European Journal of Soil Science*, **57**, 235–246.
doi:10.1111/j.1365-2389.2005.00733.x

**Companion experiment:** De Gryze, S., Six, J., Brits, C. & Merckx, R. (2005).
A quantification of short-term macroaggregate dynamics: influences of wheat
residue input and texture. *Soil Biology & Biochemistry*, **37**, 55–66.
doi:10.1016/j.soilbio.2004.07.024
```

Also still wrong in that README: "4-week (28-day)" → **21 days**; "5 soil textures
(sandy loam → clay)" → sandy loam → clay **loam**.

---

## 2. Why the SBB 2005 paper is the better replication target

**[MY ANALYSIS]** Having now read it: this paper is a substantially better match
to what your model does than the 2006 one, for four reasons.

1. **It has a 0 % residue control.** This directly resolves the open question from
   `degryze_replication_spec.md` §5.3 — see §5 below, where I quantify it. The
   2006 paper has no control at all.
2. **It measures fungal hyphal length** (m g⁻¹ soil) across a residue gradient.
   That is a direct observable for your Falconer fungal sub-model — currently you
   have no fungal data to constrain `F_i`/`F_n`/`F_m` at all.
3. **It sweeps residue input over six levels**, so aggregate formation per unit
   substrate is measured rather than inferred from one loading.
4. **It fits the exact model structure you are building** — aggregate formation
   rate proportional to respiration rate — and reports the fitted coefficients
   and turnover times (§4.4 below). That is a published benchmark to check
   against.

**The bridge between the two papers:** the SBB 2005 **1 wt % treatment is
4.43 mg C g⁻¹ soil**, which is *exactly* the 2006 loading (1.5 g / 150 g × 44.3 %).
Same site, same residue batch, same wet-sieving protocol. So the 1 % column of
SBB 2005 and the whole of EJSS 2006 are the same carbon input.

⚠ **But the soils are different.** Same meadow, different selection — see §3.1.
Do not carry parameters across without re-mapping texture.

---

## 3. Experimental spec — De Gryze et al. (2005), SBB 37, 55–66

Same site as the 2006 paper: fertilized sheep meadow, **K.U. Leuven Centre for
Animal Husbandry, Lovenjoel, Belgium**, **0–10 cm**.

### 3.1 Soils — Table 1, p. 56 (verbatim)

| Soil | % sand | % silt | % clay | Texture class | % OC |
|---|---|---|---|---|---|
| 1 | 63 | 30 | 7 | Sandy loam | 1.2 |
| 2 | 41 | 52 | 7 | Silt loam | 1.5 |
| 3 | 10 | 66 | 24 | Silty clay loam | 4.4 |

**[DERIVED] `f_clay_silt`:** soil 1 = **0.37**, soil 2 = **0.59**, soil 3 = **0.90**.

> ⚠ **These are not the 2006 soils.** The 2006 paper used five soils with 7–27 %
> clay and 1.72–3.10 % OC; here there are three with 7–24 % clay and 1.2–4.4 % OC.
> The silty clay loam here (10 % sand, 4.4 % OC) is more extreme than anything in
> the 2006 set. Note also that this table reports **% silt** honestly, unlike the
> 2006 Table 1 which labels the same column "Loam".

### 3.2 Pretreatment (p. 56, verbatim) — more aggressive than 2006

> "The soil was subsequently **crushed through a 250 and a 53 µm sieve**. This was
> done by vibrating the soil in a sieve together with large (1.5 cm diameter)
> rubber and small (0.5 cm diameter) agate balls. **This treatment minimized
> breakup of the sand material while breaking up all aggregates.** Native POM and
> clay material were removed from the sand fractions using **15 % H₂O₂ (24 h)** and
> **5 % sodium hexametaphosphate (shaking for 18 hrs)**. After rinsing on a 53 µm
> sieve, the rest of the organic material was removed by burning the fractions in
> a **muffle furnace at 550 °C for 3 h**. Ashes were removed by shaking the material
> in water and subsequent sieving. The powdery material < 53 µm was not further
> treated. After this treatment, the soil fractions (two sand fractions and the
> material < 53 µm) were reconstituted."

Abstract, p. 55: "**all structures > 53 µm were destroyed**".

**Consequence for the model IC.** This is a cleaner zero than the 2006 experiment:
the 2006 soils were crushed to 250 µm only, so their day-0 250–2000 µm mass was
reconstituted sand. Here the crush goes to **53 µm**, so day-0 macroaggregates are
genuinely ~0 — the paper reports WSA > 2000 µm starting at **3 %** and rising to
40 %. **If you want a clean "start from fully dispersed" initial condition, use
this experiment, not the 2006 one.**

### 3.3 Residue (p. 57, verbatim)

> "For the incubations, wheat residue (*Triticum aestivum* L.) was used as an
> organic substrate. **Nodes, leaves and grains were carefully removed. The
> remaining stems were then cut in pieces of 1–2 mm.** Carbon (C) and nitrogen (N)
> percentages were measured on a Vario Max total element analyzer (Elementar,
> Hanau, Germany) and were found to be **44.3 % C and 0.17 % N (C/N ratio is 261)**.
> Residue composition was measured on a Fibertec System M (Foss Tecator, Sweden)
> according to the Van Soest method. The material had **11 % soluble fraction, 3 %
> hemicellulose, 46 % cellulose, 11 % lignin, and 1 % inorganic fraction**."

✅ **Chemistry is identical to the 2006 paper** — same batch. Only the cut length
differs: **1–2 mm** here vs 0.5–2 mm in 2006.

### 3.4 Incubation 1 — residue gradient (p. 57, verbatim)

> "In a first incubation, **80 g** of each reconstituted soil were mixed with a
> range of wheat residue additions (**0, 0.5, 1, 1.5, 2 and 3 wt %**, corresponding
> with **0, 2.2, 4.4, 6.8, 8.8 and 13.3 mg C g⁻¹ soil**). Each soil × residue
> addition combination was **replicated six times**. Water was added to **65 % of
> water filled pore space**. Water uptake by the wheat residue was accounted for by
> **adding an extra gram of water for each gram of organic residue**. After mixing,
> the soil was transferred to **high-density PVC cylinders (diameter of 5 cm)** with
> a closed bottom. During transfer, a plunger was used to compress the soil
> material evenly and bring all soils to the following bulk densities: **0.97 g cm⁻³
> for the silty clay loam soil, 1.35 g cm⁻³ for the silt loam and 1.51 g cm⁻³ for
> the sandy loam.**"

On how those bulk densities were set — worth knowing, because it is iterative and
not a free choice:

> "These bulk densities were determined in an iterative process. As a first
> approximation, a **pedotransfer function was used to calculate the bulk density
> given the textural data (Saxton *et al.*, 1986)**. Using this bulk density value,
> the porosity and the amount of water corresponding to 65 % of water filled pore
> space were calculated. The soil was then put in cylinders and compressed. Soil
> volume and bulk density were measured. If the measured bulk density differed
> from bulk density used to calculate the porosity, the whole procedure was
> repeated starting from the measured bulk density."

> "The filled cylinders were then put in glass canning jars, **together with a test
> tube filled with water to maintain air humidity**. The jars were closed with lids
> having a septum for gas measurements. Soils were incubated for **22 days at
> 25 °C**. During the incubation, respiration was measured using a **PIR-2000R
> infrared gas analyzer (Horiba, Kyoto, Japan) at days 1, 2, 3, 5, 8, 12, 15, 17,
> 19 and 21.** On day 21, the soil material was removed from the cylinders by using
> a plunger to minimize disturbance. The soil was then gently crumbled along the
> natural planes of weakness, until all material passed an **8 mm sieve**. The
> crumbled soil was split into two parts: (1) 10 g of soil was put in falcon tubes
> and stored moist at 4 °C for microbial biomass measurements, (2) the remaining
> soil was put in aluminum pans to air dry for water stable aggregate separations."

> ⚠ Note the jars are **sealed with a septum**, not flushed daily as in the 2006
> experiment. CO₂ accumulates between measurements and O₂ is drawn down. Your
> constant-21 %-O₂ assumption is **less well justified here** than it was for 2006.
> The paper does not report headspace O₂ or jar volume.

### 3.5 Incubation 2 — time course (p. 57, verbatim)

> "In a second experiment, **24 cylinders** were filled with the **silt loam** soil
> and **1.5 wt %** of organic residue was added. The cylinders were prepared and
> incubated as described above. At **days 0, 1, 3, 5, 7, 9, 13, 15, 17, 19 and 21**,
> **two cylinders** were prepared as described above for aggregate separation."

**[MY ANALYSIS]** This is the one to fit trajectories against — it is the only
densely time-resolved aggregate series in either paper (11 dates), and Fig. 5
plots it with four candidate models overlaid. Single soil, single loading, so it
isolates the temporal dynamics from the texture and dose effects.

### 3.6 Aggregate separation (p. 57)

**Identical protocol to the 2006 paper** — Elliott (1986); 50 g submerged 5 min in
deionized water; three sieves → > 2000 / 250–2000 / 53–250 / < 53 µm; sieve moved
2 cm in and out for **50 cycles in 2 min**; micro fraction washed until the water
runs clear; all fractions oven-dried at **50 °C**; **six replicates sieved
immediately without incubation** for day 0.

> ⚠ **No MWD is reported in this paper.** The observable is **% water-stable large
> macroaggregates (> 2000 µm)**, in wt %. So this paper constrains the `[A%]` term
> of the 2006 MWD equation directly, and nothing else. That is arguably cleaner —
> it avoids the fixed-weight artefact entirely.

### 3.7 Fungal measurement (p. 57, verbatim) — the method behind the hyphal lengths

> "Fungal abundance and biomass was determined following the method of **Bloem
> *et al.* (1995)**. Ten grams of soil were taken and suspended in 100 ml of
> filter-sterilized (0.22 µm) water. This mixture was blended using a Waring lab
> blender and allowed to settle for 30 s. A 10 ml subsample was taken, which was
> used to make three dilutions (1:50, 1:100 and 1:200). From each sample and each
> dilution, a smear was prepared… Fungi were stained with **Calcifluor M2R**
> (fluorescent brightener)… **Fungal lengths were calculated by taking 30 pictures
> from each well** using a Zeiss Axiophot 2 epifluorescence microscope (200×
> magnification)… Lengths of fungal hyphae were measured **manually** from the
> images. Image analysis was done using… ImageJ."

**[MY ANALYSIS]** To compare your `F_i + F_n` to this you need hyphal radius and
carbon density: length = biomass_C / (π r² ρ_hyphal · f_C). Neither this paper
nor Falconer gives r or ρ. Typical literature values are r ≈ 1–2.5 µm — but pick
one, cite it, and treat the conversion as a stated assumption, because the answer
scales as r⁻².

---

## 4. Results — what to calibrate against

### 4.1 Respiration

- **Peak rate at day 2** for all three soils (p. 57). *(The 2006 paper reports
  days 3–6 — a real difference, possibly the finer crush or the sealed jars.)*
- Cumulative respiration rose linearly with residue added, P < 0.0001.

**Table 2, p. 59 — linear regressions vs residue added (g per 100 g soil):**

| Variable | Soil | Intercept ± SE | Slope ± SE |
|---|---|---|---|
| Cumulative respiration /µg C g⁻¹ soil | Sandy loam | **670 ± 100**\*\*\* | 380 ± 50\*\*\* |
| | Silt loam | **760 ± 170**\*\* | 660 ± 92\*\*\* |
| | Silty clay loam | **2300 ± 67**\*\* | 1400 ± 41\*\*\* |
| Water-stable macroaggregates /g 100 g⁻¹ soil | Sandy loam | 5.8 ± 3.8 | 10 ± 1.4\*\*\* |
| | Silt loam | −2.6 ± 3.8 | 13 ± 1.4\*\*\* |
| | Silty clay loam | 0.1 ± 2.8 | 14 ± 1.1\*\*\* |
| | **All three soils** | 1 ± 1.4 | **12 ± 0.7**\*\*\* |
| Fungal lengths /m g⁻¹ soil | Sandy loam | 12 ± 5.3\* | 20 ± 3.8\*\*\* |
| | Silt loam | 1.4 ± 2.3 | 6.3 ± 1.2\*\* |
| | Silty clay loam | 6.7 ± 1.2\*\*\* | 0.66 ± 1.7 (NS) |

\*\*\* P < 0.0001; \*\* P < 0.001; \* P < 0.05.

**The intercepts are the unamended controls.** That is the number the 2006 paper
does not give you.

### 4.2 Aggregation

- WSA > 2000 µm rose **from 3 % to 40 %** across the residue gradient (p. 59).
- **No significant texture effect on formation** (P = 0.36 across the three slopes)
  — despite clay ranging 7–24 %.
- Pooled formation rate: **12.0 ± 0.7 g aggregates g⁻¹ residue added**, equivalently
  **27.1 ± 1.6 g aggregates g⁻¹ C added** (p. 63).
- Field WSA, by contrast, *did* differ: silty clay loam 40 ± 1 % > silt loam
  25 ± 3 % > sandy loam 15 ± 3 % (p. 59).
- Fungal lengths, control → 3 % addition: sandy loam **5.2 ± 0.5 → 67.8 ± 1.7** m g⁻¹;
  silt loam **1.6 ± 0.3 → 20.0 ± 4.4**; silty clay loam **no clear increase**, 7.9 ± 0.6
  at 3 % (p. 60).

### 4.3 Correlations — Table 3, p. 61

| Soil | WSA vs cumulative respiration (r) | WSA vs fungal length (r) |
|---|---|---|
| Sandy loam | 0.91 (P < 0.0001, n = 15) | 0.90 (P < 0.0001, n = 15) |
| Silt loam | 0.89 (P < 0.0001, n = 15) | 0.85 (P = 0.0001, n = 15) |
| Silty clay loam | 0.91 (P < 0.0001, n = 25) | **0.25 (P = 0.27, NS, n = 22)** |
| All three together | 0.47 (P = 0.0005, n = 53) | 0.42 (P = 0.0026, n = 50) |

**Key conclusion (p. 65, verbatim):** "**Cumulative respiration correlated well with
water-stable macroaggregates for all three soils. Neither cumulative respiration
nor fungal lengths was a good predictor for water-stable macroaggregates when data
from all three soils were combined.**"

And on fungi specifically (p. 64): "In the sandy loam and the silt loam, aggregate
formation was much more related to fungal length than in the silty clay loam.
Fungal hyphae are able to grow more in macroporous soils, which explains the
greater fungal lengths in the sandy soil compared to the other two soils.
**Cumulative respiration was a better predictor for macroaggregation than fungal
length.**"

### 4.4 The four aggregation models — Table 4, p. 62

This is the published benchmark. Their Model 4 is structurally what your model is
trying to be.

**Model 2 (sigmoidal), eq. (1):** `M = k₁ / (1 + exp(−(t−t₀)/k₂))`

**Models 3 & 4, eq. (2):** `dM/dt = f_t·U − b_t·M`, with `M + U = 100 %` (wt %).

- **Model 3** — `f_t`, `b_t` both constant → eq. (3),
  `M = f_t/(f_t+b_t) · (1 − exp(−(f_t+b_t)t))`
- **Model 4** — eq. (4): `f_t = k₃·C(t − t₀)`, `b_t = k₄`, where **C is the
  respiration rate** and **t₀ = 4 days** is a lag between microbial activity and
  aggregate formation.

**Fits to the incubation-2 silt loam time course (Fig. 5):**

| Model | Description | Adj. R² | Parameters |
|---|---|---|---|
| 1 | Linear | 0.85 | a = 0.65 ± 0.06 % day⁻¹; y₀ = 2.04 ± 0.67 |
| 2 | Sigmoidal | 0.83 | k₁ = 13.9 ± 1.1; k₂ = 3.59 ± 0.92; T₀ = 6.7 ± 1.1 |
| 3 | Constant formation | 0.84 | f_t = 1.34 ± 0.25; b_t = 5.88 ± 2.61 |
| 4 | **Formation ∝ respiration** | **0.86** | k₃ = 1.9×10⁻³ ± 5.8×10⁻⁴; k₄ = 1.8×10⁻² ± 2.8×10⁻² |

**Model 4 fitted across the residue gradient, per soil:**

| Soil | k₃ | k₄ | Turnover time | Adj. R² |
|---|---|---|---|---|
| Sandy loam | 6.4 × 10⁻³ | 1.7 × 10⁻² | **59 d** | 0.69 |
| Silt loam | 4.2 × 10⁻³ | 2.4 × 10⁻² | **42 d** | 0.62 |
| Silty clay loam | 1.6 × 10⁻³ | 1.6 × 10⁻² | **63 d** | 0.63 |

**Their verdict on Model 3 (p. 62):** the constant-formation model gives
`b_t = 5.88 day⁻¹`, i.e. a turnover time of **0.17 days** — "unrealistically
fast". Model 4 gives a mean breakdown parameter of 1.8 × 10⁻² day⁻¹ → **56 days**,
with a lower bound of 22 days. Their conclusion (p. 64): "**assuming a constant
aggregate formation rate may lead to skewed results**… if it is assumed that
aggregate formation is proportional to microbial activity, this slowing can be
explained by a decrease in aggregate formation."

**[MY ANALYSIS] Two things to take from this.** First, `k₃` differs ~4× across
soils while `k₄` is nearly constant (1.6–2.4 × 10⁻²) — so in their fitting, texture
enters through **formation efficiency per unit respiration**, not through
breakdown. If your model has a texture-dependent binding term, that is where it
should live. Second, over 21 days breakdown is nearly irrelevant: with
k₄ ≈ 0.018 day⁻¹, less than a third of formed aggregate turns over in the window.
Their own framing, p. 62: "**during the first 3 weeks of incubation, aggregate
breakdown is minimal**." You can probably run formation-only for this comparison
and lose very little.

---

## 5. This resolves the CO₂ attribution problem

`degryze_replication_spec.md` §5.3 flagged that the 2006 respiration curves have
no unamended control, so you cannot tell how much of the signal is residue-derived.
**Table 2's intercepts and slopes answer this.**

**[DERIVED]** Evaluating the Table 2 regressions at **1 wt % residue** — the same
loading as the whole 2006 experiment:

| Soil | Native (intercept) | Residue-derived (slope × 1) | Total | **Residue share** |
|---|---|---|---|---|
| Sandy loam | 670 | 380 | 1050 | **36 %** |
| Silt loam | 760 | 660 | 1420 | **46 %** |
| Silty clay loam | 2300 | 1400 | 3700 | **38 %** |

*(all µg C g⁻¹ soil over 3 weeks)*

**Only about 36–46 % of cumulative CO₂ at this loading is residue-derived.** The
totals (1050–3700) also bracket the 2006 measurements (1448–3157) nicely, which
cross-validates both experiments.

**Implication for your calibration.** The note at `run_degryze.jl:181–186` records
a "~1.9× CO₂ overshoot" against the 2006 data. If your model respires only the
added POM, then the correct comparison target is not the full 2006 curve but
roughly **40 % of it** — which would turn a 1.9× overshoot into something closer
to agreement, or even an undershoot. **Do not tune `R_P_max` down until this is
settled.** Combined with the 20 °C → 25 °C correction (which pushes respiration
*up* ~40 %), the two corrections work in opposite directions and could easily
account for the whole discrepancy without touching any biological parameter.

There is also a mechanistic hint worth heeding, p. 65: "**Possibly, only the
decomposition of added residue (and not of native residue) will affect aggregate
formation.**" If you adopt that, aggregate formation should be driven by
POM-derived respiration while total CO₂ (for comparison with data) needs the
native SOC term added on top.

---

## 6. Two quotes that endorse your modelling approach

**On the high residue loadings — this justifies the single-POM-particle domain
tessellation in `run_degryze.jl` (p. 62, verbatim):**

> "Our residue additions may seem unreasonably high (ranging from 2.22 to
> 13.3 mg C g⁻¹ soil or 13–78 ton C ha⁻¹). However, it is found that **fresh organic
> substrate forms hot spots of microbial activity and there, a new soil aggregate
> is developed** (Guggenberger *et al.*, 1999). This can occur near a root tip,
> where exudation is ample or locally next to a piece of organic matter
> incorporated in the soil by tillage or soil fauna. **With our high additions, we
> aimed to simulate the situation near such hot spots.**"

**On where fungi live, relevant to your fungal diffusion coefficients (p. 63,
verbatim):**

> "**Whereas fungi predominantly proliferate in larger pore spaces occurring
> between microaggregates of size 53–250 µm; bacteria reside in smaller pores
> within microaggregates.**"

---

## 7. Consolidated experimental comparison

| | **SBB 2005** | **EJSS 2006** |
|---|---|---|
| Soils | 3 (sandy loam, silt loam, silty clay loam) | 5 (sandy loam → clay loam) |
| Clay range | 7–24 % | 7–27 % |
| OC range | 1.2–4.4 % | 1.72–3.10 % |
| Crushed to | **53 µm** | 250 µm |
| Soil mass | 80 g | 150 g |
| Residue | **0, 0.5, 1, 1.5, 2, 3 wt %** | 1 wt % only |
| Residue C | 0–13.3 mg C g⁻¹ | 4.43 mg C g⁻¹ |
| Residue cut | 1–2 mm | 0.5–2 mm |
| Residue chemistry | 44.3 %C, 0.17 %N, C:N 261 | **identical** |
| Water | **65 % WFPS** + 1 g per g residue | 60 % WFPS + 1.5 ml |
| Bulk density | 0.97 / 1.35 / 1.51 | 1.28–1.42 |
| Cylinder | PVC, 5 cm ⌀ | stainless steel, 5 cm ⌀ |
| Duration | **22 d** at 25 °C | 21 d at 25 °C |
| Headspace | **sealed, septum** | **flushed daily with air** |
| Respiration | IRGA, 10 dates | GC, every 3 h |
| Aggregate observable | **% WSA > 2000 µm** | **MWD** (4-class weighted) |
| Time-resolved aggregates | 11 dates (silt loam, 1.5 %) | ~8 dates × 5 soils |
| Fungal length | ✅ **yes** | ✗ |
| Unamended control | ✅ **yes** | ✗ |
| Wet sieving | Elliott (1986), identical | Elliott (1986), identical |

**[MY ANALYSIS] Suggested use of the two together:**

- **Fit** on SBB 2005 incubation 2 (silt loam, 1.5 %, 11 time points) for temporal
  dynamics, and on the six-level residue gradient for dose response.
- **Constrain the fungal sub-model** on the SBB 2005 hyphal lengths — this is the
  only fungal observable available anywhere in this project.
- **Validate** on EJSS 2006, which is genuinely out-of-sample: different soils,
  different loading regime, and a different aggregate statistic (MWD).
- Use SBB 2005 Table 4 `k₃`/`k₄` and the 40–60 day turnover times as sanity
  bounds on whatever your emergent formation/breakdown rates turn out to be.
