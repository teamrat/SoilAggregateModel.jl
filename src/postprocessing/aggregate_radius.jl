# aggregate_radius.jl
# Aggregate stability radius from binding agent profiles

"""
    wet_sieving_stress(; L_s, f_s, μ, ν_w) -> (τ_w, δ_s, v_s, Ω)

Hydrodynamics of the Eijkelkamp wet-sieving apparatus. Defaults are the standard
settings: 13 mm stroke, 34 oscillations per minute, water at 20 °C.

- `v_s = π f_s L_s` — peak sieving velocity [mm/s], ≈ 23.1
- `δ_s = √(2ν_w/Ω)` — Stokes boundary layer thickness [mm], ≈ 0.751
- `τ_w = √2 μ v_s / δ_s` — peak wall shear stress [Pa], ≈ 0.0437

The √2 is the phase difference between velocity and wall stress in the
oscillatory Stokes solution (Batchelor §4.3). `v_s` is in mm/s and `μ` in
Pa·s = kg/(m·s), so both `v_s` and `δ_s` are converted to metres for τ_w.

Extracted from `compute_r_agg` so diagnostics use these numbers rather than
restating them (CLAUDE.md §8).
"""
function wet_sieving_stress(; L_s::Real = 13.0, f_s::Real = 34.0/60.0,
                              μ::Real = 1.002e-3, ν_w::Real = 1.004)
    v_s = π * f_s * L_s
    Ω   = 2π * f_s
    δ_s = sqrt(2.0 * ν_w / Ω)
    τ_w = sqrt(2.0) * μ * (v_s * 1e-3) / (δ_s * 1e-3)
    return (τ_w = τ_w, δ_s = δ_s, v_s = v_s, Ω = Ω)
end

"""
    critical_binding(soil::SoilProperties)     -> G_ref
    critical_binding(soil::SoilProperties, r)  -> G_c(r)

Binder concentration [μg/mm³] an aggregate must carry to survive wet sieving.

    G_ref  = τ_w · d_32 / κ_b                    at r = δ_s
    G_c(r) = G_ref · (r / δ_s)^p_Gc

`r` may be a scalar or a vector; the return has the same shape.

Rises with `d_32`: a coarser soil packs fewer bonded contacts into the same
failure area and needs more binder to hold (REFERENCE.md §4.4).

# The size dependence

`p_Gc = 0` is a threshold that does not depend on aggregate size — a 2 mm
aggregate faces the same threshold as a 0.1 mm one. That is what this model did
before 2026-07-29 and it remains the package default, so `critical_binding(soil)`
and `critical_binding(soil, r)` return the same number unless a problem opts in.

`p_Gc > 0` makes a larger aggregate need more binder per unit volume, i.e. be
weaker at fixed binder concentration.

**The exponent is fitted, not derived.** Two mechanisms are known to act in this
direction and neither fixes a value:

  * *Flaw statistics.* A larger body samples more of the flaw-size distribution,
    so its strength falls with size (Weibull). For soil aggregates specifically,
    measured tensile strength falls with aggregate diameter — Dexter and Watts
    report roughly `Y ∝ d^-0.5` to `d^-1`, which is `p_Gc` between 0.5 and 1.
  * *Boundary-layer protrusion.* `δ_s ≈ 0.751 mm` sits inside the size range the
    sieve series resolves. An aggregate smaller than δ_s sits inside the Stokes
    layer and sees a velocity difference growing with its own radius; one larger
    protrudes into the free stream. The flat-wall closure that gives τ_w does not
    hold across that transition, which is why δ_s and not some fitted length is
    the pivot: at `r = δ_s` the size-dependent form and the flat-wall form agree
    exactly, whatever `p_Gc` is.

A third possibility — that the disruptive stress genuinely is size-independent,
as the Stokes-drag balance for a sphere in uniform shear gives — is `p_Gc = 0`,
which is inside the family. The parameter is not a correction bolted onto a
correct model; it spans a structural question the model does not settle.

**`p_Gc < 0` is a form that was already rejected.** Before the 2026-07 rewrite
this function tested `[F_i + w_E·E]·r ≥ G_c`, which is `p_Gc = -1`: larger
aggregates inherently more stable. It was removed because its derivation
distributed Stokes drag over the aggregate surface — the wrong failure model —
and because it made the criterion easier to satisfy the further out it was
tested, so nothing stopped the aggregate growing. `p_Gc > 0` is the opposite
sign. Negative values are reachable because the family is continuous, not
because the question is open.

# Why it was added

Measured 2026-07-29 at `p_Gc = 0`: `r_agg` jumps on day 3 and is then flat for
42 days, because the crossing sits where the binder profile is flat while the
binder accumulates 6-fold nearer the POM. A threshold rising with `r` moves the
crossing inward, into the region still accumulating — the only way this model
can produce a radius that grows for 21 days. The MATLAB precursor had the same
idea as a commented-out `strength ./ x`.
"""
function critical_binding(soil::SoilProperties)
    return wet_sieving_stress().τ_w * soil.d_32 / soil.κ_b
end

function critical_binding(soil::SoilProperties, r)
    G_ref = critical_binding(soil)
    # p_Gc == 0 short-circuits to the flat threshold WITHOUT evaluating x^0.
    # `one.(r)` carries the shape of r; multiplying by 1.0 is exact, so the
    # default path is bit-identical to the pre-2026-07-29 scalar threshold.
    soil.p_Gc == 0.0 && return G_ref .* one.(r)
    return G_ref .* (r ./ wet_sieving_stress().δ_s) .^ soil.p_Gc
end

"""
    compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)

Stable aggregate radius [mm] from oscillatory wall shear stress criterion.

The aggregate is the region where biologically produced binding agents
generate sufficient cohesive strength to resist hydrodynamic disaggregation
during wet sieving. The criterion compares the local cohesive strength τ_c(r)
against the wall shear stress τ_w from the oscillatory Stokes boundary layer
(Batchelor, 1967, §4.3, §5.14):

    τ_c(r) = (κ_b / d_32) · [F_i(r) + w_E·E(r)]  ≥  τ_w(r)

so the criterion is a concentration threshold:

    F_i(r) + w_E·E(r)  ≥  G_c(r) = (τ_w · d_32 / κ_b) · (r/δ_s)^p_Gc

With `p_Gc = 0` — the default — the threshold does not depend on position and
this reduces to the pre-2026-07-29 form `G_c = τ_w·d_32/κ_b`, exactly. See
`critical_binding` for what the exponent means and why it is fitted.

where:
- μ: dynamic viscosity of water [Pa·s]
- v_s = π·f_s·L_s: maximum sieving velocity [mm/s]
- δ_s = √(2·ν_w/Ω): Stokes boundary layer thickness [mm]
- Ω = 2π·f_s: angular frequency [rad/s]
- κ_b: specific binding strength per unit binder carbon [Pa·mm/(μg/mm³)]
- w_E: binding weight of EPS relative to insulated hyphae, per unit carbon [-]
- d_32: Sauter mean particle diameter [mm] (`sauter_from_texture`)

# Why strength scales as 1/d_32

Failure separates particles at their bonded contacts. Particles of size d tile a
failure surface, so the number of bonded contacts crossing unit area scales as
1/d². Strength is that count times the force per bond. Where the bond force
scales linearly with particle size — the geometry of a binder bridge held at a
contact, as for a capillary bridge — strength ∝ (1/d²)·d = 1/d.

At fixed binder concentration a FINER soil is therefore stronger: it packs more
contacts into the same failure area. Equivalently a coarser soil needs a higher
binder concentration to hold, so G_c ∝ d_32. This is the form derived and tested
against three sand sizes by Albalasmeh and Ghezzehei (2014), whose model makes
aggregate stability depend on soil texture together with the strength, density
and mass fraction of the binding organic matter.

The Sauter (surface-weighted) mean is the correct average here because the bond
count follows interfacial area; a geometric or arithmetic mean is dominated by
the coarse fraction, which supplies the fewest contacts per unit mass.

**κ_b and w_E are fitted.** The argument above fixes the FORM, not the
coefficients.

# Arguments
- `record::OutputRecord`: State snapshot at a specific time
- `grid::GridInfo`: Grid geometry
- `params::ParameterSet`: Biological and soil parameters

# Returns
- `r_agg::Float64`: Stable aggregate radius [mm]
  - Returns r_grid[1] (the POM surface) if the criterion is met at no node —
    the bare residue particle, which is still sieve-retained
  - Otherwise the interpolated crossing outward of the outermost passing node

# Notes
- Scans ALL nodes and takes the outermost passing one (no early break),
  because binding agent profiles can be non-monotonic.
- The flat-wall approximation is accurate for aggregates with diameter > ~3 mm
  (δ_s ≈ 0.75 mm). For smaller aggregates, curvature effects cause the model
  to overestimate stability. Because δ_s falls INSIDE the size range that
  standard sieve series resolve, this is a limitation of the closure rather
  than a negligible correction — it is not absorbed anywhere. `p_Gc` is a
  one-parameter empirical stand-in for that size dependence; it does not derive
  the curvature correction, and calling it one would be overclaiming.

# References
- Batchelor, G. K. (1967). An Introduction to Fluid Dynamics. Cambridge
  University Press. §4.3 (oscillating plate solution), §5.14 (extension to
  oscillating rigid bodies).
- Albalasmeh, A. A. and Ghezzehei, T. A. (2014). Interplay between soil drying
  and root exudation in rhizosheath development. Plant and Soil 374, 739–751.
  doi:10.1007/s11104-013-1910-y
- Sauter, J. (1926). Die Größenbestimmung der in Gemischnebeln von
  Verbrennungskraftmaschinen vorhandenen Brennstoffteilchen. VDI-Verlag, Berlin.
"""
function compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)
    # Threshold at every node. A vector even when p_Gc = 0, so there is one code
    # path and the size-dependent case cannot drift from the flat one.
    G_c = critical_binding(params.soil, grid.r_grid)

    # State variables
    F_i = record.state.F_i
    E = record.state.E

    # Compute binding concentration at all nodes. w_E weights EPS carbon against
    # insulated hyphal carbon; it is a fitted parameter, not a derived ratio.
    #
    # When domains overlap (ω > 1) this sum carries two scalings. Pre-existing
    # background binding agents are apportioned across the domains sharing that
    # soil, so they are reduced by ω. Residue-derived binding agents come from an
    # undiluted surface flux and near the residue are set by release against
    # diffusion, largely independent of domain size. Both are correct: newly
    # produced EPS adds to what the soil already held, so the pre-existing part
    # must be distributed rather than counted once per domain. G_c is a single
    # physical threshold, so the comparison is exact only where one contribution
    # dominates — the residue-derived one, once decomposition is underway. κ_b
    # carries the residual and is fitted. REFERENCE.md §4.4.
    binding = F_i .+ params.soil.w_E .* E

    # The criterion is `excess ≥ 0`. Writing it as a single field rather than
    # comparing two makes the sub-grid crossing below a plain root of `excess`,
    # which is the same arithmetic whether or not the threshold varies with r.
    excess = binding .- G_c

    # Floor at the POM surface, then take the outermost node that passes.
    #
    # r_grid[1] IS the POM surface, so an aggregate can never be smaller than
    # its own core. When no node passes the criterion the result is the bare
    # residue particle, not "nothing" — a residue fragment is still retained by
    # a sieve whether or not binding agents have built a shell around it.
    # This matches the MATLAB precursor, which sets radius = x(1) in that case
    # (aggregation_model_v8.m; see dev_notes/MATLAB_aggregation_logic.md).
    #
    # Scanning outward-in and stopping at the first pass selects the same node
    # as MATLAB's find(strength > threshold, 1, 'last'): a sub-threshold zone in
    # the interior does not terminate the aggregate, because destructive forces
    # act from the outside and a strong outer shell can hold the interior
    # together.
    # Sub-grid interpolation of the crossing, NOT the node itself.
    #
    # Returning grid.r_grid[i] quantises r_agg in units of h, which makes every
    # quantity derived from it — r_agg, the sieve class, MWD — a STEP function
    # of every parameter. Measured 2026-07-29: in a ±10 % sensitivity sweep four
    # unrelated parameters returned |E_MWD| = 0.1950 exactly (one node) and
    # twelve returned exactly 0. A step-function objective has zero gradient
    # almost everywhere and is undefined at the jumps, so no optimiser can work
    # on it. See docs/BACKLOG.md item 6.
    #
    # `binding` is a continuous field, so the threshold crossing has a
    # well-defined sub-cell position. Linear interpolation between the outermost
    # passing node and its outward neighbour makes r_agg continuous and
    # differentiable in the parameters, at the cost of assuming `binding` is
    # linear across one cell — the same assumption the spatial discretisation
    # already makes.
    r_agg = grid.r_grid[1]
    for i in grid.n:-1:1  # scan from outer to inner
        if excess[i] >= 0.0
            if i == grid.n
                # Crossing is at or beyond the domain edge; nothing to
                # interpolate against.
                r_agg = grid.r_grid[grid.n]
            else
                # excess[i] >= 0 > excess[i+1] by construction: i is the
                # OUTERMOST passing node, so i+1 failed. The root of the linear
                # segment through the two is where the aggregate ends.
                denom = excess[i] - excess[i+1]
                frac  = denom > 0.0 ? excess[i] / denom : 0.0
                frac  = clamp(frac, 0.0, 1.0)
                r_agg = grid.r_grid[i] + frac * (grid.r_grid[i+1] - grid.r_grid[i])
            end
            break  # outermost passing node
        end
    end

    return r_agg
end
