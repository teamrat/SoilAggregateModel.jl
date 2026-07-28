# aggregate_radius.jl
# Aggregate stability radius from binding agent profiles

"""
    compute_r_agg(record::OutputRecord, grid::GridInfo, params::ParameterSet)

Stable aggregate radius [mm] from oscillatory wall shear stress criterion.

The aggregate is the region where biologically produced binding agents
generate sufficient cohesive strength to resist hydrodynamic disaggregation
during wet sieving. The criterion compares the local cohesive strength τ_c(r)
against the wall shear stress τ_w from the oscillatory Stokes boundary layer
(Batchelor, 1967, §4.3, §5.14):

    τ_c(r) = (κ_b / d_32) · [F_i(r) + w_E·E(r)]  ≥  τ_w

so the criterion is a concentration threshold, independent of position r:

    F_i(r) + w_E·E(r)  ≥  G_c = τ_w · d_32 / κ_b

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
  - Otherwise r_grid[i] of the outermost node satisfying the criterion

# Notes
- Scans ALL nodes and takes the outermost passing one (no early break),
  because binding agent profiles can be non-monotonic.
- The flat-wall approximation is accurate for aggregates with diameter > ~3 mm
  (δ_s ≈ 0.75 mm). For smaller aggregates, curvature effects cause the model
  to overestimate stability. Because δ_s falls INSIDE the size range that
  standard sieve series resolve, this is a limitation of the closure rather
  than a negligible correction — it is not absorbed anywhere.

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
    # Standard wet sieving parameters (Eijkelkamp apparatus)
    L_s = 13.0          # stroke length [mm]
    f_s = 34.0 / 60.0   # frequency [Hz] (34 oscillations/min)
    μ = 1.002e-3         # dynamic viscosity of water at 20°C [Pa·s]
    ν_w = 1.004          # kinematic viscosity of water at 20°C [mm²/s] = 1.004e-6 m²/s

    # Maximum sieving velocity [mm/s]
    v_s = π * f_s * L_s  # ≈ 23.1 mm/s

    # Angular frequency [rad/s]
    Ω = 2π * f_s  # ≈ 3.56 rad/s

    # Stokes boundary layer thickness [mm]
    δ_s = sqrt(2.0 * ν_w / Ω)  # ≈ 0.75 mm

    # Wall shear stress [Pa] — peak amplitude from exact flat-plate solution
    # τ_w = √2 · μ · v_s / δ_s
    # (the √2 arises from the phase difference between velocity and stress
    #  at the wall in the oscillatory Stokes solution; Batchelor §4.3)
    #
    # NOTE: v_s is in mm/s but μ is in Pa·s = kg/(m·s).
    # Convert v_s to m/s and δ_s to m for consistent Pa units:
    τ_w = sqrt(2.0) * μ * (v_s * 1e-3) / (δ_s * 1e-3)  # [Pa]

    # Specific binding strength and texture length scale
    κ_b  = params.soil.κ_b     # [Pa·mm/(μg/mm³)]
    d_32 = params.soil.d_32    # [mm] Sauter mean particle diameter

    # Critical threshold [μg/mm³]. Rises with d_32: a coarser soil packs fewer
    # bonded contacts into the same failure area and needs more binder to hold.
    G_c = τ_w * d_32 / κ_b

    # State variables
    F_i = record.state.F_i
    E = record.state.E

    # Compute binding concentration at all nodes. w_E weights EPS carbon against
    # insulated hyphal carbon; it is a fitted parameter, not a derived ratio.
    binding = F_i .+ params.soil.w_E .* E

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
    r_agg = grid.r_grid[1]
    for i in grid.n:-1:1  # scan from outer to inner
        if binding[i] >= G_c
            r_agg = grid.r_grid[i]
            break  # outermost passing node
        end
    end

    return r_agg
end
