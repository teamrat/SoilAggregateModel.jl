"""
    maoc.jl

Mineral-associated organic carbon (MAOC) sorption/desorption kinetics.

Implements smooth switching between rate-limited sorption and desorption toward
local equilibrium determined by Langmuir-Freundlich isotherm.

# Dependencies
Requires: SoilProperties

# Critical C usage
- M_eq uses C_eq (equilibrium sorbed concentration), NOT C or C_aq
- Compute C_eq = k_d · C_aq = k_d·C / (θ + ρ_b·k_d) ONCE per node
- S_C coupling is `-J_M`, with NO conversion factor. An earlier version of this
  file insisted on `-J_M·(θ + ρ_b·k_d)/k_d`; that factor was removed by
  REFERENCE.md §26 erratum 2 in February and the claim survived here.

# Softplus regularization
- Replaces max(0, x) with smooth C∞ approximation (see math_utils.jl)
- Prevents discontinuous derivatives at M = M_eq
- Uses numerically stable form to avoid overflow
"""

"""
    sorption_capacity(θ, ρ_b, k_d) -> θ + ρ_b·k_d

Denominator of the instantaneous DOC partition. `C` is total soluble carbon per
bulk volume and splits as `C = θ·C_aq + ρ_b·k_d·C_aq`, so this factor converts
between the two.

**The only definition.** It appeared at six sites in `src/` with no owner, which
is how it came to be spelled two different ways.
"""
sorption_capacity(θ::Real, ρ_b::Real, k_d::Real) = θ + ρ_b * k_d

"""
    C_aqueous(C, θ, soil) -> C_aq

Aqueous DOC concentration [µg per mm³ of water] from total soluble carbon per
bulk volume. `C_eq = k_d · C_aq` is the sorbed concentration the isotherm sees —
NOT `C` or `C_aq` directly.
"""
C_aqueous(C::Real, θ::Real, soil::SoilProperties) =
    C / sorption_capacity(θ, soil.ρ_b, soil.k_d_eq)

"""
    M_eq_langmuir_freundlich(C_eq::Real, M_max::Real, k_L::Real, n_LF::Real)

Equilibrium MAOC concentration from Langmuir-Freundlich isotherm.

# Formula (manuscript Eq. 264)
    M_eq = M_max · (k_L·C_eq)^n_LF / [1 + (k_L·C_eq)^n_LF]

# Arguments
- `C_eq`: Equilibrium sorbed DOC concentration [μg/mm³ solid]
- `M_max`: Maximum sorption capacity [μg/mm³]
- `k_L`: Langmuir affinity constant [mm³/μg]
- `n_LF`: Freundlich heterogeneity exponent [-]

# Returns
- Equilibrium MAOC [μg/mm³]

# Notes
- **CRITICAL**: Uses C_eq = k_d·C_aq, NOT C or C_aq directly
- n_LF < 1: heterogeneous sites (typical for soil minerals)
- n_LF = 1: recovers standard Langmuir isotherm
- Pass `maoc_capacity(soil)` for `M_max`. Do not recompute it.
- As C_eq → 0: M_eq → 0
- As C_eq → ∞: M_eq → M_max
"""
function M_eq_langmuir_freundlich(C_eq::Real, M_max::Real, k_L::Real, n_LF::Real)
    k_L_C = k_L * C_eq
    k_L_C_n = k_L_C^n_LF

    M_max * k_L_C_n / (1.0 + k_L_C_n)
end

"""
    maoc_capacity(soil::SoilProperties) -> M_max

Maximum MAOC concentration the mineral phase can hold [µg-C/mm³]:

    M_max = k_ma · f_clay_silt · ρ_b

**This is the only definition of `M_max` in the package.** Everything that needs
the capacity — initialisation, the solver's isotherm, post-processing — calls
this. There is no `SoilProperties.M_max` field.

Units close because `k_ma` is DIMENSIONLESS (g-C per g of clay+silt):

    [-] · [-] · [µg/mm³] = [µg-C/mm³]

# History (REFERENCE.md §26 erratum 12)

Until 2026-07-29 there were two: `initial_conditions.jl` computed this formula
while `reactions.jl`, `api.jl` and `derived.jl` read a `SoilProperties.M_max`
field defaulting to 10.0. For De Gryze soil 3 that was 368 against 10 — so the
state was initialised above the ceiling the solver then enforced, and desorbed
from t=0 with no equilibrium to reach. The 368 was itself wrong: `k_ma = 0.48`
was Georgiou's 48 mg-C/g with the mg→g conversion dropped, so it read 10× high
in a formula that needs a mass fraction. The guard in `create_initial_state`
compared against the texture value, not the solver's, so it never fired.

Both arms of the same defect: one quantity, two implementations, and the
disagreement was absorbed into a parameter (`κ_d_ref` was lowered tenfold in
`paper/de_gryze/degryze_config.jl` to suppress the resulting desorption flux).
"""
maoc_capacity(soil::SoilProperties) = soil.k_ma * soil.f_clay_silt * soil.ρ_b

"""
    J_M(M::Real, M_eq::Real, κ_s_T::Real, κ_d_T::Real, ε_maoc::Real)

Net MAOC flux with smooth switching between sorption and desorption.

# Formula (manuscript Eq. 237)
    J_M = κ_s(T)·φ_ε(M_eq - M) - κ_d(T)·φ_ε(M - M_eq)

where φ_ε is softplus regularization.

# Arguments
- `M`: Current MAOC concentration [μg/mm³]
- `M_eq`: Equilibrium MAOC (from M_eq_langmuir_freundlich) [μg/mm³]
- `κ_s_T`: Sorption rate at current T [1/day]
- `κ_d_T`: Desorption rate at current T [1/day]
- `ε_maoc`: Softplus smoothing width [μg/mm³] (default 0.01)

# Returns
- Net MAOC flux [μg-C/mm³/day]

# Notes
- When M < M_eq: sorption dominates (J_M > 0, M increases)
- When M > M_eq: desorption dominates (J_M < 0, M decreases)
- κ_s_T, κ_d_T already include Arrhenius factors
- Softplus ensures C∞ smoothness near M = M_eq
- S_C coupling is `-J_M`, no conversion factor (REFERENCE §26 erratum 2).
"""
function J_M(M::Real, M_eq::Real, κ_s_T::Real, κ_d_T::Real, ε_maoc::Real)
    sorption = κ_s_T * softplus(M_eq - M, ε_maoc)
    desorption = κ_d_T * softplus(M - M_eq, ε_maoc)

    sorption - desorption
end
