# Julia Implementation Deviations from Falconer et al. (2005, 2008) / MATLAB

Two deviations were discovered between the Julia implementation and the 
Falconer/MATLAB formulation during debugging of fungal dynamics. Both affect 
the fungal subsystem and both are now corrected.

---

## Deviation 1: Mobile Fungi Diffusivity (D_Fm)

### Falconer (2005) Eq. 2.1 / MATLAB line 394
Translocation of mobile biomass requires an established hyphal network. The 
effective diffusivity depends on local network density:

**Falconer (2005):**
$$D_n = \begin{cases} 10^{-7} D_b & n \leq n_0 \\ D_b & n > n_0 \end{cases}$$

Mobile biomass diffusivity switches from negligible to full based on a threshold 
concentration — transport requires connected pathways.

**MATLAB (single_aggregate_beta.m, line 394):**
```matlab
D_Fm_eff = D_Fm * (u(3)+u(4))/(u(3)+u(4) + 10*(Fi_min + Fn_min));
```

Smooth Michaelis-Menten form: D_Fm scales with local sessile fungal biomass 
(F_i + F_n), reaching half-maximum when the network density equals 10× the 
minimum viable fungal concentration.

### Julia (before correction)
```julia
D_Fm_val = D_Fm0_ref * f_Arrhenius   # spatially uniform, no network dependence
```

D_Fm was constant everywhere — mobile fungi could translocate even where no 
hyphal network existed. This is unphysical: translocation is internal to hyphae, 
so where there are no hyphae, there is no translocation pathway.

### Julia (corrected)
```julia
function D_eff_fungi_mobile(D_Fm0_ref, f_Arrhenius, F_n, F_i, K_D)
    F_network = F_n + F_i
    D_Fm0_ref * f_Arrhenius * F_network / (F_network + K_D)
end
```

where K_D = 20 × F_i_min (matching MATLAB's 10 × (Fi_min + Fn_min)).

### Physical consequence
Without network dependence, F_m diffuses uniformly across the domain regardless 
of where hyphae exist, preventing the development of radial fungal gradients. 
With the correction, fungi can only spread into territory where their own network 
already extends, creating a propagating front that matches the biological picture 
of hyphal extension followed by translocation.

### Files changed
- `src/physics/effective_diffusion.jl`: D_eff_fungi_mobile() now takes F_n, F_i, K_D
- `src/physics/effective_diffusion.jl`: update_effective_diffusion!() now takes state
- `src/solver/timestepper.jl`: passes state to update_effective_diffusion!()

---

## Deviation 2: Insulation as Splitting Fraction vs. Independent Sink

### Falconer (2005, 2008) — Insulation as splitting fraction
In Falconer's formulation (Box 1, 2008 paper), ζ is a **splitting fraction** 
applied to the entire non-insulated biomass tendency. A fraction ζ of all new 
b_n production is redirected to b_i:

$$\frac{\partial b_i}{\partial t} = \zeta \frac{\partial b_n}{\partial t} + g(\alpha_i p^q - \beta_i p) b_i$$

$$\frac{\partial b_n}{\partial t} = (1-\zeta) \left[\frac{\partial}{\partial x}\left(D_b \frac{\partial b_n}{\partial x}\right) + g(\alpha_n p^q - \beta_n p) b_n \right]$$

Insulation is not a separate loss term — it is a fraction of the total F_n 
production (including diffusion) that matures into F_i. This means:
- If F_n is growing, a fraction ζ of that growth becomes F_i
- If F_n is shrinking (mobilization), a fraction ζ of that loss is also 
  redirected — insulation tracks the net tendency

### MATLAB (single_aggregate_beta.m, lines 389-390, 404-406)
```matlab
immobil_i = alpha_i * rel_Fm * Fi + zeta * alpha_n * rel_Fm * Fn;
immobil_n = (1-zeta) * alpha_n * rel_Fm * Fn;
...
Fi_rate = immobil_i - death_Fi;
Fn_rate = immobil_n - death_Fn;
Fm_rate = RF * epsilon_F - immobil_i - immobil_n;
```

ζ splits the immobilization flux from F_m → F_n: fraction ζ goes to F_i 
instead, fraction (1-ζ) stays as F_n. No separate insulation drain.

### Julia (before correction)
```julia
# In fungal_transitions():
insulation = ζ_T * F_n                               # independent drain!

# In reactions.jl:
S_Fn = trans_n - insulation                           # trans_n minus ζ·F_n
S_Fi = insulation + trans_i - R_rec_F                 # gains ζ·F_n
```

Insulation was implemented as an **independent first-order sink** on F_n with 
rate ζ·F_n [μg/mm³/day]. This is fundamentally different from Falconer's 
splitting fraction.

### Quantitative impact
With ζ = 0.2/day:
- **Julia (incorrect)**: F_n has an unconditional half-life of ~3.5 days from 
  insulation alone, regardless of growth. Even if transitions feed F_n, the 
  insulation drain (0.2 × F_n) overwhelms the supply 
  (η·β_n·Π·F_n ≈ 0.0004·F_n at t=0), causing F_n → 0 within days.
- **Falconer/MATLAB (correct)**: F_n only loses the fraction ζ of its net 
  production to F_i. If production is zero, insulation is zero. F_n persists 
  at its initial value until transitions begin.

### Julia (corrected)
The correction requires ζ to act as a splitting fraction on the entire F_n 
tendency, matching Falconer Box 1. This affects both the reaction step (source 
terms) and the diffusion step (ζ multiplies the diffusion contribution too).

**Reaction source terms:**
```julia
S_Fn = (1-ζ) * trans_n
S_Fi = ζ * trans_n + trans_i - R_rec_F
```

**Diffusion step:**
The Crank-Nicolson solve for F_n should produce the full tendency, then split:
- (1-ζ) fraction stays as F_n
- ζ fraction is added to F_i

### Files changed
- `src/physics/fungi.jl`: fungal_transitions() — remove insulation as separate term
- `src/physics/reactions.jl`: S_Fn, S_Fi — use ζ as splitting fraction
- `src/solver/diffusion_step.jl`: apply ζ splitting after F_n diffusion solve

---

## Summary

| Aspect | Falconer/MATLAB | Julia (before) | Julia (corrected) |
|--------|----------------|----------------|-------------------|
| D_Fm | Network-dependent | Spatially uniform | Network-dependent (Michaelis-Menten) |
| Insulation (ζ) | Splitting fraction of ∂b_n/∂t | Independent sink ζ·F_n | Splitting fraction of ∂F_n/∂t |

Both deviations suppressed fungal spatial gradients:
1. Uniform D_Fm prevented F_m from concentrating near the POM
2. Independent insulation drained F_n to extinction before gradients could develop

Together, they explain why the Julia model showed spatially uniform, declining 
fungal profiles while the MATLAB model produced the expected radial gradients 
with fungal accumulation near the POM surface.
